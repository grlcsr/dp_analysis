//! Traits for sampling from probability distributions.
//! 
//! //! Traits for sampling from probability distributions.
//!
//! ## Custom Randomness Source
//! This module includes functions to read randomness from binary files on disk
//! instead of system entropy sources. The randomness is consumed sequentially:
//! bytes are read from files in sorted order, with per-file offsets tracked in memory.
//! Files are not modified or deleted; reading resumes from the last offset on the next refill.
//!
//! - `get_random_bytes_from_buffer()`: Extracts random bytes from an in-memory buffer,
//!   refilling it from disk files if needed.
//! - `refill_buffer_from_bin_files()`: Reads sequential chunks from .bin files in ./bins/,
//!   advancing per-file offsets and preserving files for repeated access.
//! 
//! Required: add once_cell = "1.19" to opendp's Cargo.toml dependencies.
//! With this you can also track the total amount of randomness consumed via
//! `RANDOMNESS_CONSUMED` static variable on your main program.

#[cfg(test)]
pub(crate) mod test;

mod bernoulli;
pub use bernoulli::*;

mod cks20;
pub use cks20::*;

mod geometric;
pub use geometric::*;

mod psrn;
pub use psrn::*;

mod uniform;
pub use uniform::*;

use rand::prelude::SliceRandom;
use rand::RngCore;

use crate::error::Fallible;

use once_cell::sync::Lazy;
use std::io::{Read, Seek, Write};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Mutex;
use std::collections::HashMap;
use std::path::PathBuf;

const BUFFER_SIZE: usize = 1024 * 1024;
static RANDOMNESS_BUFFER: Lazy<Mutex<Vec<u8>>> = Lazy::new(|| Mutex::new(Vec::new()));
pub static RANDOMNESS_CONSUMED: Lazy<AtomicUsize> = Lazy::new(|| AtomicUsize::new(0));

// Track per-file read offsets so files are not truncated/removed.
static FILE_READ_OFFSETS: Lazy<Mutex<HashMap<PathBuf, u64>>> =
    Lazy::new(|| Mutex::new(HashMap::new()));

fn get_random_bytes_from_buffer(request: &mut [u8]) -> std::io::Result<()> {
    let mut buffer = RANDOMNESS_BUFFER.lock().unwrap();

    if buffer.len() < request.len() {
        refill_buffer_from_bin_files(&mut buffer)?;
    }

    if buffer.len() < request.len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            "Not enough randomness from files",
        ));
    }

    request.copy_from_slice(&buffer[..request.len()]);
    buffer.drain(..request.len());

    RANDOMNESS_CONSUMED.fetch_add(request.len(), Ordering::Relaxed);

    Ok(())
}

fn refill_buffer_from_bin_files(buffer: &mut Vec<u8>) -> std::io::Result<()> {
    let bins_path = std::path::Path::new("./bins");
    let mut bytes_needed = BUFFER_SIZE - buffer.len();

    while bytes_needed > 0 {
        let mut entries: Vec<_> = std::fs::read_dir(bins_path)?
            .filter_map(Result::ok)
            .filter(|e| {
                e.path().extension().map_or(false, |ext| ext == "bin")
                    && e.file_type().map(|ft| ft.is_file()).unwrap_or(false)
            })
            .collect();

        entries.sort_by_key(|e| e.path());

        if entries.is_empty() {
            break;
        }

        // iterate entries and read from each according to saved offsets until we
        // satisfy bytes_needed or run out of readable bytes.
        let mut offsets = FILE_READ_OFFSETS.lock().unwrap();
        let mut read_any = false;
        for entry in &entries {
            if bytes_needed == 0 {
                break;
            }
            let path = entry.path();
            let metadata = match entry.metadata() {
                Ok(m) => m,
                Err(_) => continue,
            };
            let file_size = metadata.len();
            let offset = offsets.entry(path.clone()).or_insert(0u64);

            if *offset >= file_size {
                // no new bytes available in this file
                continue;
            }

            let mut file = std::fs::OpenOptions::new().read(true).open(&path)?;
            file.seek(std::io::SeekFrom::Start(*offset))?;

            let bytes_available = (file_size - *offset) as usize;
            let bytes_to_read = bytes_needed.min(bytes_available);
            let mut temp_buffer = vec![0u8; bytes_to_read];
            let bytes_read = file.read(&mut temp_buffer)?;

            if bytes_read == 0 {
                // nothing to read currently from this file
                continue;
            }

            buffer.extend_from_slice(&temp_buffer[..bytes_read]);
            bytes_needed -= bytes_read;
            *offset += bytes_read as u64;
            read_any = true;
        }

        // If we couldn't read anything from any file, break to avoid busy loop.
        if !read_any {
            break;
        }
    }

    Ok(())
}

///
/// # Proof Definition
/// For any input `buffer`, fill the `buffer` with random bits, where each bit is an iid draw from Bernoulli(p=0.5).
/// Return `Err(e)` if there is insufficient system entropy, otherwise return `Ok(())`.
#[cfg(feature = "use-openssl")]
pub fn fill_bytes(buffer: &mut [u8]) -> Fallible<()> {
    use openssl::rand::rand_bytes;
    if let Err(e) = rand_bytes(buffer) {
        return fallible!(FailedFunction, "OpenSSL error: {:?}", e);
    } else {
        RANDOMNESS_CONSUMED.fetch_add(buffer.len(), Ordering::Relaxed);
        return Ok(());
    }
    
    if let Err(e) = get_random_bytes_from_buffer(buffer) {
        fallible!(FailedFunction, "OpenSSL error: {:?}", e)
    } else {
        Ok(())
    }
}

/// Non-securely fill a byte buffer with random bits.
///
/// Enable `use-openssl` for a secure implementation.
#[cfg(not(feature = "use-openssl"))]
pub fn fill_bytes(buffer: &mut [u8]) -> Fallible<()> {
    /*use rand::Rng;
    if let Err(e) = rand::thread_rng().try_fill(buffer) {
        return fallible!(FailedFunction, "Rand error: {:?}", e);
    } else {
        RANDOMNESS_CONSUMED.fetch_add(buffer.len(), Ordering::Relaxed);
        return Ok(());
    }*/

    if let Err(e) = get_random_bytes_from_buffer(buffer) {
        fallible!(FailedFunction, "OpenSSL error: {:?}", e)
    } else {
        Ok(())
    }
}
/// An OpenDP random number generator that implements [`rand::RngCore`].
pub(crate) struct GeneratorOpenDP {
    /// If an error happens while sampling, it is packed into this struct and thrown later.
    pub error: Fallible<()>,
}

impl GeneratorOpenDP {
    pub fn new() -> Self {
        GeneratorOpenDP { error: Ok(()) }
    }
}
impl Default for GeneratorOpenDP {
    fn default() -> Self {
        Self::new()
    }
}

impl RngCore for GeneratorOpenDP {
    fn next_u32(&mut self) -> u32 {
        let mut buffer = [0u8; 4];
        self.fill_bytes(&mut buffer);
        u32::from_ne_bytes(buffer)
    }

    fn next_u64(&mut self) -> u64 {
        let mut buffer = [0u8; 8];
        self.fill_bytes(&mut buffer);
        u64::from_ne_bytes(buffer)
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        if let Err(e) = fill_bytes(dest) {
            self.error = Err(e)
        }
    }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand::Error> {
        fill_bytes(dest).map_err(rand::Error::new)
    }
}

/// Shuffle a mutable reference to a collection.
pub trait Shuffle {
    /// # Proof Definition
    /// For any input `self` of type `Self`,
    /// mutate `self` such that the elements within are ordered randomly.
    /// Returns `Err(e)` if there is insufficient system entropy,
    /// or `Ok(())` otherwise.
    fn shuffle(&mut self) -> Fallible<()>;
}
impl<T> Shuffle for Vec<T> {
    fn shuffle(&mut self) -> Fallible<()> {
        let mut rng = GeneratorOpenDP::new();
        SliceRandom::shuffle(self.as_mut_slice(), &mut rng);
        rng.error
    }
}
