// Shared OpenSSL RNG wrapper 
// The rand crate fallback is kept commented in case you want to switch (Windows users may have trouble with OpenSSL)

pub type Fallible<T> = Result<T, String>;

pub fn fill_bytes(buffer: &mut [u8]) -> Fallible<()> {
    /*use rand::Rng;
    rand::thread_rng()
        .try_fill(buffer)
        .map_err(|e| format!("RNG error: {:?}", e))?;*/
    use openssl::rand::rand_bytes;

    rand_bytes(buffer).map_err(|e| format!("OpenSSL error: {:?}", e))?;
    Ok(())
}