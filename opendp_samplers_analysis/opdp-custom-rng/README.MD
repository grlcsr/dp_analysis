# Using Custom File-Based Randomness in OpenDP

This modified version of OpenDP uses random bytes from binary files on disk instead of system entropy sources. This is useful for reproducible experiments, testing, or using pre-generated randomness datasets.

## Setup

1. **Add dependency to Cargo.toml**

   In `opendp/Cargo.toml`, add:
   ```toml
   [dependencies]
   once_cell = "1.19"
   ```

2. **Replace the samplers module**

   Replace mod.rs with the modified version provided. This file now includes:
   - `get_random_bytes_from_buffer()`: Reads from an in-memory buffer, refilling from disk as needed
   - `refill_buffer_from_bin_files()`: Sequentially reads `.bin` files from `./bins/` directory
   - Per-file offset tracking so files are never deleted or truncated
   - `RANDOMNESS_CONSUMED` static counter to track total bytes consumed

## How It Works

- The code reads binary files from a `./bins/` directory in sorted order
- Files are read sequentially: the first file is fully consumed before moving to the second
- Each file maintains an in-memory offset, so subsequent calls resume where they left off
- Files are never modified or deleted, and offsets are lost only when the process restarts
- If you need to reset and re-read files, restart the application

## Using It

In your OpenDP-based program:

```rust
use opendp::traits::samplers::RANDOMNESS_CONSUMED;
use std::sync::atomic::Ordering;

fn main() {
    // Your OpenDP code here (it will now use randomness from ./bins/*.bin files)
    
    // Track how much randomness was consumed:
    let consumed = RANDOMNESS_CONSUMED.load(Ordering::Relaxed);
    println!("Total bytes consumed: {}", consumed);
	
	// Track randomness consumed for each iteration:
	let consumed = RANDOMNESS_CONSUMED.swap(0, Ordering::Relaxed);
}
```

Run your program in a directory that has a `./bins/` subdirectory containing `.bin` files.

## Notes

- The in-memory buffer refills at 1 MB chunks by default (modify `BUFFER_SIZE` constant if needed)
- Offsets are tracked per file path and renaming files will treat them as new files
- Not thread-safe for concurrent reads across processes; single-process use only
- To return to system entropy, rebuild without this modified module or disable the custom `fill_bytes()` implementation