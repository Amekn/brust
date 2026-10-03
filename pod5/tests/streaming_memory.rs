//! Streaming signal reads must not keep every decoded read in memory.
//!
//! This file is its own test binary so the counting allocator sees only this
//! test's allocations.

use brust_pod5::Pod5Reader;
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicIsize, Ordering};

const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/A_100.pod5");

static LIVE_BYTES: AtomicIsize = AtomicIsize::new(0);

struct CountingAllocator;

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc(layout) };
        if !pointer.is_null() {
            LIVE_BYTES.fetch_add(layout.size() as isize, Ordering::Relaxed);
        }
        pointer
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        unsafe { System.dealloc(pointer, layout) };
        LIVE_BYTES.fetch_sub(layout.size() as isize, Ordering::Relaxed);
    }
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

#[test]
fn streaming_signal_reads_do_not_retain_decoded_samples() {
    let mut reader = Pod5Reader::from_path(FIXTURE).unwrap();
    let before = LIVE_BYTES.load(Ordering::Relaxed);
    let mut decoded_bytes = 0;
    let mut largest_read = 0;
    while let Some(record) = reader.read_record().unwrap() {
        let samples = reader.signal_for_record(&record).unwrap();
        let bytes = samples.len() * size_of::<i16>();
        decoded_bytes += bytes;
        largest_read = largest_read.max(bytes);
    }
    let retained = (LIVE_BYTES.load(Ordering::Relaxed) - before) as usize;

    // About 2.25 MB decodes in total; keeping one read's samples is fine.
    assert!(decoded_bytes > 2_000_000);
    assert!(
        retained < largest_read + 256 * 1024,
        "reader kept {retained} bytes after decoding {decoded_bytes}"
    );
}
