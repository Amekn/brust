//! Writing a POD5 must not hold extra copies of the whole signal payload.
//!
//! This file is its own test binary so the counting allocator sees only this
//! test's allocations.

mod common;

use brust_pod5::Pod5Writer;
use std::alloc::{GlobalAlloc, Layout, System};
use std::io;
use std::sync::atomic::{AtomicUsize, Ordering};

static LIVE_BYTES: AtomicUsize = AtomicUsize::new(0);
static PEAK_BYTES: AtomicUsize = AtomicUsize::new(0);

struct CountingAllocator;

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc(layout) };
        if !pointer.is_null() {
            let live = LIVE_BYTES.fetch_add(layout.size(), Ordering::Relaxed) + layout.size();
            PEAK_BYTES.fetch_max(live, Ordering::Relaxed);
        }
        pointer
    }

    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        unsafe { System.dealloc(pointer, layout) };
        LIVE_BYTES.fetch_sub(layout.size(), Ordering::Relaxed);
    }
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

#[test]
fn writing_streams_signal_batches_instead_of_copying_the_payload() {
    let pod5 = common::pod5_with_reads(2_000);
    let payload: usize = pod5
        .signals
        .iter()
        .filter_map(|signal| signal.compressed_bytes())
        .map(<[u8]>::len)
        .sum();

    let before = LIVE_BYTES.load(Ordering::Relaxed);
    PEAK_BYTES.store(before, Ordering::Relaxed);
    let mut writer = Pod5Writer::from_writer(io::sink());
    writer.write_all(&pod5).unwrap();
    let extra = PEAK_BYTES.load(Ordering::Relaxed) - before;

    // The whole-file encoder held about three copies of the payload.
    assert!(payload > 15_000_000, "{payload}");
    assert!(
        extra < payload / 4,
        "writing a {payload}-byte payload used {extra} more bytes at its peak"
    );
}
