//! BGZF block framing.
//!
//! BGZF is a series of independent gzip members, each holding a limited amount
//! of uncompressed data, so a reader can jump to a block and decompress it
//! without reading the rest of the file. This module has:
//!
//! - [`BgzfReader`], which reads a BGZF stream and reports
//!   [`BgzfVirtualOffset`]s.
//! - [`BgzfWriter`], which writes any byte stream as BGZF blocks and ends it
//!   with [`EOF_BLOCK`].
//! - [`compress_block`], which frames one block of data as a BGZF block.
//! - [`EOF_BLOCK`], the empty block that ends a BGZF stream.

use super::{invalid_data, read_exact_or_eof};
use flate2::read::DeflateDecoder;
use flate2::write::DeflateEncoder;
use flate2::{Compression, Crc};
use std::io::{self, Read, Write};
use std::mem;
use std::panic::{RefUnwindSafe, UnwindSafe};
use std::sync::mpsc::{self, Receiver, SyncSender, TryRecvError};
use std::sync::{Mutex, PoisonError};
use std::thread::{self, JoinHandle};

/// Most uncompressed bytes [`compress_block`] accepts for one block (65,280).
///
/// Chosen so a framed block stays within 64 KiB even when the data does not
/// compress.
pub const MAX_BLOCK_DATA: usize = 64 * 1024 - 256;

/// Largest uncompressed size a BGZF block may hold (the SAM spec's 64 KiB).
const MAX_BLOCK_UNCOMPRESSED: usize = 64 * 1024;

/// The standard 28-byte BGZF EOF marker: an empty BGZF block that ends a stream.
pub const EOF_BLOCK: [u8; 28] = *b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00";

/// Packed BGZF virtual offset (`compressed_block_offset << 16 | block_offset`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BgzfVirtualOffset(u64);

impl BgzfVirtualOffset {
    /// Creates a virtual offset from a compressed block start and uncompressed
    /// offset within that block.
    pub fn new(compressed_offset: u64, uncompressed_offset: u16) -> Self {
        Self((compressed_offset << 16) | u64::from(uncompressed_offset))
    }

    /// Creates a virtual offset from its packed `u64` representation.
    pub fn from_raw(raw: u64) -> Self {
        Self(raw)
    }

    /// Returns the packed `u64` representation.
    pub fn raw(self) -> u64 {
        self.0
    }

    /// Returns the compressed BGZF block start offset.
    pub fn compressed_offset(self) -> u64 {
        self.0 >> 16
    }

    /// Returns the uncompressed offset inside the BGZF block.
    pub fn uncompressed_offset(self) -> u16 {
        (self.0 & 0xffff) as u16
    }
}

/// Reader for BGZF blocks with virtual-offset tracking.
pub struct BgzfReader<R: Read> {
    inner: R,
    buffer: Vec<u8>,
    position: usize,
    current_block_start: u64,
    next_block_start: u64,
    eof: bool,
    require_eof_block: bool,
    last_block_is_eof_marker: bool,
    last_read_failed: bool,
}

/// Diagnostic for a strict reader that ran out of input without the EOF block.
const MISSING_EOF_BLOCK: &str =
    "BGZF stream ended without the EOF block; the file may be truncated";

impl<R: Read> BgzfReader<R> {
    /// Creates a BGZF reader over a compressed byte stream.
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            buffer: Vec::new(),
            position: 0,
            current_block_start: 0,
            next_block_start: 0,
            eof: false,
            require_eof_block: false,
            last_block_is_eof_marker: false,
            last_read_failed: false,
        }
    }

    /// Sets whether the stream must end with the standard BGZF EOF block.
    ///
    /// The default is `false`: a stream that ends without the marker is read
    /// as complete, as before. When `true`, reaching the end of the input is an
    /// [`io::ErrorKind::InvalidData`] error unless the last block read is
    /// byte-for-byte the 28-byte [`EOF_BLOCK`]. An empty input fails too. A
    /// final empty block that differs in any byte, such as its MTIME, is not
    /// accepted.
    ///
    /// The check happens at the end of the stream and is reported once; later
    /// reads return `Ok(0)`, as they do after a lenient end. Set this before
    /// the reader reaches the end, because a lenient end is cached and is not
    /// re-checked. Empty blocks in the middle of the stream are still skipped,
    /// and only the last block has to be the marker.
    ///
    /// A block that fails to read, such as one cut short, is reported by its
    /// own error, not also as a missing EOF block.
    ///
    /// A joined stream cut exactly after an interior EOF marker ends with a
    /// valid marker, so this check cannot detect that truncation.
    pub fn set_require_eof_block(&mut self, require: bool) {
        self.require_eof_block = require;
    }

    /// Returns the current virtual offset.
    pub fn virtual_offset(&self) -> BgzfVirtualOffset {
        if self.position >= self.buffer.len() {
            BgzfVirtualOffset::new(self.next_block_start, 0)
        } else {
            BgzfVirtualOffset::new(self.current_block_start, self.position as u16)
        }
    }

    /// Consumes this reader and returns the wrapped stream.
    pub fn into_inner(self) -> R {
        self.inner
    }

    fn fill_block(&mut self) -> io::Result<bool> {
        if self.eof {
            return Ok(false);
        }

        loop {
            let block = match read_bgzf_block(&mut self.inner, self.next_block_start) {
                Ok(block) => block,
                Err(error) => {
                    // A failed block is reported once, by its own error. A
                    // transient error is not a failed block: a retry that then
                    // reaches the end must still check for the EOF block.
                    if !matches!(
                        error.kind(),
                        io::ErrorKind::WouldBlock
                            | io::ErrorKind::Interrupted
                            | io::ErrorKind::TimedOut
                    ) {
                        self.last_read_failed = true;
                    }
                    return Err(error);
                }
            };
            let Some(block) = block else {
                self.eof = true;
                self.buffer.clear();
                self.position = 0;
                if self.require_eof_block
                    && !self.last_block_is_eof_marker
                    && !self.last_read_failed
                {
                    return Err(invalid_data(MISSING_EOF_BLOCK));
                }
                return Ok(false);
            };

            self.last_block_is_eof_marker = block.is_eof_marker;
            self.last_read_failed = false;

            self.current_block_start = self.next_block_start;
            self.next_block_start = self
                .next_block_start
                .checked_add(block.compressed_size as u64)
                .ok_or_else(|| invalid_data("BGZF compressed offset overflow"))?;
            self.buffer = block.uncompressed;
            self.position = 0;

            if !self.buffer.is_empty() {
                return Ok(true);
            }
        }
    }
}

impl<R: Read> Read for BgzfReader<R> {
    fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
        if output.is_empty() {
            return Ok(0);
        }

        if self.position >= self.buffer.len() && !self.fill_block()? {
            return Ok(0);
        }

        let available = self.buffer.len() - self.position;
        let take = available.min(output.len());
        output[..take].copy_from_slice(&self.buffer[self.position..self.position + take]);
        self.position += take;
        Ok(take)
    }
}

struct BgzfBlock {
    compressed_size: usize,
    uncompressed: Vec<u8>,
    /// Whether the raw block bytes equal [`EOF_BLOCK`] exactly.
    is_eof_marker: bool,
}

fn read_bgzf_block<R: Read>(
    reader: &mut R,
    compressed_offset: u64,
) -> io::Result<Option<BgzfBlock>> {
    let mut prefix = [0u8; 12];
    if !read_exact_or_eof(reader, &mut prefix, "BGZF block header")? {
        return Ok(None);
    }

    if prefix[..4] != [0x1f, 0x8b, 0x08, 0x04] {
        return Err(invalid_data("invalid BGZF gzip header"));
    }

    let xlen = u16::from_le_bytes([prefix[10], prefix[11]]) as usize;
    if xlen < 6 {
        return Err(invalid_data("BGZF extra field is too short"));
    }

    let mut extra = vec![0u8; xlen];
    reader.read_exact(&mut extra)?;
    let bsize = bgzf_bsize(&extra)?;
    let compressed_size = usize::from(bsize) + 1;
    let header_size = 12usize
        .checked_add(xlen)
        .ok_or_else(|| invalid_data("BGZF header size overflow"))?;
    if compressed_size < header_size + 8 {
        return Err(invalid_data("BGZF block size is too small"));
    }

    let remaining_len = compressed_size - header_size;
    let mut remaining = vec![0u8; remaining_len];
    reader.read_exact(&mut remaining)?;

    let is_eof_marker = compressed_size == EOF_BLOCK.len()
        && prefix[..] == EOF_BLOCK[..12]
        && extra[..] == EOF_BLOCK[12..12 + xlen]
        && remaining[..] == EOF_BLOCK[12 + xlen..];

    let deflate_len = remaining_len - 8;
    let deflate = &remaining[..deflate_len];
    let footer = &remaining[deflate_len..];
    let expected_crc = u32::from_le_bytes(footer[..4].try_into().unwrap());
    let expected_isize = u32::from_le_bytes(footer[4..8].try_into().unwrap());

    // A BGZF block holds at most 64 KiB, which virtual offsets rely on, so a
    // larger ISIZE is invalid and decompression stops just past the limit.
    if expected_isize as usize > MAX_BLOCK_UNCOMPRESSED {
        return Err(invalid_data("BGZF ISIZE exceeds the 64 KiB block limit"));
    }
    let mut decoder = DeflateDecoder::new(deflate);
    let mut uncompressed = Vec::with_capacity(expected_isize as usize);
    (&mut decoder)
        .take(MAX_BLOCK_UNCOMPRESSED as u64 + 1)
        .read_to_end(&mut uncompressed)
        .map_err(|error| invalid_data(format!("BGZF block has corrupt deflate data: {error}")))?;
    if uncompressed.len() != expected_isize as usize {
        return Err(invalid_data(
            "BGZF ISIZE does not match decompressed length",
        ));
    }
    if decoder.total_in() != deflate.len() as u64 {
        return Err(invalid_data(
            "BGZF block has bytes after its deflate stream",
        ));
    }

    let mut crc = Crc::new();
    crc.update(&uncompressed);
    if crc.sum() != expected_crc {
        return Err(invalid_data(format!(
            "BGZF CRC mismatch at compressed offset {compressed_offset}"
        )));
    }

    Ok(Some(BgzfBlock {
        compressed_size,
        uncompressed,
        is_eof_marker,
    }))
}

fn bgzf_bsize(extra: &[u8]) -> io::Result<u16> {
    let mut offset = 0usize;
    while offset + 4 <= extra.len() {
        let si1 = extra[offset];
        let si2 = extra[offset + 1];
        let slen = u16::from_le_bytes([extra[offset + 2], extra[offset + 3]]) as usize;
        offset += 4;
        let end = offset
            .checked_add(slen)
            .ok_or_else(|| invalid_data("BGZF extra subfield length overflow"))?;
        if end > extra.len() {
            return Err(invalid_data("BGZF extra subfield is truncated"));
        }
        if si1 == b'B' && si2 == b'C' {
            if slen != 2 {
                return Err(invalid_data("BGZF BC subfield has invalid length"));
            }
            return Ok(u16::from_le_bytes([extra[offset], extra[offset + 1]]));
        }
        offset = end;
    }

    Err(invalid_data("BGZF BC subfield is missing"))
}

/// Compresses `data` as one BGZF block and appends the framed block to `out`.
///
/// Returns [`io::ErrorKind::InvalidInput`], leaving `out` unchanged, when `data`
/// is longer than [`MAX_BLOCK_DATA`] or the framed block would exceed 64 KiB.
pub fn compress_block(data: &[u8], out: &mut Vec<u8>) -> io::Result<()> {
    if data.len() > MAX_BLOCK_DATA {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "BGZF block data is {} bytes, more than the {MAX_BLOCK_DATA} byte limit",
                data.len()
            ),
        ));
    }

    let mut encoder = DeflateEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data)?;
    let compressed = encoder.finish()?;
    let block_size = 18usize
        .checked_add(compressed.len())
        .and_then(|size| size.checked_add(8))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "BGZF block size overflow"))?;

    if block_size > 64 * 1024 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "BGZF compressed block exceeds 64 KiB",
        ));
    }

    let bsize = u16::try_from(block_size - 1)
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "BGZF block size exceeds u16"))?;
    let mut crc = Crc::new();
    crc.update(data);

    out.extend_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C', 0x02,
        0x00,
    ]);
    out.extend_from_slice(&bsize.to_le_bytes());
    out.extend_from_slice(&compressed);
    out.extend_from_slice(&crc.sum().to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    Ok(())
}

/// Most worker threads [`BgzfWriter::with_threads`] starts (256).
///
/// Larger thread counts are clamped to this value. It limits threads and
/// memory only; output does not depend on it.
pub const MAX_THREADS: usize = 256;

/// Message for every call on a [`BgzfWriter`] after one of its calls failed.
const FAILED_EARLIER: &str = "BGZF writer failed earlier";

/// Message for a worker thread that stopped before returning its blocks.
const WORKER_STOPPED: &str = "BGZF worker stopped";

/// Most blocks in flight per worker, and the capacity of each worker's job and
/// result channels.
const PER_WORKER: usize = 2;

/// Writer that compresses a byte stream into BGZF blocks.
///
/// Bytes are buffered until [`MAX_BLOCK_DATA`] are held, then written to the
/// wrapped writer as one BGZF block. [`Write::flush`] writes any partial block
/// early and flushes the wrapped writer. [`BgzfWriter::finish`] writes the
/// last partial block and [`EOF_BLOCK`], and returns the wrapped writer.
///
/// [`BgzfWriter::new`] compresses on the calling thread.
/// [`BgzfWriter::with_threads`] compresses on worker threads and still writes
/// the blocks in order.
///
/// Block boundaries depend only on the bytes written and on where
/// [`Write::flush`] is called, so the same calls always produce the same
/// bytes, whatever the thread count. A block that is only partly full is
/// written only by `flush` or `finish`, so avoid calling `flush` after every
/// small write: each call produces a small, poorly compressed block.
///
/// Deflate output comes from flate2's backend. Bytes match within a build,
/// but switching flate2 to another backend (for example zlib-ng) changes them.
/// Output is still identical for every thread count.
///
/// Dropping a `BgzfWriter` without calling `finish` loses any buffered bytes
/// and leaves the stream without its EOF block, so `finish` is required for a
/// complete stream. A dropped threaded writer stops and joins its workers.
///
/// If any call returns an error, the stream is incomplete, so every later
/// [`Write::write`], [`Write::flush`] or [`BgzfWriter::finish`] returns an
/// error too. The call that hit the original error returns that error. The
/// same holds after a call that panicked, for example because the wrapped
/// writer panicked and the caller caught the panic.
pub struct BgzfWriter<W: Write> {
    inner: W,
    /// Uncompressed bytes not yet written as a block; never longer than
    /// [`MAX_BLOCK_DATA`].
    pending: Vec<u8>,
    /// Scratch buffer that holds one framed block on its way to `inner`; used
    /// only when compressing inline.
    framed: Vec<u8>,
    /// Worker threads that compress blocks; `None` compresses inline.
    pool: Option<Pool>,
    /// Set once any call has returned an error or panicked, and while a call
    /// runs.
    failed: bool,
}

impl<W: Write> BgzfWriter<W> {
    /// Creates a BGZF writer that writes compressed blocks to `inner`.
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            pending: Vec::with_capacity(MAX_BLOCK_DATA),
            framed: Vec::new(),
            pool: None,
            failed: false,
        }
    }

    /// Creates a BGZF writer that compresses blocks on `threads` worker
    /// threads and writes them to `inner` in stream order.
    ///
    /// `threads` of 0 or 1 compresses on the calling thread, exactly like
    /// [`BgzfWriter::new`]. Larger values start `threads` workers, clamped to
    /// [`MAX_THREADS`]. If the system cannot start them all, the writer uses
    /// the workers it did start, or compresses inline if it started none.
    ///
    /// Workers only compress. Every write to `inner` happens on the thread
    /// that calls `write`, `flush` or `finish`, so `inner` need not be
    /// [`Send`]. At most two blocks per worker are in flight at once, which
    /// bounds memory use.
    ///
    /// Output does not depend on the thread count: given the same sequence of
    /// `write`, `flush` and `finish` calls, every thread count gives the same
    /// bytes as [`BgzfWriter::new`]. Deflate output comes from flate2's
    /// backend, so bytes match within a build; another flate2 backend (for
    /// example zlib-ng) gives different bytes, still identical across thread
    /// counts.
    ///
    /// Call [`BgzfWriter::finish`] to complete the stream. Dropping the
    /// writer stops and joins its workers, but discards blocks not yet
    /// written and leaves the stream without its EOF block.
    pub fn with_threads(inner: W, threads: usize) -> Self {
        let mut writer = Self::new(inner);
        if threads > 1 {
            writer.pool = Pool::start(threads.min(MAX_THREADS));
        }
        writer
    }

    /// Writes the last partial block and [`EOF_BLOCK`], flushes the wrapped
    /// writer, and returns it.
    ///
    /// A threaded writer first waits for every block its workers hold, and
    /// joins the workers once the EOF block is written.
    ///
    /// Returns an error if an earlier call on this writer failed, if writing
    /// or flushing fails now, or if a worker thread stopped unexpectedly.
    pub fn finish(mut self) -> io::Result<W> {
        // `self` is consumed, so there is nothing left to poison on error.
        // Returning early drops `self`, which stops and joins any workers.
        self.ensure_not_failed()?;
        self.write_pending_partial()?;
        self.write_in_flight()?;
        self.inner.write_all(&EOF_BLOCK)?;
        self.inner.flush()?;
        if let Some(mut pool) = self.pool.take() {
            pool.shut_down()?;
        }
        Ok(self.inner)
    }

    fn ensure_not_failed(&self) -> io::Result<()> {
        if self.failed {
            Err(io::Error::other(FAILED_EARLIER))
        } else {
            Ok(())
        }
    }

    /// Runs one `write` or `flush` call, leaving the writer failed unless the
    /// call returns `Ok`.
    ///
    /// If an earlier call failed, returns the poison error without running
    /// `call`. `failed` stays set while `call` runs, so a call that panics part
    /// way also leaves the writer failed, even if the caller catches the panic:
    /// the buffered block and the blocks in flight may no longer match the
    /// stream.
    fn guarded<T>(&mut self, call: impl FnOnce(&mut Self) -> io::Result<T>) -> io::Result<T> {
        self.ensure_not_failed()?;
        self.failed = true;
        let result = call(self);
        self.failed = result.is_err();
        result
    }

    /// Copies `data` into `pending`, writing a block each time it fills up.
    fn buffer(&mut self, mut data: &[u8]) -> io::Result<()> {
        while !data.is_empty() {
            let take = (MAX_BLOCK_DATA - self.pending.len()).min(data.len());
            self.pending.extend_from_slice(&data[..take]);
            data = &data[take..];

            if self.pending.len() == MAX_BLOCK_DATA {
                self.write_pending()?;
            }
        }
        Ok(())
    }

    /// Writes `pending`, if it holds anything, as one block.
    fn write_pending_partial(&mut self) -> io::Result<()> {
        if self.pending.is_empty() {
            Ok(())
        } else {
            self.write_pending()
        }
    }

    /// Writes `pending` as one block and empties it.
    ///
    /// A threaded writer sends the block to a worker instead, then writes any
    /// finished blocks.
    fn write_pending(&mut self) -> io::Result<()> {
        match &mut self.pool {
            None => {
                self.framed.clear();
                compress_block(&self.pending, &mut self.framed)?;
                self.inner.write_all(&self.framed)?;
                self.pending.clear();
                Ok(())
            }
            Some(pool) => {
                let data = mem::replace(&mut self.pending, Vec::with_capacity(MAX_BLOCK_DATA));
                pool.send(data, &mut self.inner)
            }
        }
    }

    /// Waits for every block the workers hold and writes them to `inner`.
    fn write_in_flight(&mut self) -> io::Result<()> {
        match &mut self.pool {
            None => Ok(()),
            Some(pool) => pool.write_in_flight(&mut self.inner),
        }
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    /// Buffers all of `buf` and returns `buf.len()`.
    ///
    /// Each time [`MAX_BLOCK_DATA`] bytes have built up, they are written to the
    /// wrapped writer as one block. A threaded writer sends them to a worker
    /// instead and writes any blocks already compressed, in order. If two
    /// blocks per worker are already in flight, it waits for the oldest first.
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.guarded(|writer| writer.buffer(buf).map(|()| buf.len()))
    }

    /// Writes any buffered bytes as a block, then flushes the wrapped writer.
    ///
    /// A threaded writer waits for and writes every block in flight before it
    /// flushes. No block is added if no bytes are buffered.
    fn flush(&mut self) -> io::Result<()> {
        self.guarded(|writer| {
            writer
                .write_pending_partial()
                .and_then(|()| writer.write_in_flight())
                .and_then(|()| writer.inner.flush())
        })
    }
}

/// Worker threads that compress blocks for a threaded [`BgzfWriter`].
///
/// Block `k` goes to worker `k % workers`. Each worker compresses its jobs in
/// the order it gets them, so taking results from the workers in turn gives
/// the blocks in stream order, with no reordering.
///
/// The blocks in flight are the run of sequence numbers from `written` up to
/// `sent`, and `send` keeps that run at most `PER_WORKER × workers` long. So no
/// worker ever holds more than `PER_WORKER` blocks, which is the capacity of
/// both of its channels: a job send never blocks, and a worker never blocks
/// sending a result. The block the writer waits for is always the oldest
/// one its worker holds, so the wait always ends.
struct Pool {
    workers: Vec<Worker>,
    /// Blocks sent to workers so far.
    sent: u64,
    /// Blocks written to the inner writer so far; never more than `sent`.
    written: u64,
}

/// One compression thread and its channels.
struct Worker {
    jobs: SyncSender<Vec<u8>>,
    /// In a `Mutex` only so that `BgzfWriter` stays `Sync` (a `Receiver` is
    /// not). It is reached only through [`Mutex::get_mut`], which never locks.
    results: Mutex<Receiver<io::Result<Vec<u8>>>>,
    thread: JoinHandle<()>,
}

// `JoinHandle` is neither `UnwindSafe` nor `RefUnwindSafe`: the thread's
// result sits in an `UnsafeCell` inside it. Only the worker thread writes that
// cell, as it exits, and only `join` reads it, which consumes the handle. A
// panic on the writer's side cannot leave it half updated, and a shared
// `&Worker` cannot reach it. With these impls `BgzfWriter`, and `BamWriter`
// over it, stay as unwind safe as their inner writer, as they were before
// worker threads.
impl UnwindSafe for Worker {}
impl RefUnwindSafe for Worker {}

impl Worker {
    /// This worker's result receiver.
    fn results(&mut self) -> &Receiver<io::Result<Vec<u8>>> {
        self.results
            .get_mut()
            .unwrap_or_else(PoisonError::into_inner)
    }
}

impl Pool {
    /// Starts up to `count` workers. Returns `None` if none could start.
    fn start(count: usize) -> Option<Self> {
        let mut workers = Vec::with_capacity(count);
        for i in 0..count {
            let (jobs, job_queue) = mpsc::sync_channel(PER_WORKER);
            let (result_sender, results) = mpsc::sync_channel(PER_WORKER);
            let spawned = thread::Builder::new()
                .name(format!("bgzf-worker-{i}"))
                .spawn(move || compress_jobs(job_queue, result_sender));
            match spawned {
                Ok(thread) => workers.push(Worker {
                    jobs,
                    results: Mutex::new(results),
                    thread,
                }),
                // The workers already running give the same output.
                Err(_) => break,
            }
        }

        if workers.is_empty() {
            None
        } else {
            Some(Self {
                workers,
                sent: 0,
                written: 0,
            })
        }
    }

    /// Index of the worker that compresses block `block`.
    fn worker_for(&self, block: u64) -> usize {
        (block % self.workers.len() as u64) as usize
    }

    /// Sends `data` to the next worker as block `sent`, then writes any
    /// blocks that are ready.
    ///
    /// If the limit of blocks in flight is reached, first waits for the
    /// oldest block and writes it.
    fn send<W: Write>(&mut self, data: Vec<u8>, inner: &mut W) -> io::Result<()> {
        let limit = (PER_WORKER * self.workers.len()) as u64;
        if self.sent - self.written == limit {
            self.write_next(inner)?;
        }

        let index = self.worker_for(self.sent);
        self.workers[index]
            .jobs
            .send(data)
            .map_err(|_| io::Error::other(WORKER_STOPPED))?;
        self.sent += 1;
        self.write_ready(inner)
    }

    /// Waits for the oldest block in flight and writes it to `inner`.
    fn write_next<W: Write>(&mut self, inner: &mut W) -> io::Result<()> {
        let index = self.worker_for(self.written);
        let block = self.workers[index]
            .results()
            .recv()
            .map_err(|_| io::Error::other(WORKER_STOPPED))??;
        inner.write_all(&block)?;
        self.written += 1;
        Ok(())
    }

    /// Writes blocks to `inner`, oldest first, while the oldest is already
    /// compressed.
    fn write_ready<W: Write>(&mut self, inner: &mut W) -> io::Result<()> {
        while self.written < self.sent {
            let index = self.worker_for(self.written);
            let block = match self.workers[index].results().try_recv() {
                Ok(result) => result?,
                Err(TryRecvError::Empty) => break,
                Err(TryRecvError::Disconnected) => return Err(io::Error::other(WORKER_STOPPED)),
            };
            inner.write_all(&block)?;
            self.written += 1;
        }
        Ok(())
    }

    /// Waits for every block in flight and writes them to `inner` in order.
    fn write_in_flight<W: Write>(&mut self, inner: &mut W) -> io::Result<()> {
        while self.written < self.sent {
            self.write_next(inner)?;
        }
        Ok(())
    }

    /// Closes every worker's channels, then joins every worker thread.
    ///
    /// Returns an error if any worker panicked. Calling it again does nothing.
    fn shut_down(&mut self) -> io::Result<()> {
        // Dropping a worker's job sender ends its loop, and dropping its result
        // receiver ends it after the block it is on. All channels close before
        // the first join, so the workers stop together.
        let threads: Vec<JoinHandle<()>> =
            self.workers.drain(..).map(|worker| worker.thread).collect();

        let mut result = Ok(());
        for thread in threads {
            if thread.join().is_err() {
                result = Err(io::Error::other(WORKER_STOPPED));
            }
        }
        result
    }
}

impl Drop for Pool {
    /// Stops and joins the workers if [`BgzfWriter::finish`] did not.
    fn drop(&mut self) {
        // Nothing can report a worker panic from here, so it is ignored.
        let _ = self.shut_down();
    }
}

/// Worker loop: compresses each job and sends back the result, in order.
///
/// Returns when the job channel closes or the result receiver is dropped.
fn compress_jobs(jobs: Receiver<Vec<u8>>, results: SyncSender<io::Result<Vec<u8>>>) {
    for data in jobs {
        let mut out = Vec::new();
        let result = compress_block(&data, &mut out).map(|()| out);
        if results.send(result).is_err() {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Crc;
    use flate2::read::MultiGzDecoder;
    use std::fmt::Debug;
    use std::io::{self, Read, Write};
    use std::sync::mpsc;
    use std::time::Duration;
    use std::{panic, thread};

    fn crc32(data: &[u8]) -> u32 {
        let mut crc = Crc::new();
        crc.update(data);
        crc.sum()
    }

    fn decode_all(bytes: &[u8]) -> Vec<u8> {
        let mut decoded = Vec::new();
        MultiGzDecoder::new(bytes)
            .read_to_end(&mut decoded)
            .expect("MultiGzDecoder should decode the stream");
        decoded
    }

    /// Frames `data` as one BGZF block.
    fn block(data: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        compress_block(data, &mut out).expect("block should compress");
        out
    }

    /// Deterministic bytes that deflate cannot shrink much.
    fn noise(len: usize) -> Vec<u8> {
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (state >> 24) as u8
            })
            .collect()
    }

    /// Inner writer that accepts writes only until `fail_after` bytes are in.
    ///
    /// A write is accepted whole while `written + buf.len() <= fail_after`, and
    /// fails with `"injected failure"` otherwise.
    #[derive(Debug)]
    struct FailingWriter {
        written: usize,
        fail_after: usize,
    }

    impl Write for FailingWriter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if self.written + buf.len() <= self.fail_after {
                self.written += buf.len();
                Ok(buf.len())
            } else {
                Err(io::Error::other("injected failure"))
            }
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// Inner writer that keeps every byte written to it and records each
    /// `flush` call, so a test can tell whether, and when, it was flushed.
    #[derive(Debug, Default)]
    struct FlushRecorder {
        bytes: Vec<u8>,
        /// One entry per `flush` call: how many bytes had been written by then.
        flushed_at: Vec<usize>,
    }

    impl Write for FlushRecorder {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            self.bytes.extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            self.flushed_at.push(self.bytes.len());
            Ok(())
        }
    }

    /// Inner writer that panics with "injected panic" once a write would take
    /// it past `panic_after` bytes.
    #[derive(Debug)]
    struct PanickingWriter {
        written: usize,
        panic_after: usize,
    }

    impl Write for PanickingWriter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if self.written + buf.len() > self.panic_after {
                panic!("injected panic");
            }
            self.written += buf.len();
            Ok(buf.len())
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// A `with_threads(threads)` writer whose inner writer panicked during a
    /// `write` call that the caller caught.
    fn writer_after_a_caught_panic(threads: usize) -> BgzfWriter<PanickingWriter> {
        let inner = PanickingWriter {
            written: 0,
            panic_after: 70_000,
        };
        let mut writer = BgzfWriter::with_threads(inner, threads);
        let data = text(2_000_000);

        let payload = data
            .chunks(10_000)
            .find_map(|chunk| {
                panic::catch_unwind(panic::AssertUnwindSafe(|| {
                    let written = writer
                        .write(chunk)
                        .expect("writes should succeed until the inner writer panics");
                    assert_eq!(written, chunk.len());
                }))
                .err()
            })
            .expect("the inner writer should panic");
        assert_eq!(
            payload.downcast_ref::<&str>(),
            Some(&"injected panic"),
            "threads = {threads}"
        );
        writer
    }

    fn assert_poisoned<T: Debug>(result: io::Result<T>) {
        let err = result.expect_err("a failed writer must keep failing");
        assert_eq!(err.kind(), io::ErrorKind::Other);
        assert_eq!(err.to_string(), "BGZF writer failed earlier");
    }

    /// Runs `f` on its own thread and returns its result, panicking with
    /// "timed out" if `f` takes longer than `secs` seconds, so a deadlock fails
    /// the test instead of stalling the run.
    fn within<T, F>(secs: u64, f: F) -> T
    where
        T: Send + 'static,
        F: FnOnce() -> T + Send + 'static,
    {
        let (done, finished) = mpsc::channel();
        let handle = thread::spawn(move || {
            // The receiver is gone only if the test has already timed out.
            let _ = done.send(f());
        });
        match finished.recv_timeout(Duration::from_secs(secs)) {
            Ok(value) => value,
            Err(mpsc::RecvTimeoutError::Timeout) => panic!("timed out after {secs} s"),
            // `f` panicked: re-raise its panic so the test reports the real failure.
            Err(mpsc::RecvTimeoutError::Disconnected) => match handle.join() {
                Err(payload) => panic::resume_unwind(payload),
                Ok(()) => unreachable!("the thread only exits without sending by panicking"),
            },
        }
    }

    /// One call in a sequence of writer calls replayed by [`replay`].
    #[derive(Debug, Clone, Copy)]
    enum Op {
        /// Write the next `len` bytes of [`text`].
        Write(usize),
        Flush,
    }

    /// Deterministic DNA-like bytes: an LCG over `ACGT`, so blocks compress.
    fn text(len: usize) -> Vec<u8> {
        let mut state = 0x9e37_79b9_u32;
        (0..len)
            .map(|_| {
                state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                b"ACGT"[(state >> 30) as usize]
            })
            .collect()
    }

    /// Replays `ops` on `writer` over the bytes of [`text`], then finishes it.
    fn replay(mut writer: BgzfWriter<Vec<u8>>, ops: &[Op]) -> Vec<u8> {
        let total = ops
            .iter()
            .map(|op| match op {
                Op::Write(len) => *len,
                Op::Flush => 0,
            })
            .sum();
        let data = text(total);
        let mut position = 0;
        for op in ops {
            match *op {
                Op::Write(len) => {
                    writer
                        .write_all(&data[position..position + len])
                        .expect("write should succeed");
                    position += len;
                }
                Op::Flush => writer.flush().expect("flush should succeed"),
            }
        }
        writer.finish().expect("finish should succeed")
    }

    /// Output of an inline writer for `ops`.
    fn inline(ops: &[Op]) -> Vec<u8> {
        replay(BgzfWriter::new(Vec::new()), ops)
    }

    /// Output of a `with_threads(threads)` writer for `ops`.
    fn threaded(threads: usize, ops: &[Op]) -> Vec<u8> {
        replay(BgzfWriter::with_threads(Vec::new(), threads), ops)
    }

    /// Worker threads `writer` runs; 0 means it compresses inline.
    fn worker_count<W: Write>(writer: &BgzfWriter<W>) -> usize {
        writer.pool.as_ref().map_or(0, |pool| pool.workers.len())
    }

    /// Blocks sent to workers but not yet written; 0 for an inline writer.
    fn in_flight<W: Write>(writer: &BgzfWriter<W>) -> u64 {
        writer
            .pool
            .as_ref()
            .map_or(0, |pool| pool.sent - pool.written)
    }

    /// Number of blocks in a stream framed by [`compress_block`], counted by
    /// following each block's BSIZE field.
    fn block_count(bytes: &[u8]) -> usize {
        let mut count = 0;
        let mut offset = 0;
        while offset < bytes.len() {
            let bsize = u16::from_le_bytes([bytes[offset + 16], bytes[offset + 17]]);
            offset += usize::from(bsize) + 1;
            count += 1;
        }
        assert_eq!(offset, bytes.len(), "the last block must end the stream");
        count
    }

    /// Panics unless `actual == expected`, without printing megabytes of bytes.
    fn assert_same_bytes(actual: &[u8], expected: &[u8], context: &str) {
        if actual != expected {
            let first = actual
                .iter()
                .zip(expected)
                .position(|(a, e)| a != e)
                .unwrap_or(actual.len().min(expected.len()));
            panic!(
                "{context}: {} bytes where {} were expected, first difference at byte {first}",
                actual.len(),
                expected.len()
            );
        }
    }

    /// How a test builds its writer: `None` is [`BgzfWriter::new`] and
    /// `Some(n)` is [`BgzfWriter::with_threads`] with `n` threads.
    type Mode = Option<usize>;

    /// Inline, then two workers.
    const MODES: [Mode; 2] = [None, Some(2)];

    fn writer_in<W: Write>(mode: Mode, inner: W) -> BgzfWriter<W> {
        match mode {
            None => BgzfWriter::new(inner),
            Some(threads) => BgzfWriter::with_threads(inner, threads),
        }
    }

    /// Panics unless `out` is exactly `blocks`, each framed by
    /// [`compress_block`], followed by one [`EOF_BLOCK`]: no empty block and no
    /// second end-of-file marker.
    fn assert_blocks_then_eof(out: &[u8], blocks: &[&[u8]], context: &str) {
        let mut expected: Vec<u8> = blocks.iter().flat_map(|data| block(data)).collect();
        expected.extend_from_slice(&EOF_BLOCK);
        assert_eq!(
            block_count(out),
            blocks.len() + 1,
            "{context}: data blocks and the EOF block"
        );
        assert_same_bytes(out, &expected, context);
    }

    #[test]
    fn compress_block_frames_a_bgzf_block() {
        let data: Vec<u8> = (0..1_000).map(|i| (i % 7) as u8).collect();
        let mut out = Vec::new();
        compress_block(&data, &mut out).expect("block should compress");

        assert_eq!(
            out[..16],
            [
                0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00
            ]
        );
        assert_eq!(
            usize::from(u16::from_le_bytes([out[16], out[17]])),
            out.len() - 1
        );
        assert_eq!(
            out[out.len() - 8..out.len() - 4],
            crc32(&data).to_le_bytes()
        );
        assert_eq!(out[out.len() - 4..], 1_000u32.to_le_bytes());

        assert_eq!(decode_all(&out), data);

        let mut stream = out.clone();
        stream.extend_from_slice(&EOF_BLOCK);
        let mut read_back = Vec::new();
        BgzfReader::new(&stream[..])
            .read_to_end(&mut read_back)
            .expect("BgzfReader should read the block");
        assert_eq!(read_back, data);
    }

    #[test]
    fn compress_block_appends_to_existing_output() {
        let mut out = b"xyz".to_vec();
        compress_block(b"hello", &mut out).expect("block should compress");

        assert_eq!(&out[..3], b"xyz");
        assert_eq!(&out[3..7], &[0x1f, 0x8b, 0x08, 0x04]);
        assert_eq!(decode_all(&out[3..]), b"hello");
    }

    #[test]
    fn compress_block_rejects_more_than_max_block_data() {
        assert_eq!(MAX_BLOCK_DATA, 65_280);

        let mut out = Vec::new();
        let err = compress_block(&vec![0; MAX_BLOCK_DATA + 1], &mut out)
            .expect_err("oversized block should be rejected");
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert!(out.is_empty(), "a rejected block must not touch the output");

        compress_block(&vec![0; MAX_BLOCK_DATA], &mut out)
            .expect("a full block should be accepted");
    }

    #[test]
    fn eof_block_is_an_empty_bgzf_block() {
        assert_eq!(EOF_BLOCK.len(), 28);

        let mut read_back = Vec::new();
        let read = BgzfReader::new(&EOF_BLOCK[..])
            .read_to_end(&mut read_back)
            .expect("BgzfReader should accept the EOF block");
        assert_eq!(read, 0);

        assert!(decode_all(&EOF_BLOCK).is_empty());
    }

    #[test]
    fn writer_cuts_blocks_at_max_block_data() {
        let m = MAX_BLOCK_DATA;
        let data = noise(2 * m + 10);

        let mut writer = BgzfWriter::new(Vec::new());
        assert_eq!(
            writer.write(&data).expect("write should succeed"),
            data.len(),
            "write should accept the whole buffer"
        );
        let out = writer.finish().expect("finish should succeed");

        let expected = [
            block(&data[..m]),
            block(&data[m..2 * m]),
            block(&data[2 * m..]),
            EOF_BLOCK.to_vec(),
        ]
        .concat();
        assert_eq!(out, expected);
    }

    #[test]
    fn stream_ending_on_a_block_boundary_has_no_empty_block() {
        within(30, || {
            let m = MAX_BLOCK_DATA;
            let data = noise(2 * m);

            for mode in MODES {
                let mut one = writer_in(mode, Vec::new());
                one.write_all(&data[..m]).expect("write should succeed");
                let out = one.finish().expect("finish should succeed");
                assert_blocks_then_eof(&out, &[&data[..m]], &format!("{mode:?}, one block"));

                let mut two = writer_in(mode, Vec::new());
                two.write_all(&data).expect("write should succeed");
                let out = two.finish().expect("finish should succeed");
                assert_blocks_then_eof(
                    &out,
                    &[&data[..m], &data[m..]],
                    &format!("{mode:?}, two blocks"),
                );
            }
        });
    }

    #[test]
    fn flush_on_a_block_boundary_writes_no_empty_block() {
        within(30, || {
            let m = MAX_BLOCK_DATA;
            let data = noise(m + 10);

            for mode in MODES {
                let mut writer = writer_in(mode, Vec::new());
                writer.write_all(&data[..m]).expect("write should succeed");
                writer.flush().expect("flush should succeed");
                writer.write_all(&data[m..]).expect("write should succeed");
                let out = writer.finish().expect("finish should succeed");

                assert_blocks_then_eof(&out, &[&data[..m], &data[m..]], &format!("{mode:?}"));
            }
        });
    }

    #[test]
    fn writer_output_decodes_to_input_for_any_write_sizes() {
        let data = noise(200_000);
        let sizes = [1, 7, 1_000, 65_280, 70_000];

        let mut writer = BgzfWriter::new(Vec::new());
        let mut position = 0;
        let mut turn = 0;
        while position < data.len() {
            let end = (position + sizes[turn % sizes.len()]).min(data.len());
            writer
                .write_all(&data[position..end])
                .expect("write should succeed");
            position = end;
            turn += 1;
        }
        let out = writer.finish().expect("finish should succeed");

        assert_eq!(decode_all(&out), data);
        assert!(out.ends_with(&EOF_BLOCK));
    }

    #[test]
    fn flush_emits_the_partial_block() {
        let data = noise(20);

        let mut writer = BgzfWriter::new(Vec::new());
        writer.write_all(&data[..10]).expect("write should succeed");
        writer.flush().expect("flush should succeed");
        writer.write_all(&data[10..]).expect("write should succeed");
        let out = writer.finish().expect("finish should succeed");

        let expected = [block(&data[..10]), block(&data[10..]), EOF_BLOCK.to_vec()].concat();
        assert_eq!(out, expected);
    }

    #[test]
    fn flush_with_nothing_pending_writes_no_block() {
        let mut writer = BgzfWriter::new(Vec::new());
        writer.flush().expect("flush should succeed");
        writer.flush().expect("flush should succeed");
        let out = writer.finish().expect("finish should succeed");

        assert_eq!(out, EOF_BLOCK);
    }

    #[test]
    fn finish_on_an_empty_writer_writes_only_eof() {
        let out = BgzfWriter::new(Vec::new())
            .finish()
            .expect("finish should succeed");

        assert_eq!(out, EOF_BLOCK);
    }

    #[test]
    fn flush_flushes_the_inner_writer_after_writing_every_block() {
        within(30, || {
            let m = MAX_BLOCK_DATA;
            // Three full blocks and a partial one, so a threaded writer still
            // has blocks in flight when `flush` is called.
            let data = noise(3 * m + 10);
            let expected = [
                block(&data[..m]),
                block(&data[m..2 * m]),
                block(&data[2 * m..3 * m]),
                block(&data[3 * m..]),
            ]
            .concat();

            for mode in MODES {
                let mut writer = writer_in(mode, FlushRecorder::default());
                writer.write_all(&data).expect("write should succeed");
                writer.flush().expect("flush should succeed");

                assert_same_bytes(&writer.inner.bytes, &expected, &format!("{mode:?}"));
                assert_eq!(
                    writer.inner.flushed_at,
                    [expected.len()],
                    "{mode:?}: one inner flush, after every block was written"
                );
            }
        });
    }

    #[test]
    fn finish_flushes_the_inner_writer_after_the_eof_block() {
        within(30, || {
            let data = noise(MAX_BLOCK_DATA + 10);

            for mode in MODES {
                let mut writer = writer_in(mode, FlushRecorder::default());
                writer.write_all(&data).expect("write should succeed");
                let inner = writer.finish().expect("finish should succeed");

                assert!(inner.bytes.ends_with(&EOF_BLOCK), "{mode:?}: EOF block");
                assert_eq!(
                    inner.flushed_at,
                    [inner.bytes.len()],
                    "{mode:?}: one inner flush, after the EOF block"
                );
            }
        });
    }

    #[test]
    fn writer_is_poisoned_after_an_inner_error() {
        let data = noise(3 * MAX_BLOCK_DATA);
        let mut writer = BgzfWriter::new(FailingWriter {
            written: 0,
            fail_after: 100,
        });

        let first = writer
            .write(&data)
            .expect_err("the inner writer should fail");
        assert_eq!(first.to_string(), "injected failure");

        assert_poisoned(writer.write(b"more"));
        assert_poisoned(writer.flush());
        assert_poisoned(writer.finish());
    }

    #[test]
    fn flush_error_poisons_the_writer() {
        let mut writer = BgzfWriter::new(FailingWriter {
            written: 0,
            fail_after: 0,
        });
        writer
            .write_all(b"pending")
            .expect("buffering does not touch the inner writer");

        let first = writer.flush().expect_err("the inner writer should fail");
        assert_eq!(first.to_string(), "injected failure");

        assert_poisoned(writer.write(b"more"));
        assert_poisoned(writer.finish());
    }

    #[test]
    fn threaded_output_matches_inline_for_every_thread_count() {
        within(30, || {
            let sizes = [1, 7, 4_096, 65_280, 100_000];
            let mut ops = Vec::new();
            let mut remaining = 1_000_000;
            let mut turn = 0;
            while remaining > 0 {
                let len = sizes[turn % sizes.len()].min(remaining);
                ops.push(Op::Write(len));
                remaining -= len;
                turn += 1;
                if turn == 3 {
                    ops.push(Op::Flush);
                }
            }
            let expected = inline(&ops);

            for threads in [0, 1, 2, 3, 8] {
                let writer = BgzfWriter::with_threads(Vec::new(), threads);
                // A host that limits threads may start fewer workers, never more.
                let workers = worker_count(&writer);
                if threads < 2 {
                    assert_eq!(workers, 0, "threads = {threads} must compress inline");
                } else {
                    assert!(
                        (1..=threads).contains(&workers),
                        "threads = {threads} started {workers} workers"
                    );
                }
                assert_same_bytes(
                    &replay(writer, &ops),
                    &expected,
                    &format!("threads = {threads}"),
                );
            }
        });
    }

    #[test]
    fn threaded_writer_exceeds_the_in_flight_limit() {
        within(30, || {
            // 41 data blocks and no flush, so every count below is driven far
            // past its limit of 2 blocks per worker in flight: 4 blocks at 2
            // threads, 6 at 3 and 16 at 8.
            let data = text(40 * MAX_BLOCK_DATA + 5);
            let expected = inline(&[Op::Write(data.len())]);
            assert_eq!(
                block_count(&expected),
                42,
                "41 data blocks and the EOF block"
            );

            for threads in [2, 3, 8] {
                let mut writer = BgzfWriter::with_threads(Vec::new(), threads);
                // A host that limits threads may start fewer workers.
                let limit = 2 * worker_count(&writer) as u64;
                for chunk in data.chunks(10_000) {
                    writer.write_all(chunk).expect("write should succeed");
                    assert!(
                        in_flight(&writer) <= limit,
                        "threads = {threads}: no more than 2 blocks per worker may be in flight"
                    );
                }
                let out = writer.finish().expect("finish should succeed");

                assert_same_bytes(&out, &expected, &format!("threads = {threads}"));
            }
        });
    }

    #[test]
    fn one_huge_write_matches_inline() {
        within(30, || {
            let ops = [Op::Write(5_000_000)];
            assert_same_bytes(&threaded(2, &ops), &inline(&ops), "threads = 2");
        });
    }

    #[test]
    fn inner_error_at_the_in_flight_limit_is_returned_without_hanging() {
        within(30, || {
            let inner = FailingWriter {
                written: 0,
                fail_after: 70_000,
            };
            let mut writer = BgzfWriter::with_threads(inner, 2);
            let data = text(2_000_000);

            let first = data
                .chunks(10_000)
                .find_map(|chunk| writer.write(chunk).err())
                .expect("the inner writer should fail");
            assert_eq!(first.to_string(), "injected failure");

            assert_poisoned(writer.write(b"more"));
            assert_poisoned(writer.flush());
            assert_poisoned(writer.finish());
        });
    }

    #[test]
    fn flush_after_every_write_matches_inline() {
        within(30, || {
            let ops: Vec<Op> = (0..300).flat_map(|_| [Op::Write(97), Op::Flush]).collect();
            let expected = inline(&ops);
            assert_eq!(
                block_count(&expected),
                301,
                "one block per flush and the EOF block"
            );

            for threads in [2, 8] {
                assert_same_bytes(
                    &threaded(threads, &ops),
                    &expected,
                    &format!("threads = {threads}"),
                );
            }
        });
    }

    #[test]
    fn absurd_thread_counts_are_clamped() {
        within(30, || {
            assert_eq!(MAX_THREADS, 256);
            let writer = BgzfWriter::with_threads(Vec::new(), 100_000);
            // A host that limits threads may start fewer workers, never more.
            let workers = worker_count(&writer);
            assert!(
                (1..=MAX_THREADS).contains(&workers),
                "threads = 100_000 started {workers} workers"
            );

            let ops = [Op::Write(200_000)];
            assert_same_bytes(&replay(writer, &ops), &inline(&ops), "threads = 100_000");
        });
    }

    #[test]
    fn dropping_a_threaded_writer_without_finish_returns() {
        within(30, || {
            let mut writer = BgzfWriter::with_threads(Vec::new(), 4);
            writer
                .write_all(&text(500_000))
                .expect("write should succeed");
            drop(writer);
        });
    }

    #[test]
    fn threaded_output_reads_back() {
        within(30, || {
            let data = text(300_000);
            let mut writer = BgzfWriter::with_threads(Vec::new(), 3);
            writer
                .write_all(&data[..1_000])
                .expect("write should succeed");
            writer.flush().expect("flush should succeed");
            writer
                .write_all(&data[1_000..])
                .expect("write should succeed");
            let out = writer.finish().expect("finish should succeed");

            let mut read_back = Vec::new();
            BgzfReader::new(&out[..])
                .read_to_end(&mut read_back)
                .expect("BgzfReader should read the stream");
            assert!(read_back == data, "BgzfReader must give back the input");
            assert!(
                decode_all(&out) == data,
                "MultiGzDecoder must give back the input"
            );
        });
    }

    #[test]
    fn a_call_that_panics_poisons_the_writer() {
        within(30, || {
            // 1 compresses inline; with 4 workers the panic comes while blocks
            // are in flight.
            for threads in [1, 4] {
                let mut writer = writer_after_a_caught_panic(threads);
                assert_poisoned(writer.write(b"more"));
                assert_poisoned(writer.flush());
                assert_poisoned(writer.finish());

                // Dropping instead of finishing must also return.
                drop(writer_after_a_caught_panic(threads));
            }
        });
    }

    #[test]
    fn writer_keeps_the_auto_traits_of_the_inner_writer() {
        // BamWriter wraps BgzfWriter. Before worker threads it was Send, Sync,
        // UnwindSafe and RefUnwindSafe whenever its inner writer was.
        fn assert_auto_traits<T: Send + Sync + panic::UnwindSafe + panic::RefUnwindSafe>() {}
        assert_auto_traits::<BgzfWriter<Vec<u8>>>();
        assert_auto_traits::<crate::BamWriter<Vec<u8>>>();
    }
}
