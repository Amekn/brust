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

/// Most uncompressed bytes [`compress_block`] accepts for one block (65,280).
///
/// Chosen so a framed block stays within 64 KiB even when the data does not
/// compress.
pub const MAX_BLOCK_DATA: usize = 64 * 1024 - 256;

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
}

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
        }
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
            let Some(block) = read_bgzf_block(&mut self.inner, self.next_block_start)? else {
                self.eof = true;
                self.buffer.clear();
                self.position = 0;
                return Ok(false);
            };

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
}

fn read_bgzf_block<R: Read>(
    reader: &mut R,
    compressed_offset: u64,
) -> io::Result<Option<BgzfBlock>> {
    let mut prefix = [0u8; 12];
    if !read_exact_or_eof(reader, &mut prefix)? {
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

    let deflate_len = remaining_len - 8;
    let deflate = &remaining[..deflate_len];
    let footer = &remaining[deflate_len..];
    let expected_crc = u32::from_le_bytes(footer[..4].try_into().unwrap());
    let expected_isize = u32::from_le_bytes(footer[4..8].try_into().unwrap());

    let mut decoder = DeflateDecoder::new(deflate);
    let mut uncompressed = Vec::with_capacity(expected_isize as usize);
    decoder.read_to_end(&mut uncompressed)?;
    if uncompressed.len() != expected_isize as usize {
        return Err(invalid_data(
            "BGZF ISIZE does not match decompressed length",
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

/// Message for every call on a [`BgzfWriter`] after one of its calls failed.
const FAILED_EARLIER: &str = "BGZF writer failed earlier";

/// Writer that compresses a byte stream into BGZF blocks.
///
/// Bytes are buffered until [`MAX_BLOCK_DATA`] are held, then written to the
/// wrapped writer as one BGZF block. [`Write::flush`] writes any partial block
/// early and flushes the wrapped writer. [`BgzfWriter::finish`] writes the
/// last partial block and [`EOF_BLOCK`], and returns the wrapped writer.
///
/// Block boundaries depend only on the bytes written and on where
/// [`Write::flush`] is called, so the same calls always produce the same
/// bytes. A block that is only partly full is written only by `flush` or
/// `finish`, so avoid calling `flush` after every small write: each call
/// produces a small, poorly compressed block.
///
/// Dropping a `BgzfWriter` without calling `finish` loses any buffered bytes
/// and leaves the stream without its EOF block.
///
/// If any call returns an error, the stream is incomplete, so every later
/// [`Write::write`], [`Write::flush`] or [`BgzfWriter::finish`] returns an
/// error too. The call that hit the original error returns that error.
pub struct BgzfWriter<W: Write> {
    inner: W,
    /// Uncompressed bytes not yet written as a block; never longer than
    /// [`MAX_BLOCK_DATA`].
    pending: Vec<u8>,
    /// Scratch buffer that holds one framed block on its way to `inner`.
    framed: Vec<u8>,
    /// Set once any call has returned an error.
    failed: bool,
}

impl<W: Write> BgzfWriter<W> {
    /// Creates a BGZF writer that writes compressed blocks to `inner`.
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            pending: Vec::with_capacity(MAX_BLOCK_DATA),
            framed: Vec::new(),
            failed: false,
        }
    }

    /// Writes the last partial block and [`EOF_BLOCK`], flushes the wrapped
    /// writer, and returns it.
    ///
    /// Returns an error if an earlier call on this writer failed, or if
    /// writing or flushing fails now.
    pub fn finish(mut self) -> io::Result<W> {
        // `self` is consumed, so there is nothing left to poison on error.
        self.ensure_not_failed()?;
        self.write_pending_partial()?;
        self.inner.write_all(&EOF_BLOCK)?;
        self.inner.flush()?;
        Ok(self.inner)
    }

    fn ensure_not_failed(&self) -> io::Result<()> {
        if self.failed {
            Err(io::Error::other(FAILED_EARLIER))
        } else {
            Ok(())
        }
    }

    /// Passes `result` through, marking the writer failed if it is an error.
    fn track<T>(&mut self, result: io::Result<T>) -> io::Result<T> {
        if result.is_err() {
            self.failed = true;
        }
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
    fn write_pending(&mut self) -> io::Result<()> {
        self.framed.clear();
        compress_block(&self.pending, &mut self.framed)?;
        self.inner.write_all(&self.framed)?;
        self.pending.clear();
        Ok(())
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    /// Buffers all of `buf` and returns `buf.len()`.
    ///
    /// Each time [`MAX_BLOCK_DATA`] bytes have built up, they are written to the
    /// wrapped writer as one block.
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.ensure_not_failed()?;
        let result = self.buffer(buf).map(|()| buf.len());
        self.track(result)
    }

    /// Writes any buffered bytes as a block, then flushes the wrapped writer.
    ///
    /// Nothing is written if no bytes are buffered.
    fn flush(&mut self) -> io::Result<()> {
        self.ensure_not_failed()?;
        let result = self
            .write_pending_partial()
            .and_then(|()| self.inner.flush());
        self.track(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Crc;
    use flate2::read::MultiGzDecoder;
    use std::fmt::Debug;
    use std::io::{self, Read, Write};

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

    fn assert_poisoned<T: Debug>(result: io::Result<T>) {
        let err = result.expect_err("a failed writer must keep failing");
        assert_eq!(err.kind(), io::ErrorKind::Other);
        assert_eq!(err.to_string(), "BGZF writer failed earlier");
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
}
