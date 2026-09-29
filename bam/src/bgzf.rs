//! BGZF block framing.
//!
//! BGZF is a series of independent gzip members, each holding a limited amount
//! of uncompressed data, so a reader can jump to a block and decompress it
//! without reading the rest of the file. This module has:
//!
//! - [`BgzfReader`], which reads a BGZF stream and reports
//!   [`BgzfVirtualOffset`]s.
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

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Crc;
    use flate2::read::MultiGzDecoder;
    use std::io::{self, Read};

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
}
