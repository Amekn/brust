//! Opt-in strict BGZF end-of-file check.
//!
//! A strict reader must see the 28-byte EOF marker as the last BGZF block.
//! Lenient readers (the default) must keep accepting streams without it.

use brust_bam::bgzf::{EOF_BLOCK, compress_block};
use brust_bam::{BamReader, BgzfReader};
use std::fs;
use std::io::{self, Read};

const ALIGNED_BAM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/aligned.bam");
const MESSAGE: &str = "BGZF stream ended without the EOF block; the file may be truncated";

fn bytes() -> Vec<u8> {
    fs::read(ALIGNED_BAM).unwrap()
}

/// Record count of the fixture, read leniently.
fn full() -> usize {
    BamReader::from_path(ALIGNED_BAM)
        .unwrap()
        .records()
        .collect::<io::Result<Vec<_>>>()
        .unwrap()
        .len()
}

fn strict_reader(data: &[u8]) -> BamReader<&[u8]> {
    let mut reader = BamReader::from_reader(data).unwrap();
    reader.set_require_eof_block(true);
    reader
}

fn block(data: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    compress_block(data, &mut out).unwrap();
    out
}

fn read_to_end(stream: &[u8], strict: bool) -> io::Result<Vec<u8>> {
    let mut reader = BgzfReader::new(stream);
    reader.set_require_eof_block(strict);
    let mut out = Vec::new();
    reader.read_to_end(&mut out)?;
    Ok(out)
}

#[test]
fn strict_reader_accepts_complete_bam() {
    let full = full();
    let mut reader = BamReader::from_path(ALIGNED_BAM).unwrap();
    reader.set_require_eof_block(true);
    for _ in 0..full {
        reader.read_record().unwrap().unwrap();
    }
    assert!(reader.read_record().unwrap().is_none());
}

#[test]
fn missing_eof_block_is_lenient_by_default_and_rejected_when_strict() {
    let bytes = bytes();
    assert!(bytes.ends_with(&EOF_BLOCK));
    let no_eof = &bytes[..bytes.len() - EOF_BLOCK.len()];
    let full = full();

    let mut lenient = BamReader::from_reader(no_eof).unwrap();
    for _ in 0..full {
        lenient.read_record().unwrap().unwrap();
    }
    assert!(lenient.read_record().unwrap().is_none());

    let mut strict = strict_reader(no_eof);
    for _ in 0..full {
        strict.read_record().unwrap().unwrap();
    }
    let error = strict.read_record().unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    assert_eq!(error.to_string(), format!("invalid BAM: {MESSAGE}"));
}

#[test]
fn strict_error_surfaces_through_records_and_read_all() {
    let bytes = bytes();
    let no_eof = &bytes[..bytes.len() - EOF_BLOCK.len()];

    let collected = strict_reader(no_eof)
        .records()
        .collect::<io::Result<Vec<_>>>();
    assert!(collected.is_err());

    assert!(strict_reader(no_eof).read_all().is_err());
}

#[test]
fn strict_error_is_reported_once() {
    let bytes = bytes();
    let no_eof = &bytes[..bytes.len() - EOF_BLOCK.len()];
    let full = full();

    let mut strict = strict_reader(no_eof);
    let items: Vec<_> = strict.records().take(full + 3).collect();
    assert_eq!(items.len(), full + 1);
    assert!(items[..full].iter().all(|item| item.is_ok()));
    assert!(items[full].is_err());
    assert!(strict.read_record().unwrap().is_none());
}

#[test]
fn cut_inside_a_block_fails_in_both_modes() {
    let bytes = bytes();
    let cut = &bytes[..bytes.len() - 38];

    assert!(BamReader::from_reader(cut).unwrap().read_all().is_err());
    assert!(strict_reader(cut).read_all().is_err());
}

#[test]
fn strict_bgzf_accepts_joined_streams() {
    let stream = [
        block(b"abc"),
        EOF_BLOCK.to_vec(),
        block(b"def"),
        EOF_BLOCK.to_vec(),
    ]
    .concat();
    assert_eq!(read_to_end(&stream, true).unwrap(), b"abcdef");
}

#[test]
fn strict_bgzf_rejects_joined_streams_without_final_marker() {
    let stream = [block(b"abc"), EOF_BLOCK.to_vec(), block(b"def")].concat();
    assert_eq!(read_to_end(&stream, false).unwrap(), b"abcdef");
    let error = read_to_end(&stream, true).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

#[test]
fn strict_bgzf_accepts_non_canonical_middle_empty_block() {
    let mut odd = EOF_BLOCK;
    odd[4] = 1;
    let stream = [
        block(b"abc"),
        odd.to_vec(),
        block(b"def"),
        EOF_BLOCK.to_vec(),
    ]
    .concat();
    assert_eq!(read_to_end(&stream, true).unwrap(), b"abcdef");
}

#[test]
fn strict_bgzf_rejects_empty_input() {
    assert_eq!(read_to_end(b"", false).unwrap(), b"");
    let error = read_to_end(b"", true).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

#[test]
fn strict_bgzf_rejects_stream_without_eof_block() {
    let stream = block(b"abc");
    assert_eq!(read_to_end(&stream, false).unwrap(), b"abc");
    let error = read_to_end(&stream, true).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

#[test]
fn strict_bgzf_rejects_final_empty_block_that_is_not_the_marker() {
    let mut marker = EOF_BLOCK;
    marker[4] = 1; // MTIME byte: still a valid empty block, but not byte-identical
    let stream = [block(b"abc"), marker.to_vec()].concat();
    assert_eq!(read_to_end(&stream, false).unwrap(), b"abc");
    let error = read_to_end(&stream, true).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}
