//! Lengths and counts read from a BAM file must not size allocations before
//! the data behind them has been read, and BGZF blocks must stay within spec.

use brust_bam::bgzf::{EOF_BLOCK, compress_block};
use brust_bam::{BamAuxArray, BamAuxValue, BamReader, BamRecord, BamRecordAuxiliary, BgzfReader};
use flate2::Compression;
use flate2::write::DeflateEncoder;
use std::io::{self, Read, Write};

/// Frames raw deflate data as one BGZF block with the given CRC and ISIZE.
fn raw_block(deflate: &[u8], crc: u32, isize: u32) -> Vec<u8> {
    let bsize = (12 + 6 + deflate.len() + 8 - 1) as u16;
    let mut block = vec![
        0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, b'B', b'C', 2, 0,
    ];
    block.extend_from_slice(&bsize.to_le_bytes());
    block.extend_from_slice(deflate);
    block.extend_from_slice(&crc.to_le_bytes());
    block.extend_from_slice(&isize.to_le_bytes());
    block
}

fn deflate(data: &[u8]) -> Vec<u8> {
    let mut encoder = DeflateEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data).unwrap();
    encoder.finish().unwrap()
}

fn crc(data: &[u8]) -> u32 {
    let mut crc = flate2::Crc::new();
    crc.update(data);
    crc.sum()
}

fn read_all(stream: &[u8]) -> io::Result<Vec<u8>> {
    let mut data = Vec::new();
    BgzfReader::new(stream).read_to_end(&mut data)?;
    Ok(data)
}

/// A BGZF stream holding `data` in one block, then the EOF block.
fn bgzf(data: &[u8]) -> Vec<u8> {
    let mut stream = Vec::new();
    compress_block(data, &mut stream).unwrap();
    stream.extend_from_slice(&EOF_BLOCK);
    stream
}

fn assert_invalid(result: io::Result<impl std::fmt::Debug>) {
    let error = result.expect_err("malformed input should be rejected");
    assert!(
        matches!(
            error.kind(),
            io::ErrorKind::InvalidData | io::ErrorKind::UnexpectedEof
        ),
        "{error}"
    );
}

#[test]
fn bgzf_block_larger_than_64_kib_is_rejected() {
    // Virtual offsets keep 16 bits for the in-block position.
    let data = vec![b'a'; 70_000];
    let stream = raw_block(&deflate(&data), crc(&data), data.len() as u32);

    assert_invalid(read_all(&stream));
}

#[test]
fn bgzf_isize_claiming_4_gib_is_rejected() {
    let data = b"hello";
    let stream = raw_block(&deflate(data), crc(data), u32::MAX);

    assert_invalid(read_all(&stream));
}

#[test]
fn bytes_after_a_blocks_deflate_stream_are_rejected() {
    let data = b"hello";
    let mut payload = deflate(data);
    payload.extend_from_slice(b"junk");
    let stream = raw_block(&payload, crc(data), data.len() as u32);

    assert_invalid(read_all(&stream));
}

#[test]
fn huge_reference_count_is_an_error_not_an_abort() {
    let mut header = b"BAM\x01".to_vec();
    header.extend_from_slice(&0u32.to_le_bytes());
    header.extend_from_slice(&0x7fff_ffffu32.to_le_bytes());

    assert_invalid(BamReader::from_reader(&bgzf(&header)[..]).map(|_| ()));
}

#[test]
fn huge_header_text_and_name_lengths_are_errors() {
    let mut text = b"BAM\x01".to_vec();
    text.extend_from_slice(&u32::MAX.to_le_bytes());
    assert_invalid(BamReader::from_reader(&bgzf(&text)[..]).map(|_| ()));

    let mut name = b"BAM\x01".to_vec();
    name.extend_from_slice(&0u32.to_le_bytes());
    name.extend_from_slice(&1u32.to_le_bytes());
    name.extend_from_slice(&u32::MAX.to_le_bytes());
    assert_invalid(BamReader::from_reader(&bgzf(&name)[..]).map(|_| ()));
}

#[test]
fn huge_record_block_size_is_an_error() {
    let mut data = b"BAM\x01".to_vec();
    data.extend_from_slice(&0u32.to_le_bytes());
    data.extend_from_slice(&0u32.to_le_bytes());
    data.extend_from_slice(&u32::MAX.to_le_bytes());
    data.extend_from_slice(&[0; 40]);
    let stream = bgzf(&data);
    let mut reader = BamReader::from_reader(&stream[..]).unwrap();

    assert_invalid(reader.read_record());
}

#[test]
fn b_array_count_beyond_its_record_is_an_error() {
    let record = BamRecord {
        fixed: brust_bam::BamRecordFixed {
            block_size: 0,
            ref_id: -1,
            pos: -1,
            l_read_name: 0,
            mapq: 0,
            bin: 4680,
            n_cigar_op: 0,
            flag: 4,
            l_seq: 0,
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
        },
        variable: brust_bam::BamRecordVariable {
            read_name: "r".into(),
            cigar: Vec::new(),
            seq: Vec::new(),
            qual: Vec::new(),
        },
        auxiliary: vec![BamRecordAuxiliary {
            tag: "XB".into(),
            value: BamAuxValue::B(BamAuxArray::f(vec![1.0])),
        }],
    };
    let mut encoded = Vec::new();
    record.encode(&mut encoded).unwrap();
    // The payload ends with the array's u32 count and its one f32 value.
    let count = encoded.len() - 8;
    encoded[count..count + 4].copy_from_slice(&u32::MAX.to_le_bytes());
    let block_size = u32::from_le_bytes(encoded[..4].try_into().unwrap());

    assert_invalid(BamRecord::new(block_size, encoded[4..].to_vec()));
}
