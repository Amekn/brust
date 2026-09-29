//! `BamWriter` output must not depend on the number of compression threads.
//!
//! Each test writes the same header and records at several thread counts and
//! compares the bytes with the single-threaded output, which in turn must equal
//! the plain `BamWriter::from_writer` output.

use brust_bam::bgzf::{EOF_BLOCK, MAX_BLOCK_DATA};
use brust_bam::{
    BamAuxArray, BamAuxValue, BamRecord, BamRecordAuxiliary, BamRecordFixed, BamRecordVariable,
    BamWriter, BgzfReader, SamToBamConverter,
};
use sam::{SamHeader, SamHeaderField, SamHeaderRecord};
use std::fs;
use std::io::Read;

/// Thread counts compared with one thread. 0 and 1 are both inline.
const THREADS: [usize; 4] = [0, 2, 3, 8];

/// `@HD VN:1.6 SO:unsorted` and `@SQ SN:fc LN:687`, then `extra_contigs` more
/// `@SQ` lines named `contig-{i}`.
fn sam_header(extra_contigs: usize) -> SamHeader {
    let sq = |name: String, length: usize| {
        SamHeaderRecord::new(
            "SQ".into(),
            vec![
                SamHeaderField::new("SN".into(), name),
                SamHeaderField::new("LN".into(), length.to_string()),
            ],
        )
    };
    let mut records = vec![
        SamHeaderRecord::new(
            "HD".into(),
            vec![
                SamHeaderField::new("VN".into(), "1.6".into()),
                SamHeaderField::new("SO".into(), "unsorted".into()),
            ],
        ),
        sq("fc".into(), 687),
    ];
    records.extend((0..extra_contigs).map(|i| sq(format!("contig-{i}"), 687)));
    SamHeader { records }
}

fn converter(extra_contigs: usize) -> SamToBamConverter {
    SamToBamConverter::new(&sam_header(extra_contigs)).unwrap()
}

fn aux(tag: &str, value: BamAuxValue) -> BamRecordAuxiliary {
    BamRecordAuxiliary {
        tag: tag.into(),
        value,
    }
}

fn record(name: &str, length: usize, auxiliary: Vec<BamRecordAuxiliary>) -> BamRecord {
    BamRecord {
        fixed: BamRecordFixed {
            block_size: 0,
            ref_id: 0,
            pos: 10,
            l_read_name: 0,
            mapq: 60,
            bin: 4680,
            n_cigar_op: 1,
            flag: 0,
            l_seq: length as u32,
            next_ref_id: -1,
            next_pos: -1,
            tlen: 0,
        },
        variable: BamRecordVariable {
            read_name: name.into(),
            cigar: vec![(length as u32) << 4],
            seq: (0..length.div_ceil(2)).map(|i| (i % 251) as u8).collect(),
            qual: (0..length).map(|i| (i % 41) as u8).collect(),
        },
        auxiliary,
    }
}

/// Length of `record` in BAM binary form, including `block_size`.
fn encoded_len(record: &BamRecord) -> usize {
    let mut bytes = Vec::new();
    record.encode(&mut bytes).unwrap();
    bytes.len()
}

/// A record whose serialised length is exactly `target` bytes (at least 43).
/// With a one-character name and one CIGAR operation, a record of length l
/// takes 42 + ceil(l / 2) + l bytes, which reaches every total except 1
/// (mod 3); a two-character name covers those.
fn record_of_length(target: usize) -> BamRecord {
    let (name, rest) = match (target - 42) % 3 {
        1 => ("nn", target - 43),
        _ => ("n", target - 42),
    };
    let length = rest / 3 * 2 + usize::from(rest % 3 == 2);
    let candidate = record(name, length, vec![]);
    assert_eq!(encoded_len(&candidate), target);
    candidate
}

/// The BAM stream from `BamWriter::from_writer_with_threads`.
fn write(threads: usize, converter: &SamToBamConverter, records: &[BamRecord]) -> Vec<u8> {
    let mut writer = BamWriter::from_writer_with_threads(Vec::new(), threads);
    writer
        .write_header(converter.header(), converter.refs())
        .unwrap();
    for record in records {
        writer.write_record(record).unwrap();
    }
    writer.finish().unwrap()
}

/// The BAM stream from the plain `BamWriter::from_writer`.
fn write_plain(converter: &SamToBamConverter, records: &[BamRecord]) -> Vec<u8> {
    let mut writer = BamWriter::from_writer(Vec::new());
    writer
        .write_header(converter.header(), converter.refs())
        .unwrap();
    for record in records {
        writer.write_record(record).unwrap();
    }
    writer.finish().unwrap()
}

/// Checks that every thread count writes the same bytes as `from_writer`, and
/// returns them.
fn assert_matches_for_every_thread_count(
    converter: &SamToBamConverter,
    records: &[BamRecord],
) -> Vec<u8> {
    let expected = write_plain(converter, records);
    assert_eq!(write(1, converter, records), expected, "threads = 1");
    for threads in THREADS {
        assert_eq!(
            write(threads, converter, records),
            expected,
            "threads = {threads}"
        );
    }
    expected
}

/// The uncompressed BAM stream inside a BGZF stream.
fn decompress(bgzf: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    BgzfReader::new(bgzf).read_to_end(&mut out).unwrap();
    out
}

/// Number of blocks in a BGZF stream, counted by following each block's BSIZE
/// field.
fn block_count(bgzf: &[u8]) -> usize {
    let mut count = 0;
    let mut offset = 0;
    while offset < bgzf.len() {
        let bsize = u16::from_le_bytes([bgzf[offset + 16], bgzf[offset + 17]]);
        offset += usize::from(bsize) + 1;
        count += 1;
    }
    assert_eq!(offset, bgzf.len(), "the last block must end the stream");
    count
}

#[test]
fn header_only_stream_matches_for_every_thread_count() {
    let converter = converter(0);
    assert_matches_for_every_thread_count(&converter, &[]);
}

#[test]
fn header_longer_than_a_block_matches_for_every_thread_count() {
    // About 3,000 reference lines of about 30 bytes: several blocks of header,
    // no records.
    let converter = converter(3_000);
    let bytes = assert_matches_for_every_thread_count(&converter, &[]);
    assert!(
        decompress(&bytes).len() > 2 * MAX_BLOCK_DATA,
        "the header must span several blocks"
    );
}

#[test]
fn stream_ending_exactly_on_a_block_boundary_matches() {
    // Header plus records fill exactly one block, then exactly two. The stream
    // must be those full blocks and one end-of-file marker, with no empty block
    // between or after them. The unit tests in `bgzf.rs` check the same for the
    // writer on its own.
    let converter = converter(0);
    let header_len = decompress(&write(1, &converter, &[])).len();
    let first = record_of_length(MAX_BLOCK_DATA - header_len);
    let second = record_of_length(MAX_BLOCK_DATA);
    for (records, blocks) in [(vec![first.clone()], 1), (vec![first, second], 2)] {
        let bytes = assert_matches_for_every_thread_count(&converter, &records);
        assert_eq!(decompress(&bytes).len(), blocks * MAX_BLOCK_DATA);
        assert!(bytes.ends_with(&EOF_BLOCK), "{blocks} block(s): EOF block");
        assert_eq!(
            block_count(&bytes),
            blocks + 1,
            "{blocks} block(s): data blocks and one EOF block"
        );
    }
}

#[test]
fn incompressible_records_match() {
    // Pseudo-random sequence and quality bytes, so deflate stores the blocks
    // nearly raw, close to the 64 KiB BGZF limit.
    let mut state = 12345_u64;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (state >> 33) as u8
    };
    let mut noisy = record("noisy", 2 * MAX_BLOCK_DATA, vec![]);
    noisy.variable.seq.iter_mut().for_each(|b| *b = next());
    noisy.variable.qual.iter_mut().for_each(|b| *b = next());
    let records = vec![noisy, record("after", 40, vec![])];
    assert_matches_for_every_thread_count(&converter(0), &records);
}

#[test]
fn every_aux_type_and_a_record_larger_than_a_block_match() {
    let every_type = vec![
        aux("XA", BamAuxValue::A(b'x')),
        aux("Xc", BamAuxValue::c(-5)),
        aux("XC", BamAuxValue::C(200)),
        aux("Xs", BamAuxValue::s(-300)),
        aux("XS", BamAuxValue::S(60_000)),
        aux("Xi", BamAuxValue::i(-70_000)),
        aux("XI", BamAuxValue::I(4_000_000_000)),
        aux("Xf", BamAuxValue::f(0.25)),
        aux("XZ", BamAuxValue::Z("text".into())),
        aux("XH", BamAuxValue::H("1AE3".into())),
        aux("Ba", BamAuxValue::B(BamAuxArray::c(vec![-1, 2]))),
        aux("Bb", BamAuxValue::B(BamAuxArray::C(vec![1, 2, 3]))),
        aux("Bc", BamAuxValue::B(BamAuxArray::s(vec![-1, 2]))),
        aux("Bd", BamAuxValue::B(BamAuxArray::S(vec![1, 60_000]))),
        aux("Be", BamAuxValue::B(BamAuxArray::i(vec![-1, 70_000]))),
        aux("Bf", BamAuxValue::B(BamAuxArray::I(vec![1, 4_000_000_000]))),
        aux("Bg", BamAuxValue::B(BamAuxArray::f(vec![0.5, -1.5]))),
    ];
    let records = vec![
        record("short", 20, every_type.clone()),
        record("long", 3 * MAX_BLOCK_DATA, every_type),
        record("after", 30, vec![]),
    ];
    assert_matches_for_every_thread_count(&converter(0), &records);
}

#[test]
fn many_records_across_many_blocks_match() {
    let records = (0..600)
        .map(|i| {
            record(
                &format!("read-{i}"),
                150 + i % 7 * 90,
                vec![
                    aux("NM", BamAuxValue::i(i as i32)),
                    aux("MD", BamAuxValue::Z(format!("{}A{}", i % 50, i % 30))),
                ],
            )
        })
        .collect::<Vec<_>>();
    let stream = records.iter().map(encoded_len).sum::<usize>();
    assert!(
        stream > 3 * MAX_BLOCK_DATA,
        "the records must span several blocks"
    );
    assert_matches_for_every_thread_count(&converter(0), &records);
}

#[test]
fn from_path_with_threads_writes_the_same_file() {
    let converter = converter(0);
    // About 200 KiB of records, so the file holds several blocks.
    let records = (0..300)
        .map(|i| record(&format!("read-{i}"), 400, vec![]))
        .collect::<Vec<_>>();
    let expected = write_plain(&converter, &records);

    let path = std::env::temp_dir().join(format!("brust-bam-threads-{}.bam", std::process::id()));
    let mut writer = BamWriter::from_path_with_threads(&path, 4).unwrap();
    writer
        .write_header(converter.header(), converter.refs())
        .unwrap();
    for record in &records {
        writer.write_record(record).unwrap();
    }
    writer.finish().unwrap();
    let written = fs::read(&path);
    fs::remove_file(&path).unwrap();

    assert_eq!(written.unwrap(), expected);
}
