use brust_bam::{Bam, BamWriter, SamToBamConverter};
use flate2::Crc;
use sam::SamReader;

const ALIGNED_BAM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/aligned.bam");
const UNALIGNED_BAM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/unaligned.bam");
const ALIGNED_SAM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../sam/aligned.sam");

fn fingerprint(bytes: &[u8]) -> (usize, u32) {
    let mut crc = Crc::new();
    crc.update(bytes);
    (bytes.len(), crc.sum())
}

#[test]
fn aligned_bam_rewrite_bytes_are_pinned() {
    let bam = Bam::from_path(ALIGNED_BAM).unwrap();
    let mut out = Vec::new();
    bam.to_writer(&mut out).unwrap();
    assert_eq!(fingerprint(&out), (55475, 1371713897));
}

#[test]
fn unaligned_bam_rewrite_bytes_are_pinned() {
    let bam = Bam::from_path(UNALIGNED_BAM).unwrap();
    let mut out = Vec::new();
    bam.to_writer(&mut out).unwrap();
    assert_eq!(fingerprint(&out), (129915, 112315190));
}

#[test]
fn sam_to_bam_stream_bytes_are_pinned() {
    let mut reader = SamReader::from_path(ALIGNED_SAM).unwrap();
    let converter = SamToBamConverter::new(&reader.header).unwrap();
    let mut writer = BamWriter::from_writer(Vec::new());
    writer
        .write_header(converter.header(), converter.refs())
        .unwrap();
    while let Some(sam_record) = reader.read_record().unwrap() {
        let bam_record = converter.convert_record(&sam_record).unwrap();
        writer.write_record(&bam_record).unwrap();
    }
    let out = writer.finish().unwrap();
    assert_eq!(fingerprint(&out), (55711, 4017805689));
}

#[test]
fn flush_after_header_bytes_are_pinned() {
    let bam = Bam::from_path(ALIGNED_BAM).unwrap();
    let mut out = Vec::new();
    let mut writer = BamWriter::from_writer(&mut out);
    writer.write_header(&bam.header, &bam.refs).unwrap();
    writer.flush().unwrap();
    for record in &bam.records {
        writer.write_record(record).unwrap();
    }
    writer.finish().unwrap();
    assert_eq!(fingerprint(&out), (55425, 2406701901));
}
