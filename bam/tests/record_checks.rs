//! Record-level checks a BAM reader applies, matching htslib's `bam_read1`,
//! and the `bin` values htslib computes.

use brust_bam::{BamReader, BamRecord, BamWriter, BgzfReader, SamToBamConverter};
use sam::Sam;
use std::io::{self, Read};

const HEADER: &[u8] = b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:2000000000\n";

fn converter() -> SamToBamConverter {
    SamToBamConverter::new(&Sam::from_reader(HEADER).unwrap().header).unwrap()
}

fn record(converter: &SamToBamConverter, line: &str) -> BamRecord {
    let mut text = HEADER.to_vec();
    text.extend_from_slice(line.as_bytes());
    text.push(b'\n');
    converter
        .convert_record(&Sam::from_reader(&text[..]).unwrap().records.remove(0))
        .unwrap()
}

/// Writes `record` to an in-memory BAM and reads it back.
fn read_back(converter: &SamToBamConverter, mut record: BamRecord) -> io::Result<BamRecord> {
    record.fixed.block_size = 0;
    let mut writer = BamWriter::from_writer(Vec::new());
    writer.write_header(converter.header(), converter.refs())?;
    writer.write_record(&record)?;
    let bytes = writer.finish()?;
    let mut reader = BamReader::from_reader(&bytes[..])?;
    Ok(reader.read_record()?.expect("one record"))
}

#[test]
fn bins_match_samtools() {
    // Values samtools 1.22.1 writes for the same SAM records.
    let converter = converter();
    for (line, bin) in [
        ("unplaced\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII", 4680),
        ("insonly\t0\tref\t1\t60\t4I\t*\t0\t0\tACGT\tIIII", 4681),
        ("placedunmapped\t4\tref\t1\t0\t2M\t*\t0\t0\tAC\tII", 4681),
        ("far\t0\tref\t1000000001\t60\t4M\t*\t0\t0\tACGT\tIIII", 180),
        ("normal\t0\tref\t100\t60\t4M\t*\t0\t0\tACGT\tIIII", 4681),
    ] {
        assert_eq!(record(&converter, line).fixed.bin, bin, "{line}");
    }
}

/// Writes `record` validly, then sets its reference IDs in the raw stream,
/// which a BamWriter would refuse to write, and reads it back.
fn read_back_with_ids(
    converter: &SamToBamConverter,
    record: &BamRecord,
    ref_id: i32,
    next_ref_id: i32,
) -> io::Result<BamRecord> {
    use std::io::Write;
    let mut writer = BamWriter::from_writer(Vec::new());
    writer.write_header(converter.header(), converter.refs())?;
    writer.write_record(record)?;
    let compressed = writer.finish()?;
    let mut raw = Vec::new();
    BgzfReader::new(&compressed[..]).read_to_end(&mut raw)?;
    let mut encoded = Vec::new();
    record.encode(&mut encoded)?;
    // After block_size: refID at +4, next_refID at +24.
    let start = raw
        .windows(encoded.len())
        .position(|window| window == encoded)
        .expect("record bytes");
    raw[start + 4..start + 8].copy_from_slice(&ref_id.to_le_bytes());
    raw[start + 24..start + 28].copy_from_slice(&next_ref_id.to_le_bytes());
    let mut bgzf = brust_bam::BgzfWriter::new(Vec::new());
    bgzf.write_all(&raw)?;
    let patched = bgzf.finish()?;
    let mut reader = BamReader::from_reader(&patched[..])?;
    Ok(reader.read_record()?.expect("one record"))
}

#[test]
fn reference_ids_outside_the_dictionary_are_rejected_on_read() {
    // samtools: "Numerical result out of range".
    let converter = converter();
    let mapped = record(&converter, "r1\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\tIIII");
    for (ref_id, next_ref_id) in [(1, -1), (0, 1), (-2, -1), (0, -2)] {
        let error = read_back_with_ids(&converter, &mapped, ref_id, next_ref_id).unwrap_err();
        assert_eq!(
            error.kind(),
            io::ErrorKind::InvalidData,
            "{ref_id} {next_ref_id}"
        );
    }
    assert!(read_back_with_ids(&converter, &mapped, 0, -1).is_ok());
}

#[test]
fn cigar_and_sequence_lengths_must_agree_for_mapped_records() {
    // samtools: "CIGAR and query sequence lengths differ".
    let converter = converter();
    let mut mapped = record(&converter, "r1\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\tIIII");
    mapped.variable.cigar = vec![5 << 4];
    let error = read_back(&converter, mapped).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);

    // Unmapped records and records without SEQ are exempt, as in htslib.
    let mut unmapped = record(&converter, "r1\t4\tref\t1\t0\t4M\t*\t0\t0\tACGT\tIIII");
    unmapped.variable.cigar = vec![5 << 4];
    assert!(read_back(&converter, unmapped).is_ok());
    let no_seq = record(&converter, "r1\t0\tref\t1\t60\t5M\t*\t0\t0\t*\t*");
    assert!(read_back(&converter, no_seq).is_ok());
}

#[test]
fn zero_length_cigar_operations_survive_sam_to_bam_and_back() {
    let converter = converter();
    let bam = record(&converter, "r1\t0\tref\t1\t60\t0M1M0I\t*\t0\t0\tA\tI");

    let back = read_back(&converter, bam).unwrap();
    assert_eq!(
        back.to_sam_record(converter.refs()).unwrap().cigar,
        "0M1M0I"
    );
}

#[test]
fn sequence_string_of_an_inconsistent_record_does_not_panic() {
    // The fields are public, so a caller can leave l_seq above the packed data.
    let converter = converter();
    let mut bam = record(&converter, "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII");
    bam.variable.seq.clear();

    assert!(bam.sequence_string().len() <= 4);
}

#[test]
fn truncated_bgzf_header_is_reported_as_such() {
    let mut data = Vec::new();
    let error = BgzfReader::new(&[0x1f, 0x8b][..])
        .read_to_end(&mut data)
        .unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
    assert!(error.to_string().contains("BGZF"), "{error}");
}

/// A mapped 4-base record stored as `stored` with `CG:B:I` holding `real`.
fn with_cg(converter: &SamToBamConverter, stored: Vec<u32>, real: Vec<u32>) -> BamRecord {
    let mut bam = record(converter, "r1\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\tIIII");
    bam.fixed.n_cigar_op = stored.len() as u16;
    bam.variable.cigar = stored;
    bam.auxiliary.push(brust_bam::BamRecordAuxiliary {
        tag: "CG".into(),
        value: brust_bam::BamAuxValue::B(brust_bam::BamAuxArray::I(real)),
    });
    bam
}

#[test]
fn query_length_is_checked_against_the_cg_cigar() {
    // htslib moves the CG CIGAR in before comparing lengths.
    let converter = converter();
    let placeholder = vec![4 << 4 | 4, 5 << 4 | 3];
    let bad = with_cg(&converter, placeholder, vec![3 << 4, 2 << 4]);
    let error = read_back(&converter, bad).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);

    let good = with_cg(&converter, vec![4 << 4 | 4, 1 << 4], vec![2 << 4, 2 << 4]);
    let back = read_back(&converter, good).unwrap();
    assert_eq!(back.to_sam_record(converter.refs()).unwrap().cigar, "2M2M");
}
