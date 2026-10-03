//! SAM spec 4.2.2: an alignment with more than 65,535 CIGAR operations is
//! stored in BAM as a `kSmN` placeholder (k = l_seq, m = reference length)
//! with the real CIGAR in a `CG:B:I` tag.

use brust_bam::{BamAuxArray, BamAuxValue, BamRecord, SamToBamConverter};
use sam::Sam;
use std::io;

const OPS: usize = 70_000;

fn header() -> &'static [u8] {
    b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:1000000\n"
}

fn converter() -> SamToBamConverter {
    SamToBamConverter::new(&Sam::from_reader(header()).unwrap().header).unwrap()
}

/// `OPS` alternating 1M/1D operations over an `OPS / 2`-base read.
fn long_cigar() -> String {
    "1M1D".repeat(OPS / 2)
}

fn long_sam_record(flag: u16, rname: &str, pos: u32, extra: &str) -> sam::SamRecord {
    let bases = OPS / 2;
    let line = format!(
        "r1\t{flag}\t{rname}\t{pos}\t60\t{}\t*\t0\t0\t{}\t{}{extra}\n",
        long_cigar(),
        "A".repeat(bases),
        "I".repeat(bases)
    );
    let mut text = header().to_vec();
    text.extend_from_slice(line.as_bytes());
    Sam::from_reader(&text[..]).unwrap().records.remove(0)
}

fn round_trip(record: &BamRecord) -> BamRecord {
    let mut encoded = Vec::new();
    record.encode(&mut encoded).unwrap();
    let block_size = u32::from_le_bytes(encoded[..4].try_into().unwrap());
    BamRecord::new(block_size, encoded[4..].to_vec()).unwrap()
}

#[test]
fn sam_with_more_than_65535_operations_round_trips_through_bam() {
    let converter = converter();
    let bam = converter
        .convert_record(&long_sam_record(0, "ref", 1, ""))
        .unwrap();

    // Placeholder: l_seq soft-clipped, then the reference length skipped.
    let reference_length = (OPS as u32 / 2) * 2;
    assert_eq!(bam.fixed.n_cigar_op, 2);
    assert_eq!(
        bam.variable.cigar,
        [(OPS as u32 / 2) << 4 | 4, reference_length << 4 | 3]
    );
    assert!(matches!(
        &bam.auxiliary.last().unwrap().value,
        BamAuxValue::B(BamAuxArray::I(ops)) if ops.len() == OPS
    ));

    let sam = round_trip(&bam).to_sam_record(converter.refs()).unwrap();
    assert_eq!(sam.cigar, long_cigar());
    assert!(sam.optional.iter().all(|field| field.tag != "CG"));
    assert_eq!(round_trip(&bam).cigar_string(), long_cigar());
}

#[test]
fn sam_with_a_long_cigar_and_its_own_cg_tag_is_refused() {
    let record = long_sam_record(0, "ref", 1, "\tCG:B:I,1");

    let error = converter().convert_record(&record).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

#[test]
fn unplaced_records_keep_their_stored_cigar() {
    // Like htslib, only placed records have their CG CIGAR moved back.
    let converter = converter();
    let bam = converter
        .convert_record(&long_sam_record(4, "*", 0, ""))
        .unwrap();

    let sam = round_trip(&bam).to_sam_record(converter.refs()).unwrap();
    assert_eq!(sam.cigar, format!("{}S{}N", OPS / 2, OPS));
    assert!(sam.optional.iter().any(|field| field.tag == "CG"));
}

#[test]
fn unplaced_placeholders_with_zero_lengths_still_convert_to_sam_text() {
    // No SEQ gives `0S`, no reference-consuming operation gives `0N`; samtools
    // prints both, and SAM's CIGAR grammar allows zero-length operations.
    let converter = converter();
    for cigar in ["1M1D", "1I1P"] {
        let line = format!("r1\t4\t*\t0\t0\t{}\t*\t0\t0\t*\t*\n", cigar.repeat(OPS / 2));
        let mut text = header().to_vec();
        text.extend_from_slice(line.as_bytes());
        let sam = Sam::from_reader(&text[..]).unwrap().records.remove(0);
        let bam = round_trip(&converter.convert_record(&sam).unwrap());

        let back = bam.to_sam_record(converter.refs()).unwrap();
        let mut writer = sam::SamWriter::from_writer(Vec::new());
        writer.write_record(&back).unwrap();
        let reread = Sam::from_reader(&writer.into_inner()[..]).unwrap();
        assert_eq!(reread.records[0].cigar, back.cigar, "{cigar}");
    }
}
