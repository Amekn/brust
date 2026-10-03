use brust_bam::{BamAuxValue, BamRecord, BamRecordAuxiliary, SamToBamConverter};
use sam::Sam;
use std::io;

const HEADER: &[u8] = b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:1000000000\n";

fn converter() -> SamToBamConverter {
    SamToBamConverter::new(&Sam::from_reader(HEADER).unwrap().header).unwrap()
}

fn sam_record(fields: &str) -> sam::SamRecord {
    let mut text = HEADER.to_vec();
    text.extend_from_slice(fields.as_bytes());
    text.push(b'\n');
    Sam::from_reader(&text[..]).unwrap().records.remove(0)
}

#[test]
fn cigar_lengths_beyond_28_bits_are_errors_not_truncated() {
    let converter = converter();
    let largest = sam_record("r1\t0\tref\t1\t60\t2M268435455N2M\t*\t0\t0\tACGT\tIIII");
    let bam = converter.convert_record(&largest).unwrap();
    assert_eq!(
        bam.to_sam_record(converter.refs()).unwrap().cigar,
        "2M268435455N2M"
    );

    let too_long = sam_record("r1\t0\tref\t1\t60\t2M268435456N2M\t*\t0\t0\tACGT\tIIII");
    let error = converter.convert_record(&too_long).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

#[test]
fn bases_bam_has_no_code_for_are_stored_as_n_and_u_as_t() {
    // SAM spec 4.2.3 maps other characters to N; htslib also stores U as T.
    let converter = converter();
    let record = sam_record("r1\t4\t*\t0\t0\t*\t*\t0\t0\t.UXacgu\tIIIIIII");
    let bam = converter.convert_record(&record).unwrap();

    assert_eq!(bam.to_sam_record(converter.refs()).unwrap().seq, "NTNACGT");
}

#[test]
fn non_hex_bam_hex_values_are_errors_when_converted_to_sam() {
    let converter = converter();
    let mut bam = converter
        .convert_record(&sam_record("r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"))
        .unwrap();
    // A BAM H value with bytes FF 61 reads back through from_utf8_lossy.
    bam.auxiliary.push(BamRecordAuxiliary {
        tag: "XH".into(),
        value: BamAuxValue::H("\u{FFFD}a".into()),
    });

    let error = bam.to_sam_record(converter.refs()).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

fn unmapped_with_name(name: String) -> BamRecord {
    let mut record = converter()
        .convert_record(&sam_record("r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"))
        .unwrap();
    record.fixed.block_size = 0;
    record.fixed.l_read_name = 0;
    record.variable.read_name = name;
    record
}

#[test]
fn read_names_longer_than_254_bytes_are_errors_when_encoding() {
    let mut longest = Vec::new();
    unmapped_with_name("r".repeat(254))
        .encode(&mut longest)
        .unwrap();
    let parsed = BamRecord::new(
        u32::from_le_bytes(longest[..4].try_into().unwrap()),
        longest[4..].to_vec(),
    )
    .unwrap();
    assert_eq!(parsed.variable.read_name, "r".repeat(254));

    for length in [255, 256, 300] {
        let mut output = Vec::new();
        let error = unmapped_with_name("r".repeat(length))
            .encode(&mut output)
            .unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput, "{length}");
        assert!(output.is_empty());
    }
}
