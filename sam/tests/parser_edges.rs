use brust_sam::{Sam, SamRecord, flags};
use std::io;

#[test]
fn malformed_alignment_reports_line_number() {
    let error = Sam::from_reader(
        &b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:10\nr1\t0\tref\t1\t60\t1M\t*\t0\t0\tA\n"[..],
    )
    .unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    assert_eq!(
        error.to_string(),
        "invalid SAM at line 3: SAM alignment line has fewer than 11 fields"
    );
}

#[test]
fn header_after_alignment_reports_line_number() {
    let error = Sam::from_reader(
        &b"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:10\nr1\t0\tref\t1\t60\t1M\t*\t0\t0\tA\tI\n@CO\tlate\n"[..],
    )
    .unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    assert_eq!(
        error.to_string(),
        "invalid SAM at line 4: header record found after alignment section"
    );
}

#[test]
fn flag_helpers_cover_standard_sam_bits() {
    let record = SamRecord::new(
        "r1".to_string(),
        flags::MULTIPLE_SEGMENTS | flags::FIRST_SEGMENT | flags::DUPLICATE,
        "*".to_string(),
        0,
        0,
        "*".to_string(),
        "*".to_string(),
        0,
        0,
        "*".to_string(),
        "*".to_string(),
        Vec::new(),
    );

    assert!(record.has_multiple_segments());
    assert!(record.is_first_segment());
    assert!(record.is_duplicate());
    assert!(!record.is_last_segment());
    assert!(!record.is_unmapped());
}

#[test]
fn hex_values_must_be_hex_digit_pairs() {
    // Non-ASCII once panicked slicing the string; '+' passed u8::from_str_radix.
    for value in ["A€BC", "+A+b", "0G", "ABC"] {
        let line = format!("r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tXH:H:{value}\n");
        let error = Sam::from_reader(line.as_bytes()).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData, "{value}");
    }
}

#[test]
fn hex_values_accept_either_case() {
    let sam =
        Sam::from_reader(&b"r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tXH:H:0aFf\n"[..]).unwrap();

    assert_eq!(
        sam.records[0].optional[0].value,
        brust_sam::SamOptionalValue::Hex(vec![0x0a, 0xff])
    );
}

const UTF8_HEADER: &str = "@HD\tVN:1.6\n\
    @SQ\tSN:ref\tLN:10\tDS:café\n\
    @RG\tID:g1\tDS:übersicht\n\
    @PG\tID:p1\tCL:align --ref /data/réf.fa\tDS:Ausrichtung\n\
    @CO\tcafé\twith a tab\n";

#[test]
fn utf8_header_values_are_read_where_the_spec_allows_them() {
    // SAMv1 allows UTF-8 in @SQ DS, @RG DS, @PG CL, @PG DS and @CO.
    let sam = Sam::from_reader(UTF8_HEADER.as_bytes()).unwrap();
    let mut writer = brust_sam::SamWriter::from_writer(Vec::new());
    writer.write_header(&sam.header).unwrap();

    assert_eq!(String::from_utf8(writer.into_inner()).unwrap(), UTF8_HEADER);
}

#[test]
fn utf8_and_control_characters_are_rejected_elsewhere_in_headers() {
    for header in [
        "@SQ\tSN:réf\tLN:10\n",
        "@RG\tID:grüppe\n",
        "@SQ\tSN:ref\tLN:10\tDS:bell\u{7}\n",
    ] {
        let error = Sam::from_reader(header.as_bytes()).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData, "{header:?}");
    }
}

#[test]
fn writer_refuses_records_that_would_read_back_differently() {
    let record = |qual: &str| {
        let line = "r1\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n";
        let mut record = Sam::from_reader(line.as_bytes()).unwrap().records.remove(0);
        record.qual = qual.to_string();
        record
    };
    let mut writer = brust_sam::SamWriter::from_writer(Vec::new());

    // A tab inside a field starts a new field: here an extra XX tag.
    let error = writer
        .write_record(&record("I\tXX:Z:injected"))
        .unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    assert!(writer.into_inner().is_empty());
}

#[test]
fn b_arrays_with_empty_elements_are_rejected() {
    let line = |array: &str| format!("r1\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\tXB:B:{array}\n");

    for array in ["i,1,,2", "i,", "i,,1", "f,1.5,"] {
        assert!(Sam::from_reader(line(array).as_bytes()).is_err(), "{array}");
    }
    for array in ["i", "i,1,2", "f,1.5,-2"] {
        assert!(Sam::from_reader(line(array).as_bytes()).is_ok(), "{array}");
    }
}

#[test]
fn non_finite_float_tags_still_round_trip() {
    // BAM files can hold NaN in `f` tags; htslib writes and reads them as text.
    let mut sam = Sam::from_reader(&b"r1\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\tXF:f:1.5\n"[..]).unwrap();
    sam.records[0].optional[0].value = brust_sam::SamOptionalValue::Float(f32::NAN);
    let mut writer = brust_sam::SamWriter::from_writer(Vec::new());
    writer.write_record(&sam.records[0]).unwrap();

    let read_back = Sam::from_reader(&writer.into_inner()[..]).unwrap();
    assert!(matches!(
        read_back.records[0].optional[0].value,
        brust_sam::SamOptionalValue::Float(value) if value.is_nan()
    ));
}

#[test]
fn to_sam_line_refuses_fields_that_would_split() {
    let mut sam = Sam::from_reader(&b"r1\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI\n"[..]).unwrap();
    sam.records[0].qual = "I\tXX:Z:injected".to_string();

    let error = sam.records[0].to_sam_line().unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

#[test]
fn comments_with_nul_are_refused() {
    // samtools turns a NUL in @CO into a line break, so the text would not
    // read back.
    let error = Sam::from_reader(&b"@CO\tbefore\0after\n"[..]).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);

    let mut sam = Sam::from_reader(&b"@CO\tplain\n"[..]).unwrap();
    sam.header.records[0].comment = Some("before\0after".to_string());
    let mut writer = brust_sam::SamWriter::from_writer(Vec::new());
    assert!(writer.write_header(&sam.header).is_err());
}

#[test]
fn zero_length_cigar_operations_are_read() {
    // SAMv1's CIGAR grammar allows them, and samtools accepts them.
    let sam = Sam::from_reader(&b"r1\t0\t*\t0\t0\t0M1M0I\t*\t0\t0\tA\tI\n"[..]).unwrap();

    assert_eq!(sam.records[0].cigar, "0M1M0I");
}

#[test]
fn hd_version_is_as_lenient_as_samtools() {
    // samtools 1.22 reads and writes @HD without VN, or with a VN outside
    // SAMv1's /^[0-9]+\.[0-9]+$/, so these files must keep working here too.
    for header in [
        "@HD\tSO:unsorted\n",
        "@HD\tVN:1\n",
        "@HD\tVN:1.6.0\n",
        "@HD\tVN:v1.6\n",
        "@HD\tVN:1.6\n",
    ] {
        let sam = Sam::from_reader(header.as_bytes()).unwrap();
        let mut writer = brust_sam::SamWriter::from_writer(Vec::new());
        writer.write_header(&sam.header).unwrap();

        assert_eq!(writer.into_inner(), header.as_bytes());
    }
}
