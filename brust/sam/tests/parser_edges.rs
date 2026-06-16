use sam::{Sam, SamRecord, flags};
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
