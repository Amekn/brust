use brust_fasta::{Fasta, FastaRecord};
use std::io;

#[test]
fn blank_lines_are_ignored_around_records() {
    let fasta = Fasta::from_reader(&b"\n>seq1 description\nAC\n\nGT\n\n>seq2\nTTAA\n"[..])
        .expect("FASTA should parse with blank spacer lines");

    assert_eq!(
        fasta.records,
        vec![
            FastaRecord::new(
                "seq1".to_string(),
                Some("description".to_string()),
                "ACGT".to_string()
            ),
            FastaRecord::new("seq2".to_string(), None, "TTAA".to_string()),
        ]
    );
}

#[test]
fn empty_id_is_rejected_on_read() {
    let error = Fasta::from_reader(&b">\nACGT\n"[..]).unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    assert_eq!(
        error.to_string(),
        "invalid FASTA at line 1: FASTA record ID must be non-empty"
    );
}

#[test]
fn non_header_junk_before_first_record_is_rejected() {
    let error = Fasta::from_reader(&b"ACGT\n>seq1\nTTAA\n"[..]).unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    assert_eq!(
        error.to_string(),
        "invalid FASTA at line 1: FASTA data before first header line"
    );
}

#[test]
fn cr_only_line_endings_are_rejected() {
    // Read with LF splitting, the whole file would become one header line.
    let error = Fasta::from_reader(&b">seq1\rACGT\r>seq2\rTTAA\r"[..]).unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

#[test]
fn lone_carriage_return_inside_a_sequence_line_is_rejected() {
    let error = Fasta::from_reader(&b">seq1\nAC\rGT\n"[..]).unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

#[test]
fn crlf_line_endings_are_read() {
    let fasta =
        Fasta::from_reader(&b">seq1 description\r\nAC\r\nGT\r\n>seq2\r\nTTAA\r\n"[..]).unwrap();

    assert_eq!(fasta.records[0].sequence, "ACGT");
    assert_eq!(fasta.records[0].description.as_deref(), Some("description"));
    assert_eq!(fasta.records[1].sequence, "TTAA");
}

#[test]
fn carriage_return_ending_the_input_is_read() {
    let fasta = Fasta::from_reader(&b">seq1\r\nAC\r\nGT\r"[..]).unwrap();

    assert_eq!(fasta.records[0].sequence, "ACGT");
}
