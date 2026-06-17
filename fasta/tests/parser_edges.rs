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
