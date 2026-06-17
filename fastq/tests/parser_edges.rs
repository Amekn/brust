use brust_fastq::Fastq;
use std::io;

#[test]
fn overlong_quality_reports_line_and_lengths() {
    let error = Fastq::from_reader(&b"@seq1\nACGT\n+\nIIIII\n"[..]).unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    assert_eq!(
        error.to_string(),
        "invalid FASTQ at line 4: FASTQ quality length exceeds sequence length (5 > 4)"
    );
}

#[test]
fn trailing_junk_after_complete_record_is_rejected() {
    let error = Fastq::from_reader(&b"@seq1\nACGT\n+\nIIII\njunk\n"[..]).unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    assert_eq!(
        error.to_string(),
        "invalid FASTQ at line 5: FASTQ header line must start with @"
    );
}

#[test]
fn fastq_to_fasta_drops_quality_scores() {
    let fastq = Fastq::from_reader(&b"@seq1 description\nACGT\n+\nIIII\n"[..]).unwrap();
    let fasta = fastq.to_fasta();

    assert_eq!(fasta.records.len(), 1);
    assert_eq!(fasta.records[0].id, "seq1");
    assert_eq!(fasta.records[0].description.as_deref(), Some("description"));
    assert_eq!(fasta.records[0].sequence, "ACGT");
}
