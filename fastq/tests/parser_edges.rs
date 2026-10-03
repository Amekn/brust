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

#[test]
fn zero_length_reads_are_read() {
    // cutadapt writes reads trimmed to nothing as an empty SEQ and QUAL pair.
    let fastq = Fastq::from_reader(&b"@empty\n\n+\n\n@seq1\nACGT\n+\nIIII\n"[..]).unwrap();

    assert_eq!(fastq.records.len(), 2);
    assert_eq!(fastq.records[0].id, "empty");
    assert_eq!(fastq.records[0].sequence, "");
    assert_eq!(fastq.records[0].quality, "");
    assert_eq!(fastq.records[1].sequence, "ACGT");
}

#[test]
fn zero_length_read_needs_its_quality_line() {
    let error = Fastq::from_reader(&b"@empty\n\n+\n"[..]).unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

#[test]
fn zero_length_read_with_quality_is_rejected() {
    let error = Fastq::from_reader(&b"@empty\n\n+\nI\n"[..]).unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

#[test]
fn blank_line_inside_a_sequence_is_still_rejected() {
    let error = Fastq::from_reader(&b"@seq1\nAC\n\nGT\n+\nIIII\n"[..]).unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

#[test]
fn records_iterator_ends_after_an_error() {
    // A gzip member whose deflate data is corrupt: every read fails the same way.
    let mut data = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
    data.extend_from_slice(&[0xff; 32]);
    let mut reader = brust_fastq::FastqReader::from_reader(&data[..]).unwrap();

    let results = reader.records().take(10).collect::<Vec<_>>();

    assert_eq!(results.len(), 1);
    assert!(results[0].is_err());
}

#[test]
fn separator_line_naming_another_record_is_rejected() {
    let error = Fastq::from_reader(&b"@seq1\nACGT\n+seq2\nIIII\n"[..]).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);

    for separator in ["+", "+seq1", "+seq1 description"] {
        let input = format!("@seq1 description\nACGT\n{separator}\nIIII\n");
        assert!(Fastq::from_reader(input.as_bytes()).is_ok(), "{separator}");
    }
}

#[test]
fn quality_outside_printable_ascii_is_rejected() {
    let error = Fastq::from_reader(&b"@seq1\nACGT\n+\nII\x01I\n"[..]).unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

#[test]
fn trailing_whitespace_or_control_characters_in_quality_are_rejected() {
    for quality in ["IIII\x0b", "IIII\t", "IIII\x0c", "IIII ", "III\u{2003}"] {
        let input = format!("@seq1\nACGT\n+\n{quality}\n");
        let error = Fastq::from_reader(input.as_bytes()).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData, "{quality:?}");
    }
    for quality in ["\x0b", " ", "\t"] {
        let input = format!("@empty\n\n+\n{quality}\n");
        let error = Fastq::from_reader(input.as_bytes()).unwrap_err();
        assert_eq!(
            error.kind(),
            io::ErrorKind::InvalidData,
            "empty read {quality:?}"
        );
    }
}

#[test]
fn crlf_line_endings_are_still_read() {
    let fastq =
        Fastq::from_reader(&b"@seq1\r\nACGT\r\n+\r\nIIII\r\n@empty\r\n\r\n+\r\n\r\n"[..]).unwrap();

    assert_eq!(fastq.records[0].quality, "IIII");
    assert_eq!(fastq.records[1].sequence, "");
}
