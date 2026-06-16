use brust::{Error, Format, bam, fasta, fastq, pod5, sam};

#[test]
fn facade_reexports_format_crates_and_shared_error() {
    let _ = fasta::FastaRecord::new("seq1".to_string(), None, "ACGT".to_string());
    let _ = fastq::FastqRecord::new(
        "seq1".to_string(),
        None,
        "ACGT".to_string(),
        "IIII".to_string(),
    );
    let _ = sam::SamHeader::default();
    let _ = bam::BgzfVirtualOffset::new(10, 4);
    let _ = pod5::Pod5SectionKind::Reads;

    let error = Error::invalid(Format::Fasta, "empty ID");
    assert_eq!(error.to_string(), "invalid FASTA: empty ID");
}
