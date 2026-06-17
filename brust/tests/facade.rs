use brust::{Conversion, Error, Format, Stats, bam, fasta, fastq, pod5, sam};

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
    let _ = Conversion::FastqToFasta;
    let _ = Stats::Fasta(brust::stats::FastaStats {
        file_size_bytes: 0,
        records: 0,
        sequence_lengths: brust::stats::LengthStats {
            count: 0,
            total: 0,
            min: None,
            max: None,
            mean: None,
            n50: None,
            n90: None,
        },
        bases: brust::stats::BaseComposition::default(),
        records_with_description: 0,
        empty_records: 0,
        unique_ids: 0,
        duplicate_id_records: 0,
    });

    let error = Error::invalid(Format::Fasta, "empty ID");
    assert_eq!(error.to_string(), "invalid FASTA: empty ID");
}
