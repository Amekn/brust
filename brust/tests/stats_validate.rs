mod common;

use brust::{Format, Stats, stats, validate};
use std::fs;

#[test]
fn validation_dispatcher_accepts_all_fixture_formats() {
    validate::validate(Format::Fasta, common::fixture("fasta/ace2_fragments.fasta")).unwrap();
    validate::validate(Format::Fastq, common::fixture("fastq/UDP0057_sub100.fastq")).unwrap();
    validate::validate(Format::Sam, common::fixture("sam/aligned.sam")).unwrap();
    validate::validate(Format::Bam, common::fixture("bam/aligned.bam")).unwrap();
    validate::validate(Format::Pod5, common::fixture("pod5/A_100.pod5")).unwrap();
}

#[test]
fn validation_preserves_structured_parser_errors() {
    let temp = common::TempDir::new("validate-error");
    let input = temp.join("bad.fastq");
    fs::write(&input, "@read1\nAC\n+\nIII\n").unwrap();

    let error = validate::validate(Format::Fastq, &input).unwrap_err();

    assert_eq!(error.format(), Some(Format::Fastq));
    assert_eq!(
        error.diagnostic().and_then(|diagnostic| diagnostic.line),
        Some(4)
    );
    assert!(
        error
            .to_string()
            .contains("quality length exceeds sequence length")
    );
}

#[test]
fn stats_dispatcher_returns_typed_variants_with_fixture_rollups() {
    let fasta = stats::stats(Format::Fasta, common::fixture("fasta/ace2_fragments.fasta")).unwrap();
    let fastq = stats::stats(Format::Fastq, common::fixture("fastq/UDP0057_sub100.fastq")).unwrap();
    let sam = stats::stats(Format::Sam, common::fixture("sam/aligned.sam")).unwrap();
    let bam = stats::stats(Format::Bam, common::fixture("bam/aligned.bam")).unwrap();
    let pod5 = stats::stats(Format::Pod5, common::fixture("pod5/A_100.pod5")).unwrap();

    match fasta {
        Stats::Fasta(stats) => {
            assert_eq!(stats.records, 6);
            assert_eq!(stats.sequence_lengths.total, 13_500);
        }
        _ => panic!("expected FASTA stats"),
    }

    match fastq {
        Stats::Fastq(stats) => {
            assert_eq!(stats.reads, 100);
            assert_eq!(stats.qualities.q30_bases, 36_286);
        }
        _ => panic!("expected FASTQ stats"),
    }

    match (sam, bam) {
        (Stats::Sam(sam), Stats::Bam(bam)) => {
            assert_eq!(sam.alignments.records, 100);
            assert_eq!(sam.alignments.query_lengths, bam.alignments.query_lengths);
            assert_eq!(sam.alignments.cigar_ops, bam.alignments.cigar_ops);
        }
        _ => panic!("expected SAM and BAM stats"),
    }

    match pod5 {
        Stats::Pod5(stats) => {
            assert_eq!(stats.read_count, 100);
            assert_eq!(stats.total_samples, 1_126_116);
        }
        _ => panic!("expected POD5 stats"),
    }
}
