mod common;

use brust::{Format, Stats, stats, validate};
use std::fs;

#[test]
fn validation_dispatcher_accepts_all_fixture_formats() {
    validate::validate(Format::Fasta, common::fixture("fasta/ace2_fragments.fasta")).unwrap();
    validate::validate(Format::Fastq, common::fixture("fastq/UDP0057_sub100.fastq")).unwrap();
    // Gzip input follows the same parser-backed validation path.
    validate::validate(
        Format::Fastq,
        common::fixture("fastq/UDP0057_sub100.fastq.gz"),
    )
    .unwrap();
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

#[test]
fn compressed_fastq_stats_match_plain_fastq_stats() {
    // Biological rollups are computed after streaming decompression.
    let plain = stats::fastq_stats(common::fixture("fastq/UDP0057_sub100.fastq")).unwrap();
    let gzip = stats::fastq_stats(common::fixture("fastq/UDP0057_sub100.fastq.gz")).unwrap();

    assert_eq!(gzip.reads, plain.reads);
    assert_eq!(gzip.read_lengths, plain.read_lengths);
    assert_eq!(gzip.qualities, plain.qualities);
}

fn fastq_qscore_stats(
    temp: &common::TempDir,
    name: &str,
    contents: &str,
) -> brust::stats::FastqStats {
    let input = temp.join(name);
    fs::write(&input, contents).unwrap();
    match stats::stats(Format::Fastq, &input).unwrap() {
        Stats::Fastq(stats) => stats,
        _ => panic!("expected FASTQ stats"),
    }
}

#[test]
fn per_read_qscore_of_uniform_reads_is_their_phred() {
    let temp = common::TempDir::new("qscore-uniform");
    let stats = fastq_qscore_stats(
        &temp,
        "uniform.fastq",
        "@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n@r3\nACGT\n+\nIIII\n",
    );

    let qscore = stats.qualities.per_read_qscore;
    assert_eq!(qscore.count, 3);
    for value in [qscore.min, qscore.max, qscore.mean] {
        assert!((value.unwrap() - 40.0).abs() < 1e-9);
    }
}

#[test]
fn per_read_qscore_averages_error_probabilities() {
    let temp = common::TempDir::new("qscore-mixed");
    let stats = fastq_qscore_stats(&temp, "mixed.fastq", "@r1\nAC\n+\n!I\n");

    let qscore = stats.qualities.per_read_qscore;
    assert_eq!(qscore.count, 1);
    assert!((qscore.mean.unwrap() - 3.0098656839).abs() < 1e-9);
    // The arithmetic mean of Q0 and Q40 is 20, far above the qscore.
    assert_eq!(stats.qualities.per_read_mean_phred.mean, Some(20.0));
}

#[test]
fn per_read_qscore_matches_independent_values_on_the_fixture() {
    let stats = stats::fastq_stats(common::fixture("fastq/UDP0057_sub100.fastq")).unwrap();

    let qscore = &stats.qualities.per_read_qscore;
    let arithmetic = &stats.qualities.per_read_mean_phred;
    assert_eq!(qscore.count, 100);
    assert_eq!(qscore.non_finite_count, 0);
    // Expected values were computed independently in Python, not with brust-seq.
    assert!((qscore.min.unwrap() - 6.37049031904244).abs() < 1e-9);
    assert!((qscore.max.unwrap() - 7.616430289220246).abs() < 1e-9);
    assert!((qscore.mean.unwrap() - 7.229083912490575).abs() < 1e-9);
    // A qscore never exceeds the arithmetic mean Phred.
    assert!(qscore.min.unwrap() <= arithmetic.min.unwrap());
    assert!(qscore.max.unwrap() <= arithmetic.max.unwrap());
    assert!(qscore.mean.unwrap() <= arithmetic.mean.unwrap());
}

#[test]
fn every_fixture_read_has_qscore_at_most_its_arithmetic_mean_phred() {
    let mut reader =
        brust::fastq::FastqReader::from_path(common::fixture("fastq/UDP0057_sub100.fastq"))
            .unwrap();

    let mut reads = 0;
    while let Some(record) = reader.read_record().unwrap() {
        let quality = record.quality.as_bytes();
        let qscore = brust::seq::read_mean_phred(quality).unwrap();
        let arithmetic = quality
            .iter()
            .map(|byte| f64::from(byte.saturating_sub(33)))
            .sum::<f64>()
            / quality.len() as f64;
        assert!(
            qscore <= arithmetic + 1e-9,
            "read {}: qscore {qscore} above arithmetic mean {arithmetic}",
            record.id
        );
        reads += 1;
    }
    assert_eq!(reads, 100);
}

#[test]
fn empty_fastq_has_empty_qscore_summary() {
    let temp = common::TempDir::new("qscore-empty");
    let stats = fastq_qscore_stats(&temp, "empty.fastq", "");

    let qscore = stats.qualities.per_read_qscore;
    assert_eq!(qscore.count, 0);
    assert_eq!(qscore.min, None);
    assert_eq!(qscore.max, None);
    assert_eq!(qscore.mean, None);
}

#[test]
fn validate_and_stats_reject_missing_eof_block() {
    let temp = common::TempDir::new("stats-validate-missing-eof");
    let mut bytes = fs::read(common::fixture("bam/aligned.bam")).unwrap();
    bytes.truncate(bytes.len() - 28);
    let input = temp.join("no_eof.bam");
    fs::write(&input, bytes).unwrap();

    let error = validate::validate_bam(&input).unwrap_err();
    assert!(error.to_string().contains("EOF block"), "{error}");
    let error = stats::bam_stats(&input).unwrap_err();
    assert!(error.to_string().contains("EOF block"), "{error}");
}
