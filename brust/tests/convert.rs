mod common;

use brust::{Conversion, ConvertOptions, Format, bam, convert, fasta, fastq, sam};
use std::fs;

#[test]
fn convert_dispatcher_reports_formats_and_converts_fastq_to_fasta() {
    let temp = common::TempDir::new("fastq-to-fasta");
    let output = temp.join("reads.fasta");
    let conversion = Conversion::FastqToFasta;

    assert_eq!(conversion.name(), "fastq-to-fasta");
    assert_eq!(conversion.input_format(), Format::Fastq);
    assert_eq!(conversion.output_format(), Format::Fasta);

    convert::convert(
        conversion,
        common::fixture("fastq/UDP0057_sub100.fastq"),
        &output,
    )
    .unwrap();

    let fasta = fasta::Fasta::from_path(&output).unwrap();
    assert_eq!(fasta.records.len(), 100);
    assert_eq!(fasta.records[0].sequence.len(), 687);
}

#[test]
fn fastq_alignment_conversions_stream_unmapped_records() {
    let temp = common::TempDir::new("fastq-alignment-conversions");
    let sam_output = temp.join("reads.sam");
    let bam_output = temp.join("reads.bam");
    let input = common::fixture("fastq/UDP0057_sub100.fastq");

    convert::fastq_to_sam(&input, &sam_output).unwrap();
    let sam = sam::Sam::from_path(&sam_output).unwrap();
    assert_eq!(sam.records.len(), 100);
    assert!(sam.header.records.is_empty());
    assert!(sam.records.iter().all(sam::SamRecord::is_unmapped));
    assert!(sam.records.iter().all(|record| record.rname == "*"));

    convert::fastq_to_bam(&input, &bam_output).unwrap();
    let bam = bam::Bam::from_path(&bam_output).unwrap();
    assert_eq!(bam.header.n_ref, 0);
    assert!(bam.refs.is_empty());
    assert_eq!(bam.records.len(), 100);
    assert!(bam.records.iter().all(bam::BamRecord::is_unmapped));
}

#[test]
fn compressed_fastq_inputs_and_outputs_stream_through_conversions() {
    let temp = common::TempDir::new("compressed-fastq-conversions");
    let fasta_output = temp.join("reads.fasta");
    let fastq_output = temp.join("aligned.fq.gz");

    // Compressed FASTQ input is decoded record-by-record during conversion.
    convert::fastq_to_fasta(
        common::fixture("fastq/UDP0057_sub100.fastq.gz"),
        &fasta_output,
    )
    .unwrap();
    assert_eq!(
        fasta::Fasta::from_path(fasta_output).unwrap().records.len(),
        100
    );

    // A gzip FASTQ destination remains compressed through the atomic temp path.
    convert::sam_to_fastq(common::fixture("sam/aligned.sam"), &fastq_output).unwrap();
    assert_eq!(&fs::read(&fastq_output).unwrap()[..2], &[0x1f, 0x8b]);
    assert_eq!(
        fastq::Fastq::from_path(fastq_output).unwrap().records.len(),
        100
    );
}

#[test]
fn supported_alignment_conversions_round_trip_through_real_fixtures() {
    let temp = common::TempDir::new("alignment-conversions");
    let bam_output = temp.join("aligned.bam");
    let sam_output = temp.join("aligned.sam");
    let fastq_output = temp.join("aligned.fastq");
    let sam_fastq_output = temp.join("aligned-from-sam.fastq");

    convert::sam_to_bam(common::fixture("sam/aligned.sam"), &bam_output).unwrap();
    let bam = bam::Bam::from_path(&bam_output).unwrap();
    assert_eq!(bam.refs.len(), 1);
    assert_eq!(bam.records.len(), 100);

    convert::bam_to_sam(&bam_output, &sam_output).unwrap();
    let sam = sam::Sam::from_path(&sam_output).unwrap();
    assert_eq!(sam.header.sequence_names(), vec!["fc_reference"]);
    assert_eq!(sam.records.len(), 100);

    convert::bam_to_fastq(&bam_output, &fastq_output).unwrap();
    let fastq = fastq::Fastq::from_path(&fastq_output).unwrap();
    assert_eq!(fastq.records.len(), 100);
    assert_eq!(
        fastq.records[0].sequence.len(),
        fastq.records[0].quality.len()
    );

    convert::sam_to_fastq(common::fixture("sam/aligned.sam"), &sam_fastq_output).unwrap();
    let fastq = fastq::Fastq::from_path(&sam_fastq_output).unwrap();
    assert_eq!(fastq.records.len(), 100);
    assert_eq!(
        fastq.records[0].sequence.len(),
        fastq.records[0].quality.len()
    );
}

#[test]
fn failed_conversion_preserves_existing_output_file() {
    let temp = common::TempDir::new("conversion-failure");
    let input = temp.join("missing-sequence.sam");
    let output = temp.join("reads.fastq");
    fs::write(
        &input,
        "@HD\tVN:1.6\tSO:unknown\nread1\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n",
    )
    .unwrap();
    fs::write(&output, "existing output\n").unwrap();

    let error = convert::sam_to_fastq(&input, &output).unwrap_err();

    assert_eq!(error.format(), Some(Format::Sam));
    assert!(
        error
            .to_string()
            .contains("cannot convert record without SEQ")
    );
    assert_eq!(fs::read_to_string(&output).unwrap(), "existing output\n");
}

/// Writes a four-record SAM with one primary forward read, one primary reverse
/// read and one secondary and one supplementary record.
fn write_strand_sam(temp: &common::TempDir) -> std::path::PathBuf {
    let input = temp.join("strands.sam");
    fs::write(
        &input,
        concat!(
            "@HD\tVN:1.6\tSO:unsorted\n",
            "@SQ\tSN:ref\tLN:100\n",
            "fwd\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t!#%'\n",
            "rev\t16\tref\t1\t60\t4M\t*\t0\t0\tAACG\t!#%'\n",
            "sec\t256\tref\t1\t60\t4M\t*\t0\t0\t*\t*\n",
            "sup\t2048\tref\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
        ),
    )
    .unwrap();
    input
}

/// Reads a FASTQ file back as `(id, sequence, quality)` triples.
fn fastq_triples(path: &std::path::Path) -> Vec<(String, String, String)> {
    fastq::Fastq::from_path(path)
        .unwrap()
        .records
        .into_iter()
        .map(|record| (record.id, record.sequence, record.quality))
        .collect()
}

fn triple(id: &str, sequence: &str, quality: &str) -> (String, String, String) {
    (id.to_string(), sequence.to_string(), quality.to_string())
}

#[test]
fn sam_to_fastq_writes_primary_reads_in_original_orientation() {
    let temp = common::TempDir::new("sam-to-fastq-strands");
    let input = write_strand_sam(&temp);
    let output = temp.join("reads.fastq");

    convert::sam_to_fastq(&input, &output).unwrap();

    // The reverse read is restored; the secondary and supplementary are skipped.
    assert_eq!(
        fastq_triples(&output),
        [triple("fwd", "ACGT", "!#%'"), triple("rev", "CGTT", "'%#!")]
    );
}

#[test]
fn bam_to_fastq_writes_primary_reads_in_original_orientation() {
    let temp = common::TempDir::new("bam-to-fastq-strands");
    let input = write_strand_sam(&temp);
    let bam_output = temp.join("strands.bam");
    let output = temp.join("reads.fastq");

    convert::sam_to_bam(&input, &bam_output).unwrap();
    convert::bam_to_fastq(&bam_output, &output).unwrap();

    assert_eq!(
        fastq_triples(&output),
        [triple("fwd", "ACGT", "!#%'"), triple("rev", "CGTT", "'%#!")]
    );
}

#[test]
fn fixture_reverse_strand_reads_come_out_reverse_complemented() {
    let temp = common::TempDir::new("fixture-reverse-strand");
    let output = temp.join("reads.fastq");

    // Independent of brust::seq: the fixture holds only A, C, G, T and N.
    fn reverse_complement(sequence: &str) -> String {
        sequence
            .chars()
            .rev()
            .map(|base| match base {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                'N' => 'N',
                other => panic!("unexpected base {other:?} in fixture"),
            })
            .collect()
    }

    let sam = sam::Sam::from_path(common::fixture("sam/aligned.sam")).unwrap();
    convert::sam_to_fastq(common::fixture("sam/aligned.sam"), &output).unwrap();
    let fastq = fastq::Fastq::from_path(&output).unwrap();

    assert_eq!(fastq.records.len(), 100);
    let mut reverse = 0;
    for (record, read) in sam.records.iter().zip(&fastq.records) {
        assert_eq!(read.id, record.qname);
        if record.flag == 16 {
            reverse += 1;
            assert_eq!(read.sequence, reverse_complement(&record.seq));
            assert_eq!(read.quality, record.qual.chars().rev().collect::<String>());
        } else {
            assert_eq!(record.flag, 0);
            assert_eq!(read.sequence, record.seq);
            assert_eq!(read.quality, record.qual);
        }
    }
    assert_eq!(reverse, 32);
}

#[test]
fn soft_masked_iupac_reverse_read_keeps_case() {
    let temp = common::TempDir::new("reverse-soft-masked");
    let input = temp.join("masked.sam");
    let output = temp.join("reads.fastq");
    fs::write(
        &input,
        concat!(
            "@SQ\tSN:ref\tLN:100\n",
            "masked\t16\tref\t1\t60\t5M\t*\t0\t0\tacgRN\t!!!!!\n",
        ),
    )
    .unwrap();

    convert::sam_to_fastq(&input, &output).unwrap();

    assert_eq!(fastq_triples(&output), [triple("masked", "NYcgt", "!!!!!")]);
}

#[test]
fn bam_to_fastq_rejects_primary_reverse_read_without_qualities() {
    let temp = common::TempDir::new("bam-reverse-without-qual");
    let input = temp.join("no-qual.sam");
    let bam_output = temp.join("no-qual.bam");
    let output = temp.join("reads.fastq");
    fs::write(
        &input,
        concat!(
            "@SQ\tSN:ref\tLN:100\n",
            "noqual\t16\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*\n",
        ),
    )
    .unwrap();

    convert::sam_to_bam(&input, &bam_output).unwrap();
    let error = convert::bam_to_fastq(&bam_output, &output).unwrap_err();

    assert_eq!(error.format(), Some(Format::Bam));
    assert!(
        error
            .to_string()
            .contains("cannot convert record without QUAL to FASTQ"),
        "{error}"
    );
    assert!(!output.exists());
}

#[test]
fn convert_options_default_to_one_thread() {
    assert_eq!(ConvertOptions::default().threads, 1);
    assert_eq!(ConvertOptions::default().threads(8).threads, 8);
}

#[test]
fn threaded_bam_conversions_match_default_bytes() {
    let temp = common::TempDir::new("threaded-bam-conversions");
    let cases = [
        (Conversion::SamToBam, "sam/aligned.sam"),
        (Conversion::SamToBam, "sam/unaligned.sam"),
        (Conversion::FastqToBam, "fastq/UDP0057_sub100.fastq"),
    ];

    for (index, (conversion, fixture)) in cases.into_iter().enumerate() {
        let input = common::fixture(fixture);
        let expected_output = temp.join(&format!("expected-{index}.bam"));
        convert::convert(conversion, &input, &expected_output).unwrap();
        let expected = fs::read(&expected_output).unwrap();

        // 0 and 1 compress inline; 4 uses worker threads. All must match.
        for threads in [0, 4] {
            let output = temp.join(&format!("threads-{threads}-{index}.bam"));
            let options = ConvertOptions::default().threads(threads);
            convert::convert_with(conversion, &input, &output, &options).unwrap();
            assert_eq!(
                fs::read(&output).unwrap(),
                expected,
                "{} on {fixture} with {threads} threads",
                conversion.name()
            );
        }
    }
}

#[test]
fn threads_do_not_change_non_bam_outputs() {
    let temp = common::TempDir::new("threads-non-bam");
    let input = common::fixture("fastq/UDP0057_sub100.fastq");
    let expected_output = temp.join("expected.fasta");
    let threaded_output = temp.join("threaded.fasta");

    convert::convert(Conversion::FastqToFasta, &input, &expected_output).unwrap();
    convert::convert_with(
        Conversion::FastqToFasta,
        &input,
        &threaded_output,
        &ConvertOptions::default().threads(4),
    )
    .unwrap();

    assert_eq!(
        fs::read(&threaded_output).unwrap(),
        fs::read(&expected_output).unwrap()
    );
}

#[test]
fn failed_threaded_conversion_preserves_existing_output() {
    let temp = common::TempDir::new("threaded-conversion-failure");
    let input = temp.join("bad-final-line.sam");
    let output = temp.join("reads.bam");

    // Valid header and records, then a final line that is not a SAM record.
    let mut sam = fs::read_to_string(common::fixture("sam/aligned.sam")).unwrap();
    if !sam.ends_with('\n') {
        sam.push('\n');
    }
    sam.push_str("bad\tline\n");
    fs::write(&input, sam).unwrap();
    fs::write(&output, b"sentinel").unwrap();

    let result = convert::convert_with(
        Conversion::SamToBam,
        &input,
        &output,
        &ConvertOptions::default().threads(4),
    );

    let error = result.unwrap_err();
    assert_eq!(error.format(), Some(Format::Sam));
    assert_eq!(fs::read(&output).unwrap(), b"sentinel");
    // No temporary output is left beside the input and the output.
    let mut names: Vec<_> = fs::read_dir(temp.join("."))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    names.sort();
    assert_eq!(names, ["bad-final-line.sam", "reads.bam"]);
}
