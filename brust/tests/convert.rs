mod common;

use brust::{Conversion, Format, bam, convert, fasta, fastq, sam};
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
