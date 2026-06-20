mod common;

use brust::fasta;
use std::fs;
use std::process::Command;

fn brust() -> Command {
    Command::new(env!("CARGO_BIN_EXE_brust"))
}

#[test]
fn stats_cli_prints_human_readable_summary() {
    let output = brust()
        .args([
            "stats",
            "fasta",
            common::fixture("fasta/ace2_fragments.fasta")
                .to_str()
                .unwrap(),
        ])
        .output()
        .unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("FASTA statistics"));
    assert!(stdout.contains("records: 6"));
    assert!(String::from_utf8(output.stderr).unwrap().is_empty());
}

#[test]
fn validate_cli_reports_structured_errors_and_nonzero_exit() {
    let temp = common::TempDir::new("cli-validate");
    let input = temp.join("bad.fastq");
    fs::write(&input, "@read1\nAC\n+\nIII\n").unwrap();

    let output = brust()
        .args(["validate", "fastq", input.to_str().unwrap()])
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("Validation failed"));
    assert!(stderr.contains("invalid FASTQ at line 4"));
    assert!(stderr.contains("quality length exceeds sequence length"));
}

#[test]
fn convert_cli_writes_expected_output() {
    let temp = common::TempDir::new("cli-convert");
    let output_path = temp.join("reads.fasta");
    let input_path = common::fixture("fastq/UDP0057_sub100.fastq");

    let output = brust()
        .args([
            "convert",
            "fastq-to-fasta",
            input_path.to_str().unwrap(),
            output_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .contains("Conversion completed")
    );
    assert!(String::from_utf8(output.stderr).unwrap().is_empty());

    let fasta = fasta::Fasta::from_path(&output_path).unwrap();
    assert_eq!(fasta.records.len(), 100);
}

#[test]
fn cli_validates_and_converts_compressed_fastq() {
    let input = common::fixture("fastq/UDP0057_sub100.fastq.gz");
    let temp = common::TempDir::new("cli-compressed-fastq");
    let output_path = temp.join("reads.fasta");

    // CLI validation and conversion share the gzip-aware facade readers.
    let validation = brust()
        .args(["validate", "fastq", input.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(validation.status.success());

    let conversion = brust()
        .args([
            "convert",
            "fastq-to-fasta",
            input.to_str().unwrap(),
            output_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(conversion.status.success());
    assert_eq!(
        fasta::Fasta::from_path(output_path).unwrap().records.len(),
        100
    );
}
