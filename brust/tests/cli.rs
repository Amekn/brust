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
fn stats_cli_prints_per_read_qscore() {
    let output = brust()
        .args([
            "stats",
            "fastq",
            common::fixture("fastq/UDP0057_sub100.fastq")
                .to_str()
                .unwrap(),
        ])
        .output()
        .unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("per_read_qscore:"));

    // An empty FASTQ also succeeds, with an empty qscore summary.
    let temp = common::TempDir::new("cli-qscore-empty");
    let empty = temp.join("empty.fastq");
    fs::write(&empty, "").unwrap();
    let output = brust()
        .args(["stats", "fastq", empty.to_str().unwrap()])
        .output()
        .unwrap();

    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    let block = stdout
        .split("  per_read_qscore:\n")
        .nth(1)
        .expect("per_read_qscore block");
    let block: Vec<&str> = block.lines().take(5).collect();
    assert_eq!(
        block,
        [
            "    count: 0",
            "    non_finite_count: 0",
            "    min: -",
            "    max: -",
            "    mean: -",
        ]
    );
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

#[test]
fn convert_cli_threads_match_default_bytes() {
    let temp = common::TempDir::new("cli-threads-bytes");
    let sam_input = common::fixture("sam/aligned.sam");
    let fastq_input = common::fixture("fastq/UDP0057_sub100.fastq");

    let sam_default = temp.join("sam-default.bam");
    let sam_threaded = temp.join("sam-threaded.bam");
    let fastq_default = temp.join("fastq-default.bam");
    let fastq_threaded = temp.join("fastq-threaded.bam");

    let run = |args: &[&str]| {
        let output = brust().args(args).output().unwrap();
        assert!(output.status.success());
        assert!(
            String::from_utf8(output.stdout)
                .unwrap()
                .contains("Conversion completed")
        );
        assert!(String::from_utf8(output.stderr).unwrap().is_empty());
    };

    run(&[
        "convert",
        "sam-to-bam",
        sam_input.to_str().unwrap(),
        sam_default.to_str().unwrap(),
    ]);
    run(&[
        "convert",
        "sam-to-bam",
        sam_input.to_str().unwrap(),
        sam_threaded.to_str().unwrap(),
        "--threads",
        "4",
    ]);
    assert_eq!(
        fs::read(&sam_default).unwrap(),
        fs::read(&sam_threaded).unwrap()
    );

    run(&[
        "convert",
        "fastq-to-bam",
        fastq_input.to_str().unwrap(),
        fastq_default.to_str().unwrap(),
    ]);
    run(&[
        "convert",
        "fastq-to-bam",
        fastq_input.to_str().unwrap(),
        fastq_threaded.to_str().unwrap(),
        "-t",
        "4",
    ]);
    assert_eq!(
        fs::read(&fastq_default).unwrap(),
        fs::read(&fastq_threaded).unwrap()
    );
}

#[test]
fn convert_cli_rejects_zero_threads() {
    let temp = common::TempDir::new("cli-threads-zero");
    let input = common::fixture("sam/aligned.sam");
    let output_path = temp.join("aligned.bam");

    let output = brust()
        .args([
            "convert",
            "sam-to-bam",
            input.to_str().unwrap(),
            output_path.to_str().unwrap(),
            "--threads",
            "0",
        ])
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("invalid value '0' for '--threads"));
    assert!(!output_path.exists());
}

#[test]
fn convert_cli_threads_only_on_bam_outputs() {
    let temp = common::TempDir::new("cli-threads-fasta");
    let input = common::fixture("fastq/UDP0057_sub100.fastq");
    let output_path = temp.join("reads.fasta");

    let output = brust()
        .args([
            "convert",
            "fastq-to-fasta",
            input.to_str().unwrap(),
            output_path.to_str().unwrap(),
            "--threads",
            "2",
        ])
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("unexpected argument '--threads'"));
    assert!(!output_path.exists());
}

#[test]
fn convert_cli_help_describes_threads() {
    let output = brust()
        .args(["convert", "sam-to-bam", "--help"])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert!(
        String::from_utf8(output.stdout)
            .unwrap()
            .contains("BGZF compression threads")
    );
}

fn names(dir: &std::path::Path) -> Vec<String> {
    let mut names: Vec<_> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect();
    names.sort();
    names
}

#[test]
fn cli_rejects_bam_without_eof_block() {
    let temp = common::TempDir::new("cli-missing-eof");
    let mut bytes = fs::read(common::fixture("bam/aligned.bam")).unwrap();
    bytes.truncate(bytes.len() - 28);
    let input = temp.join("no_eof.bam");
    fs::write(&input, bytes).unwrap();
    let input = input.to_str().unwrap();
    let out_sam = temp.join("out.sam");
    let out_fastq = temp.join("out.fastq");

    let commands: [Vec<&str>; 4] = [
        vec!["validate", "bam", input],
        vec!["stats", "bam", input],
        vec!["convert", "bam-to-sam", input, out_sam.to_str().unwrap()],
        vec![
            "convert",
            "bam-to-fastq",
            input,
            out_fastq.to_str().unwrap(),
        ],
    ];
    for args in commands {
        let output = brust().args(&args).output().unwrap();
        assert!(!output.status.success(), "{args:?}");
        let stderr = String::from_utf8(output.stderr).unwrap();
        assert!(stderr.contains("EOF block"), "{args:?}: {stderr}");
    }
    assert_eq!(names(&temp.join(".")), ["no_eof.bam"]);
}

#[test]
fn cli_convert_accepts_bare_relative_output() {
    let temp = common::TempDir::new("cli-relative-output");
    let input = std::path::absolute(common::fixture("fastq/UDP0057_sub100.fastq")).unwrap();

    let output = brust()
        .current_dir(temp.join("."))
        .args([
            "convert",
            "fastq-to-fasta",
            input.to_str().unwrap(),
            "out.fasta",
        ])
        .output()
        .unwrap();

    assert!(output.status.success(), "{output:?}");
    assert_eq!(names(&temp.join(".")), ["out.fasta"]);
}
