//! Validation helpers for the `brust` facade.
//!
//! Validation is intentionally parser-backed: each function streams through the
//! target format reader and relies on the lower-level crate to perform the
//! format checks it already owns. This avoids duplicating validation logic while
//! still surfacing structured `brust::Error` diagnostics.

use crate::Format;
use crate::Result;
use std::path::Path;

/// Validates a file for a format selected at runtime.
pub fn validate<P: AsRef<Path>>(format: Format, input: P) -> Result<()> {
    match format {
        Format::Fasta => validate_fasta(input),
        Format::Fastq => validate_fastq(input),
        Format::Sam => validate_sam(input),
        Format::Bam => validate_bam(input),
        Format::Pod5 => validate_pod5(input),
    }
}

/// Validates FASTQ by streaming every record through the parser.
pub fn validate_fastq<P: AsRef<Path>>(input: P) -> Result<()> {
    let mut reader = fastq::FastqReader::from_path(input.as_ref())?;
    while let Some(_record) = reader.read_record()? {
        // Reading is validation; parser errors carry line/context where available.
    }
    Ok(())
}

/// Validates FASTA by streaming every record through the parser.
pub fn validate_fasta<P: AsRef<Path>>(input: P) -> Result<()> {
    let mut reader = fasta::FastaReader::from_path(input.as_ref())?;
    while let Some(_record) = reader.read_record()? {
        // Reading is validation; parser errors carry line/context where available.
    }
    Ok(())
}

/// Validates SAM headers and records by streaming through the parser.
pub fn validate_sam<P: AsRef<Path>>(input: P) -> Result<()> {
    let mut reader = sam::SamReader::from_path(input.as_ref())?;
    while let Some(_record) = reader.read_record()? {
        // Reading is validation; parser errors carry line/context where available.
    }
    Ok(())
}

/// Validates BAM headers, references, BGZF blocks, and records by streaming.
pub fn validate_bam<P: AsRef<Path>>(input: P) -> Result<()> {
    let mut reader = bam::BamReader::from_path(input.as_ref())?;
    while let Some(_record) = reader.read_record()? {
        // Reading is validation; parser errors carry context where available.
    }
    Ok(())
}

/// Validates POD5 wrapper metadata and read rows by streaming through the parser.
pub fn validate_pod5<P: AsRef<Path>>(input: P) -> Result<()> {
    let mut reader = pod5::Pod5Reader::from_path(input.as_ref())?;
    while let Some(_record) = reader.read_record()? {
        // Reading is validation; parser errors carry context where available.
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEST_DIR_COUNTER: AtomicU64 = AtomicU64::new(0);

    struct TestDir {
        path: PathBuf,
    }

    impl TestDir {
        fn new() -> Self {
            let counter = TEST_DIR_COUNTER.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir().join(format!(
                "brust-validate-test-{}-{counter}",
                std::process::id()
            ));
            fs::create_dir(&path).unwrap();
            Self { path }
        }

        fn path(&self, name: &str) -> PathBuf {
            self.path.join(name)
        }
    }

    impl Drop for TestDir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.path);
        }
    }

    #[test]
    fn validate_dispatch_accepts_valid_text_formats() {
        // Dispatch should route to the same parser-backed validation as direct calls.
        let dir = TestDir::new();
        let fasta = dir.path("seqs.fasta");
        let fastq = dir.path("reads.fastq");
        let sam = dir.path("reads.sam");
        fs::write(&fasta, ">seq1 description\nACGT\n").unwrap();
        fs::write(&fastq, "@read1\nACGT\n+\nIIII\n").unwrap();
        fs::write(&sam, "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n").unwrap();

        validate(Format::Fasta, &fasta).unwrap();
        validate(Format::Fastq, &fastq).unwrap();
        validate(Format::Sam, &sam).unwrap();
    }

    #[test]
    fn validate_fastq_reports_parser_errors() {
        // The lower-level FASTQ parser should surface malformed records as FASTQ errors.
        let dir = TestDir::new();
        let input = dir.path("broken.fastq");
        fs::write(&input, "@read1\nACGT\n+\nIII\n").unwrap();

        let error = validate_fastq(&input).unwrap_err();

        assert_eq!(error.format(), Some(Format::Fastq));
    }

    #[test]
    fn validate_sam_reports_parser_errors() {
        // Invalid mandatory fields are rejected while streaming alignment records.
        let dir = TestDir::new();
        let input = dir.path("broken.sam");
        fs::write(
            &input,
            "read1\tnot-a-flag\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
        )
        .unwrap();

        let error = validate_sam(&input).unwrap_err();

        assert_eq!(error.format(), Some(Format::Sam));
    }
}
