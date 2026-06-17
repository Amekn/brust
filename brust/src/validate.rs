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
