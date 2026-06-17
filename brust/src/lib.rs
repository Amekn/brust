//! Public facade for the Brust bioinformatics format crates.
//!
//! This crate keeps the individual format implementations available under one
//! namespace, so consumers can depend on `brust` and use paths such as
//! `brust::fasta::FastaReader` or `brust::pod5::Pod5Reader`.
/// Format conversion helpers.
pub mod convert;
/// Structured summary statistics helpers.
pub mod stats;
/// Parser-backed validation helpers.
pub mod validate;

/// Re-export of the BAM format crate.
pub use bam;
/// Shared Brust diagnostic, error, format, and result types.
pub use brust_core::{Diagnostic, Error, Format, Result};
/// Supported conversion paths.
pub use convert::Conversion;
/// Re-export of the FASTA format crate.
pub use fasta;
/// Re-export of the FASTQ format crate.
pub use fastq;
/// Re-export of the POD5 format crate.
pub use pod5;
/// Re-export of the SAM format crate.
pub use sam;
/// Runtime-dispatched statistics value.
pub use stats::Stats;

#[cfg(test)]
mod tests {
    #[test]
    fn facade_exports_format_crates() {
        let _ = crate::fasta::FastaRecord {
            id: "seq1".to_string(),
            description: None,
            sequence: "ACGT".to_string(),
        };
        let _ = crate::fastq::FastqRecord {
            id: "seq1".to_string(),
            description: None,
            sequence: "ACGT".to_string(),
            quality: "IIII".to_string(),
        };
        let _ = crate::sam::SamHeader::default();
    }
}
