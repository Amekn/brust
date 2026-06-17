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
        // Compile-time construction verifies the facade re-exports the format crates.
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

    #[test]
    fn facade_exports_runtime_helpers() {
        // Public helper types should be reachable from the top-level crate namespace.
        assert_eq!(crate::Conversion::FastqToFasta.name(), "fastq-to-fasta");
        assert_eq!(
            crate::Stats::Fasta(crate::stats::FastaStats {
                file_size_bytes: 0,
                records: 0,
                sequence_lengths: crate::stats::LengthStats {
                    count: 0,
                    total: 0,
                    min: None,
                    max: None,
                    mean: None,
                    n50: None,
                    n90: None,
                },
                bases: crate::stats::BaseComposition::default(),
                records_with_description: 0,
                empty_records: 0,
                unique_ids: 0,
                duplicate_id_records: 0,
            })
            .format(),
            crate::Format::Fasta
        );
    }
}
