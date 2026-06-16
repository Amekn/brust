//! Public facade for the Brust bioinformatics format crates.
//!
//! This crate keeps the individual format implementations available under one
//! namespace, so consumers can depend on `brust` and use paths such as
//! `brust::fasta::FastaReader` or `brust::pod5::Pod5Reader`.

pub use bam;
pub use fasta;
pub use fastq;
pub use pod5;
pub use sam;

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
