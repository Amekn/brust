//! Sequence primitives for the Brust bioinformatics crates.
//!
//! Every function works on bytes, so callers pass `record.sequence.as_bytes()`.
//! Every byte these functions return is ASCII, so a result converts to a
//! `String` without loss.
//!
//! - [`complement`], [`reverse_complement`]: IUPAC nucleotide complement, keeping case.
//! - [`iupac_matches`]: whether a concrete base fits an IUPAC code.
//! - [`translate_codon`], [`translate`]: standard genetic code (NCBI table 1), frame 0.
//! - [`PhredMean`], [`read_mean_phred`]: Phred mean taken over error probabilities, the
//!   read qscore.
//!
//! Most applications reach this crate as `brust::seq`.
#![forbid(unsafe_code)]

mod codon;
mod nucleotide;
mod quality;

pub use codon::{translate, translate_codon};
pub use nucleotide::{complement, iupac_matches, reverse_complement};
pub use quality::{PhredMean, read_mean_phred};
