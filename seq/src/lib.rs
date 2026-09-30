//! Sequence primitives for the Brust bioinformatics crates.
//!
//! Every function works on bytes, so callers pass `record.sequence.as_bytes()`.
//! Every byte these functions return is ASCII, so a result converts to a
//! `String` without loss.
//!
//! - [`complement`], [`reverse_complement`]: IUPAC nucleotide complement, keeping case.
//! - [`iupac_matches`]: whether a concrete base fits an IUPAC code.
//! - [`translate_codon`], [`translate`]: standard genetic code (NCBI table 1), frame 0.
//!
//! Most applications reach this crate as `brust::seq`.
#![forbid(unsafe_code)]

mod codon;
mod nucleotide;

pub use codon::{translate, translate_codon};
pub use nucleotide::{complement, iupac_matches, reverse_complement};
