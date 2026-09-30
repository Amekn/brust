# Sequence helpers, read qscore and strand-aware FASTQ conversion — design

Date: 30/09/2026
Status: approved in conversation (including the 0.2.0 bump), awaiting written-spec review
Branch: `feat/seq-helpers`

## Background

This is the second item from the review of NPTune's Rust code (`~/NPTune/src`) for code worth
backporting into brust. The first item, the parallel BGZF writer, is merged.

NPTune carries small, general sequence primitives that brust lacks:

- IUPAC reverse complement (`src/sequence.rs`). NPTune has two different versions: one that
  upper-cases and turns unknown symbols into `N`, and one in `move_profile.rs` that keeps case
  and passes unknown symbols through unchanged.
- IUPAC base matching (`src/assay.rs`).
- Standard codon translation (`src/sequence.rs`).
- An error-probability-space Phred mean (`src/statistics.rs`, `PhredMean`).

brust's FASTQ stats report `per_read_mean_phred` as the arithmetic mean of Phred values, which
is not how Nanopore tools usually report read quality.

While reading the code, one related gap turned up: `sam-to-fastq` and `bam-to-fastq`
(`brust/src/convert.rs`, `sam_record_to_fastq`) copy SEQ and QUAL exactly as stored. SAM/BAM
stores reverse-strand reads (flag 0x10) reverse-complemented, so brust writes those reads in
the wrong orientation. It also writes secondary and supplementary records, so one read can
appear more than once. This spec fixes that as well.

## Goals

1. A zero-dependency `brust-seq` crate with the sequence and quality primitives above,
   re-exported as `brust::seq`.
2. A per-read qscore (the Phred value of the mean base error probability) in FASTQ stats, with
   the existing fields kept.
3. SAM/BAM record helpers that return SEQ and QUAL in their original sequencing orientation.
4. `sam-to-fastq` and `bam-to-fastq` write reads in their original orientation and skip
   secondary and supplementary records.
5. Every workspace crate moves to version 0.2.0, which the user approved on 30/09/2026. The new
   `brust-seq` crate starts at 0.2.0.

## Non-goals

- Genetic codes other than the standard code (NCBI table 1), reading frames other than 0, and
  six-frame translation.
- Handling RNA and DNA as separate alphabets.
- The rest of NPTune's `statistics.rs` (Jensen–Shannon divergence, `mean`).
- Ignoring a read's first N bases when computing its qscore, as Dorado can. Its exact default
  was not confirmed, so it is not copied.
- A median in the stats, and quality stats for SAM/BAM.
- Paired-read naming (`/1`, `/2`) or interleaving in FASTQ output, `samtools fastq`'s other
  filters and options, and a switch to restore the old conversion behaviour.
- Publishing to crates.io, and changes to NPTune.

## Design

### Crate `brust-seq`

A new workspace member in `seq/`, package `brust-seq`, imported as `brust_seq` by users who
depend on it directly. It has no dependencies and uses the shared workspace metadata.

Inside the workspace, it is depended on the way the facade names the other format crates:
`seq = { package = "brust-seq", path = "../seq", version = "0.2.0" }`.

- The facade depends on it and re-exports it as `pub use seq;`, giving `brust::seq`.
- `brust-sam` and `brust-bam` depend on it for the record helpers.

Layout:

- `seq/src/lib.rs`: crate docs, and re-exports of the functions below at the crate root.
- `seq/src/nucleotide.rs`: `complement`, `reverse_complement`, `iupac_matches`.
- `seq/src/codon.rs`: `translate_codon`, `translate`.
- `seq/src/quality.rs`: `PhredMean`, `read_mean_phred`.

Every function works on bytes (`u8` / `&[u8]`), so callers pass `record.sequence.as_bytes()`.
Every output byte is ASCII, so it converts to `String` without loss.

### Nucleotides

```rust
pub fn complement(base: u8) -> u8;
pub fn reverse_complement(seq: &[u8]) -> Vec<u8>;
pub fn iupac_matches(code: u8, base: u8) -> bool;
```

- `complement` handles every IUPAC code and keeps case:
  - A↔T, C↔G, R↔Y, K↔M, B↔V, D↔H; S, W and N map to themselves.
  - U → A.
  - Lowercase input gives lowercase output (a → t, r → y, and so on).
  - The gap symbols `-` and `.` are returned unchanged.
  - Any other byte becomes `N`.
  - This differs from NPTune, which upper-cases. Callers that want uppercase apply it themselves.
- `reverse_complement` returns the reversed sequence with `complement` applied to each byte.
- `iupac_matches(code, base)` reports whether a concrete base fits an IUPAC code.
  - Both arguments are case-insensitive, and U in either counts as T.
  - `base` must be a concrete A, C, G or T. An ambiguous or unknown `base` never matches.
  - `N` as the `code` matches any concrete base. A byte that is not an IUPAC code matches
    nothing.
  - Otherwise this is NPTune's `iupac_matches`.

### Codons

```rust
pub fn translate_codon(codon: &[u8]) -> Option<u8>;
pub fn translate(seq: &[u8]) -> Vec<u8>;
```

- `translate_codon` uses the standard genetic code (NCBI table 1).
  - Input is case-insensitive, and U counts as T.
  - It returns the one-letter amino acid, `*` for TAA, TAG and TGA, or `None` unless given
    exactly 3 bases that are all A, C, G or T.
- `translate` reads frame 0. Each codon becomes `translate_codon(...)` or `X` when that returns
  `None`. A trailing 1 or 2 bases are ignored, and translation does not stop at `*`.

### Quality

```rust
#[derive(Debug, Clone, Default)]
pub struct PhredMean { /* count, error-probability sum */ }
impl PhredMean {
    pub fn add(&mut self, phred: f64);
    pub fn mean(&self) -> Option<f64>;
}
pub fn read_mean_phred(phred33: &[u8]) -> Option<f64>;
```

- Both give the Phred value of the mean error probability, −10·log₁₀(mean(10^(−Q/10))),
  clamped at 0.0 so the result is never −0.0.
  - An empty input gives `None`.
  - `PhredMean::add` takes finite, non-negative Phred values, not Phred+33 bytes.
- `read_mean_phred` decodes Phred+33 bytes with `saturating_sub(33)`, the same way brust's
  stats code does today.
  - It looks up error probabilities in a 256-entry table built once (`std::sync::LazyLock`), not
    `powf` per base, so FASTQ stats stay fast on large files.
  - Its result must equal `PhredMean` fed the same decoded values, within 1e-12.

### FASTQ stats (`brust/src/stats.rs`)

- `QualityStats` becomes `#[non_exhaustive]`. It is a one-time break for code outside the crate
  that builds it as a struct literal; later fields will not break anything.
- It gains one field after `per_read_mean_phred`:

```rust
/// Summary of per-read quality scores: the Phred value of each read's mean base error
/// probability, −10·log10(mean(10^(−Q/10))). Never above the read's arithmetic mean Phred,
/// and noticeably lower when its base qualities vary.
pub per_read_qscore: FloatStats,
```

- The rustdoc of `per_read_mean_phred` is clarified to say it is the arithmetic mean of each
  read's Phred values. Its behaviour does not change.
- `QualityAccumulator::push` calls `brust_seq::read_mean_phred(quality.as_bytes())` once per
  read and adds the result to a new `FloatAccumulator`.
  - Reads with no quality bases are skipped, as `per_read_mean_phred` does now.
- The text display (`brust stats fastq`) prints a `per_read_qscore` block, in the same format as
  `per_read_mean_phred`, straight after it. Nothing else in the output changes.

### Record helpers

In `brust-sam`:

```rust
impl SamRecord {
    /// SEQ in its original sequencing orientation: reverse-complemented when the record is on
    /// the reverse strand (flag 0x10), unchanged otherwise. `*` stays `*`.
    pub fn original_seq(&self) -> String;
    /// QUAL in its original sequencing orientation: reversed when the record is on the reverse
    /// strand, unchanged otherwise. `*` stays `*`.
    pub fn original_qual(&self) -> String;
}
```

In `brust-bam`, with the same rules:

```rust
impl BamRecord {
    pub fn original_sequence_string(&self) -> String;  // based on sequence_string()
    pub fn original_quality_string(&self) -> String;   // based on quality_string()
}
```

### Strand-aware FASTQ conversion (`brust/src/convert.rs`)

`sam_record_to_fastq` is used by both `sam_to_fastq` and `bam_to_fastq`. It changes to:

- Skip records with flag 0x100 (secondary) or 0x800 (supplementary). A skipped record never
  causes an error, even when its SEQ or QUAL is `*`.
- For every other record, including unmapped, QC-fail (0x200) and duplicate (0x400) records,
  build the FASTQ record from `original_seq()` and `original_qual()`.
  - The read name (QNAME) is unchanged, with no `/1` or `/2` suffix.
  - The description stays `None`.
- Keep today's error messages for a non-skipped record with SEQ or QUAL `*`.

Atomic output, gzip selection, and the functions' signatures do not change. The rustdoc on
`sam_to_fastq`, `bam_to_fastq` and `Conversion::{SamToFastq, BamToFastq}` states the new rules.

This changes output for any input with reverse-strand, secondary or supplementary records. The
README's Convert section says so in one line.

### Version 0.2.0

The `#[non_exhaustive]` change to `QualityStats` and the new conversion output are breaking
changes for a 0.x crate, so the release is a minor bump:

- `[workspace.package] version` becomes `"0.2.0"`.
- `brust`, `brust-core` and `brust-fastq`, which set their own versions today, move to
  `version.workspace = true`, so every crate is 0.2.0.
- Every internal path dependency's `version = "..."` becomes `"0.2.0"`, including the new
  `brust-seq` dependencies.
- Version numbers in the docs follow: the root README's "still early at version `0.1.1`" and
  `brust/README.md`'s `brust-fasta = "0.1.1"` example.
- `Cargo.lock` is updated by the build.

## Testing

Written test-first.

1. **`brust-seq` unit tests.**
   - `complement` and `reverse_complement`:
     - every IUPAC pair in both cases;
     - U → A;
     - gaps unchanged;
     - an unknown byte becomes N;
     - reverse complement applied twice gives the original sequence when the input has no
       unknown bytes or U.
   - `iupac_matches`: NPTune's NC-IUB table test (`~/NPTune/src/assay.rs`), all 15 codes against
     `ACGTN` in both cases, plus U counted as T.
   - `translate_codon`:
     - NPTune's check of all 64 codons against NCBI table 1, and its check of the 20 amino-acid
       identities (`~/NPTune/src/sequence.rs`);
     - lowercase and U input;
     - `None` for ambiguous or wrong-length input.
   - `translate`: frame 0, `X` for ambiguous codons, trailing bases ignored, continues past `*`.
   - `PhredMean` and `read_mean_phred`:
     - NPTune's worked values: empty → None; 10 → 10; 5 × 30 → 30; 93, 93 → 93; 0, 0 → 0;
       0 and 40 → 3.0098656839; 10 and 30 → 12.9670862188;
     - `read_mean_phred` equals `PhredMean` on the same decoded values, for every Phred 0–93;
     - bytes below 33 decode as Phred 0.
2. **FASTQ stats.**
   - A FASTQ whose reads are all Q40 gives `per_read_qscore` min = max = mean = 40.
   - A read of one `!` (Q0) and one `I` (Q40) gives 3.0098656839.
   - For every read in `fastq/UDP0057_sub100.fastq`, the qscore is ≤ the arithmetic mean Phred.
   - The fixture's `per_read_qscore` min, max and mean equal values computed independently
     (Python, recorded in the plan) within 1e-9.
   - The CLI output contains a `per_read_qscore` block.
   - The existing stats tests still pass.
3. **Record helpers.**
   - `SamRecord` and `BamRecord` original-orientation methods, for flag 0 and 16, `*` SEQ and
     QUAL, and lowercase and IUPAC bases.
   - The BAM method agrees with the SAM method for the same record after `SamToBamConverter`.
4. **Conversion.**
   - A synthetic SAM with one forward, one reverse, one secondary (SEQ `*`) and one
     supplementary record. Run through `sam_to_fastq`, and through `sam_to_bam` then
     `bam_to_fastq`. Both must give exactly the forward record unchanged, then the reverse
     record reverse-complemented with QUAL reversed.
   - The existing fixture tests still give 100 records. Also, each of the 32 reverse-strand
     records in `sam/aligned.sam` comes out as the reverse complement of its SEQ, and the 68
     forward records come out unchanged.
   - A primary record with SEQ `*` still gives today's error.

### Verification before completion

`cargo fmt --check`, `cargo clippy --workspace --all-targets` (no warnings),
`cargo test --workspace`, `cargo doc --workspace --no-deps` (no warnings), and
`cargo package --list -p brust-seq` (the new crate packages), and
`cargo metadata --no-deps --format-version 1` showing every workspace package at 0.2.0.

## Documentation

- `seq/README.md`: purpose, API list and a short example, in the style of the other crates'
  READMEs.
- Root `README.md`:
  - `brust-seq` in the package table and the `use brust::{...}` line;
  - the `per_read_qscore` value in the Statistics list;
  - the new FASTQ conversion rules in the Convert section.
- `sam/README.md` and `bam/README.md`: the original-orientation methods.
- Rustdoc for every new public item.

## Follow-up (not in this work)

- NPTune can replace its `sequence.rs` helpers, `iupac_matches` and `PhredMean` with
  `brust::seq` once a brust release includes them. It needs to upper-case where it relies on
  that today. The paper's run stays pinned.
- Publishing 0.2.0 to crates.io, crate by crate in dependency order, is the user's to run.
