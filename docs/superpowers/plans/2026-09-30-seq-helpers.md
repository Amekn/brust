# Sequence Helpers, Read Qscore and Strand-Aware FASTQ Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a zero-dependency `brust-seq` crate (IUPAC reverse complement and matching, standard codon translation, error-space Phred means), a `per_read_qscore` FASTQ stat, original-orientation SAM/BAM record helpers, and strand-aware `sam-to-fastq`/`bam-to-fastq`, all released as 0.2.0.

**Architecture:** The workspace moves to 0.2.0 first, so new path dependencies can say `version = "0.2.0"`. `seq/` is a new crate with three focused modules. The stats code, the `sam`/`bam` record types, and the facade's conversion code consume it. The conversion change is local to `sam_record_to_fastq` in `brust/src/convert.rs`.

**Tech Stack:** Rust 2024 (1.96), std only (`std::sync::LazyLock`), and the existing workspace crates.

**Spec:** `docs/superpowers/specs/2026-09-30-seq-helpers-design.md`

## Global Constraints

- `brust-seq` has no dependencies. No other crate gains new third-party dependencies.
- Every workspace crate is version 0.2.0, and every internal path dependency says `version = "0.2.0"`.
- `brust-seq` lives in `seq/`. Other crates depend on it as `seq = { package = "brust-seq", path = "../seq", version = "0.2.0" }`, and the facade re-exports it as `brust::seq`.
- Existing public signatures do not change. The only intended breaking change is `#[non_exhaustive]` on `QualityStats`.
- Every sequence output byte is ASCII.
- `complement` keeps case, maps U → A, returns gap symbols `-` and `.` unchanged, and maps every other non-IUPAC byte to `N`.
- Standard genetic code only (NCBI table 1). `translate` reads frame 0, gives `X` for untranslatable codons, ignores a trailing 1 or 2 bases, and continues past `*`.
- Qscore is −10·log₁₀(mean(10^(−Q/10))), clamped at 0.0, and `None` for an empty input. Phred+33 decoding uses `saturating_sub(33)`.
- Conversion skips flags 0x100 and 0x800 only. The read name is unchanged, and today's error messages stay for a non-skipped `*` SEQ or QUAL.
- Rustdoc and README text match the style of the file being edited. Every new public item has rustdoc.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. An empty FASTQ (0 reads): `per_read_qscore` has count 0, and min, max and mean are `None`. The CLI prints `-` without panicking. Test in Task 3.
2. A FASTQ with a zero-length record: it is left out of `per_read_qscore` (count = reads with bases), and nothing is NaN. Test in Task 3.
3. A 1,000,000-base Q20 read: its qscore is 20 within 1e-9, so summing a long read loses no precision. Test in Task 3.
4. A reverse-strand BAM record with missing qualities (all 0xff): `original_quality_string()` is `"*"`, and `bam-to-fastq` rejects it with the same "without QUAL" error as SAM. Tests in Tasks 4 and 5.
5. A soft-masked or IUPAC reverse-strand SEQ through `sam-to-fastq`: the FASTQ keeps case and complements IUPAC codes (`acgRN` with flag 16 → `NYcgt`). Test in Task 5.

---

### Task 1: Version 0.2.0

**Files:**
- Modify: `Cargo.toml`: `[workspace.package] version = "0.2.0"`.
- Modify: `brust-core/Cargo.toml`, `brust/Cargo.toml` and `fastq/Cargo.toml`: `version = "0.1.2"` becomes `version.workspace = true`.
- Modify: every internal path dependency in `brust/`, `fasta/`, `fastq/`, `sam/`, `bam/` and `pod5/` `Cargo.toml`: `version = "0.2.0"`.
- Modify: `README.md` (line 357, "still early at version `0.1.1`" → `0.2.0`) and `brust/README.md` (line 168, `brust-fasta = "0.1.1"` → `"0.2.0"`).
- Modify: `Cargo.lock` (updated by the build).

**Interfaces:**
- Produces: all packages at 0.2.0. Later tasks write `version = "0.2.0"` on new path dependencies.

- [ ] **Step 1: Record the current state.** Run `cargo metadata --no-deps --format-version 1 | python3 -c "import json,sys; print(sorted((p['name'],p['version']) for p in json.load(sys.stdin)['packages']))"`. Expected: every package at 0.1.2.
- [ ] **Step 2: Make the edits listed above.** Leave third-party versions alone.
- [ ] **Step 3: Verify.** Run the Step 1 command again: every package reads `0.2.0`. `cargo run -q -p brust -- --version` prints `brust 0.2.0`. Run `cargo test --workspace` and `cargo clippy --workspace --all-targets`: everything passes (170 tests) with no warnings. `grep -rn '"0\.1\.' --include=Cargo.toml .` finds no internal crate versions.
- [ ] **Step 4: Commit.**

```bash
git add Cargo.toml Cargo.lock */Cargo.toml README.md brust/README.md
git commit -m "chore: version 0.2.0 across the workspace"
```

---

### Task 2: `brust-seq` crate with nucleotides and codons

**Files:**
- Create: `seq/Cargo.toml`. Package `brust-seq`, with every metadata field set to `.workspace = true` (as in `sam/Cargo.toml`), `readme = "README.md"`, the description `"Sequence and quality primitives for Brust: IUPAC reverse complement and matching, codon translation, and Phred means."`, and no `[dependencies]`.
- Create: `seq/src/lib.rs`: crate docs, `mod nucleotide; mod codon;` (Task 3 adds `quality`), and re-exports of their public functions at the root.
- Create: `seq/src/nucleotide.rs` and `seq/src/codon.rs`, each with a unit-test module.
- Create: `seq/README.md`, in the style of `sam/README.md`.
- Modify: root `Cargo.toml`: add `"seq"` to `members`.
- Modify: `brust/Cargo.toml` (the `seq` dependency) and `brust/src/lib.rs` (`/// Re-export of the sequence helper crate.` then `pub use seq;`).
- Modify: `brust/tests/facade.rs`, and `README.md` (a `brust-seq` row in the package table, and `seq` in the `use brust::{...}` line).

**Interfaces:**
- Produces: `brust_seq::{complement(base: u8) -> u8, reverse_complement(seq: &[u8]) -> Vec<u8>, iupac_matches(code: u8, base: u8) -> bool, translate_codon(codon: &[u8]) -> Option<u8>, translate(seq: &[u8]) -> Vec<u8>}`, reachable as `brust::seq::...`.

- [ ] **Step 1: Write the failing tests.**

In `nucleotide.rs`:

```rust
#[test] fn complement_pairs_every_iupac_code_in_both_cases()
// b"ACGTRYKMBVDHSWN" maps bytewise to b"TGCAYRMKVBHDSWN"; the lowercase input maps to the lowercase output
#[test] fn complement_maps_u_to_a_and_keeps_gaps()          // U→A, u→a, '-'→'-', '.'→'.'
#[test] fn complement_turns_other_bytes_into_n()           // b'X', b'*', b'=', b'1', 0xC3 → b'N'
#[test] fn reverse_complement_keeps_case_and_reverses()    // b"acgtRykm" → b"kmrYacgt"; b"" → b""
#[test] fn reverse_complement_twice_is_identity()          // for b"ACGTRYKMBVDHSWNacgtrykmbvdhswn-."
#[test] fn iupac_memberships_match_nc_iub()
// Port NPTune's table test (~/NPTune/src/assay.rs:180-205): all 15 codes against b"ACGTN" in both cases
#[test] fn iupac_matches_treats_u_as_t()
// (b'U', b'T'), (b'T', b'u') and (b'W', b'U') → true; (b'N', b'N') and (b'X', b'A') → false
```

In `codon.rs`:

```rust
#[test] fn standard_code_all_64_codons_and_20_amino_acid_identities()
// Port NPTune's test (~/NPTune/src/sequence.rs:118-176) unchanged, including its reference comments
#[test] fn translate_codon_accepts_lowercase_and_u()      // b"atg"→Some(b'M'); b"AUG"→Some(b'M'); b"uaa"→Some(b'*')
#[test] fn translate_codon_rejects_ambiguous_and_wrong_length()   // b"NNK", b"AT", b"ATGA", b"" → None
#[test] fn translate_reads_frame_zero()
// b"ATGGCCTAAGGN" → b"MA*X"; b"ATGGC" → b"M"; b"" → b""
```

In `brust/tests/facade.rs`:

```rust
#[test] fn facade_reexports_seq_helpers()   // brust::seq::reverse_complement(b"AACG") == b"CGTT"; brust::seq::translate(b"ATG") == b"M"
```

- [ ] **Step 2: Run the tests and confirm they fail.** Run `cargo test -p brust-seq; cargo test -p brust --test facade`. Expected: they don't compile (the crate and functions do not exist yet).
- [ ] **Step 3: Implement the functions.** Use `match` tables; the byte-level rules are in the spec's Nucleotides and Codons sections. The whole crate carries `#![forbid(unsafe_code)]`.
- [ ] **Step 4: Run the checks.** Run `cargo test -p brust-seq && cargo test -p brust && cargo clippy --workspace --all-targets && cargo doc -p brust-seq --no-deps && cargo package --list -p brust-seq`. Expected: every test passes, there are no warnings, and the package list includes `src/nucleotide.rs`, `src/codon.rs` and `README.md`.
- [ ] **Step 5: Commit.**

```bash
git add Cargo.toml Cargo.lock seq brust/Cargo.toml brust/src/lib.rs brust/tests/facade.rs README.md
git commit -m "feat(seq): brust-seq crate with IUPAC complement, matching and codon translation"
```

---

### Task 3: Phred means and `per_read_qscore`

**Files:**
- Create: `seq/src/quality.rs` (with tests). Modify `seq/src/lib.rs` (`mod quality;` plus its re-exports) and `seq/README.md`.
- Modify: `brust/src/stats.rs`:
  - `QualityStats` (around lines 185-207): add `#[non_exhaustive]` and the new field, and clarify the `per_read_mean_phred` doc;
  - `QualityAccumulator` (around lines 876-931);
  - `write_quality_stats` (around lines 1440-1482).
- Test: `brust/tests/stats_validate.rs`, `brust/tests/cli.rs`.
- Modify: `README.md`, the Statistics list, which mentions the per-read qscore.

**Interfaces:**
- Consumes: the `seq` dependency of `brust` from Task 2.
- Produces:
  - `brust_seq::PhredMean` (`#[derive(Debug, Clone, Default)]`) with `add(&mut self, phred: f64)` and `mean(&self) -> Option<f64>`;
  - `brust_seq::read_mean_phred(phred33: &[u8]) -> Option<f64>`;
  - `brust::stats::QualityStats::per_read_qscore: FloatStats`.

- [ ] **Step 1: Write the failing tests.**

In `quality.rs`:

```rust
#[test] fn phred_mean_matches_worked_values()
// Port NPTune's cases (~/NPTune/src/statistics.rs:65-83): [] → None; [10] → 10; [30;5] → 30; [93;2] → 93;
// [0;2] → 0; [0, 40] → 3.0098656839; [10, 30] → 12.9670862188 (tolerance 1e-9)
#[test] fn read_mean_phred_equals_phred_mean_for_every_score()
// for p in 0..=93: read_mean_phred(&[33 + p, 33 + (93 - p)]) equals PhredMean{add(p), add(93 - p)}.mean() within 1e-12
#[test] fn read_mean_phred_decodes_like_the_stats_code()
// b"" → None; b"\x20!" → Some(0.0) (bytes below 33 decode as 0); the b"!!" result is exactly 0.0 and positive
#[test] fn a_long_read_keeps_precision()                  // Review Focus 3: vec![b'5'; 1_000_000] (Q20) → 20 within 1e-9
```

In `brust/tests/stats_validate.rs`, each test writes a FASTQ into `common::TempDir` and uses `brust::stats::stats(Format::Fastq, path)`:

```rust
#[test] fn per_read_qscore_of_uniform_reads_is_their_phred()   // 3 reads of "IIII" (Q40): min = max = mean = 40
#[test] fn per_read_qscore_averages_error_probabilities()      // one read, quality "!I" → 3.0098656839 (1e-9)
#[test] fn per_read_qscore_matches_independent_values_on_the_fixture()
// fastq/UDP0057_sub100.fastq: count 100; min 6.37049031904244, max 7.616430289220246,
// mean 7.229083912490575 (1e-9; computed in Python, not with brust-seq); every
// per_read_qscore summary value ≤ its per_read_mean_phred counterpart
#[test] fn empty_fastq_has_empty_qscore_summary()          // Review Focus 1: count 0, min/max/mean None
#[test] fn zero_length_records_are_left_out_of_qscore()     // Review Focus 2: records "" and "IIII" → count 1, mean 40, no NaN
```

In `brust/tests/cli.rs`:

```rust
#[test] fn stats_cli_prints_per_read_qscore()
// `stats fastq fastq/UDP0057_sub100.fastq` stdout contains "per_read_qscore:"; an empty FASTQ also succeeds
```

- [ ] **Step 2: Run the tests and confirm they fail.** Run `cargo test -p brust-seq quality; cargo test -p brust --test stats_validate --test cli`. Expected: they don't compile (`quality`, `PhredMean` and `per_read_qscore` do not exist yet).
- [ ] **Step 3: Implement.** `read_mean_phred` sums `ERROR[p]` from a `static ERROR: LazyLock<[f64; 256]>` holding `10^(-p/10)`. `PhredMean::add` uses `powf`. Both return `(-10 * log10(sum / n)).max(0.0)`. In the stats, `QualityAccumulator::push` adds `seq::read_mean_phred(quality.as_bytes())` (the facade's alias for `brust-seq`) to a new `per_read_qscore: FloatAccumulator` only when the read has bases. The display block goes straight after `per_read_mean_phred`, through `write_float_stats`.
- [ ] **Step 4: Run the checks.** Run `cargo test --workspace && cargo clippy --workspace --all-targets && cargo doc --workspace --no-deps`. Expected: everything passes, including the existing `stats_validate.rs` tests, with no warnings.
- [ ] **Step 5: Commit.**

```bash
git add seq brust/src/stats.rs brust/tests/stats_validate.rs brust/tests/cli.rs README.md
git commit -m "feat(stats): per_read_qscore from error-space Phred means"
```

---

### Task 4: Original-orientation record helpers

**Files:**
- Modify: `sam/Cargo.toml` and `bam/Cargo.toml` (the `seq` dependency).
- Modify: `sam/src/lib.rs`: in the `impl SamRecord` block at line 171, next to `is_reverse_complemented`, add the two methods and their tests in the existing test module.
- Modify: `bam/src/lib.rs`: in the `impl BamRecord` block, next to `sequence_string`/`quality_string` (around line 855), add the two methods and their tests in `bam_tests`.
- Modify: `sam/README.md` and `bam/README.md` (API list).

**Interfaces:**
- Consumes: `seq::reverse_complement` from Task 2.
- Produces:
  - `SamRecord::original_seq(&self) -> String` and `SamRecord::original_qual(&self) -> String`;
  - `BamRecord::original_sequence_string(&self) -> String` and `BamRecord::original_quality_string(&self) -> String`.
  - Flag 0x10: SEQ is reverse-complemented and QUAL reversed. Otherwise both are unchanged, and `*` stays `*`.

- [ ] **Step 1: Write the failing tests.**

In `sam/src/lib.rs`, build records with `SamRecord::new` or by parsing lines:

```rust
#[test] fn original_orientation_leaves_forward_records_unchanged()   // flag 0, "ACgtN"/"!#%'+" → unchanged
#[test] fn original_orientation_reverses_reverse_strand_records()    // flag 16, "AACGr"/"!#%')" → "yCGTT"/")'%#!"
#[test] fn original_orientation_keeps_star()                          // flag 16, seq "*", qual "*" → "*", "*"
```

In `bam/src/lib.rs` (`bam_tests`), build records through `SamToBamConverter` from the same SAM records:

```rust
#[test] fn bam_original_orientation_matches_sam()
// for flags 0 and 16: original_sequence_string/original_quality_string == the SAM record's original_seq/original_qual
#[test] fn bam_reverse_record_without_qualities_reports_star()
// Review Focus 4: flag 16, QUAL "*" (stored as 0xff) → original_quality_string() == "*"
```

- [ ] **Step 2: Run the tests and confirm they fail.** Run `cargo test -p brust-sam original; cargo test -p brust-bam original`. Expected: they don't compile (the methods do not exist yet).
- [ ] **Step 3: Implement.** Each method checks `is_reverse_complemented()` and the `*` case, then calls `seq::reverse_complement` on the bytes. Build the `String` with `String::from_utf8(..).expect("reverse_complement output is ASCII")`. The BAM methods start from `sequence_string()` and `quality_string()`.
- [ ] **Step 4: Run the checks.** Run `cargo test --workspace && cargo clippy --workspace --all-targets`. Expected: everything passes, including `bam/tests/golden_bytes.rs`, with no warnings.
- [ ] **Step 5: Commit.**

```bash
git add Cargo.lock sam bam
git commit -m "feat(sam,bam): original-orientation sequence and quality helpers"
```

---

### Task 5: Strand-aware `sam-to-fastq` and `bam-to-fastq`

**Files:**
- Modify: `brust/src/convert.rs`:
  - `sam_record_to_fastq` (around line 300) returns `Result<Option<fastq::FastqRecord>>`;
  - `sam_to_fastq` and `bam_to_fastq` write only `Some` records;
  - rustdoc on both functions and on `Conversion::{SamToFastq, BamToFastq}`.
- Test: `brust/tests/convert.rs`.
- Modify: `README.md`, the Convert section, with one line on the new rules and that output differs from 0.1.x for reverse-strand, secondary and supplementary records.

**Interfaces:**
- Consumes: `sam::SamRecord::{original_seq, original_qual}` from Task 4. `bam_to_fastq` already converts each BAM record with `to_sam_record`, so the SAM methods cover both paths.
- Produces: no new public API.

- [ ] **Step 1: Write the failing tests.** The shared input for the synthetic tests is a SAM with an `@SQ SN:ref LN:100` header and four records:
  - `fwd`, flag 0, SEQ `ACGT`, QUAL `!#%'`;
  - `rev`, flag 16, SEQ `AACG`, QUAL `!#%'`;
  - `sec`, flag 256, SEQ `*`, QUAL `*`;
  - `sup`, flag 2048, SEQ `ACGT`, QUAL `IIII`.

  Mapped records use RNAME `ref`, POS 1 and CIGAR `4M` (for example `5M` for a 5-base read), so they pass SAM validation and `SamToBamConverter`.

```rust
#[test] fn sam_to_fastq_writes_primary_reads_in_original_orientation()
// output records exactly: ("fwd", "ACGT", "!#%'"), ("rev", "CGTT", "'%#!")
#[test] fn bam_to_fastq_writes_primary_reads_in_original_orientation()   // same input via sam_to_bam then bam_to_fastq; same two records
#[test] fn fixture_reverse_strand_reads_come_out_reverse_complemented()
// sam/aligned.sam through sam_to_fastq: 100 records. Each flag-16 record equals the SAM SEQ reverse-complemented with a
// local ACGTN-only table in the test (not brust::seq), and QUAL reversed. The 68 forward records are unchanged.
#[test] fn soft_masked_iupac_reverse_read_keeps_case()       // Review Focus 5: flag 16, SEQ "acgRN", QUAL "!!!!!" → SEQ "NYcgt"
#[test] fn bam_to_fastq_rejects_primary_reverse_read_without_qualities()
// Review Focus 4: flag 16, QUAL "*" via sam_to_bam then bam_to_fastq → Err whose message contains
// "cannot convert record without QUAL to FASTQ"
```

The existing tests for a `*` SEQ or QUAL (`sam_to_fastq_rejects_records_without_quality_scores` in `convert.rs`, and `failed_conversion_preserves_existing_output_file` in `brust/tests/convert.rs`) and the fixture tests must still pass unchanged.
- [ ] **Step 2: Run the tests and confirm they fail.** Run `cargo test -p brust --test convert`. Expected: FAIL. `rev` still comes out as `AACG`, and `sec` fails with "cannot convert record without SEQ to FASTQ".
- [ ] **Step 3: Implement.** The skip check is `record.flag & 0x900 != 0` and comes first. Then keep today's `*` checks, then call `FastqRecord::new(record.qname.clone(), None, record.original_seq(), record.original_qual())`.
- [ ] **Step 4: Run the checks.** Run `cargo fmt --check && cargo clippy --workspace --all-targets && cargo test --workspace && cargo doc --workspace --no-deps`. Expected: clean, with every test passing.
- [ ] **Step 5: Commit.**

```bash
git add brust/src/convert.rs brust/tests/convert.rs README.md
git commit -m "fix(convert): write SAM/BAM reads to FASTQ in original orientation; skip secondary and supplementary"
```

---

## Final verification (after Task 5)

- [ ] Run `cargo fmt --check && cargo clippy --workspace --all-targets && cargo test --workspace && cargo doc --workspace --no-deps && cargo package --list -p brust-seq`. Expected: clean, with every test passing.
- [ ] Run the Task 1 `cargo metadata` command. Expected: every package, including `brust-seq`, is at `0.2.0`.
- [ ] Run `cargo run -q -p brust -- stats fastq fastq/UDP0057_sub100.fastq | grep -A5 per_read_qscore` and check that the mean shown is about 7.229.
