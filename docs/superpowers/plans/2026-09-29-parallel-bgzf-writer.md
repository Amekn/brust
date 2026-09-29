# Parallel BGZF Writer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Multi-threaded BGZF compression for brust's BAM output, in the library and the CLI, with output byte-identical to today's `BamWriter` at any thread count.

**Architecture:** BGZF code moves into a new `bam/src/bgzf.rs`. It gains a public `compress_block` and a `BgzfWriter` that compresses inline or on writer-owned std worker threads, with blocks assigned round-robin and written back on the caller's thread. `BamWriter` wraps `BgzfWriter` and gains `*_with_threads` constructors. The facade passes a `ConvertOptions { threads }` through to the BAM-writing conversions, and the CLI exposes `-t/--threads`.

**Tech Stack:** Rust 2024 (1.96), std threads and `std::sync::mpsc`, flate2 (already a dependency), clap 4.

**Spec:** `docs/superpowers/specs/2026-09-29-parallel-bgzf-writer-design.md`

## Global Constraints

- No new dependencies or dev-dependencies in any crate.
- Existing public signatures do not change. `brust_bam::BgzfReader` and `brust_bam::BgzfVirtualOffset` stay reachable at the crate root.
- Format crates return `std::io::Result`; the `brust` facade returns `brust::Result`.
- Bytes written are identical for every thread count, and identical to the pre-change `BamWriter`, given the same sequence of `write`, `flush` and `finish` calls.
- `MAX_BLOCK_DATA = 65_280` (64 × 1024 − 256).
- `threads` of 0 or 1 means inline compression. Workers are capped at `MAX_THREADS = 256`, a plan decision that makes absurd values safe (it does not affect output).
- At most `2 × workers` blocks are in flight. Each worker's job and result channels are `sync_channel(2)`.
- The CLI `-t/--threads` flag exists only on `sam-to-bam` and `fastq-to-bam`: default 1, 0 rejected, help text `BGZF compression threads`.
- Rustdoc and README text match the style of the file being edited. Every new public item has rustdoc.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. A single `write()` of several MB in threaded mode: identical bytes to inline, and no deadlock at the in-flight limit. Test in Task 4.
2. The inner writer failing while the writer waits at the in-flight limit: an error comes back, the writer stays poisoned, and nothing hangs. Test in Task 4.
3. `flush()` after every small write: one block per flush, identical across thread counts. Test in Task 4.
4. An absurd thread count (100_000): clamped to 256 workers, correct output, prompt finish. Test in Task 4.
5. A threaded SAM-to-BAM conversion failing part-way (malformed SAM): an error comes back, the existing output file is untouched, no temp file is left, and nothing hangs. Test in Task 6.

---

### Task 1: Pin today's BAM bytes

A characterisation test. It must pass on the unmodified code, and it guards every later task.

**Files:**
- Create: `bam/tests/golden_bytes.rs`

**Interfaces:**
- Consumes: the current public API (`Bam`, `BamWriter`, `SamToBamConverter`, `sam::SamReader`).
- Produces: nothing new. Later tasks must keep this test passing unchanged.

- [ ] **Step 1: Write the test**

A helper `fn fingerprint(bytes: &[u8]) -> (usize, u32)` returns the length and the CRC32 from `flate2::Crc`. Fixture paths come from `env!("CARGO_MANIFEST_DIR")`. Four tests:

```rust
#[test] fn aligned_bam_rewrite_bytes_are_pinned()      // Bam::from_path("aligned.bam"), bam.to_writer(&mut out)
    { assert_eq!(fingerprint(&out), (55475, 1371713897)); }
#[test] fn unaligned_bam_rewrite_bytes_are_pinned()    // same with "unaligned.bam"
    { assert_eq!(fingerprint(&out), (129915, 112315190)); }
#[test] fn sam_to_bam_stream_bytes_are_pinned()        // SamReader "../sam/aligned.sam" → SamToBamConverter →
    { assert_eq!(fingerprint(&out), (55711, 4017805689)); }  // BamWriter::from_writer(Vec::new()), header, records, finish
#[test] fn flush_after_header_bytes_are_pinned()       // aligned.bam: write_header, flush(), every record, finish
    { assert_eq!(fingerprint(&out), (55425, 2406701901)); }
```

- [ ] **Step 2: Run it on the unmodified code**

Run: `cargo test -p brust-bam --test golden_bytes`
Expected: 4 passed. (These values were captured from commit `492acbd` and matched on two runs. If they do not match, stop and report; do not change the values.)

- [ ] **Step 3: Commit**

```bash
git add bam/tests/golden_bytes.rs
git commit -m "test(bam): pin current BamWriter output bytes"
```

---

### Task 2: `bgzf` module, `compress_block` and `BamRecord::encode`

Moves the existing BGZF code into its own module without changing behaviour, and makes the two encoders public.

**Files:**
- Create: `bam/src/bgzf.rs`
- Modify: `bam/src/lib.rs`
  - move out: `BgzfVirtualOffset` (lines 80–110), the `BgzfReader` struct (121–129), its `impl` and `Read` impl plus `BgzfBlock`, `read_bgzf_block` and `bgzf_bsize` (1307–1470), `BGZF_MAX_UNCOMPRESSED_BLOCK`, `BGZF_EOF_BLOCK` and `write_bgzf_block` (2014–2046)
  - keep in `lib.rs`: `read_exact_or_eof` and `invalid_data`. `bgzf.rs` uses them through `super::`
  - add: `pub mod bgzf;` and `pub use bgzf::{BgzfReader, BgzfVirtualOffset};`
- Test: unit tests in `bam/src/bgzf.rs`, plus the `bam_tests` module in `bam/src/lib.rs`

**Interfaces:**
- Produces:
  - `brust_bam::bgzf::MAX_BLOCK_DATA: usize = 65_280`
  - `brust_bam::bgzf::EOF_BLOCK: [u8; 28]` (the current `BGZF_EOF_BLOCK` bytes)
  - `brust_bam::bgzf::compress_block(data: &[u8], out: &mut Vec<u8>) -> io::Result<()>`: appends one framed block, and returns `ErrorKind::InvalidInput` when `data.len() > MAX_BLOCK_DATA` or when the framed block would exceed 64 KiB
  - `BamRecord::encode(&self, out: &mut Vec<u8>) -> io::Result<()>`: appends `block_size` (u32 LE) then the payload, with the same validation as today's `write_bam_record`

- [ ] **Step 1: Write the failing tests**

In `bgzf.rs`:

```rust
#[test] fn compress_block_frames_a_bgzf_block()
// data: 1_000 bytes of (i % 7) as u8. After compress_block(&data, &mut out):
//   out[..16] == [0x1f,0x8b,0x08,0x04,0,0,0,0,0,0xff,0x06,0x00,b'B',b'C',0x02,0x00]
//   u16::from_le_bytes(out[16..18]) as usize == out.len() - 1
//   out[out.len()-8..out.len()-4] == CRC32(data).to_le_bytes()
//   out[out.len()-4..] == 1_000u32.to_le_bytes()
//   MultiGzDecoder over out decodes to data; BgzfReader over (out ++ EOF_BLOCK) reads data
#[test] fn compress_block_appends_to_existing_output()   // out starts as b"xyz" → still starts with b"xyz"
#[test] fn compress_block_rejects_more_than_max_block_data()
// vec![0; MAX_BLOCK_DATA + 1] → Err with kind InvalidInput; vec![0; MAX_BLOCK_DATA] → Ok
#[test] fn eof_block_is_an_empty_bgzf_block()
// EOF_BLOCK.len() == 28; BgzfReader over EOF_BLOCK reads 0 bytes; MultiGzDecoder decodes it to empty
```

In `lib.rs` `bam_tests`:

```rust
#[test] fn record_encode_matches_writer_bytes()
// For every record of aligned.bam: the concatenated encode() outputs equal the decompressed
// BamWriter stream (header + records) minus its header prefix. The prefix length comes from
// decompressing a header-only BamWriter output.
#[test] fn record_encode_rejects_what_write_record_rejects()
// Start from a valid aligned.bam record and make one change per case: l_read_name = 1 with a
// longer name; n_cigar_op + 1; l_seq + 2; qual.pop(); aux tag "X"; fixed.block_size = 1.
// For each, encode() fails, write_record() fails, and their to_string() values are equal.
```

- [ ] **Step 2: Run to confirm they fail**

Run: `cargo test -p brust-bam bgzf:: ; cargo test -p brust-bam record_encode`
Expected: compile errors (`bgzf` module and `encode` do not exist).

- [ ] **Step 3: Do the move and add the two public encoders**

- Move the listed items into `bgzf.rs`, along with the `flate2` imports they need. Rename the constants to `MAX_BLOCK_DATA` and `EOF_BLOCK` and make them `pub`.
- Turn `write_bgzf_block` into `compress_block`, adding the `MAX_BLOCK_DATA` check. `BamWriter::flush_pending` calls `compress_block` into a scratch `Vec`, then `write_all`s it.
- Turn `write_bam_record` into `BamRecord::encode`. `BamWriter::write_record` calls `record.encode(&mut data)`.
- Add a rustdoc module comment to `bgzf.rs`.

- [ ] **Step 4: Run the whole crate**

Run: `cargo test -p brust-bam && cargo clippy -p brust-bam --all-targets`
Expected: every test passes, including `golden_bytes` and the existing `conversions.rs` (which imports `BgzfVirtualOffset` from the root). No new clippy warnings.

- [ ] **Step 5: Commit**

```bash
git add bam/src/bgzf.rs bam/src/lib.rs
git commit -m "refactor(bam): move BGZF into its own module; public compress_block and BamRecord::encode"
```

---

### Task 3: Inline `BgzfWriter`, and `BamWriter` built on it

**Files:**
- Modify: `bam/src/bgzf.rs` (add `BgzfWriter` and its tests)
- Modify: `bam/src/lib.rs`: the `BamWriter` struct (lines 55–63 before Task 2) becomes `{ writer: BgzfWriter<W>, header_written: bool }`, and the `pending`, `write_uncompressed` and `flush_pending` logic goes
- Modify: `bam/src/lib.rs`: `pub use bgzf::{BgzfReader, BgzfVirtualOffset, BgzfWriter};`

**Interfaces:**
- Consumes: `compress_block`, `MAX_BLOCK_DATA` and `EOF_BLOCK` from Task 2.
- Produces:
  - `pub struct BgzfWriter<W: Write>`
  - `BgzfWriter::new(inner: W) -> Self`
  - `impl<W: Write> Write for BgzfWriter<W>`: `write` fills chunks and emits a block for each full `MAX_BLOCK_DATA` chunk; `flush` emits a non-empty partial chunk as a block, then flushes `inner`
  - `BgzfWriter::finish(self) -> io::Result<W>`: flush, write `EOF_BLOCK`, return `inner`
  - Poisoning: after any call returns `Err`, every later `write`, `flush` or `finish` returns `Err(io::Error::other("BGZF writer failed earlier"))`
  - A test helper `FailingWriter { written: usize, fail_after: usize }` in the `bgzf` test module. It accepts a write only while `written + buf.len() <= fail_after`, and otherwise returns `Err(io::Error::other("injected failure"))`. Task 4 reuses it.

- [ ] **Step 1: Write the failing tests** (in `bgzf.rs`; `block(data)` is a test helper wrapping `compress_block`)

```rust
#[test] fn writer_cuts_blocks_at_max_block_data()
// write 2 * MAX_BLOCK_DATA + 10 bytes in one call, then finish →
// == block(&d[..M]) ++ block(&d[M..2M]) ++ block(&d[2M..]) ++ EOF_BLOCK
#[test] fn writer_output_decodes_to_input_for_any_write_sizes()
// 200_000 bytes written in chunks cycling [1, 7, 1_000, 65_280, 70_000] → MultiGzDecoder == input
#[test] fn flush_emits_the_partial_block()
// write 10 bytes, flush, write 10 more, finish → block(first) ++ block(second) ++ EOF_BLOCK
#[test] fn flush_with_nothing_pending_writes_no_block()   // flush, flush, finish → == EOF_BLOCK
#[test] fn finish_on_an_empty_writer_writes_only_eof()   // finish → == EOF_BLOCK
#[test] fn writer_is_poisoned_after_an_inner_error()
// FailingWriter { fail_after: 100 }: write 3 * MAX_BLOCK_DATA → Err; a second write → Err;
// flush → Err; finish → Err
```

- [ ] **Step 2: Run to confirm they fail**

Run: `cargo test -p brust-bam bgzf::`
Expected: compile error (`BgzfWriter` does not exist).

- [ ] **Step 3: Implement `BgzfWriter` (inline mode) and move `BamWriter` onto it**

`BamWriter::from_writer` builds `BgzfWriter::new(writer)`. `write_header` and `write_record` call `self.writer.write_all(&data)`. `flush` and `finish` delegate. The public `BamWriter` API and behaviour are unchanged.

- [ ] **Step 4: Run the crate**

Run: `cargo test -p brust-bam && cargo clippy -p brust-bam --all-targets`
Expected: every test passes, including all four `golden_bytes` tests.

- [ ] **Step 5: Commit**

```bash
git add bam/src/bgzf.rs bam/src/lib.rs
git commit -m "feat(bam): public BgzfWriter; BamWriter writes through it"
```

---

### Task 4: Threaded `BgzfWriter`

**Files:**
- Modify: `bam/src/bgzf.rs`

**Interfaces:**
- Consumes: Task 3's `BgzfWriter` and `FailingWriter`.
- Produces:
  - `BgzfWriter::with_threads(inner: W, threads: usize) -> Self`: 0 or 1 is identical to `new`; otherwise `min(threads, MAX_THREADS)` workers
  - `pub const MAX_THREADS: usize = 256` in `bgzf`
  - `W` gains no new bounds

- [ ] **Step 1: Write the failing tests**

`within(secs, f)` is a test helper that runs `f` on a spawned thread and panics with "timed out" if `recv_timeout(Duration::from_secs(secs))` expires. `inline(ops)` and `threaded(n, ops)` replay the same list of `Write(len)` and `Flush` operations over deterministic data (an LCG over a 4-letter alphabet, so blocks compress). They return the finished bytes.

```rust
#[test] fn threaded_output_matches_inline_for_every_thread_count()
// ops: 1_000_000 bytes in writes cycling [1, 7, 4_096, 65_280, 100_000], one Flush after the
// third write; for n in [0, 1, 2, 3, 8]: threaded(n, ops) == inline(ops)
#[test] fn threaded_writer_exceeds_the_in_flight_limit()     // n = 2, 40 * MAX_BLOCK_DATA + 5 bytes
#[test] fn one_huge_write_matches_inline()                   // Review Focus 1: n = 2, one 5_000_000-byte write, within(30)
#[test] fn inner_error_at_the_in_flight_limit_is_returned_without_hanging()
// Review Focus 2: within(30): n = 2, FailingWriter { fail_after: 70_000 }, writes of 10_000 up to
// 2_000_000 bytes → some write returns Err; the next write → Err; finish → Err
#[test] fn flush_after_every_write_matches_inline()          // Review Focus 3: 300 × (Write(97), Flush), n in [2, 8]
#[test] fn absurd_thread_counts_are_clamped()                // Review Focus 4: within(30): n = 100_000, 200_000 bytes == inline
#[test] fn dropping_a_threaded_writer_without_finish_returns() // within(30): n = 4, write 500_000, drop
#[test] fn threaded_output_reads_back()                      // n = 3: BgzfReader and MultiGzDecoder both give the input
```

- [ ] **Step 2: Run to confirm they fail**

Run: `cargo test -p brust-bam bgzf::`
Expected: compile error (`with_threads` does not exist).

- [ ] **Step 3: Implement the worker mode**

The mode is private, for example `enum Mode { Inline, Threaded(Pool) }`. Algorithm:

- Each worker holds a job `SyncSender<Vec<u8>>` and a result `Receiver<io::Result<Vec<u8>>>`, both `sync_channel(2)`. Spawn with `std::thread::Builder` named `bgzf-worker-{i}`. If a spawn fails, keep the workers already started; with none started, fall back to inline.
- Counters `sent` and `written` (u64). Block `k` goes to worker `k % workers`. The next block to write always comes from worker `written % workers`.
- After each send, `try_recv` from the head worker and write every ready block in order. Before a send, if `sent - written == 2 * workers`, do a blocking `recv` of the head block first.
- `flush`: send the partial chunk if non-empty, block until `written == sent`, then flush `inner`. `finish`: flush, write `EOF_BLOCK`, close the channels, join the workers (a join error is `io::Error::other("BGZF worker stopped")`), return `inner`.
- A disconnected `recv` is `io::Error::other("BGZF worker stopped")`. Every error poisons the writer (Task 3's rule).
- `Drop` (when `finish` did not run): drop each worker's job sender *and* result receiver, then join it, ignoring join errors.
- Worker loop: `for data in jobs { let mut out = Vec::new(); let r = compress_block(&data, &mut out).map(|()| out); if results.send(r).is_err() { break } }`.

Add rustdoc on `with_threads`: what the thread count means, the clamp, determinism (including the flate2-backend caveat from the spec's Guarantees section), and that `finish()` is required.

- [ ] **Step 4: Run the crate**

Run: `cargo test -p brust-bam && cargo clippy -p brust-bam --all-targets`
Expected: all pass, and no test takes more than a few seconds.

- [ ] **Step 5: Commit**

```bash
git add bam/src/bgzf.rs
git commit -m "feat(bam): multi-threaded BgzfWriter with in-order output"
```

---

### Task 5: Threaded `BamWriter` constructors and NPTune's BAM cases

**Files:**
- Modify: `bam/src/lib.rs` (the `BamWriter` impls)
- Create: `bam/tests/threaded_writer.rs`
- Modify: `bam/README.md` (API list, and a "Parallel compression" example); `bam/src/lib.rs` crate docs (mention `bgzf` and threads)

**Interfaces:**
- Consumes: `BgzfWriter::with_threads` from Task 4, and `BamRecord::encode` from Task 2.
- Produces:
  - `impl BamWriter<File> { pub fn from_path_with_threads<P: AsRef<Path>>(path: P, threads: usize) -> io::Result<Self> }`
  - `impl<W: Write> BamWriter<W> { pub fn from_writer_with_threads(writer: W, threads: usize) -> Self }`

- [ ] **Step 1: Write the failing tests**

Port the helpers `aux`, `record` and `record_of_length` from `~/NPTune/src/bam_blocks.rs:231-340`. The header comes from `SamToBamConverter::new(&header)`, where `header` has `@HD VN:1.6 SO:unsorted` and `@SQ SN:fc LN:687` (as in `~/NPTune/src/alignment.rs:773`, without the `@PG`). `write(threads, header, records)` uses `from_writer_with_threads` and returns the bytes. Each test asserts `write(n, ..) == write(1, ..)` for n in `[0, 2, 3, 8]`, and `write(1, ..)` equals the plain `BamWriter::from_writer` output.

```rust
#[test] fn header_only_stream_matches_for_every_thread_count()
#[test] fn header_longer_than_a_block_matches_for_every_thread_count()      // plus 3_000 @SQ lines "contig-{i}"
#[test] fn stream_ending_exactly_on_a_block_boundary_matches()
// records record_of_length(M - header_len) then record_of_length(M), where M = MAX_BLOCK_DATA and
// header_len is the decompressed header-only length; the decompressed stream length == 2 * M
#[test] fn incompressible_records_match()                 // NPTune's noisy 2 * M record, then "after"
#[test] fn every_aux_type_and_a_record_larger_than_a_block_match()  // NPTune's 17-tag list; 20, 3 * M, 30 bases
#[test] fn many_records_across_many_blocks_match()        // 600 records with NM and MD tags, > 3 blocks
#[test] fn from_path_with_threads_writes_the_same_file()
// std::env::temp_dir()/"brust-bam-threads-{pid}.bam", threads 4; fs::read == from_writer bytes; remove file
```

- [ ] **Step 2: Run to confirm they fail**

Run: `cargo test -p brust-bam --test threaded_writer`
Expected: compile error (`from_writer_with_threads` does not exist).

- [ ] **Step 3: Add the two constructors** as thin wrappers over `BgzfWriter::with_threads`. Rustdoc on both: what `threads` means, and that the output is identical to `from_writer` / `from_path`. Update `bam/README.md` and the crate docs.

- [ ] **Step 4: Run the crate and its docs**

Run: `cargo test -p brust-bam && cargo doc -p brust-bam --no-deps`
Expected: all pass, and the docs build with no warnings.

- [ ] **Step 5: Commit**

```bash
git add bam/src/lib.rs bam/tests/threaded_writer.rs bam/README.md
git commit -m "feat(bam): threaded BamWriter constructors"
```

---

### Task 6: `ConvertOptions` and `convert_with` in the facade

**Files:**
- Modify: `brust/src/convert.rs`
- Modify: `brust/src/lib.rs`: `pub use convert::{Conversion, ConvertOptions};`
- Test: `brust/tests/convert.rs`
- Modify: `README.md`, Library API section (a `ConvertOptions` example)

**Interfaces:**
- Consumes: `bam::BamWriter::from_path_with_threads(path, threads)` from Task 5.
- Produces:
  - `#[derive(Debug, Clone, PartialEq, Eq)] #[non_exhaustive] pub struct ConvertOptions { pub threads: usize }` with `impl Default` giving `threads: 1`, and `pub fn threads(self, threads: usize) -> Self`
  - `pub fn convert_with<I: AsRef<Path>, O: AsRef<Path>>(conversion: Conversion, input: I, output: O, options: &ConvertOptions) -> Result<()>`
  - `convert` delegates with `&ConvertOptions::default()`. The public `fastq_to_bam` and `sam_to_bam` keep their signatures and delegate to private option-taking versions.

- [ ] **Step 1: Write the failing tests** (in `brust/tests/convert.rs`, using `common::TempDir` and `common::fixture`)

```rust
#[test] fn convert_options_default_to_one_thread()
// ConvertOptions::default().threads == 1; ConvertOptions::default().threads(8).threads == 8
#[test] fn threaded_bam_conversions_match_default_bytes()
// (SamToBam, "sam/aligned.sam"), (SamToBam, "sam/unaligned.sam"), (FastqToBam, "fastq/UDP0057_sub100.fastq"):
// fs::read(convert output) == fs::read(convert_with output) for threads 0, 4
#[test] fn threads_do_not_change_non_bam_outputs()          // FastqToFasta with threads 4 == default bytes
#[test] fn failed_threaded_conversion_preserves_existing_output()
// Review Focus 5: input = header and records of sam/aligned.sam plus a final line "bad\tline";
// output pre-written as b"sentinel"; convert_with(SamToBam, threads 4) → Err; output == b"sentinel";
// the temp dir holds only the input and the output
```

- [ ] **Step 2: Run to confirm they fail**

Run: `cargo test -p brust --test convert`
Expected: compile error (`ConvertOptions` and `convert_with` do not exist).

- [ ] **Step 3: Implement it**, with rustdoc stating that `threads` affects only BAM outputs.

- [ ] **Step 4: Run the package**

Run: `cargo test -p brust && cargo clippy -p brust --all-targets`
Expected: all pass, including the existing `facade.rs` and `convert.rs` tests.

- [ ] **Step 5: Commit**

```bash
git add brust/src/convert.rs brust/src/lib.rs brust/tests/convert.rs README.md
git commit -m "feat(brust): ConvertOptions and convert_with for threaded BAM output"
```

---

### Task 7: CLI `-t/--threads`

**Files:**
- Modify: `brust/src/main.rs`
- Test: `brust/tests/cli.rs`
- Modify: `README.md`, CLI Convert section (a `--threads` example, and the note that only `*-to-bam` take it)

**Interfaces:**
- Consumes: `brust::convert::convert_with` and `brust::ConvertOptions` from Task 6.
- Produces: `ConvertCommands::FastqToBam` and `ConvertCommands::SamToBam` gain a `threads: u32` field:
  `#[arg(short = 't', long, default_value_t = 1, value_parser = clap::value_parser!(u32).range(1..))]` with the doc comment `/// BGZF compression threads`. `into_parts` returns `(Conversion, PathBuf, PathBuf, ConvertOptions)`, and `run` calls `convert_with`.

- [ ] **Step 1: Write the failing tests** (in `brust/tests/cli.rs`)

```rust
#[test] fn convert_cli_threads_match_default_bytes()
// sam-to-bam sam/aligned.sam with "--threads", "4" vs no flag → success, identical file bytes;
// fastq-to-bam fastq/UDP0057_sub100.fastq with "-t", "4" vs no flag → identical
#[test] fn convert_cli_rejects_zero_threads()
// sam-to-bam ... --threads 0 → !status.success(), stderr non-empty, output file does not exist
#[test] fn convert_cli_threads_only_on_bam_outputs()        // fastq-to-fasta ... --threads 2 → !success
#[test] fn convert_cli_help_describes_threads()             // "convert sam-to-bam --help" stdout contains "BGZF compression threads"
```

- [ ] **Step 2: Run to confirm they fail**

Run: `cargo test -p brust --test cli`
Expected: FAIL (`--threads` is an unexpected argument).

- [ ] **Step 3: Implement it**, and update `README.md`.

- [ ] **Step 4: Run the package**

Run: `cargo test -p brust`
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add brust/src/main.rs brust/tests/cli.rs README.md
git commit -m "feat(cli): --threads for sam-to-bam and fastq-to-bam"
```

---

## Final verification (after Task 7)

- [ ] `cargo fmt --check && cargo clippy --workspace --all-targets && cargo test --workspace && cargo doc --workspace --no-deps`. Expected: clean, all tests pass, no new warnings.
- [ ] Timing, on a release build, in the scratchpad (not committed):

```bash
cargo build --release -p brust
S=<scratchpad>; grep '^@' sam/aligned.sam > $S/big.sam; grep -v '^@' sam/aligned.sam > $S/body.sam
for i in $(seq 2000); do cat $S/body.sam; done >> $S/big.sam      # about 300 MB
time target/release/brust convert sam-to-bam $S/big.sam $S/t1.bam
time target/release/brust convert sam-to-bam $S/big.sam $S/t8.bam --threads 8
cmp $S/t1.bam $S/t8.bam && echo identical
```

Expected: `identical`, and the 8-thread wall time clearly lower. Report both times; don't just claim a speed-up.
