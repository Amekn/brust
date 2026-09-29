# Parallel BGZF writer — design

Date: 29/09/2026
Status: approved in conversation, awaiting written-spec review
Branch: `feat/parallel-bgzf-writer`

## Background

A review of NPTune's Rust code (`~/NPTune/src`) for code worth backporting into brust found one
clear gap in `brust-bam`. NPTune's `src/bam_blocks.rs` compresses BAM BGZF blocks in parallel
and produces output byte-identical to `brust::bam::BamWriter`. To do so it had to copy
brust-bam's private serialisers (`write_bam_record`, `encode_bam_auxiliary`,
`write_bgzf_block`), because brust exposes neither record encoding nor BGZF block framing, and
`BamWriter` compresses every block on the calling thread.

This spec covers that one item. The other accepted candidate from the review (sequence helpers
and error-space read mean Q in FASTQ stats) gets its own spec. Paper-specific NPTune code
(Fc alignment, evidence, UMI families, figures, Bonito training helpers, SPOA) is out of scope
for brust.

## Goals

1. Faster BAM writing through multi-threaded BGZF compression, in the library and the CLI.
2. Output byte-identical for every thread count, and identical to the current `BamWriter`.
3. Public building blocks (record encoding, block framing, a general BGZF writer) so callers
   such as NPTune no longer copy private code.
4. No new dependencies and no breaking changes to existing public API.

## Non-goals

- Public BAM header encoding.
- A compression-level setting.
- BGZF output for FASTQ `.gz` (it stays plain gzip).
- Parallel BAM or BGZF reading.
- Version bumps, publishing, or changes to NPTune.

## Design

### Module layout

A new module `bam/src/bgzf.rs` holds all BGZF code. The existing reader side moves there
unchanged: `BgzfReader`, `BgzfVirtualOffset`, block parsing, `bgzf_bsize` and the constants.
The writer side is new. The crate root re-exports every item that is public today, so paths
such as `brust::bam::BgzfReader` and `brust_bam::BgzfVirtualOffset` keep working. `lib.rs`
loses roughly 250 lines.

### Public API

```rust
pub mod bgzf {
    /// Uncompressed bytes per full BGZF block (the current BGZF_MAX_UNCOMPRESSED_BLOCK).
    pub const MAX_BLOCK_DATA: usize = 64 * 1024 - 256; // 65_280

    /// The 28-byte BGZF end-of-file marker block (the current BGZF_EOF_BLOCK bytes).
    pub const EOF_BLOCK: [u8; 28] = *b"\x1f\x8b\x08\x04...";

    /// Deflate `data` (at most MAX_BLOCK_DATA bytes) as one BGZF block and append it to `out`,
    /// framed exactly as BamWriter frames blocks today.
    pub fn compress_block(data: &[u8], out: &mut Vec<u8>) -> io::Result<()>;

    pub struct BgzfWriter<W: Write> { /* ... */ }

    impl<W: Write> BgzfWriter<W> {
        /// Compresses on the calling thread.
        pub fn new(inner: W) -> Self;
        /// `threads` compression workers; 0 or 1 compresses on the calling thread.
        pub fn with_threads(inner: W, threads: usize) -> Self;
        /// Writes any partial block and the EOF marker, joins workers, returns `inner`.
        pub fn finish(self) -> io::Result<W>;
    }

    impl<W: Write> Write for BgzfWriter<W> { /* write, flush */ }
}
```

`BgzfWriter` is re-exported at the crate root beside `BgzfReader`.

`BamRecord` gains:

```rust
/// Append this record as BamWriter writes it: `block_size` (u32 LE), then the payload.
/// Applies the same length and tag validation as `BamWriter::write_record`.
pub fn encode(&self, out: &mut Vec<u8>) -> io::Result<()>;
```

`BamWriter<W: Write>` wraps a `BgzfWriter<W>` instead of its own `pending` buffer and gains:

```rust
impl BamWriter<File> {
    pub fn from_path_with_threads<P: AsRef<Path>>(path: P, threads: usize) -> io::Result<Self>;
}
impl<W: Write> BamWriter<W> {
    pub fn from_writer_with_threads(writer: W, threads: usize) -> Self;
}
```

`from_path`, `from_writer`, `write_header`, `write_record`, `write`, `write_all`, `flush` and
`finish` keep their signatures and behaviour. `W` gains no new bounds (no `Send`, no `'static`).

### Threading model (`threads >= 2`)

- The writer owns `threads` std worker threads. Each has a bounded job channel
  (`std::sync::mpsc::sync_channel`) and a result channel.
- Blocks get sequence numbers. Block `k` goes to worker `k % threads`, and its compressed
  result comes back on that worker's result channel. Workers handle jobs first in, first out,
  so reading results round-robin from worker `0, 1, 2, …` gives stream order with no reorder
  buffer and no shared mutex.
- Compressed blocks are written to `W` on the caller's thread only, during `write`, `flush` and
  `finish`. That is why `W` needs no `Send` bound.
- At most `2 × threads` blocks are in flight in total, and each worker's job channel has
  capacity 2. Because of round-robin assignment, one worker never holds more than 2 blocks, so
  a send never blocks on a full channel. When the limit is reached, the writer waits for the
  oldest block and writes it before sending another. Memory is bounded at roughly
  `2 × threads × 64 KiB` of uncompressed data plus the compressed results.
- Jobs and results are owned `Vec<u8>` buffers. Reusing buffers is a possible later
  optimisation, not part of this design.

### Data flow

- `write(buf)`: bytes fill the current chunk up to `MAX_BLOCK_DATA`. Each full chunk is sent
  (or compressed inline when single-threaded). After each send, results already finished at
  the head of the queue are written to `W` without blocking.
- `flush()`: a non-empty partial chunk is sent as its own block (the current `BamWriter::flush`
  behaviour). The writer waits for and writes every in-flight block, then calls `W::flush()`.
- `finish()`: `flush()`, then `EOF_BLOCK` is written. Workers are shut down and joined, and `W`
  is returned.

### Errors

- Compression errors (for example a compressed block over 64 KiB) are returned through the
  result channel. `W` errors happen on the caller's thread. Both surface as `io::Error` from
  the `write`, `flush` or `finish` call that meets them.
- The first error poisons the writer. That call and every later call return an error, so
  output with a missing block can never be finished.
- If a worker panics, its channels disconnect, and the caller gets
  `io::Error("BGZF worker stopped")` rather than hanging.
- Format crates keep returning `io::Result`.

### Drop

Dropping a `BgzfWriter` without calling `finish()` closes the job channels and joins the
workers, so no threads outlive it. Unwritten data is discarded, as it is with the current
`BamWriter`. The rustdoc states that `finish()` is required for a complete stream.

### Guarantees

- For the same sequence of `write`, `flush` and `finish` calls, output bytes are identical for
  every thread count, and identical to the current `BamWriter`.
- Block boundaries depend only on the uncompressed stream and the `flush` calls, never on
  thread timing.
- Caveat: deflate output depends on flate2's backend. Bytes match within a build. Switching
  flate2 to another backend (for example zlib-ng) changes the bytes, but not the thread-count
  guarantee. This is documented on `BgzfWriter`.

## Facade and CLI

### Library (`brust::convert`)

```rust
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct ConvertOptions {
    /// BGZF compression threads for BAM outputs (default 1).
    pub threads: usize,
}
impl Default for ConvertOptions { /* threads: 1 */ }
impl ConvertOptions {
    pub fn threads(self, threads: usize) -> Self;
}

pub fn convert_with<I: AsRef<Path>, O: AsRef<Path>>(
    conversion: Conversion,
    input: I,
    output: O,
    options: &ConvertOptions,
) -> Result<()>;
```

- `convert()` calls `convert_with` with the default options. The public per-conversion
  functions (`sam_to_bam` and the rest) keep their signatures and defaults.
- Inside, the FASTQ-to-BAM and SAM-to-BAM paths take the options and open
  `BamWriter::from_path_with_threads`.
- `threads` has no effect on conversions that don't write BAM. This is documented.
- Atomic temp-file-then-rename output is unchanged.
- `ConvertOptions` is re-exported from the `brust` crate root beside `Conversion`.

### CLI

```text
brust convert sam-to-bam   aligned.sam aligned.bam --threads 8
brust convert fastq-to-bam reads.fastq reads.bam  -t 8
```

- `-t/--threads <N>` appears on the `sam-to-bam` and `fastq-to-bam` subcommands only.
- Default 1. clap rejects 0 (range `1..`). Help text: "BGZF compression threads".

## Testing

Written test-first, in this order:

1. **Pin current bytes (before any refactor).** A test in `bam` re-encodes `bam/aligned.bam`,
   `bam/unaligned.bam` and `sam/aligned.sam` (through `SamToBamConverter`) with the current
   `BamWriter`, and asserts each output's length and CRC32 (`flate2::Crc`). The expected values
   are captured from the unmodified code, and the test must pass on it. It then guards both
   the module move and the new writer.
2. **Thread-count invariance (from NPTune `bam_blocks.rs`).** Output of
   `BgzfWriter::with_threads(n)` and `BamWriter::from_writer_with_threads(n)` for
   n ∈ {0, 1, 2, 3, 8} equals the inline output for each of:
   - a header-only stream
   - a header longer than one block (thousands of `@SQ` lines)
   - a stream ending exactly on a block boundary (no empty block, one EOF marker)
   - incompressible pseudo-random sequence and quality bytes near the 64 KiB limit
   - a record larger than one block
   - every auxiliary scalar and array type
   - enough blocks to exceed the `2 × threads` in-flight limit
3. **Writer behaviour.**
   - `flush()` mid-stream gives identical bytes for every thread count.
   - An inner writer that fails after N bytes: the error is returned, later calls also fail,
     and no call hangs.
   - Drop without `finish()` returns promptly with the workers joined.
   - Output decodes through both `BgzfReader` and `flate2::read::MultiGzDecoder`.
4. **`BamRecord::encode`.** Its bytes equal the decompressed record bytes `BamWriter` writes,
   for every auxiliary type, and it gives the same validation errors (name length, `n_cigar_op`,
   `l_seq`, quality length, tag syntax, `block_size` mismatch).
5. **Facade and CLI.**
   - `convert_with(.., &ConvertOptions::default().threads(4))` and `--threads 4` produce files
     byte-identical to the defaults for `sam-to-bam` and `fastq-to-bam`.
   - `--threads 0` exits non-zero.

### Verification before completion

- `cargo fmt --check`
- `cargo clippy --workspace --all-targets`
- `cargo test --workspace`
- `cargo doc --workspace --no-deps`
- A timing of `brust convert sam-to-bam` at `--threads 1` and `--threads 8` on a synthetic SAM
  of a few hundred MB, to confirm the speed-up. The result goes in the final report, not in the
  test suite.

## Documentation

- `bam/README.md`: `BgzfWriter`, `compress_block`, `BamRecord::encode` and the threaded
  constructors, with a short example.
- Root `README.md`: `--threads` under Convert, and `ConvertOptions` under Library API.
- Rustdoc for every new public item, including the determinism caveat and the
  `finish()` requirement.

## Follow-up (not in this work)

- Spec 2: sequence helpers (IUPAC reverse complement and matching, standard codon translation)
  and an error-space per-read mean Q field in FASTQ stats.
- NPTune: once a brust release includes this, `src/bam_blocks.rs` can be replaced with
  `BamWriter::from_writer_with_threads` in a later NPTune commit. The paper's recorded run stays
  pinned to its existing commit.
