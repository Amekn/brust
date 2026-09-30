# Atomic writes and strict BAM reading — design

Date: 01/10/2026
Status: approved by the user 01/10/2026, then amended after Codex reviewed the plan
(01/10/2026): close and folder-sync semantics, the strict error's format, and error-once
Branch: `feat/atomic-io`

## Background

Today only the CLI's `brust convert` writes atomically. Its private helper
(`brust/src/convert.rs`, `write_atomic`) writes to a temp file beside the target and renames it
over the target. It has three gaps:

- it never calls `fsync`, so a power cut or OS crash can leave an empty or partial file after
  the rename;
- it only cleans up the temp file when the conversion returns an error, not on a panic;
- library users can't reach it.

The five library writers (`FastaWriter`, `FastqWriter`, `SamWriter`, `BamWriter`,
`Pod5Writer`) open files with `File::create`, which truncates the target and writes into it in
place. Each in-memory type's `to_path` (`Fasta`, `Fastq`, `Sam`, `Bam`, `Pod5`) does the same.

On the reading side, the readers stream the file with no locking. Checked on 01/10/2026:

- gzip FASTQ already rejects truncated input. flate2 1.1.9 returns `UnexpectedEof` even when
  only the last byte is missing.
- POD5 already rejects truncated input. The reader checks the trailing magic and footer
  (`pod5/src/lib.rs`, `inspect_pod5`).
- BAM doesn't require the BGZF EOF block. A BAM cut exactly between two blocks reads as a
  complete, shorter file (`bam/src/bgzf.rs`, `fill_block`).
- No reader can tell a plain FASTA, FASTQ or SAM file cut at a record boundary from a complete
  one.

The user wants protection against four failures:

1. a crash leaves a partial file that looks finished;
2. a reader sees a file while it is being written;
3. power loss or an OS crash after a write returns;
4. truncated input read as complete.

## Decision: one shared atomic file type, not separate writer types

Three approaches were compared.

1. **Separate atomic types per format** (`AtomicBamWriter` and so on, plus a
   `StrictBamReader`). Rejected: five new public types that each repeat every writer method,
   either forwarded by hand (they drift) or through `Deref` (hides what commit does), and a
   reader type that copies the whole BAM reader API to change one check.
2. **Atomic state inside each existing writer.** Rejected: the logic is written five times, and
   deleting the temp file on drop needs a `Drop` impl on existing public generic types. That
   changes drop-check and stops `into_inner` and `finish` moving fields out, which is an API
   break.
3. **Chosen: one `AtomicFile` type in brust-core, plus small helpers on each writer.** Every
   writer is already generic over `W: Write`, so `AtomicFile` works with all of them through
   `from_writer`, and no existing type changes. The helpers add `from_path_atomic` and `commit`
   for convenience. Readers get an opt-in check, not new types.

## Goals

1. `brust_core::AtomicFile`: a buffered `Write` sink that commits with fsync, rename and a
   folder fsync, and deletes its temp file if dropped without commit.
2. `from_path_atomic` and `commit` on all five writers, `from_path_atomic_with_threads` on
   `BamWriter`, and `to_path_atomic` on all five in-memory types.
3. An opt-in strict EOF-block check on `BgzfReader` and `BamReader`.
4. `brust::convert`, `brust::validate` and `brust::stats` use the new APIs: atomic output for
   every conversion, and the strict check on every BAM input. The private `write_atomic`,
   `create_temp_output_path` and `TEMP_COUNTER` in `convert.rs` are removed.
5. Regression tests pinning today's truncation errors for gzip FASTQ and POD5.
6. Every workspace crate moves to version 0.3.0, which the user approved on 01/10/2026.

## Non-goals

- Changing the behaviour of any existing constructor or reader default in the format crates
  (`from_path`, `new`, `from_writer`, `to_path`, `from_reader`). Everything new there is
  opt-in. Only the facade's `convert`, `validate` and `stats` change behaviour (see Facade and
  CLI).
- File locking or any other coordination between processes.
- Checking for the EOF block at open by seeking to the end.
- A way to turn fsync off.
- Keeping the old file's permissions, owner, symlink or hard links.
- Cleaning up temp files after Ctrl-C or a kill (no signal handler).
- Detecting truncation in plain FASTA, FASTQ or SAM.
- Re-exporting `AtomicFile` from the format crates.
- `new_atomic` aliases, and a FASTQ `from_path_atomic_with_compression`.
- Publishing to crates.io, and changes to NPTune.

## Design

### `brust_core::AtomicFile`

Added to `brust-core/src/lib.rs`, which keeps it dependency-free. Every format crate and the
facade already depend on brust-core.

```rust
#[derive(Debug)]
pub struct AtomicFile { /* BufWriter<File>, temp path, target path, committed flag */ }

impl AtomicFile {
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<Self>;
    pub fn commit(self) -> io::Result<()>;
}

impl Write for AtomicFile { /* write and flush go to the temp file through the buffer */ }
impl Drop for AtomicFile { /* not committed: delete the temp file, ignore errors */ }
```

`AtomicFile` is `Send + Sync`, so it works with `BamWriter::from_writer_with_threads`.

**`create(path)`**

- The temp file goes in the target's parent folder, or `.` when the path has no parent part. Its
  name is `.{pid}.{counter}.tmp.{file_name}`, the same scheme `convert` uses today. `counter`
  is a process-wide `AtomicU64` in brust-core. The target's full file name comes last, so
  suffix checks such as `.gz` and `.bam` still work on the temp path.
- The temp file is opened with `create_new(true)`, so it never overwrites anything. If the name
  exists, the counter moves on, up to 100 tries; after that `create` returns an error.
- A path with no file name returns an `InvalidInput` error. A missing parent folder returns the
  open error (`NotFound`). Either way nothing is created.
- The target is not opened or touched.
- Writes go through a `BufWriter<File>` with the std default capacity. `Write::flush` pushes the
  buffer to the temp file only; it does not fsync.

**`commit(self)`**, in this order:

1. Flush the buffer and take back the `File`.
2. `sync_all` the temp file.
3. Close the temp file. This is unchecked: std ignores close errors, and the data is already
   synced.
4. `fs::rename` the temp file over the target. The new file is published here.
5. On Unix, open the parent folder and `sync_all` it. Other platforms skip this step.

If step 1, 2 or 4 fails, the temp file is deleted (best effort), the target is left as it was,
and the original error is returned. A failed flush never retries writing the buffer. If step 5
fails, whether opening or syncing the folder, the new file is already in place. The error keeps
the original `io::ErrorKind`, and its message says the file was renamed but the folder sync
failed.

**Drop without commit** covers returning an error, leaving early with `?`, and a panic that
unwinds. The temp file is deleted and errors are ignored. Nothing appears at the target, and an
existing target is unchanged. A kill, Ctrl-C (the default handler doesn't unwind) or power loss
before commit leaves the temp file behind. The target is untouched, and the `.tmp.` in the name
makes the leftover easy to spot.

**Known limits** (documented in rustdoc and the brust-core README):

- The rename replaces a symlink at the target with a plain file.
- The new file gets default permissions (umask), not the old file's.
- Hard links to the old file keep the old contents.
- Paths are kept as given, so a relative path is resolved against the current folder at each
  step. Don't change the current folder between `create` and `commit`.
- The rename is atomic only within one filesystem. Putting the temp file beside the target
  guarantees that.

### Writer helpers

Each writer crate adds an `impl Writer<AtomicFile>` block. The existing `impl Writer<File>` and
`impl<W: Write> Writer<W>` blocks are unchanged, so every existing method (`write_record`,
`write_header`, `set_line_width`, `flush`, `into_inner`, `finish` and so on) works on an atomic
writer as it is.

| Writer | New constructors | `commit(self) -> io::Result<()>` |
|---|---|---|
| `FastaWriter` | `from_path_atomic(path)` | `self.into_inner().commit()` |
| `SamWriter` | `from_path_atomic(path)` | `self.into_inner().commit()` |
| `Pod5Writer` | `from_path_atomic(path)` | `self.into_inner().commit()` |
| `FastqWriter` | `from_path_atomic(path)` | `self.finish()?.commit()` |
| `BamWriter` | `from_path_atomic(path)`, `from_path_atomic_with_threads(path, threads)` | `self.finish()?.commit()` |

- `FastqWriter::from_path_atomic` chooses compression from the **target** path
  (`Compression::from_path(path)`), then calls
  `from_writer_with_compression(AtomicFile::create(path)?, compression)`.
- `BamWriter::from_path_atomic_with_threads` follows `from_path_with_threads`: 0 or 1 compresses
  on the calling thread, and larger values are capped at `bgzf::MAX_THREADS`.
- `commit` for FASTQ writes the gzip trailer, and for BAM writes the EOF block and joins the
  workers, before committing the file.
- The rustdoc on each `from_path_atomic` states that nothing appears at `path` until `commit()`
  renames the finished file into place, and that dropping the writer discards the output. It
  points to `AtomicFile` for the one error `commit` can return after the rename. The compiler can't force a
  `commit` call, so the docs must be clear. Each existing `from_path` doc gains one line
  pointing to `from_path_atomic`.

**In-memory types.** `Fasta`, `Fastq`, `Sam`, `Bam` and `Pod5` each gain
`to_path_atomic(&self, path) -> io::Result<()>`. It writes the same bytes as `to_path`, through
`from_path_atomic`, then calls `commit`. `Bam`'s private `write_with` finishes the writer and
drops the stream, which would drop the `AtomicFile` uncommitted. So `Bam::to_path_atomic` must
not reuse `write_with` as it stands.

**Re-exports.** `AtomicFile` joins the facade's existing
`pub use brust_core::{Compression, Diagnostic, Error, Format, Result};` line. The format crates
don't re-export it, which matches how `Compression` is handled today.

### Strict BAM reading

```rust
impl<R: Read> BgzfReader<R> { pub fn set_require_eof_block(&mut self, require: bool); }
impl<R: Read> BamReader<R>  { pub fn set_require_eof_block(&mut self, require: bool); }
```

- The default is `false`, which keeps today's behaviour. The setter style matches
  `FastaWriter::set_line_width`. `BamReader` forwards to its `BgzfReader`. Callers set it
  before reading to the end of the stream.
- **The rule.** When it is `true` and the underlying stream ends cleanly (no bytes left at a
  block boundary), the last BGZF block read must be byte-for-byte `EOF_BLOCK`, the SAM spec's
  28-byte marker and the same check htslib makes. A stream with no blocks at all fails the rule.
- **On failure,** the read that would have reported the end of the stream (`read` returning
  `Ok(0)`, or `read_record` returning `Ok(None)`) returns an `io::ErrorKind::InvalidData` error.
  - It carries a BAM diagnostic, built with the crate's `invalid_data`, so the facade sees
    `error.format() == Some(Format::Bam)`.
  - The diagnostic message is "BGZF stream ended without the EOF block; the file may be
    truncated", and `to_string()` is "invalid BAM: " followed by that message.
  - The error comes from `BgzfReader`, and `BamReader` passes it through unchanged.
- **The error is reported once.** Later reads behave as at end of stream, as they already do
  after other truncation errors, so `records()` ends and doesn't repeat the error forever.
- Empty blocks in the middle of the stream are still skipped, so joined BGZF streams pass. Only
  the final block counts. As a result, joined streams cut exactly after an interior EOF marker
  can't be detected. This is documented.
- A cut in the middle of a block already fails in both modes, and that doesn't change.
- The check runs at the end of the stream, so it also covers stdin and pipes. Callers only find
  out about truncation once they have read every record, as with truncated gzip today.
- `Bam::from_path` and `Bam::from_reader` stay lenient. A strict in-memory read is
  `BamReader::from_path`, then `set_require_eof_block(true)`, then `read_all()`.

### Facade and CLI

- `brust::convert`: every conversion writes through `from_path_atomic` and `commit`, or through
  `AtomicFile` directly. `bam_to_sam` writes the header text straight to its output, so it
  creates an `AtomicFile`, writes the header to it, passes it to `SamWriter::from_writer`, and
  ends with `writer.into_inner().commit()`. The `ConvertOptions` threads setting reaches
  `BamWriter::from_path_atomic_with_threads`. `write_atomic`, `create_temp_output_path` and
  `TEMP_COUNTER` are removed.
- `bam_to_sam`, `bam_to_fastq`, `validate::validate_bam` and the BAM path in `stats` turn on
  `set_require_eof_block(true)`.
- These are public functions of the `brust` crate as well as the CLI. They now reject BAMs they
  used to accept, which is why the release is 0.3.0 and not 0.2.1.

### Version 0.3.0

- `[workspace.package] version` becomes `"0.3.0"`.
- Every internal path dependency's `version = "..."` becomes `"0.3.0"`.
- The docs follow: the root README's "still early at version `0.2.0`" and `brust/README.md`'s
  `brust-fasta = "0.2.0"` example.
- `Cargo.lock` is updated by the build.

## Testing

Written test-first.

1. **`AtomicFile`** (brust-core unit tests, in a fresh temp folder per test):
   - After `create` and writes, the target doesn't exist when it didn't exist before, and is
     unchanged when it did.
   - `commit` makes the target hold exactly the bytes written and leaves no temp file in the
     folder.
   - Dropping without commit deletes the temp file, whether or not the target existed. A panic
     inside `catch_unwind` after some writes does the same.
   - The temp file's name starts with `.`, contains `.tmp.`, and ends with the target's file
     name.
   - A path with no file name gives `InvalidInput`. A missing parent folder gives `NotFound`.
     Neither creates anything.
   - A failed rename (target is a non-empty folder) returns an error, deletes the temp file and
     leaves the folder and its contents unchanged.
   - Two `AtomicFile`s for the same target have different temp paths. Committing both leaves
     the target holding exactly one of the two payloads, whichever committed last.
   - On Unix, a folder that can't be opened for reading makes `commit` return the "into place"
     error, with the new file in place and no temp file left. The test is skipped when running
     as root.
   - A compile-time check that `AtomicFile` is `Send + Sync`.
2. **Writer helpers** (each format crate):
   - For each of the five writers, `from_path_atomic` + `commit` gives byte-for-byte the same
     file as `from_path` for the same calls.
   - For each writer, dropping it after some writes without `commit` leaves no target and no
     temp file.
   - A `FastqWriter` with a `.fq.gz` target writes gzip (`1f 8b`) that reads back to the same
     records. A `.fq` target writes plain text.
   - `BamWriter::from_path_atomic_with_threads(path, 4)` gives the same bytes as
     `from_path_atomic`.
   - For each in-memory type, `to_path_atomic` gives the same bytes as `to_path`, and leaves no
     temp file.
3. **Strict BAM** (`bam` crate):
   - A valid fixture BAM reads fully in strict mode.
   - The same BAM without its last 28 bytes: lenient mode returns every record and then
     `Ok(None)`, pinning today's behaviour. Strict mode returns every record and then
     `InvalidData` with the message above.
   - A BAM cut in the middle of a block fails in both modes.
   - In strict mode, `records()` gives every record, then one error, then ends. A read after
     that returns `Ok(None)`.
   - `BgzfReader`, strict:
     - two complete BGZF streams joined together read fully;
     - data, then the marker, then more data with no final marker, fails;
     - an altered empty block in the middle, followed by a proper final marker, passes;
     - an empty input fails;
     - a stream whose last block is empty but differs from `EOF_BLOCK` in one header byte
       fails, which pins the byte-for-byte rule.
4. **Regression pins:**
   - A gzip FASTQ with its last byte removed makes `FastqReader` return an error.
   - A POD5 fixture with its last byte removed makes `Pod5Reader::from_path` return an error.
5. **Facade and CLI** (`brust/tests`):
   - A BAM without its EOF block makes `convert` (bam-to-sam and bam-to-fastq), `validate` and
     `stats` exit non-zero with a message containing "EOF block". `convert` leaves no output and
     no temp file.
   - The existing tests pass unchanged, including the no-temp-file checks at
     `brust/tests/convert.rs:131` and `:421`.
6. **Not tested:** that data survives a real power cut. Tests confirm that `commit` runs and
   succeeds. Review checks the order of fsync, close, rename and folder fsync against this spec.

### Verification before completion

`cargo fmt --check`, `cargo clippy --workspace --all-targets` (no warnings),
`cargo test --workspace`, `cargo doc --workspace --no-deps` (no warnings), and
`cargo metadata --no-deps --format-version 1` showing every workspace package at 0.3.0.

## Documentation

- `brust-core/README.md`: an `AtomicFile` section with a short example and the known limits.
- Each format crate's README: `from_path_atomic` and `commit`, and `to_path_atomic`.
  `bam/README.md` also covers `set_require_eof_block`.
- `fasta`, `fastq` and `sam` crate docs: a note that a plain file cut at a record boundary
  can't be detected when read, with a pointer to `from_path_atomic` for files you write.
- Root `README.md` (the conversion paragraph at lines 149–150) and `brust/README.md` (line 84):
  conversions now fsync the file before the rename, and on Unix the folder after it. A
  folder-sync error after the rename is reported, with the new file already in place. BAM input
  without the EOF block is rejected.
- Rustdoc for every new public item.

## Follow-up (not in this work)

- Publishing 0.3.0 to crates.io, crate by crate in dependency order, is the user's to run.
- Possible later: a seek-at-open EOF check for file-backed readers, a fsync opt-out, keeping
  permissions, and Ctrl-C cleanup in the CLI.
