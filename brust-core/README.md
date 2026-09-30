# brust-core

`brust-core` contains shared diagnostics and result types used by the Brust
format crates and the `brust` facade.

Most applications should depend on `brust` instead. Use `brust-core` directly
when implementing a Brust-compatible format crate or when you need to inspect
domain diagnostics carried inside `std::io::Error` values from a lower-level
format crate.

## Installation

```bash
cargo add brust-core
```

## API

The crate provides:

- `Format`: the supported format enum: FASTA, FASTQ, SAM, BAM, and POD5.
- `Compression`: shared uncompressed/gzip selection for path-based format I/O.
- `Diagnostic`: message plus optional line and field context.
- `Error`: cloneable Brust error type with format-specific invalid-data
  variants and an I/O variant.
- `Result<T>`: alias for `std::result::Result<T, Error>`.

`Compression::from_path` recognizes a final `.gz` suffix case-insensitively,
including the conventional `.fq.gz` and `.fastq.gz` FASTQ names.

## Example

```rust
use brust_core::{Error, Format};

fn main() {
    let error = Error::invalid(Format::Fastq, "quality length mismatch")
        .with_line(4)
        .with_field("QUAL");

    assert_eq!(error.format(), Some(Format::Fastq));
    assert_eq!(error.diagnostic().unwrap().line, Some(4));
    println!("{error}");
}
```

## I/O Interoperability

The format crates expose `std::io::Result` APIs. Parser failures are represented
as `InvalidData` I/O errors whose inner error can be recovered as
`brust_core::Error`:

```rust
use brust_core::Error;

fn inspect(error: std::io::Error) {
    if let Some(domain) = error.get_ref().and_then(|inner| inner.downcast_ref::<Error>()) {
        eprintln!("format={:?} diagnostic={:?}", domain.format(), domain.diagnostic());
    }
}
```

## Atomic File Output

`AtomicFile` is a buffered file sink that publishes its contents at the target
path only when `commit` succeeds. It writes to a hidden temporary file beside
the target, then flushes, syncs, renames and (on Unix) syncs the folder.
Dropping it without `commit` removes the temporary file and leaves any existing
target unchanged.

```rust
use brust_core::AtomicFile;
use std::io::Write;

fn main() -> std::io::Result<()> {
    let mut file = AtomicFile::create("out.txt")?;
    file.write_all(b"hello\n")?;
    file.commit()
}
```

`commit` can return an error after the rename, if syncing the folder fails. The
new file is already in place in that case.

Known limits:

- The rename replaces a symlink at the target with a plain file.
- The new file gets default permissions (umask), not the old file's.
- Hard links to the old file keep the old contents.
- A read-only existing target is replaced, as `mv` does.
- Paths are kept as given, so do not change the current folder between `create`
  and `commit`.
- The rename is atomic only within one filesystem. The temporary file is
  created beside the target to guarantee that.
- Only Unix syncs the folder after the rename. Elsewhere a power cut soon after
  `commit` can undo the rename, leaving the old file (or no file if the target
  was new), but never a partial one.
