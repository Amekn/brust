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
- `Diagnostic`: message plus optional line and field context.
- `Error`: cloneable Brust error type with format-specific invalid-data
  variants and an I/O variant.
- `Result<T>`: alias for `std::result::Result<T, Error>`.

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
