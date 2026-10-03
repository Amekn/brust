# brust-sam

`brust-sam` provides SAM header and alignment parsing, writing, flag helpers,
CIGAR helpers, and optional-field handling for the Brust workspace.

Most multi-format applications should depend on `brust` and use `brust::sam`.
Use `brust-sam` directly for SAM-only parsing or writing.

## Installation

```bash
cargo add brust-sam
```

```rust
use brust_sam::SamReader;
```

## API

- `SamReader`: validates the header, then streams alignment records.
- `SamWriter`: writes headers and records, validating output before emission.
  `SamWriter::from_path` buffers its output; call `finish` to write the rest
  and see any error.
- `Sam`: materialized SAM payload with header and records.
- `SamHeader`, `SamHeaderRecord`, `SamHeaderField`: parsed header data.
- `SamRecord`: mandatory fields plus parsed optional fields.
  `SamRecord::validate` checks a record against the rules the reader applies,
  so it can be written and read back unchanged.
  `SamRecord::original_seq` and `SamRecord::original_qual` return SEQ and QUAL
  in the original sequencing orientation: reverse-complemented and reversed for
  reverse-strand records (flag `0x10`), unchanged otherwise.
- `SamOptionalValue` and `SamOptionalArray`: optional `TAG:TYPE:VALUE` data.
- `flags`: standard SAM bit flags.

## Streaming Alignments

```rust
use brust_sam::SamReader;

fn main() -> std::io::Result<()> {
    let mut reader = SamReader::from_path("aligned.sam")?;

    println!("references={:?}", reader.header.sequence_names());

    while let Some(record) = reader.read_record()? {
        println!(
            "{} mapped={} reverse={}",
            record.qname,
            !record.is_unmapped(),
            record.is_reverse_complemented(),
        );
    }

    Ok(())
}
```

## Inspect CIGAR and Optional Fields

```rust
use brust_sam::{SamOptionalValue, SamReader};

fn main() -> std::io::Result<()> {
    let mut reader = SamReader::from_path("aligned.sam")?;

    while let Some(record) = reader.read_record()? {
        let query_len = record.query_len_from_cigar()?;
        let reference_len = record.reference_len_from_cigar()?;

        if let Some(SamOptionalValue::Integer(nm)) = record.aux("NM") {
            println!("{} query={query_len} reference={reference_len} NM={nm}", record.qname);
        }
    }

    Ok(())
}
```

## Atomic Write

```rust
use brust_sam::SamWriter;

fn main() -> std::io::Result<()> {
    let mut writer = SamWriter::from_path_atomic("out.sam")?;
    // ... write to `writer` ...
    writer.commit()
}
```

`from_path_atomic` writes to a hidden temporary file beside the target and
renames it into place on `commit`, so a crash or early return never leaves a
half-written SAM file. Dropping the writer without committing discards the
output. `Sam::to_path_atomic` does the same for an in-memory `Sam`.

## Materialized Round Trip

```rust
use brust_sam::Sam;

fn main() -> std::io::Result<()> {
    let sam = Sam::from_path("aligned.sam")?;
    sam.to_path("copy.sam")
}
```

Malformed headers, alignment records, CIGAR strings, and optional fields are
reported as `InvalidData` I/O errors with structured Brust diagnostics when
available.
