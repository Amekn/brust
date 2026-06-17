# brust-fasta

`brust-fasta` provides FASTA reader and writer primitives for the Brust
workspace. It supports streaming records from any `Read`, writing to any
`Write`, and materializing a full FASTA payload when that is more convenient.

Most multi-format applications should depend on `brust` and use
`brust::fasta`. Use `brust-fasta` directly for a small FASTA-only dependency.

## Installation

```bash
cargo add brust-fasta
```

```rust
use brust_fasta::FastaReader;
```

## API

- `FastaReader`: streaming parser over files or arbitrary readers.
- `FastaWriter`: streaming writer with optional sequence line wrapping.
- `Fasta`: materialized FASTA collection with `from_path` and `to_path`.
- `FastaRecord`: owned record with `id`, optional `description`, and
  `sequence`.

## Streaming Read

```rust
use brust_fasta::FastaReader;

fn main() -> std::io::Result<()> {
    let mut reader = FastaReader::from_path("seqs.fasta")?;

    while let Some(record) = reader.read_record()? {
        println!("{} {}", record.id, record.sequence.len());
    }

    Ok(())
}
```

## Write Records

```rust
use brust_fasta::{FastaRecord, FastaWriter};

fn main() -> std::io::Result<()> {
    let mut writer = FastaWriter::from_path("seqs.fasta")?;
    writer.set_line_width(80);

    writer.write_record(&FastaRecord::new(
        "seq1".to_string(),
        Some("example sequence".to_string()),
        "ACGTACGT".to_string(),
    ))?;

    writer.flush()
}
```

## Materialized Workflow

```rust
use brust_fasta::Fasta;

fn main() -> std::io::Result<()> {
    let fasta = Fasta::from_path("input.fasta")?;
    println!("records={}", fasta.records.len());
    fasta.to_path("copy.fasta")
}
```

Parser failures are returned as `std::io::ErrorKind::InvalidData` with a
`brust_core::Error` inner value when structured context is available.
