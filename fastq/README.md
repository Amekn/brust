# brust-fastq

`brust-fastq` provides FASTQ reader and writer primitives for the Brust
workspace, plus helpers for converting FASTQ records to FASTA records.

Most multi-format applications should depend on `brust` and use
`brust::fastq`. Use `brust-fastq` directly for a small FASTQ-focused
dependency.

## Installation

```bash
cargo add brust-fastq
```

```rust
use brust_fastq::FastqReader;
```

## API

- `FastqReader`: streaming parser over files or arbitrary readers.
- `FastqWriter`: streaming writer over files or arbitrary writers.
- `Fastq`: materialized FASTQ collection with `from_path`, `to_path`, and
  `to_fasta`.
- `FastqRecord`: owned read record with `id`, optional `description`,
  `sequence`, and `quality`.

## Streaming Read

```rust
use brust_fastq::FastqReader;

fn main() -> std::io::Result<()> {
    let mut reader = FastqReader::from_path("reads.fastq")?;

    while let Some(record) = reader.read_record()? {
        println!("{} len={}", record.id, record.sequence.len());
    }

    Ok(())
}
```

## Write Records

```rust
use brust_fastq::{FastqRecord, FastqWriter};

fn main() -> std::io::Result<()> {
    let mut writer = FastqWriter::from_path("reads.fastq")?;

    writer.write_record(&FastqRecord::new(
        "read1".to_string(),
        None,
        "ACGT".to_string(),
        "IIII".to_string(),
    ))?;

    writer.flush()
}
```

## Convert to FASTA

```rust
use brust_fastq::Fastq;

fn main() -> std::io::Result<()> {
    let fastq = Fastq::from_path("reads.fastq")?;
    let fasta = fastq.to_fasta();
    fasta.to_path("reads.fasta")
}
```

The parser checks record structure and sequence/quality length agreement.
Malformed input is returned as `std::io::ErrorKind::InvalidData` with structured
Brust diagnostics when available.
