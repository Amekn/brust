# brust-fastq

`brust-fastq` provides FASTQ reader and writer primitives for the Brust
workspace, plus helpers for converting FASTQ records to FASTA records. Plain
`.fq`/`.fastq` and gzip-compressed `.fq.gz`/`.fastq.gz` files are supported.

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

- `FastqReader`: streaming parser over plain or gzip files and arbitrary
  readers; gzip streams are recognized by their magic bytes.
- `FastqWriter`: streaming plain/gzip writer; path APIs select gzip for `.gz`
  and arbitrary writers can use an explicit `Compression` value.
- `Fastq`: materialized FASTQ collection with `from_path`, `to_path`, and
  `to_fasta`.
- `FastqRecord`: owned read record with `id`, optional `description`,
  `sequence`, and `quality`.

## Streaming Read

```rust
use brust_fastq::FastqReader;

fn main() -> std::io::Result<()> {
    let mut reader = FastqReader::from_path("reads.fastq.gz")?;

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
    let mut writer = FastqWriter::from_path("reads.fastq.gz")?;

    writer.write_record(&FastqRecord::new(
        "read1".to_string(),
        None,
        "ACGT".to_string(),
        "IIII".to_string(),
    ))?;

    writer.finish().map(drop)
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

Compression and decompression stay streaming: only parser and codec buffers plus
the current record are held in memory. `Fastq::from_path` remains the explicitly
materialized convenience API.
