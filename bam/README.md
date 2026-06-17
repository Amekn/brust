# brust-bam

`brust-bam` provides BAM and BGZF reader/writer primitives for the Brust
workspace, including virtual offsets and conversion helpers between SAM records
and BAM records.

Most multi-format applications should depend on `brust` and use `brust::bam`.
Use `brust-bam` directly for BAM-only workflows.

## Installation

```bash
cargo add brust-bam
```

```rust
use brust_bam::BamReader;
```

SAM conversion APIs also use types from `brust-sam`, so add it when calling
`Bam::from_sam` directly:

```bash
cargo add brust-sam
```

## API

- `BamReader`: reads BAM headers, reference dictionaries, and records.
- `BamWriter`: writes BAM headers and records.
- `Bam`: materialized BAM payload with `from_path`, `to_path`, and
  `from_sam`.
- `BamRecord`: decoded fixed, variable, and auxiliary record fields.
- `BamAuxValue` and `BamAuxArray`: parsed BAM auxiliary tags.
- `BgzfVirtualOffset`: compressed/uncompressed BGZF virtual offset.
- `SamToBamConverter`: converts validated SAM records into BAM records.

## Streaming BAM Records

```rust
use brust_bam::BamReader;

fn main() -> std::io::Result<()> {
    let mut reader = BamReader::from_path("aligned.bam")?;

    while let Some(record) = reader.read_record()? {
        println!(
            "{} mapped={} cigar={}",
            record.read_name(),
            !record.is_unmapped(),
            record.cigar_string(),
        );
    }

    Ok(())
}
```

## Virtual Offsets

```rust
use brust_bam::BamReader;

fn main() -> std::io::Result<()> {
    let mut reader = BamReader::from_path("aligned.bam")?;

    if let Some(positioned) = reader.read_record_with_virtual_offset()? {
        println!("raw offset={}", positioned.virtual_offset.raw());
        println!("read={}", positioned.record.read_name());
    }

    Ok(())
}
```

## Convert SAM to BAM

```rust
use brust_bam::Bam;
use brust_sam::Sam;

fn main() -> std::io::Result<()> {
    let sam = Sam::from_path("aligned.sam")?;
    let bam = Bam::from_sam(&sam)?;
    bam.to_path("aligned.bam")
}
```

The `brust` facade re-exports both crates as `brust::bam` and `brust::sam` when
you prefer one dependency for all supported formats.

Malformed BAM headers, BGZF blocks, record layouts, or auxiliary fields are
reported as `InvalidData` I/O errors with structured Brust diagnostics when
available.
