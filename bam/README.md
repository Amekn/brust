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
- `BamWriter`: writes BAM headers and records. `BamWriter::from_path` and
  `BamWriter::from_writer` compress on the calling thread;
  `BamWriter::from_path_with_threads` and `BamWriter::from_writer_with_threads`
  compress on worker threads.
- `Bam`: materialized BAM payload with `from_path`, `to_path`, and
  `from_sam`.
- `BamRecord`: decoded fixed, variable, and auxiliary record fields.
  `BamRecord::encode` appends the record in BAM binary form.
  `BamRecord::original_sequence_string` and `BamRecord::original_quality_string`
  return the sequence and qualities in the original sequencing orientation:
  reverse-complemented and reversed for reverse-strand records (flag `0x10`),
  unchanged otherwise.
- `BamAuxValue` and `BamAuxArray`: parsed BAM auxiliary tags.
- `BgzfVirtualOffset`: compressed/uncompressed BGZF virtual offset.
- `BgzfWriter`: writes any byte stream as BGZF blocks, on the calling thread or
  on worker threads (`BgzfWriter::with_threads`).
- `bgzf::compress_block`: frames one block of data as a BGZF block.
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

## Parallel Compression

Build the writer with `BamWriter::from_path_with_threads` to compress BGZF blocks
on several threads:

```rust
use brust_bam::{BamReader, BamWriter};

fn main() -> std::io::Result<()> {
    let mut reader = BamReader::from_path("aligned.bam")?;
    let mut writer = BamWriter::from_path_with_threads("copy.bam", 4)?;

    writer.write_header(&reader.header, &reader.refs)?;
    while let Some(record) = reader.read_record()? {
        writer.write_record(&record)?;
    }
    writer.finish()?;

    Ok(())
}
```

A `threads` of 0 or 1 compresses on the calling thread, and larger values are
capped at `brust_bam::bgzf::MAX_THREADS`. The file is byte-for-byte the same
for every thread count, and the same as `BamWriter::from_path` writes. Use
`BamWriter::from_writer_with_threads` to write to any other `Write`, and call
`finish` to write the last block and the BGZF EOF block.

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
