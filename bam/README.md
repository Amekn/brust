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

## Strict End-of-File Check

A complete BAM ends with the 28-byte BGZF EOF block. Readers don't require it by
default. To treat a missing marker as truncation, opt in before reading:

```rust
let mut reader = brust_bam::BamReader::from_path("aligned.bam")?;
reader.set_require_eof_block(true);
let bam = reader.read_all()?; // Err(InvalidData) if the EOF block is missing
```

The last block must match the marker byte for byte. The error is reported once,
then the reader behaves as at end of stream. `BgzfReader` has the same setter.
Joined streams cut exactly after an interior EOF marker can't be detected.

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

To publish the file atomically, use `BamWriter::from_path_atomic_with_threads`
(or `from_path_atomic` for a single thread) and call `commit` instead of
`finish`:

```rust
use brust_bam::BamWriter;

fn main() -> std::io::Result<()> {
    let mut writer = BamWriter::from_path_atomic_with_threads("copy.bam", 4)?;
    // ... write the header and records to `writer` ...
    writer.commit()
}
```

The output is written to a hidden temporary file beside the target and renamed
into place on `commit`, so a crash or early return never leaves a half-written
BAM file. Dropping the writer without committing discards the output.
`Bam::to_path_atomic` does the same for an in-memory `Bam`.

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
