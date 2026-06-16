# Brust

Brust is a Rust bioinformatics toolkit for reading, writing, validating, and
converting common sequencing formats. The public crate is `brust`, which
re-exports the format crates under one namespace:

```rust
use brust::{bam, fasta, fastq, pod5, sam};
```

## Format Examples

Read FASTA records:

```rust
use brust::fasta::FastaReader;

fn main() -> std::io::Result<()> {
    let mut reader = FastaReader::from_path("reads.fasta")?;
    while let Some(record) = reader.read_record()? {
        println!("{} {}", record.id, record.sequence.len());
    }
    Ok(())
}
```

Convert FASTQ to FASTA:

```rust
use brust::fastq::Fastq;

fn main() -> std::io::Result<()> {
    let fastq = Fastq::from_path("reads.fastq")?;
    let fasta = fastq.to_fasta();
    fasta.to_path("reads.fasta")
}
```

Read SAM alignments:

```rust
use brust::sam::SamReader;

fn main() -> std::io::Result<()> {
    let mut reader = SamReader::from_path("aligned.sam")?;
    println!("references: {:?}", reader.header.sequence_names());
    while let Some(record) = reader.read_record()? {
        println!("{} mapped={}", record.qname, !record.is_unmapped());
    }
    Ok(())
}
```

Read BAM records and format one as SAM:

```rust
use brust::bam::BamReader;

fn main() -> std::io::Result<()> {
    let mut reader = BamReader::from_path("aligned.bam")?;
    if let Some(positioned) = reader.read_record_with_virtual_offset()? {
        println!("virtual offset: {}", positioned.virtual_offset.raw());
        println!("{}", positioned.record.to_sam_line(&reader.refs)?);
    }
    Ok(())
}
```

Convert supported SAM payloads to BAM:

```rust
use brust::{bam::Bam, sam::Sam};

fn main() -> std::io::Result<()> {
    let sam = Sam::from_path("aligned.sam")?;
    let bam = Bam::from_sam(&sam)?;
    bam.to_path("aligned.bam")
}
```

Inspect POD5 reads and signal:

```rust
use brust::pod5::Pod5;

fn main() -> std::io::Result<()> {
    let pod5 = Pod5::from_path("reads.pod5")?;
    let summary = pod5.summary();
    println!("reads={} samples={}", summary.read_count, summary.total_samples);

    if let Some(record) = pod5.read_by_id("1cadb1e9-592f-4e22-9285-4626f2b7da9f") {
        let cache = pod5.signal_cache();
        let signal = cache.signal_for_record(record)?;
        println!("{} samples", signal.len());
    }
    Ok(())
}
```

Malformed format data is reported as `InvalidData` I/O errors whose inner error
is `brust::Error`, with format, line, and field context where the parser can
identify it.

## Crates

- `brust`: public facade and CLI crate.
- `brust-core`: shared error and result types.
- `fasta`, `fastq`, `sam`, `bam`, `pod5`: format-specific readers, writers,
  validators, and helpers.

## License

Source-available for viewing only. See [LICENSE](LICENSE). The project may be
released under a different license in the future.
