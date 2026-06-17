# brust

`brust` is the public facade and command-line crate for the Brust
bioinformatics workspace. It re-exports the format crates under one namespace
and provides parser-backed validation, structured statistics, and safe format
conversion helpers.

The crate is designed for applications that want one dependency for FASTA,
FASTQ, SAM, BAM, and POD5 handling:

```rust
use brust::{bam, fasta, fastq, pod5, sam};
```

## Installation

Use the library from crates.io:

```bash
cargo add brust
```

Install the CLI binary:

```bash
cargo install brust
```

Use the current repository checkout:

```bash
git clone https://github.com/Amekn/brust.git
cd brust
cargo test --workspace
cargo install --path brust
```

## CLI

Validate files with the same parsers used by the library:

```bash
brust validate fasta reads.fasta
brust validate fastq reads.fastq
brust validate sam aligned.sam
brust validate bam aligned.bam
brust validate pod5 reads.pod5
```

Print human-readable statistics:

```bash
brust stats fastq reads.fastq
brust stats bam aligned.bam
brust stats pod5 reads.pod5
```

Convert between supported formats:

```bash
brust convert fastq-to-fasta reads.fastq reads.fasta
brust convert fastq-to-sam reads.fastq reads.sam
brust convert fastq-to-bam reads.fastq reads.bam
brust convert sam-to-bam aligned.sam aligned.bam
brust convert bam-to-sam aligned.bam aligned.sam
brust convert sam-to-fastq aligned.sam reads.fastq
brust convert bam-to-fastq aligned.bam reads.fastq
```

Conversions stream records and write through a temporary output path before
renaming, so an existing output file is not replaced by a partial file when
parsing or writing fails.

## Public API

Runtime-dispatched validation, statistics, and conversion live in the facade:

```rust
use brust::{Format, Stats};
use brust::{convert, stats, validate};
use brust::convert::Conversion;

fn main() -> brust::Result<()> {
    validate::validate(Format::Fastq, "reads.fastq")?;

    convert::convert(
        Conversion::FastqToFasta,
        "reads.fastq",
        "reads.fasta",
    )?;

    match stats::stats(Format::Fasta, "reads.fasta")? {
        Stats::Fasta(summary) => {
            println!("records={}", summary.records);
            println!("gc={:?}", summary.bases.gc_fraction());
        }
        _ => unreachable!("requested FASTA stats"),
    }

    Ok(())
}
```

Each format crate is re-exported for lower-level streaming and materialized
workflows:

```rust
use brust::fasta::FastaReader;
use brust::fastq::Fastq;
use brust::sam::SamReader;
use brust::bam::BamReader;
use brust::pod5::Pod5;

fn main() -> brust::Result<()> {
    let mut fasta = FastaReader::from_path("seqs.fasta")?;
    while let Some(record) = fasta.read_record()? {
        println!("{} {}", record.id, record.sequence.len());
    }

    let fastq = Fastq::from_path("reads.fastq")?;
    fastq.to_fasta().to_path("reads.fasta")?;

    let mut sam = SamReader::from_path("aligned.sam")?;
    while let Some(record) = sam.read_record()? {
        println!("{} mapped={}", record.qname, !record.is_unmapped());
    }

    let mut bam = BamReader::from_path("aligned.bam")?;
    if let Some(positioned) = bam.read_record_with_virtual_offset()? {
        println!("offset={}", positioned.virtual_offset.raw());
        println!("{}", positioned.record.to_sam_line(&bam.refs)?);
    }

    let pod5 = Pod5::from_path("reads.pod5")?;
    println!("samples={}", pod5.total_samples());

    Ok(())
}
```

Malformed format data is surfaced as `brust::Error` diagnostics when using the
facade. The lower-level format crates expose `std::io::Result` and store the
same domain error as the inner `InvalidData` error when possible.

## Workspace Crates

- `brust`: facade library and CLI.
- `brust-core`: shared `Error`, `Diagnostic`, `Format`, and `Result` types.
- `brust-fasta`: FASTA reader, writer, and materialized records.
- `brust-fastq`: FASTQ reader, writer, and FASTA conversion helpers.
- `brust-sam`: SAM reader, writer, flags, CIGAR, and optional fields.
- `brust-bam`: BAM/BGZF reader, writer, virtual offsets, and SAM conversion.
- `brust-pod5`: POD5 metadata, reads, signal rows, and VBZ helpers.

Direct subcrate dependencies use the published package names, while Rust import
paths use underscores:

```toml
[dependencies]
brust-fasta = "0.1.1"
```

```rust
use brust_fasta::FastaReader;
```

Prefer `brust` when building an application that needs multiple formats.
