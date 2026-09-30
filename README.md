# brust

`brust` is a Rust bioinformatics toolkit for reading, writing, validating,
summarizing, and converting common sequencing formats.

The workspace publishes a top-level facade crate named `brust` plus
format-specific crates:

| Package | Rust import path | Purpose |
| --- | --- | --- |
| `brust` | `brust` | Facade library and CLI binary |
| `brust-core` | `brust_core` | Shared error, diagnostic, format, and result types |
| `brust-fasta` | `brust_fasta` | FASTA reader, writer, and materialized records |
| `brust-fastq` | `brust_fastq` | FASTQ reader, writer, and FASTA conversion helpers |
| `brust-sam` | `brust_sam` | SAM reader, writer, flags, CIGAR, and optional fields |
| `brust-bam` | `brust_bam` | BAM/BGZF reader, writer, virtual offsets, and SAM conversion |
| `brust-pod5` | `brust_pod5` | POD5 metadata, reads, signal rows, and VBZ helpers |
| `brust-seq` | `brust_seq` | IUPAC reverse complement and matching, codon translation, and Phred means |

Most applications should depend on `brust`, which re-exports the format crates
under one namespace:

```rust
use brust::{bam, fasta, fastq, pod5, sam, seq};
```

## Installation

Use the library from crates.io:

```bash
cargo add brust
```

Install the CLI binary from crates.io:

```bash
cargo install brust
```

Build and install from the GitHub source checkout:

```bash
git clone https://github.com/Amekn/brust.git
cd brust
cargo test --workspace
cargo install --path brust
```

Use a format crate directly when you want the smallest possible dependency:

```bash
cargo add brust-fasta
cargo add brust-fastq
cargo add brust-sam
cargo add brust-bam
cargo add brust-pod5
```

## CLI Usage

The `brust` binary exposes validation, statistics, and conversion commands.

```bash
brust --help
```

```text
Usage: brust <COMMAND>

Commands:
  stats
  validate
  convert
```

### Validate

Validation is parser-backed. Each command streams through the file and reports
format-specific diagnostics from the same parser used by the library.

```bash
brust validate fasta reads.fasta
brust validate fastq reads.fastq
brust validate sam aligned.sam
brust validate bam aligned.bam
brust validate pod5 reads.pod5
```

On malformed input the command exits non-zero and prints structured context when
the parser can identify it, such as line and field details.

### Statistics

Statistics are streamed where possible and returned by the library as typed
structures.

```bash
brust stats fasta reference.fasta
brust stats fastq reads.fastq
brust stats sam aligned.sam
brust stats bam aligned.bam
brust stats pod5 reads.pod5
```

Examples of reported values include:

- FASTA and FASTQ record counts, sequence length distributions, N50/N90, base
  composition, duplicate ID counts, and GC fraction.
- FASTQ quality summaries including Phred min/max/mean, Q20/Q30 fractions, and
  per-read qscore (the Phred value of each read's mean base error probability,
  reported next to the per-read arithmetic mean Phred).
- SAM and BAM header/reference summaries, alignment flags, MAPQ, template
  lengths, CIGAR operation totals, optional tag counts, and records by
  reference.
- BAM record block-size summaries and unavailable quality counts.
- POD5 read counts, signal-row counts, run/channel summaries, sample totals,
  duration estimates, pore/end-reason counts, and scaling/calibration summaries.

### Convert

Supported conversions:

```bash
brust convert fastq-to-fasta reads.fastq reads.fasta
brust convert fastq-to-sam reads.fastq reads.sam
brust convert fastq-to-bam reads.fastq reads.bam
brust convert sam-to-bam aligned.sam aligned.bam
brust convert bam-to-sam aligned.bam aligned.sam
brust convert sam-to-fastq aligned.sam reads.fastq
brust convert bam-to-fastq aligned.bam reads.fastq
```

`sam-to-fastq` and `bam-to-fastq` write reverse-strand reads in their original
orientation, skip secondary and supplementary records, and leave read names
unchanged. Unmapped, QC-fail and duplicate records are still written. Output
differs from 0.1.x for inputs with reverse-strand, secondary or supplementary
records.

Add `-t`/`--threads` to compress BAM output on several threads:

```bash
brust convert sam-to-bam aligned.sam aligned.bam --threads 4
```

Only the `*-to-bam` conversions take `--threads`. The default is 1, and the
output is identical at any thread count.

Conversions stream records into a temporary file beside the target path. On success the file
is synced to disk and renamed over the target, so readers and crashes never see a partial
output: the target holds either the old file or the complete new one. On Unix the folder is
also synced before the conversion reports success, so the rename survives a power cut. Other
platforms skip that step, so after a power cut the target may still hold the old file, or be
missing if it was new, but never a partial one. If a conversion fails before the rename, the
temporary file is deleted and an existing output is left untouched. A folder-sync error after
the rename is reported, but the new file is already in place. BAM input must end with the BGZF
EOF block; a BAM without it is rejected as possibly truncated.

## Library API

The `brust` crate is the recommended public API for multi-format applications.

### Validate and Convert

```rust
use brust::{Format, convert, validate};
use brust::convert::Conversion;

fn main() -> brust::Result<()> {
    validate::validate(Format::Fastq, "reads.fastq")?;

    convert::convert(
        Conversion::FastqToFasta,
        "reads.fastq",
        "reads.fasta",
    )?;

    Ok(())
}
```

### Convert with Options

`ConvertOptions::threads` sets the number of BGZF compression threads. It affects
only conversions that write BAM (`FastqToBam` and `SamToBam`); other conversions
ignore it. The output is byte-identical for any thread count.

```rust
use brust::{Conversion, ConvertOptions, convert};

fn main() -> brust::Result<()> {
    let options = ConvertOptions::default().threads(4);

    convert::convert_with(
        Conversion::SamToBam,
        "aligned.sam",
        "aligned.bam",
        &options,
    )?;

    Ok(())
}
```

### Typed Statistics

```rust
use brust::{Format, Stats, stats};

fn main() -> brust::Result<()> {
    let summary = stats::stats(Format::Fastq, "reads.fastq")?;

    match summary {
        Stats::Fastq(stats) => {
            println!("reads={}", stats.reads);
            println!("bases={}", stats.bases.total);
            println!("q30={:?}", stats.qualities.q30_fraction);
        }
        _ => unreachable!("requested FASTQ stats"),
    }

    Ok(())
}
```

### FASTA

```rust
use brust::fasta::{FastaReader, FastaRecord, FastaWriter};

fn main() -> brust::Result<()> {
    let mut reader = FastaReader::from_path("seqs.fasta")?;
    while let Some(record) = reader.read_record()? {
        println!("{} {}", record.id, record.sequence.len());
    }

    let mut writer = FastaWriter::from_path("filtered.fasta")?;
    writer.write_record(&FastaRecord::new(
        "seq1".to_string(),
        Some("example".to_string()),
        "ACGTACGT".to_string(),
    ))?;
    writer.flush()?;

    Ok(())
}
```

### FASTQ

```rust
use brust::fastq::{Fastq, FastqReader};

fn main() -> brust::Result<()> {
    let mut reader = FastqReader::from_path("reads.fastq")?;
    while let Some(record) = reader.read_record()? {
        println!("{} len={}", record.id, record.sequence.len());
    }

    let fastq = Fastq::from_path("reads.fastq")?;
    fastq.to_fasta().to_path("reads.fasta")?;

    Ok(())
}
```

### SAM

```rust
use brust::sam::{SamOptionalValue, SamReader};

fn main() -> brust::Result<()> {
    let mut reader = SamReader::from_path("aligned.sam")?;
    println!("references={:?}", reader.header.sequence_names());

    while let Some(record) = reader.read_record()? {
        let query_len = record.query_len_from_cigar()?;
        let reference_len = record.reference_len_from_cigar()?;

        if let Some(SamOptionalValue::Integer(nm)) = record.aux("NM") {
            println!("{} q={query_len} r={reference_len} NM={nm}", record.qname);
        }
    }

    Ok(())
}
```

### BAM

```rust
use brust::bam::{Bam, BamReader};
use brust::sam::Sam;

fn main() -> brust::Result<()> {
    let mut reader = BamReader::from_path("aligned.bam")?;

    if let Some(positioned) = reader.read_record_with_virtual_offset()? {
        println!("offset={}", positioned.virtual_offset.raw());
        println!("{}", positioned.record.to_sam_line(&reader.refs)?);
    }

    let sam = Sam::from_path("aligned.sam")?;
    let bam = Bam::from_sam(&sam)?;
    bam.to_path("aligned.bam")?;

    Ok(())
}
```

### POD5

```rust
use brust::pod5::{Pod5, Pod5Reader};

fn main() -> brust::Result<()> {
    let pod5 = Pod5::from_path("reads.pod5")?;
    let summary = pod5.summary();
    println!("reads={} samples={}", summary.read_count, summary.total_samples);

    if let Some(record) = pod5.read_by_id("read-id") {
        let cache = pod5.signal_cache();
        let signal = cache.signal_for_record(record)?;
        println!("{} samples", signal.len());
    }

    let mut reader = Pod5Reader::from_path("reads.pod5")?;
    while let Some(record) = reader.read_record()? {
        println!("{} {}", record.read_id, record.num_samples);
    }

    Ok(())
}
```

## Error Handling

The facade exposes `brust::Error`, `brust::Diagnostic`, `brust::Format`, and
`brust::Result`.

```rust
use brust::{Format, validate};

fn main() {
    let error = validate::validate(Format::Fastq, "bad.fastq").unwrap_err();

    if let Some(diagnostic) = error.diagnostic() {
        eprintln!("message={}", diagnostic.message);
        eprintln!("line={:?}", diagnostic.line);
        eprintln!("field={:?}", diagnostic.field);
    }
}
```

The lower-level format crates expose `std::io::Result` APIs. Parser failures are
`InvalidData` I/O errors whose inner error can contain the same Brust domain
diagnostic.

## Reliability Notes

Brust is structured around a few product-grade reliability choices:

- Parser-backed validation avoids duplicated validation logic between CLI and
  library entry points.
- Streaming readers and writers avoid unnecessary whole-file materialization for
  validation, conversion, and many analysis workflows.
- Atomic conversion output protects existing output files from partial writes.
- Typed statistics avoid forcing CLI-oriented text parsing on application code.
- The workspace has unit and integration tests across parsing, writing,
  validation, conversion, stats, and CLI behavior.

The project is still early at version `0.3.0`, so users should validate behavior
against their production data and report edge cases. The intended direction is a
robust, reliable bioinformatics toolkit that can serve both command-line and
Rust application workflows.

## Development

Common local checks:

```bash
cargo fmt --check
cargo test --workspace
cargo doc --workspace --no-deps
cargo package --list -p brust
```

The workspace uses Rust edition 2024 and shared package metadata from the root
`Cargo.toml`.

## License

Licensed under either of:

- [Apache License, Version 2.0](LICENSE-APACHE)
- [MIT License](LICENSE-MIT)
