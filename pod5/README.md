# brust-pod5

`brust-pod5` provides POD5 reader and writer primitives for the Brust workspace.
It focuses on metadata, read rows, signal-row references, and VBZ signal
decompression/compression helpers.

Most multi-format applications should depend on `brust` and use `brust::pod5`.
Use `brust-pod5` directly for POD5-focused workflows.

## Installation

```bash
cargo add brust-pod5
```

```rust
use brust_pod5::Pod5;
```

## API

- `Pod5Reader`: streams POD5 read rows and can load referenced signal rows.
- `Pod5Writer`: writes materialized POD5 payloads.
- `Pod5`: materialized container with reads, signals, run metadata, summaries,
  and signal lookup helpers.
- `Pod5Record`: one read row plus metadata and signal-row references.
- `Pod5Signal`: signal row with VBZ compression helpers.
- `Pod5SignalCache`: caches decompressed signal rows for repeated lookup.
- `Pod5Summary`, `Pod5ChannelSummary`, and `Pod5RunInfoSummary`: metadata
  summaries for analysis/reporting.
- `compress_vbz_signal` and `decompress_vbz_signal`: standalone VBZ helpers.

## Materialized Summary

```rust
use brust_pod5::Pod5;

fn main() -> std::io::Result<()> {
    let pod5 = Pod5::from_path("reads.pod5")?;
    let summary = pod5.summary();

    println!("reads={}", summary.read_count);
    println!("samples={}", summary.total_samples);

    for channel in summary.channels {
        println!("channel={} reads={}", channel.channel, channel.read_count);
    }

    Ok(())
}
```

## Signal Lookup

```rust
use brust_pod5::Pod5;

fn main() -> std::io::Result<()> {
    let pod5 = Pod5::from_path("reads.pod5")?;
    let cache = pod5.signal_cache();

    if let Some(record) = pod5.read_by_id("read-id") {
        let samples = cache.signal_for_record(record)?;
        println!("{} samples", samples.len());
    }

    Ok(())
}
```

## Streaming Reads

```rust
use brust_pod5::Pod5Reader;

fn main() -> std::io::Result<()> {
    let mut reader = Pod5Reader::from_path("reads.pod5")?;

    while let Some(record) = reader.read_record()? {
        println!("{} {} samples", record.read_id, record.num_samples);
    }

    Ok(())
}
```

Malformed POD5 wrapper metadata, section data, signal references, or signal
payloads are reported as `InvalidData` I/O errors with structured Brust
diagnostics when available.
