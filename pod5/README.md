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
  Each signal row is read from the Arrow batch that holds it, so reads can be
  asked for in any order without reading the batches before them.
- `Pod5Writer`: writes materialized POD5 payloads. The payload is checked
  before anything is written (sample counts, signal-row ownership, version and
  file identifier formats, end reasons). The Signal table is then streamed to
  the output in batches of 100 rows, as official POD5 writers do, so writing
  needs little memory beyond the payload. Files that declare POD5 version
  0.3.30 or later (the default is 0.3.34) include the `open_pore_level` column
  that official pod5 and dorado require.
- `Pod5`: materialized container with reads, signals, run metadata, summaries,
  and signal lookup helpers.
- `Pod5Record`: one read row plus metadata and signal-row references.
- `Pod5RunInfo`: one Run Info row with all 20 columns.
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

## Atomic Write

```rust
use brust_pod5::Pod5Writer;

fn main() -> std::io::Result<()> {
    let mut writer = Pod5Writer::from_path_atomic("out.pod5")?;
    // ... write to `writer` ...
    writer.commit()
}
```

`from_path_atomic` writes to a hidden temporary file beside the target and
renames it into place on `commit`, so a crash or early return never leaves a
half-written POD5 file. Dropping the writer without committing discards the
output. `Pod5::to_path_atomic` does the same for an in-memory `Pod5`.

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
