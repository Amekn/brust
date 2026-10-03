# Changelog

Notable changes to the brust workspace. All crates in the workspace share one
version number.

## 0.4.0 (03/10/2026)

This release fixes the bugs found in a full code review, and makes several
commands faster or use less memory. A few public types change. Files that
samtools, htslib or the official POD5 tools read are still read. Some
malformed files that brust used to accept or write are now rejected with an
error instead.

### Breaking changes

- **brust-bam:** `BamHeader.text` is now `Vec<u8>`, kept byte for byte. It used
  to be a `String` decoded lossily, so header text that wasn't UTF-8 changed
  when read and written back. Use `BamHeader::text_str()` to get the text as a
  `&str`. To build a header from a `String`, use `into_bytes()`.
- **brust-pod5:** `Pod5Record` has a new field, `open_pore_level: Option<f32>`.
  `Pod5RunInfo` has 13 new fields, so all 20 Run Info columns are kept:
  - `acquisition_start_time`, `protocol_start_time`
  - `adc_max`, `adc_min`
  - `context_tags`, `tracking_id`
  - `flow_cell_product_code`, `protocol_name`, `protocol_run_id`
  - `sequencer_position`, `sequencer_position_type`
  - `system_name`, `system_type`

  Code that builds these structs with struct literals must set the new fields.
- **brust-fasta, brust-fastq, brust-sam:** these constructors now return writers
  over `BufWriter<File>` instead of `File`:
  - `FastaWriter::from_path`
  - `FastqWriter::from_path` and `FastqWriter::from_path_with_compression`
  - `SamWriter::from_path`
  - their `new` aliases

  Write errors may now be reported later, when the buffer is written out,
  rather than on the record that caused them. Call `finish()` (new on
  `FastaWriter` and `SamWriter`) when you are done, to write the rest and see
  any error. For gzip FASTQ output, `finish()` also writes the gzip trailer,
  which `flush()` does not. Dropping a writer still writes the rest, but
  ignores errors.

### Added

- `SamRecord::validate()` checks a record against the rules the SAM reader
  applies, without formatting it as text.
- `FastaWriter::finish()` and `SamWriter::finish()` flush the output and return
  the inner writer.
- `BamHeader::text_str()`.
- `BamRecord::cigar_ops()` and `BamRecord::stores_cigar_in_cg_tag()`, for
  records that use the long-CIGAR `CG` tag (see below).
- POD5 `open_pore_level` and the full Run Info table are read and written.

### Fixed

#### POD5 (brust-pod5)

- Files written by brust said they were POD5 version 0.3.34 but lacked the
  `open_pore_level` column, so official pod5 and dorado couldn't open them. The
  column is now written for version 0.3.30 and later.
- Rewriting a POD5 file dropped 13 of the 20 Run Info columns, after which
  dorado failed with "No supported chemistry". All columns are now kept.
- `Pod5Writer` refuses a second payload, and refuses `commit` before a payload
  has been written. `Pod5::to_path` checks the payload before it replaces an
  existing file.
- Before writing anything, the writer checks:
  - each signal row's sample count
  - that a read's signal rows carry its read ID
  - the `pod5_version` and `file_identifier` formats
  - end reasons
- The writer accepts empty tables and reads with no samples. Too many distinct
  dictionary values is now an error instead of a panic.
- The reader no longer crashes, panics or runs out of memory on malformed files;
  it reports them as invalid:
  - Embedded Arrow blocks are checked before arrow-ipc reads them.
  - No memory is reserved from counts stored in the file.
  - arrow-ipc panics are turned into errors.
  - Arrow bodies compressed inside the file, which POD5 writers don't produce,
    are rejected.
- `Pod5::from_reader` uses the same parser as `Pod5::from_path`.
- Streaming signal reads no longer keep every decoded read in memory. Rows
  after an unreadable Signal batch keep their numbers.
- Sample and event totals can no longer overflow.

#### FASTQ (brust-fastq)

- Zero-length reads, as cutadapt writes them, are read.
- The record iterator stops after its first error instead of repeating it
  forever.
- `finish` on a gzip writer flushes the inner writer, so write errors are
  reported.
- The reader checks the text after `+` against the header, and only accepts
  quality characters from `!` to `~`.
- The writer refuses records that would read back differently: a sequence
  starting with `+`, or a sequence or quality ending in whitespace.

#### FASTA (brust-fasta)

- The writer refuses records that would read back differently: a written line
  starting with `>` or ending in whitespace. Line wrapping counts characters,
  not bytes.
- A carriage return is only accepted as part of a CRLF line ending, or at the
  very end of the file. Files with carriage-return-only line endings are now
  rejected instead of misread.

#### SAM (brust-sam)

- UTF-8 is accepted where SAMv1 allows it: `@SQ DS`, `@RG DS`, `@PG CL`,
  `@PG DS` and `@CO`.
- `@CO` lines containing tabs now round-trip. `@CO` lines containing NUL are
  refused.
- The writer refuses fields that contain a tab, which would split the line.
- Hex (`H`) values with non-ASCII characters are rejected instead of causing a
  panic.
- `B` arrays with empty elements are rejected.
- Zero-length CIGAR operations are accepted, as SAMv1 and htslib allow.

#### BAM and BGZF (brust-bam)

- Records with more than 65,535 CIGAR operations use the SAMv1 long-CIGAR
  convention (a `CG` tag) when read, written and summarised in statistics.
- Memory is only reserved for data actually present, not from lengths stored in
  the file. A small malformed file can no longer abort the process.
- BGZF blocks over 64 KiB, and bytes after a block's deflate data, are rejected.
- As in samtools, the reader rejects:
  - reference IDs outside the reference list
  - mapped records that have both a sequence and a CIGAR, when their lengths
    disagree
- `bin` values are computed as htslib computes them, including for unmapped
  records and very long references.
- Encoding fixes:
  - CIGAR operation lengths over 28 bits are rejected instead of truncated.
  - Unrecognised sequence characters are stored as `N`, and `U` as `T`, as
    htslib does.
  - Read names of 255 bytes or more are rejected.
- `BamWriter`:
  - writes one reference list per stream
  - refuses records whose reference IDs are outside that list
  - refuses `commit` before a header has been written
- Reference names that aren't UTF-8 are rejected.
- BAM errors carry Brust diagnostics. Faults in the SAM input to a SAM-to-BAM
  conversion are reported as SAM errors.
- The strict BGZF end-of-file check still runs after a temporary read error.

#### Core and CLI (brust-core, brust)

- Atomic output works for target names close to the 255-byte limit.
- Atomic output refuses targets that are folders, FIFOs, sockets or devices.
- The CLI no longer panics when its output pipe closes early (for example
  `brust stats ... | head -1`).
- FASTA and SAM conversions refuse an output name ending in `.gz`, instead of
  writing plain text into it.
- `bam-to-sam` builds the SAM header from the BAM header text and checks it.
  When the text has no `@SQ` lines, it adds them from the binary reference
  list. Header text that isn't UTF-8 is refused.
- FASTQ-to-BAM and BAM-to-FASTQ keep a single base of quality 9. That quality
  is the same character as SAM's `*`, so it used to be lost.
- Empty FASTQ reads convert to SAM and BAM with SEQ and QUAL set to `*`.
- `validate sam` checks RNAME and RNEXT against the `@SQ` lines.
- `validate bam` checks the header text, and that every record converts to a
  valid SAM line. A BAM that validates now converts with `bam-to-sam`.
- Statistics:
  - SAM and BAM query lengths agree.
  - The MAPQ summary leaves out MAPQ 255 and unmapped records.
  - GC fraction is counted over A, C, G, T and U only.

### Performance

Measured on release builds:

- Writers created with `from_path` are buffered. Writing 1 million reads:

  | Format | Before | After |
  | --- | --- | --- |
  | FASTQ | 1.05 s | 0.24 s |
  | FASTA | 1.2 s | 0.12 s |
  | SAM | 0.8 s | 0.44 s |

- `stats` uses much less memory:
  - It keeps a 128-bit fingerprint of each read ID instead of the ID itself.
    Two different IDs share a fingerprint with a probability below 10⁻²⁰ for a
    billion reads.
  - It keeps a histogram of lengths instead of every length. N50 and N90 stay
    exact.
  - Peak memory for 1 million FASTQ reads drops from 116 MB to 57 MB. For a
    BAM of 1 million reads it drops from 19 MB to 4 MB.
- POD5 signal reads jump straight to the Signal batch that holds each row, so
  reading in any order costs about the same as reading in file order. Reading
  4,000 reads in reverse order takes 0.075 s instead of 6.5 s.
- The POD5 writer streams the Signal table to the output in batches of 100
  rows, as official writers do, instead of building the whole file in memory.
  Writing a 75 MB payload needs 5.5 MB of extra memory instead of 380 MB.
- SAM records are checked directly instead of being formatted and parsed again.
  Output is unchanged. On 500,000 records:
  - `bam-to-sam` is 1.8 times faster.
  - `validate bam` is 2.5 times faster.
  - `sam-to-bam` with 4 threads is 2 times faster.

### Changed

- When a SAM alignment line has more than one error, a different error may be
  reported. Errors in numbers and optional-field syntax now come first.
- Files written by `Pod5Writer` hold several Signal batches of 100 rows instead
  of one. Official pod5 and dorado read them.
