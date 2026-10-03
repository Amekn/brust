//! POD5 reader, writer, and signal primitives.
//!
//! The crate exposes both cursor-style and materialized APIs:
//!
//! - [`Pod5Reader`] reads POD5 and Arrow metadata during construction, then
//!   returns one [`Pod5Record`] at a time from the Reads table without
//!   materializing the whole POD5 file.
//! - [`Pod5`] owns the parsed header, Run Info rows, Signal rows, and Reads
//!   rows and can be cloned or written back out.
//!
//! POD5 files are a container around Apache Arrow IPC/Feather V2 tables.
//! This crate validates wrapper magic, section markers, zero padding, required
//! table presence, and Arrow schema shapes, then parses the embedded Reads,
//! Signal, and Run Info tables through Arrow's IPC reader. Seekable readers use
//! the POD5 footer to open embedded Arrow sections on demand. The writer emits
//! a complete POD5 payload from the materialized representation and compresses
//! uncompressed signal rows to VBZ.

use arrow_array::builder::{
    FixedSizeBinaryBuilder, LargeBinaryBuilder, ListBuilder, MapBuilder, StringBuilder,
    StringDictionaryBuilder, UInt64Builder,
};
use arrow_array::types::Int16Type;
use arrow_array::{
    Array, ArrayRef, BooleanArray, DictionaryArray, FixedSizeBinaryArray, Float32Array, Int16Array,
    LargeBinaryArray, LargeListArray, ListArray, MapArray, RecordBatch, StringArray,
    TimestampMillisecondArray, UInt8Array, UInt16Array, UInt32Array, UInt64Array,
};
use arrow_ipc::reader::FileReader;
use arrow_ipc::writer::FileWriter;
use arrow_ipc::{MessageHeader, root_as_footer, root_as_message};
use arrow_schema::{ArrowError, DataType, Field, Schema, SchemaRef};
use brust_core::{AtomicFile, Error, Format};
use std::cell::RefCell;
use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{self, BufWriter, Cursor, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::rc::Rc;
use std::sync::Arc;

/// POD5 file signature stored at the start and end of each POD5 file.
pub const POD5_MAGIC: &[u8; 8] = b"\x8bPOD\r\n\x1a\n";
/// Marker at the start of the footer section.
pub const POD5_FOOTER_MAGIC: &[u8; 8] = b"FOOTER\0\0";
/// Arrow IPC/Feather magic used by embedded table sections.
pub const ARROW_MAGIC: &[u8; 6] = b"ARROW1";

/// Field names expected in the POD5 Reads table.
pub const READS_TABLE_FIELDS: &[&str] = &[
    "read_id",
    "signal",
    "read_number",
    "start",
    "median_before",
    "num_minknow_events",
    "tracked_scaling_scale",
    "tracked_scaling_shift",
    "predicted_scaling_scale",
    "predicted_scaling_shift",
    "num_reads_since_mux_change",
    "time_since_mux_change",
    "num_samples",
    "channel",
    "well",
    "pore_type",
    "calibration_offset",
    "calibration_scale",
    "end_reason",
    "end_reason_forced",
    "run_info",
];

/// Field names expected in the POD5 Signal table.
pub const SIGNAL_TABLE_FIELDS: &[&str] = &["read_id", "signal", "samples"];

/// Field names expected in the POD5 Run Info table.
pub const RUN_INFO_TABLE_FIELDS: &[&str] = &[
    "acquisition_id",
    "acquisition_start_time",
    "adc_max",
    "adc_min",
    "context_tags",
    "experiment_name",
    "flow_cell_id",
    "flow_cell_product_code",
    "protocol_name",
    "protocol_run_id",
    "protocol_start_time",
    "sample_id",
    "sample_rate",
    "sequencing_kit",
    "sequencer_position",
    "sequencer_position_type",
    "software",
    "system_name",
    "system_type",
    "tracking_id",
];

/// A fully materialized POD5 payload.
///
/// `Pod5` owns its header, run-info rows, and all parsed read records, so
/// cloning this type performs a deep copy of the parsed POD5 data.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5 {
    /// POD5 wrapper and table metadata.
    pub header: Pod5Header,
    /// Run Info table rows.
    pub run_infos: Vec<Pod5RunInfo>,
    /// Signal table rows.
    pub signals: Vec<Pod5Signal>,
    /// Reads table rows.
    pub records: Vec<Pod5Record>,
}

/// Cursor over POD5 read records.
///
/// POD5 stores reads, signals, and run metadata as embedded Arrow tables with a
/// footer. This reader parses metadata during construction and then exposes the
/// same `read_record`/`records` surface as the streaming crates. The underlying
/// byte stream must support seeking so the reader can use footer offsets to
/// open embedded Arrow sections on demand.
pub struct Pod5Reader<R: Read + Seek = File> {
    /// POD5 wrapper and table metadata.
    pub header: Pod5Header,
    /// Run Info table rows parsed during construction.
    pub run_infos: Vec<Pod5RunInfo>,
    shared: SharedReader<R>,
    read_sections: Vec<Pod5Section>,
    signal_sections: Vec<Pod5Section>,
    read_section_index: usize,
    read_reader: Option<FileReader<SectionReader<R>>>,
    read_buffer: VecDeque<Pod5Record>,
    signal_cursor: Pod5SignalCursor<R>,
    // Samples of the last read served by `signal_for_record`, so asking for the
    // same read again doesn't read its batch again. One read, not the file.
    last_signal: Option<LastSignal>,
}

struct LastSignal {
    read_id: String,
    signal_rows: Vec<u64>,
    samples: Vec<i16>,
}

/// POD5 file writer over any writable byte stream.
///
/// POD5 is a container of Arrow IPC files plus a footer, so this writer emits
/// one complete materialized [`Pod5`] payload rather than streaming individual
/// reads. The Signal table is written a batch at a time, not built in memory
/// first. Missing writer metadata is filled with deterministic defaults.
pub struct Pod5Writer<W: Write = File> {
    writer: W,
    // A POD5 file holds exactly one payload; a second one would corrupt it.
    payload_written: bool,
}

/// High-level counts and rollups for a POD5 payload.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pod5Summary {
    /// Number of Reads table rows.
    pub read_count: usize,
    /// Number of Signal table rows.
    pub signal_count: usize,
    /// Number of Run Info rows.
    pub run_info_count: usize,
    /// Sum of `num_samples` across all reads.
    pub total_samples: u64,
    /// Per-channel read and sample counts.
    pub channels: Vec<Pod5ChannelSummary>,
    /// Per-run read and sample counts.
    pub run_infos: Vec<Pod5RunInfoSummary>,
}

/// Per-channel POD5 summary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pod5ChannelSummary {
    /// One-indexed channel.
    pub channel: u16,
    /// Number of reads observed on this channel.
    pub read_count: usize,
    /// Sum of `num_samples` for reads on this channel.
    pub sample_count: u64,
}

/// Per-run POD5 summary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pod5RunInfoSummary {
    /// Acquisition/run identifier.
    pub acquisition_id: String,
    /// User-supplied sample identifier.
    pub sample_id: String,
    /// User-supplied experiment name.
    pub experiment_name: String,
    /// Flow-cell identifier.
    pub flow_cell_id: String,
    /// Sequencing kit name.
    pub sequencing_kit: String,
    /// Samples per second.
    pub sample_rate: u16,
    /// MinKNOW/software string from the row.
    pub software: String,
    /// Number of reads referencing this run info row.
    pub read_count: usize,
    /// Sum of `num_samples` for reads referencing this run info row.
    pub sample_count: u64,
}

/// Lazily decompresses and caches signal rows for a materialized [`Pod5`].
pub struct Pod5SignalCache<'a> {
    signals: &'a [Pod5Signal],
    decoded_cache: RefCell<HashMap<u64, Vec<i16>>>,
}

/// POD5 wrapper and table metadata.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5Header {
    /// Leading POD5 magic bytes.
    pub magic: [u8; 8],
    /// Per-file section marker used between embedded sections.
    pub section_marker: [u8; 16],
    /// Embedded Arrow and footer sections in file order.
    pub sections: Vec<Pod5Section>,
    /// `MINKNOW:file_identifier` Arrow schema metadata, when present.
    pub file_identifier: Option<String>,
    /// `MINKNOW:software` Arrow schema metadata, when present.
    pub software: Option<String>,
    /// `MINKNOW:pod5_version` Arrow schema metadata, when present.
    ///
    /// Writers stamp this version (`0.3.34` when `None`) and write the Reads
    /// columns it requires: `open_pore_level` from 0.3.30 on.
    pub pod5_version: Option<String>,
}

/// One POD5 wrapper section.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5Section {
    /// Section kind inferred from the embedded Arrow schema or footer magic.
    pub kind: Pod5SectionKind,
    /// Absolute byte offset just after the preceding section marker.
    pub offset: u64,
    /// Arrow payload length, or FlatBuffers footer payload length for footer sections.
    pub length: u64,
    /// Padded section length up to the next section marker.
    pub padded_length: u64,
    /// Number of Arrow table rows in this section, or zero for the footer.
    pub row_count: usize,
}

/// POD5 wrapper section kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pod5SectionKind {
    /// Reads table.
    Reads,
    /// Signal table.
    Signal,
    /// Run Info table.
    RunInfo,
    /// FlatBuffers footer.
    Footer,
    /// Unknown Arrow table.
    Unknown,
}

/// A POD5 Run Info table row.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5RunInfo {
    /// Acquisition/run identifier.
    pub acquisition_id: String,
    /// User-supplied sample identifier.
    pub sample_id: String,
    /// User-supplied experiment name.
    pub experiment_name: String,
    /// Flow-cell identifier.
    pub flow_cell_id: String,
    /// Sequencing kit name.
    pub sequencing_kit: String,
    /// Samples per second.
    pub sample_rate: u16,
    /// MinKNOW/software string from the row.
    pub software: String,
    /// Acquisition start, in milliseconds since the Unix epoch (UTC).
    pub acquisition_start_time: i64,
    /// Maximum ADC value.
    pub adc_max: i16,
    /// Minimum ADC value.
    pub adc_min: i16,
    /// MinKNOW context tags, as key/value pairs in file order.
    pub context_tags: Vec<(String, String)>,
    /// Flow-cell product code, such as `FLO-MIN114`.
    pub flow_cell_product_code: String,
    /// Sequencing protocol name.
    pub protocol_name: String,
    /// Protocol run identifier.
    pub protocol_run_id: String,
    /// Protocol start, in milliseconds since the Unix epoch (UTC).
    pub protocol_start_time: i64,
    /// Sequencer position, such as a device or flow-cell slot name.
    pub sequencer_position: String,
    /// Sequencer position type, such as `MinION Mk1B`.
    pub sequencer_position_type: String,
    /// Host system name.
    pub system_name: String,
    /// Host system type.
    pub system_type: String,
    /// MinKNOW tracking ID, as key/value pairs in file order.
    pub tracking_id: Vec<(String, String)>,
}

/// A POD5 Reads table row.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5Record {
    /// Read UUID.
    pub read_id: String,
    /// Zero-based Signal table row indices referenced by this read.
    pub signal_rows: Vec<u64>,
    /// Read number.
    pub read_number: u32,
    /// Sample offset on this channel at which the read starts.
    pub start_sample: u64,
    /// Current level in this well before the read.
    pub median_before: f32,
    /// Number of MinKNOW events for this read.
    pub num_minknow_events: u64,
    /// Tracked scaling scale.
    pub tracked_scaling_scale: f32,
    /// Tracked scaling shift.
    pub tracked_scaling_shift: f32,
    /// Predicted scaling scale.
    pub predicted_scaling_scale: f32,
    /// Predicted scaling shift.
    pub predicted_scaling_shift: f32,
    /// Number of selected reads since the last mux change on this channel.
    pub num_reads_since_mux_change: u32,
    /// Seconds since the last mux change on this channel.
    pub time_since_mux_change: f32,
    /// Number of signal samples in the read.
    pub num_samples: u64,
    /// One-indexed channel.
    pub channel: u16,
    /// One-indexed well/mux.
    pub well: u8,
    /// Pore type string.
    pub pore_type: String,
    /// Calibration offset.
    pub calibration_offset: f32,
    /// Calibration scale.
    pub calibration_scale: f32,
    /// End reason string.
    pub end_reason: String,
    /// Whether the end reason was forced.
    pub end_reason_forced: bool,
    /// Run Info acquisition identifier referenced by this read.
    pub run_info: String,
    /// Open pore level, or `None` when unknown. Files older than POD5 0.3.30
    /// have no such column; a stored NaN also reads as `None`.
    pub open_pore_level: Option<f32>,
}

/// A POD5 Signal table row.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5Signal {
    /// Read UUID associated with this signal row.
    pub read_id: String,
    /// Number of decoded ADC samples in this row.
    pub samples: u32,
    /// Signal payload, compressed or uncompressed depending on the parsed source.
    pub payload: Pod5SignalPayload,
}

/// Signal payload representation.
#[derive(Debug, Clone, PartialEq)]
pub enum Pod5SignalPayload {
    /// VBZ-compressed signal bytes from the POD5 Signal table.
    Vbz(Vec<u8>),
    /// Uncompressed int16 ADC samples. The writer compresses these rows to VBZ.
    Uncompressed(Vec<i16>),
}

impl Pod5 {
    /// Opens a POD5 file and materializes all records into memory.
    ///
    /// This uses the seek-based [`Pod5Reader`] internally, so the whole file is
    /// not buffered before parsing. The returned [`Pod5`] still owns all parsed
    /// reads and signal rows.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Pod5Reader::from_path(path)?.read_all()
    }

    /// Materializes a POD5 byte stream into memory.
    ///
    /// Non-seekable streams cannot use the POD5 footer for section offsets, so
    /// this compatibility API buffers the input, then parses it exactly as
    /// [`Pod5::from_path`] does. Use [`Pod5Reader::from_reader`] for seekable
    /// streaming reads.
    pub fn from_reader<R: Read>(reader: R) -> io::Result<Self> {
        let mut data = Vec::new();
        let mut reader = reader;
        reader.read_to_end(&mut data)?;
        Pod5Reader::from_reader(Cursor::new(data))?.read_all()
    }

    /// Writes this POD5 payload to a filesystem path.
    ///
    /// The payload is validated before the file is opened, so an invalid
    /// payload leaves an existing file untouched. A failed write can still
    /// leave a partial file; use [`Pod5::to_path_atomic`] to avoid that.
    pub fn to_path<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let prepared = prepare_pod5(self)?;
        prepared.write_to(File::create(path)?)
    }

    /// Writes this POD5 payload to a filesystem path atomically.
    ///
    /// Nothing appears at `path` unless the whole payload is written and
    /// committed. See [`Pod5Writer::from_path_atomic`].
    pub fn to_path_atomic<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut writer = Pod5Writer::from_path_atomic(path)?;
        writer.write_all(self)?;
        writer.commit()
    }

    /// Writes this POD5 payload to a writable byte stream.
    pub fn to_writer<W: Write>(&self, writer: W) -> io::Result<()> {
        let mut writer = Pod5Writer::from_writer(writer);
        writer.write_all(self)?;
        writer.flush()
    }

    /// Decompresses and concatenates the signal rows referenced by `record`.
    ///
    /// Returns an error if any referenced Signal row index is out of bounds or a
    /// compressed signal payload is malformed.
    pub fn signal_for_record(&self, record: &Pod5Record) -> io::Result<Vec<i16>> {
        record.signal(&self.signals)
    }

    /// Finds a read by ID and returns its decompressed signal, if present.
    pub fn signal_by_read_id(&self, read_id: &str) -> io::Result<Option<Vec<i16>>> {
        self.read_by_id(read_id)
            .map(|record| self.signal_for_record(record))
            .transpose()
    }

    /// Returns a read by ID, if present.
    pub fn read_by_id(&self, read_id: &str) -> Option<&Pod5Record> {
        self.records.iter().find(|record| record.read_id == read_id)
    }

    /// Builds a lookup map from read ID to read record.
    pub fn read_lookup(&self) -> HashMap<&str, &Pod5Record> {
        self.records
            .iter()
            .map(|record| (record.read_id.as_str(), record))
            .collect()
    }

    /// Creates a reusable lazy signal decompression cache over this payload's
    /// Signal table rows.
    pub fn signal_cache(&self) -> Pod5SignalCache<'_> {
        Pod5SignalCache::new(&self.signals)
    }

    /// Returns the sum of `num_samples` across all reads.
    ///
    /// The counts come from the file, so the sum saturates at `u64::MAX`
    /// rather than overflowing; so do the per-channel and per-run sums.
    pub fn total_samples(&self) -> u64 {
        self.records.iter().fold(0u64, |total, record| {
            total.saturating_add(record.num_samples)
        })
    }

    /// Returns per-channel read and sample summaries sorted by channel.
    pub fn channel_summaries(&self) -> Vec<Pod5ChannelSummary> {
        let mut counts = HashMap::<u16, (usize, u64)>::new();
        for record in &self.records {
            let entry = counts.entry(record.channel).or_default();
            entry.0 += 1;
            entry.1 = entry.1.saturating_add(record.num_samples);
        }

        let mut summaries = counts
            .into_iter()
            .map(|(channel, (read_count, sample_count))| Pod5ChannelSummary {
                channel,
                read_count,
                sample_count,
            })
            .collect::<Vec<_>>();
        summaries.sort_by_key(|summary| summary.channel);
        summaries
    }

    /// Returns per-run read and sample summaries in Run Info table order.
    pub fn run_info_summaries(&self) -> Vec<Pod5RunInfoSummary> {
        let mut counts = HashMap::<&str, (usize, u64)>::new();
        for record in &self.records {
            let entry = counts.entry(record.run_info.as_str()).or_default();
            entry.0 += 1;
            entry.1 = entry.1.saturating_add(record.num_samples);
        }

        self.run_infos
            .iter()
            .map(|run_info| {
                let (read_count, sample_count) = counts
                    .get(run_info.acquisition_id.as_str())
                    .copied()
                    .unwrap_or_default();
                Pod5RunInfoSummary {
                    acquisition_id: run_info.acquisition_id.clone(),
                    sample_id: run_info.sample_id.clone(),
                    experiment_name: run_info.experiment_name.clone(),
                    flow_cell_id: run_info.flow_cell_id.clone(),
                    sequencing_kit: run_info.sequencing_kit.clone(),
                    sample_rate: run_info.sample_rate,
                    software: run_info.software.clone(),
                    read_count,
                    sample_count,
                }
            })
            .collect()
    }

    /// Returns a high-level summary of reads, signals, runs, channels, and sample counts.
    pub fn summary(&self) -> Pod5Summary {
        Pod5Summary {
            read_count: self.records.len(),
            signal_count: self.signals.len(),
            run_info_count: self.run_infos.len(),
            total_samples: self.total_samples(),
            channels: self.channel_summaries(),
            run_infos: self.run_info_summaries(),
        }
    }
}

impl Pod5Reader<File> {
    /// Opens a POD5 file from a filesystem path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        Self::from_reader(file)
    }

    /// Opens a POD5 file from a filesystem path.
    ///
    /// This is a convenience alias for [`Pod5Reader::from_path`].
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_path(path)
    }
}

impl<R: Read + Seek> Pod5Reader<R> {
    /// Creates a streaming POD5 reader from a seekable byte stream.
    ///
    /// POD5 stores table offsets in a footer at the end of the file, so true
    /// streaming requires [`Seek`]. Construction validates wrapper metadata,
    /// opens embedded Arrow sections, and parses Run Info rows, but Reads and
    /// Signal rows are loaded batch-by-batch on demand.
    pub fn from_reader(reader: R) -> io::Result<Self> {
        let shared = Rc::new(RefCell::new(reader));
        let (header, run_infos, section_batch_rows) = inspect_pod5(shared.clone())?;
        let read_sections = header
            .sections
            .iter()
            .filter(|section| section.kind == Pod5SectionKind::Reads)
            .cloned()
            .collect::<Vec<_>>();
        let signal_sections = header
            .sections
            .iter()
            .filter(|section| section.kind == Pod5SectionKind::Signal)
            .cloned()
            .collect::<Vec<_>>();
        let signal_cursor = Pod5SignalCursor::new(
            shared.clone(),
            header
                .sections
                .iter()
                .zip(section_batch_rows)
                .filter(|(section, _)| section.kind == Pod5SectionKind::Signal)
                .map(|(section, batch_rows)| (section.clone(), batch_rows))
                .collect(),
        );

        Ok(Self {
            header,
            run_infos,
            shared,
            read_sections,
            signal_sections,
            read_section_index: 0,
            read_reader: None,
            read_buffer: VecDeque::new(),
            signal_cursor,
            last_signal: None,
        })
    }

    /// Reads the next POD5 record.
    ///
    /// Returns `Ok(None)` when all parsed reads have been returned.
    pub fn read_record(&mut self) -> io::Result<Option<Pod5Record>> {
        loop {
            if let Some(record) = self.read_buffer.pop_front() {
                return Ok(Some(record));
            }

            let Some(batch) = self.next_reads_batch()? else {
                return Ok(None);
            };
            self.read_buffer.extend(parse_reads_batch(&batch)?);
        }
    }

    /// Reads the next POD5 record.
    ///
    /// This is a compatibility alias for [`Pod5Reader::read_record`].
    pub fn read(&mut self) -> io::Result<Option<Pod5Record>> {
        self.read_record()
    }

    /// Returns an iterator over the remaining read records.
    pub fn records(&mut self) -> Pod5Records<'_, R> {
        Pod5Records { reader: self }
    }

    /// Decompresses and concatenates the signal rows referenced by `record`.
    ///
    /// Signal rows are read lazily from the Signal table, one Arrow batch at a
    /// time: a row is found from the batch row counts in the file's Arrow
    /// metadata, so reads can be asked for in any order without reading the
    /// batches before them. Only the last read's samples are kept, so memory
    /// stays bounded while streaming and asking for the same read twice doesn't
    /// re-read the file.
    pub fn signal_for_record(&mut self, record: &Pod5Record) -> io::Result<Vec<i16>> {
        let total = usize::try_from(record.num_samples)
            .map_err(|_| invalid_data("POD5 read sample count exceeds usize"))?;
        let samples = match &self.last_signal {
            Some(last)
                if last.read_id == record.read_id && last.signal_rows == record.signal_rows =>
            {
                last.samples.clone()
            }
            _ => {
                let mut samples = Vec::new();
                for &index in &record.signal_rows {
                    let signal = self.signal_cursor.signal_row_at(index)?;
                    if signal.read_id != record.read_id {
                        return Err(invalid_data("POD5 signal row read_id does not match read"));
                    }
                    samples.extend(signal.decompress()?);
                }
                self.last_signal = Some(LastSignal {
                    read_id: record.read_id.clone(),
                    signal_rows: record.signal_rows.clone(),
                    samples: samples.clone(),
                });
                samples
            }
        };

        if samples.len() != total {
            return Err(invalid_data(
                "POD5 read num_samples does not match referenced signal rows",
            ));
        }

        Ok(samples)
    }

    /// Reads one Signal table row by zero-based row number.
    pub fn signal_row(&mut self, row: u64) -> io::Result<Pod5Signal> {
        self.signal_cursor.signal_row_at(row)
    }

    /// Consumes this reader and materializes the remaining POD5 records.
    pub fn read_all(mut self) -> io::Result<Pod5> {
        let signals = self.read_all_signals_from_start()?;
        let mut records = Vec::new();
        while let Some(record) = self.read_record()? {
            records.push(record);
        }

        Ok(Pod5 {
            header: self.header,
            run_infos: self.run_infos,
            signals,
            records,
        })
    }

    fn next_reads_batch(&mut self) -> io::Result<Option<RecordBatch>> {
        loop {
            if let Some(reader) = &mut self.read_reader {
                match arrow_call(|| reader.next().transpose())? {
                    Some(batch) => return Ok(Some(batch)),
                    None => {
                        self.read_reader = None;
                        self.read_section_index += 1;
                    }
                }
            } else if let Some(section) = self.read_sections.get(self.read_section_index) {
                self.read_reader = Some(open_arrow_reader(self.shared.clone(), section)?);
            } else {
                return Ok(None);
            }
        }
    }

    fn read_all_signals_from_start(&self) -> io::Result<Vec<Pod5Signal>> {
        // The declared count comes from the file, so let the vector grow.
        let mut signals = Vec::new();

        for section in &self.signal_sections {
            let mut reader = open_arrow_reader(self.shared.clone(), section)?;
            while let Some(batch) = arrow_call(|| reader.next().transpose())? {
                signals.extend(parse_signal_batch(&batch)?);
            }
        }

        Ok(signals)
    }
}

impl Pod5Writer<File> {
    /// Creates or truncates a POD5 file at a filesystem path.
    ///
    /// Use [`Pod5Writer::from_path_atomic`] to keep an existing file intact until
    /// the new one is complete.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self::from_writer(file))
    }

    /// Creates or truncates a POD5 file at a filesystem path.
    ///
    /// This is a convenience alias for [`Pod5Writer::from_path`].
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_path(path)
    }
}

impl Pod5Writer<AtomicFile> {
    /// Creates a POD5 writer that publishes its output atomically.
    ///
    /// Nothing appears at `path` until [`commit`](Self::commit) renames the
    /// finished file into place; dropping the writer discards the output. See
    /// [`AtomicFile`] for the one error `commit` can return after the rename.
    pub fn from_path_atomic<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self::from_writer(AtomicFile::create(path)?))
    }

    /// Flushes, syncs and renames the output into place.
    ///
    /// Returns an [`io::ErrorKind::InvalidInput`] error and publishes nothing
    /// if no payload was written, since an empty file is not a valid POD5. See
    /// [`AtomicFile::commit`] for the exact steps and errors.
    pub fn commit(self) -> io::Result<()> {
        if !self.payload_written {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "POD5 writer has no payload to commit",
            ));
        }
        self.into_inner().commit()
    }
}

impl<W: Write> Pod5Writer<W> {
    /// Creates a POD5 writer from a writable byte stream.
    pub fn from_writer(writer: W) -> Self {
        Self {
            writer,
            payload_written: false,
        }
    }

    /// Writes a complete materialized POD5 payload.
    ///
    /// Before writing anything, the writer checks UUIDs and header metadata,
    /// that every signal row decodes to its declared `samples` (VBZ rows are
    /// decompressed to check), that each read's signal rows exist and carry
    /// the read's ID, read sample counts, end reasons, and run-info
    /// references. Signal rows are then written in Arrow batches of 100 rows,
    /// as official POD5 writers do. A POD5 file holds one payload, so a second
    /// call returns an [`io::ErrorKind::InvalidInput`] error without writing
    /// anything.
    pub fn write_all(&mut self, pod5: &Pod5) -> io::Result<()> {
        if self.payload_written {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "POD5 writer already wrote a payload; a POD5 file holds one payload",
            ));
        }
        let prepared = prepare_pod5(pod5)?;
        // Set before writing: after a failed write the output is already
        // partial, and appending another payload would not repair it.
        self.payload_written = true;
        prepared.write_to(&mut self.writer)
    }

    /// Writes a complete materialized POD5 payload.
    ///
    /// This is a compatibility alias for [`Pod5Writer::write_all`].
    pub fn write(&mut self, pod5: &Pod5) -> io::Result<()> {
        self.write_all(pod5)
    }

    /// Flushes the underlying writer.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }

    /// Consumes this writer and returns the wrapped byte stream.
    pub fn into_inner(self) -> W {
        self.writer
    }
}

impl Pod5Record {
    /// Decompresses and concatenates the signal rows referenced by this read.
    ///
    /// `signal_rows` are interpreted as indices into the supplied `signals`
    /// slice.
    pub fn signal(&self, signals: &[Pod5Signal]) -> io::Result<Vec<i16>> {
        let total = usize::try_from(self.num_samples)
            .map_err(|_| invalid_data("POD5 read sample count exceeds usize"))?;
        let mut samples = Vec::new();

        for &index in &self.signal_rows {
            let signal = signals
                .get(index as usize)
                .ok_or_else(|| invalid_data("POD5 signal row index is out of bounds"))?;
            if signal.read_id != self.read_id {
                return Err(invalid_data("POD5 signal row read_id does not match read"));
            }
            samples.extend(signal.decompress()?);
        }

        if samples.len() != total {
            return Err(invalid_data(
                "POD5 read num_samples does not match referenced signal rows",
            ));
        }

        Ok(samples)
    }
}

impl Pod5Signal {
    /// Returns `true` when this row stores VBZ-compressed signal bytes.
    pub fn is_vbz_compressed(&self) -> bool {
        matches!(self.payload, Pod5SignalPayload::Vbz(_))
    }

    /// Returns compressed VBZ bytes when this row is compressed.
    pub fn compressed_bytes(&self) -> Option<&[u8]> {
        match &self.payload {
            Pod5SignalPayload::Vbz(data) => Some(data),
            Pod5SignalPayload::Uncompressed(_) => None,
        }
    }

    /// Decompresses this signal row to raw int16 ADC samples.
    pub fn decompress(&self) -> io::Result<Vec<i16>> {
        match &self.payload {
            Pod5SignalPayload::Vbz(data) => decompress_vbz_signal(data, self.samples as usize),
            Pod5SignalPayload::Uncompressed(samples) => Ok(samples.clone()),
        }
    }

    /// Encodes this row's samples as a VBZ-compressed signal blob.
    ///
    /// Existing VBZ payloads are decompressed before recompression.
    pub fn compress(&self) -> io::Result<Vec<u8>> {
        compress_vbz_signal(&self.decompress()?)
    }
}

impl<'a> Pod5SignalCache<'a> {
    /// Creates a cache over Signal table rows.
    pub fn new(signals: &'a [Pod5Signal]) -> Self {
        Self {
            signals,
            decoded_cache: RefCell::new(HashMap::new()),
        }
    }

    /// Decompresses one Signal row by zero-based row index, reusing cached
    /// samples on subsequent calls.
    pub fn signal_row(&self, row: u64) -> io::Result<Vec<i16>> {
        if let Some(samples) = self.decoded_cache.borrow().get(&row) {
            return Ok(samples.clone());
        }

        let signal = self
            .signals
            .get(row as usize)
            .ok_or_else(|| invalid_data("POD5 signal row index is out of bounds"))?;
        let samples = signal.decompress()?;
        self.decoded_cache.borrow_mut().insert(row, samples.clone());
        Ok(samples)
    }

    /// Decompresses and concatenates signal rows referenced by a read record.
    pub fn signal_for_record(&self, record: &Pod5Record) -> io::Result<Vec<i16>> {
        let total = usize::try_from(record.num_samples)
            .map_err(|_| invalid_data("POD5 read sample count exceeds usize"))?;
        let mut samples = Vec::new();

        for &row in &record.signal_rows {
            let signal = self
                .signals
                .get(row as usize)
                .ok_or_else(|| invalid_data("POD5 signal row index is out of bounds"))?;
            if signal.read_id != record.read_id {
                return Err(invalid_data("POD5 signal row read_id does not match read"));
            }
            samples.extend(self.signal_row(row)?);
        }

        if samples.len() != total {
            return Err(invalid_data(
                "POD5 read num_samples does not match referenced signal rows",
            ));
        }

        Ok(samples)
    }

    /// Returns the number of signal rows currently cached.
    pub fn cached_row_count(&self) -> usize {
        self.decoded_cache.borrow().len()
    }
}

/// Iterator over read records from a [`Pod5Reader`].
pub struct Pod5Records<'a, R: Read + Seek> {
    reader: &'a mut Pod5Reader<R>,
}

impl<R: Read + Seek> Iterator for Pod5Records<'_, R> {
    type Item = io::Result<Pod5Record>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(error) => Some(Err(error)),
        }
    }
}

impl Pod5Header {
    /// Returns the number of parsed Reads table rows.
    pub fn read_count(&self) -> usize {
        self.row_count(Pod5SectionKind::Reads)
    }

    /// Returns the number of parsed Signal table rows.
    pub fn signal_count(&self) -> usize {
        self.row_count(Pod5SectionKind::Signal)
    }

    /// Returns the number of parsed Run Info table rows.
    pub fn run_info_count(&self) -> usize {
        self.row_count(Pod5SectionKind::RunInfo)
    }

    fn row_count(&self, kind: Pod5SectionKind) -> usize {
        self.sections
            .iter()
            .filter(|section| section.kind == kind)
            .map(|section| section.row_count)
            .sum()
    }
}

type SharedReader<R> = Rc<RefCell<R>>;

#[derive(Clone)]
struct SectionReader<R: Read + Seek> {
    reader: SharedReader<R>,
    offset: u64,
    length: u64,
    position: u64,
}

impl<R: Read + Seek> SectionReader<R> {
    fn new(reader: SharedReader<R>, section: &Pod5Section) -> Self {
        Self {
            reader,
            offset: section.offset,
            length: section.length,
            position: 0,
        }
    }
}

impl<R: Read + Seek> Read for SectionReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() || self.position >= self.length {
            return Ok(0);
        }

        let remaining = usize::try_from(self.length - self.position)
            .unwrap_or(usize::MAX)
            .min(buf.len());
        let mut reader = self.reader.borrow_mut();
        reader.seek(SeekFrom::Start(self.offset + self.position))?;
        let bytes = reader.read(&mut buf[..remaining])?;
        self.position += bytes as u64;
        Ok(bytes)
    }
}

impl<R: Read + Seek> Seek for SectionReader<R> {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        let base = match pos {
            SeekFrom::Start(position) => position as i128,
            SeekFrom::End(offset) => self.length as i128 + offset as i128,
            SeekFrom::Current(offset) => self.position as i128 + offset as i128,
        };

        if base < 0 || base > self.length as i128 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "POD5 section seek is out of bounds",
            ));
        }

        self.position = base as u64;
        Ok(self.position)
    }
}

/// Random access to Signal table rows.
///
/// Every Signal batch's row count, read from the Arrow block metadata when the
/// file was opened, maps a row number to its batch, so the cursor reads only
/// that batch. Rows are handed out once: the batch read last keeps the rows
/// not yet handed out, so reading on through it, or back to one of those rows,
/// doesn't read the file again, and it is dropped once all are handed out.
struct Pod5SignalCursor<R: Read + Seek> {
    shared: SharedReader<R>,
    sections: Vec<Pod5Section>,
    /// Every Signal batch, in row order.
    batches: Vec<SignalBatch>,
    /// An Arrow reader and the index of the section it reads.
    reader: Option<(usize, FileReader<SectionReader<R>>)>,
    loaded: Option<LoadedSignalBatch>,
}

/// The rows of the batch read last that haven't been handed out yet.
struct LoadedSignalBatch {
    /// The batch's index in `Pod5SignalCursor::batches`.
    number: usize,
    rows: Vec<Option<Pod5Signal>>,
    /// How many of `rows` are still `Some`.
    remaining: usize,
}

#[derive(Debug, Clone, Copy)]
struct SignalBatch {
    /// Index of the batch's section in `Pod5SignalCursor::sections`.
    section: usize,
    /// Index of the batch within its section.
    index: usize,
    first_row: u64,
    rows: u64,
}

impl<R: Read + Seek> Pod5SignalCursor<R> {
    /// Takes the Signal sections in file order with their batch row counts.
    fn new(shared: SharedReader<R>, sections: Vec<(Pod5Section, Vec<usize>)>) -> Self {
        let mut batches = Vec::new();
        let mut first_row = 0u64;
        for (section, (_, batch_rows)) in sections.iter().enumerate() {
            for (index, &rows) in batch_rows.iter().enumerate() {
                let rows = rows as u64;
                batches.push(SignalBatch {
                    section,
                    index,
                    first_row,
                    rows,
                });
                // Row counts are bounded by the bytes in the file.
                first_row = first_row.saturating_add(rows);
            }
        }

        Self {
            shared,
            sections: sections.into_iter().map(|(section, _)| section).collect(),
            batches,
            reader: None,
            loaded: None,
        }
    }

    /// Returns Signal row `row`, reading only the batch that holds it.
    ///
    /// Asking for a row in a batch whose rows fail to parse returns the parse
    /// error; rows in other batches are unaffected.
    fn signal_row_at(&mut self, row: u64) -> io::Result<Pod5Signal> {
        // The last batch starting at or before `row`; empty batches share the
        // next batch's first row, so they are never picked for a row inside it.
        let number = self
            .batches
            .partition_point(|batch| batch.first_row <= row)
            .checked_sub(1)
            .filter(|&number| row - self.batches[number].first_row < self.batches[number].rows)
            .ok_or_else(|| invalid_data("POD5 signal row index is out of bounds"))?;
        let batch = self.batches[number];
        let offset = (row - batch.first_row) as usize;

        let kept = self
            .loaded
            .as_mut()
            .filter(|loaded| loaded.number == number)
            .and_then(|loaded| loaded.rows[offset].take());
        let signal = match kept {
            Some(signal) => {
                if let Some(loaded) = &mut self.loaded {
                    loaded.remaining -= 1;
                }
                signal
            }
            None => {
                self.loaded = None;
                let mut rows = self
                    .read_batch(batch)?
                    .into_iter()
                    .map(Some)
                    .collect::<Vec<_>>();
                let signal = rows[offset].take().expect("the row is in its batch");
                self.loaded = Some(LoadedSignalBatch {
                    number,
                    remaining: rows.len() - 1,
                    rows,
                });
                signal
            }
        };
        if self
            .loaded
            .as_ref()
            .is_some_and(|loaded| loaded.remaining == 0)
        {
            self.loaded = None;
        }
        Ok(signal)
    }

    fn read_batch(&mut self, batch: SignalBatch) -> io::Result<Vec<Pod5Signal>> {
        if self
            .reader
            .as_ref()
            .is_none_or(|(section, _)| *section != batch.section)
        {
            self.reader = None;
            let reader = open_arrow_reader(self.shared.clone(), &self.sections[batch.section])?;
            self.reader = Some((batch.section, reader));
        }
        let (_, reader) = self.reader.as_mut().expect("the reader was just opened");
        let record_batch = match arrow_call(|| {
            reader.set_index(batch.index)?;
            reader.next().transpose()
        }) {
            Ok(record_batch) => record_batch,
            Err(error) => {
                // The reader may be part-way through a failed read; reopen it.
                self.reader = None;
                return Err(error);
            }
        };
        let record_batch =
            record_batch.ok_or_else(|| invalid_data("POD5 Signal batch is missing"))?;
        let signals = parse_signal_batch(&record_batch)?;
        if signals.len() as u64 != batch.rows {
            return Err(invalid_data(
                "POD5 Signal batch row count disagrees with its block metadata",
            ));
        }
        Ok(signals)
    }
}

fn open_arrow_reader<R: Read + Seek>(
    shared: SharedReader<R>,
    section: &Pod5Section,
) -> io::Result<FileReader<SectionReader<R>>> {
    arrow_call(|| FileReader::try_new(SectionReader::new(shared, section), None))
}

/// Runs an arrow-ipc read on file data, turning a panic into `InvalidData`.
///
/// arrow-ipc panics instead of returning an error on some malformed schemas
/// and buffers. [`validate_arrow_section`] rejects what would make it allocate
/// without bound, since an aborted allocation can't be caught, and explains
/// the common cases; this catches the rest.
fn arrow_call<T>(read: impl FnOnce() -> Result<T, ArrowError>) -> io::Result<T> {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(read))
        .map_err(|_| invalid_data("embedded Arrow data is malformed"))?
        .map_err(arrow_error)
}

/// Reads and checks a POD5 file's wrapper, footer and embedded Arrow
/// sections. Also returns each section's record batch row counts, in the
/// order of `Pod5Header::sections` (none for the footer).
fn inspect_pod5<R: Read + Seek>(
    shared: SharedReader<R>,
) -> io::Result<(Pod5Header, Vec<Pod5RunInfo>, Vec<Vec<usize>>)> {
    let file_len = shared.borrow_mut().seek(SeekFrom::End(0))?;
    if file_len < 8 + 16 + 8 + 16 + 8 {
        return Err(invalid_data("POD5 file is too short"));
    }

    let magic = slice_to_array::<8>(&read_exact_at(shared.clone(), 0, 8)?);
    if &magic != POD5_MAGIC {
        return Err(invalid_data("invalid POD5 leading magic"));
    }

    let trailing_magic = read_exact_at(shared.clone(), file_len - 8, 8)?;
    if trailing_magic.as_slice() != POD5_MAGIC {
        return Err(invalid_data("invalid POD5 trailing magic"));
    }

    let section_marker = slice_to_array::<16>(&read_exact_at(shared.clone(), 8, 16)?);
    let final_marker_offset = file_len - 8 - 16;
    let final_marker = read_exact_at(shared.clone(), final_marker_offset, 16)?;
    if final_marker.as_slice() != section_marker {
        return Err(invalid_data("POD5 final section marker is missing"));
    }

    let footer_len_offset = final_marker_offset - 8;
    let footer_len = i64::from_le_bytes(slice_to_array::<8>(&read_exact_at(
        shared.clone(),
        footer_len_offset,
        8,
    )?));
    if footer_len < 0 {
        return Err(invalid_data("POD5 footer length is negative"));
    }
    let footer_len = footer_len as u64;
    let footer_payload_start = footer_len_offset
        .checked_sub(footer_len)
        .ok_or_else(|| invalid_data("POD5 footer length exceeds file"))?;
    let footer_magic_start = footer_payload_start
        .checked_sub(POD5_FOOTER_MAGIC.len() as u64)
        .ok_or_else(|| invalid_data("POD5 footer magic offset is malformed"))?;

    let footer_magic = read_exact_at(shared.clone(), footer_magic_start, POD5_FOOTER_MAGIC.len())?;
    if footer_magic.as_slice() != POD5_FOOTER_MAGIC {
        return Err(invalid_data("POD5 footer magic is missing"));
    }

    let footer_padding = read_exact_at(shared.clone(), footer_payload_start, footer_len as usize)?;
    let footer = parse_pod5_footer(&footer_padding)?;
    let mut sections = Vec::with_capacity(footer.entries.len() + 1);
    let mut section_batch_rows = Vec::with_capacity(footer.entries.len() + 1);
    let mut run_infos = Vec::new();
    let mut metadata = Pod5Metadata {
        file_identifier: footer.file_identifier.clone(),
        software: footer.software.clone(),
        pod5_version: footer.pod5_version.clone(),
    };

    for entry in &footer.entries {
        let offset = u64::try_from(entry.offset)
            .map_err(|_| invalid_data("POD5 embedded file offset is negative"))?;
        let length = u64::try_from(entry.length)
            .map_err(|_| invalid_data("POD5 embedded file length is negative"))?;
        if length == 0 {
            return Err(invalid_data("POD5 embedded file is empty"));
        }
        let padded_length = checked_padded_len(length)?;
        let preceding_marker = offset
            .checked_sub(16)
            .ok_or_else(|| invalid_data("POD5 embedded file offset is malformed"))?;
        let following_marker = offset
            .checked_add(padded_length)
            .ok_or_else(|| invalid_data("POD5 embedded file length overflow"))?;

        if read_exact_at(shared.clone(), preceding_marker, 16)?.as_slice() != section_marker {
            return Err(invalid_data("POD5 section marker is malformed"));
        }
        if read_exact_at(shared.clone(), following_marker, 16)?.as_slice() != section_marker {
            return Err(invalid_data("POD5 section marker is malformed"));
        }
        if following_marker > footer_magic_start {
            return Err(invalid_data("POD5 embedded file overlaps footer"));
        }
        if padded_length > length {
            let padding = read_exact_at(
                shared.clone(),
                offset + length,
                usize::try_from(padded_length - length)
                    .map_err(|_| invalid_data("POD5 section padding exceeds usize"))?,
            )?;
            if !padding.iter().all(|byte| *byte == 0) {
                return Err(invalid_data("POD5 Arrow section padding is not zeroed"));
            }
        }

        // arrow-ipc trusts the block table and panics or aborts on bad values,
        // so check it before FileReader reads dictionary blocks.
        let batch_rows = validate_arrow_section(shared.clone(), offset, length)?;
        // validate_arrow_section checked that this sum fits.
        let row_count = batch_rows.iter().sum();
        let mut arrow_reader = arrow_call(|| {
            FileReader::try_new(
                SectionReader {
                    reader: shared.clone(),
                    offset,
                    length,
                    position: 0,
                },
                None,
            )
        })?;
        update_metadata(arrow_reader.custom_metadata(), &mut metadata)?;
        let schema = arrow_reader.schema();
        let kind = infer_section_kind(&schema);
        let footer_kind = section_kind_from_footer_content_type(entry.content_type);
        if footer_kind != Pod5SectionKind::Unknown && kind != footer_kind {
            return Err(invalid_data(
                "POD5 footer content type does not match Arrow schema",
            ));
        }

        if kind == Pod5SectionKind::RunInfo {
            while let Some(batch) = arrow_call(|| arrow_reader.next().transpose())? {
                run_infos.extend(parse_run_info_batch(&batch)?);
            }
        }

        sections.push(Pod5Section {
            kind,
            offset,
            length,
            padded_length,
            row_count,
        });
        section_batch_rows.push(batch_rows);
    }

    sections.push(Pod5Section {
        kind: Pod5SectionKind::Footer,
        offset: footer_magic_start,
        length: footer_len,
        padded_length: final_marker_offset - footer_magic_start,
        row_count: 0,
    });
    section_batch_rows.push(Vec::new());

    if sections
        .iter()
        .all(|section| section.kind != Pod5SectionKind::Reads)
    {
        return Err(invalid_data("POD5 Reads table is missing"));
    }
    if sections
        .iter()
        .all(|section| section.kind != Pod5SectionKind::Signal)
    {
        return Err(invalid_data("POD5 Signal table is missing"));
    }
    if sections
        .iter()
        .all(|section| section.kind != Pod5SectionKind::RunInfo)
    {
        return Err(invalid_data("POD5 Run Info table is missing"));
    }

    let header = Pod5Header {
        magic,
        section_marker,
        sections,
        file_identifier: metadata.file_identifier,
        software: metadata.software,
        pod5_version: metadata.pod5_version,
    };

    Ok((header, run_infos, section_batch_rows))
}

#[derive(Debug)]
struct ParsedPod5Footer {
    file_identifier: Option<String>,
    software: Option<String>,
    pod5_version: Option<String>,
    entries: Vec<Pod5FooterEntry>,
}

fn parse_pod5_footer(data: &[u8]) -> io::Result<ParsedPod5Footer> {
    let table = fb_root_table(data)?;
    let file_identifier = fb_string_field(data, table, 4)?;
    let software = fb_string_field(data, table, 6)?;
    let pod5_version = fb_string_field(data, table, 8)?;
    let mut entries = Vec::new();

    for entry_table in fb_table_vector_field(data, table, 10)? {
        entries.push(Pod5FooterEntry {
            offset: fb_i64_field(data, entry_table, 4)?
                .ok_or_else(|| invalid_data("POD5 footer embedded file offset is missing"))?,
            length: fb_i64_field(data, entry_table, 6)?
                .ok_or_else(|| invalid_data("POD5 footer embedded file length is missing"))?,
            content_type: fb_i16_field(data, entry_table, 10)?.unwrap_or_default(),
        });
    }

    Ok(ParsedPod5Footer {
        file_identifier,
        software,
        pod5_version,
        entries,
    })
}

fn section_kind_from_footer_content_type(content_type: i16) -> Pod5SectionKind {
    match content_type {
        0 => Pod5SectionKind::Reads,
        1 => Pod5SectionKind::Signal,
        4 => Pod5SectionKind::RunInfo,
        _ => Pod5SectionKind::Unknown,
    }
}

/// Checks an embedded Arrow file's footer blocks and the batch messages they
/// point to, then returns each record batch's row count in footer order, the
/// order arrow-ipc's `FileReader` reads them. Their sum fits in `usize`.
///
/// Every block must lie inside the section before the Arrow footer, every
/// buffer inside its block's body, and node and row counts must be
/// non-negative, with a record batch's first column as long as the batch.
/// Buffers must be long enough for their columns' lengths (see
/// [`ArrowBatchLayout`]), which bounds row counts by the bytes present; for
/// that, the first column must have a layout whose buffers grow with its
/// length (see [`column_bounds_its_length`]), and bodies must not be
/// compressed.
fn validate_arrow_section<R: Read + Seek>(
    shared: SharedReader<R>,
    offset: u64,
    length: u64,
) -> io::Result<Vec<usize>> {
    let footer = read_arrow_footer(shared.clone(), offset, length)?;
    // `read_arrow_footer` checked that the footer and its 10-byte trailer fit.
    let data_end = length - 10 - footer.len() as u64;
    let footer = root_as_footer(&footer)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))?;
    let schema_fields = footer
        .schema()
        .and_then(|schema| schema.fields())
        .ok_or_else(|| invalid_data("Arrow footer is missing its schema"))?;
    let dictionaries = footer.dictionaries().into_iter().flatten();
    let record_batches = footer.recordBatches().into_iter().flatten();

    let mut total = 0usize;
    let mut batch_rows = Vec::new();
    for (block, is_record_batch) in dictionaries
        .map(|block| (block, false))
        .chain(record_batches.map(|block| (block, true)))
    {
        let block_offset = u64::try_from(block.offset())
            .map_err(|_| invalid_data("Arrow block offset is negative"))?;
        let metadata_len = u64::try_from(block.metaDataLength())
            .map_err(|_| invalid_data("Arrow block metadata length is negative"))?;
        let body_len = u64::try_from(block.bodyLength())
            .map_err(|_| invalid_data("Arrow block body length is negative"))?;
        let block_end = block_offset
            .checked_add(metadata_len)
            .and_then(|end| end.checked_add(body_len))
            .ok_or_else(|| invalid_data("Arrow block length overflow"))?;
        if block_end > data_end {
            return Err(invalid_data("Arrow block extends past its section"));
        }

        let metadata = read_exact_at(shared.clone(), offset + block_offset, metadata_len as usize)?;
        let message = parse_arrow_message(&metadata)?;
        let (batch, dictionary_field) = match message.header_type() {
            MessageHeader::RecordBatch if is_record_batch => {
                (message.header_as_record_batch(), None)
            }
            MessageHeader::DictionaryBatch if !is_record_batch => {
                let dictionary = message
                    .header_as_dictionary_batch()
                    .ok_or_else(|| invalid_data("Arrow batch message is missing its header"))?;
                let field = find_dictionary_field(schema_fields.iter(), dictionary.id())
                    .ok_or_else(|| invalid_data("Arrow dictionary batch has no schema field"))?;
                (dictionary.data(), Some(field))
            }
            _ => return Err(invalid_data("Arrow block holds an unexpected message")),
        };
        let batch =
            batch.ok_or_else(|| invalid_data("Arrow batch message is missing its header"))?;

        for buffer in batch.buffers().into_iter().flatten() {
            let buffer_offset = u64::try_from(buffer.offset())
                .map_err(|_| invalid_data("Arrow buffer offset is negative"))?;
            let buffer_len = u64::try_from(buffer.length())
                .map_err(|_| invalid_data("Arrow buffer length is negative"))?;
            if buffer_offset
                .checked_add(buffer_len)
                .is_none_or(|end| end > body_len)
            {
                return Err(invalid_data("Arrow buffer extends past its block body"));
            }
        }
        let rows = usize::try_from(batch.length())
            .map_err(|_| invalid_data("Arrow record batch length is negative"))?;
        let nodes = batch.nodes().into_iter().flatten().collect::<Vec<_>>();
        if nodes.iter().any(|node| {
            node.length() < 0 || node.null_count() < 0 || node.null_count() > node.length()
        }) {
            return Err(invalid_data("Arrow field node counts are malformed"));
        }
        if nodes
            .first()
            .is_some_and(|node| node.length() != batch.length())
        {
            return Err(invalid_data(
                "Arrow record batch length disagrees with its columns",
            ));
        }
        // POD5 writers don't compress Arrow bodies, brust's arrow-ipc build can't
        // decode them, and compressed buffer lengths can't bound row counts.
        if batch.compression().is_some() {
            return Err(invalid_data(
                "compressed Arrow bodies are not supported in POD5",
            ));
        }
        let mut layout = ArrowBatchLayout {
            nodes: &nodes,
            buffers: batch.buffers().into_iter().flatten().collect(),
            next_node: 0,
            next_buffer: 0,
        };
        match dictionary_field {
            Some(field) => {
                layout.check_field(field, true)?;
            }
            None => {
                // The first column's length is the batch's row count, so its
                // buffers must bound that length by the bytes present; later
                // columns of types the checks don't model are left to arrow-ipc.
                let mut fields = schema_fields.iter();
                let Some(first) = fields.next() else {
                    return Err(invalid_data("Arrow section has no columns"));
                };
                if !column_bounds_its_length(&first)? || !layout.check_field(first, false)? {
                    return Err(invalid_data(
                        "first Arrow column cannot bound the batch's row count",
                    ));
                }
                for field in fields {
                    if !layout.check_field(field, false)? {
                        break;
                    }
                }
            }
        }
        if is_record_batch {
            total = total
                .checked_add(rows)
                .ok_or_else(|| invalid_data("Arrow record batch row count overflow"))?;
            batch_rows.push(rows);
        }
    }

    Ok(batch_rows)
}

fn find_dictionary_field<'a>(
    fields: impl Iterator<Item = arrow_ipc::Field<'a>>,
    id: i64,
) -> Option<arrow_ipc::Field<'a>> {
    for field in fields {
        if field
            .dictionary()
            .is_some_and(|encoding| encoding.id() == id)
        {
            return Some(field);
        }
        if let Some(found) = field
            .children()
            .and_then(|children| find_dictionary_field(children.iter(), id))
        {
            return Some(found);
        }
    }
    None
}

/// Walks an Arrow batch's field nodes and buffers in the order arrow-ipc
/// consumes them, checking each buffer is long enough for its node's length.
///
/// arrow-ipc builds arrays from these buffers before validating them and
/// panics when one is too short, so they are checked first. The checks also
/// bound every row count by the bytes actually present.
struct ArrowBatchLayout<'a> {
    nodes: &'a [&'a arrow_ipc::FieldNode],
    buffers: Vec<&'a arrow_ipc::Buffer>,
    next_node: usize,
    next_buffer: usize,
}

/// How a field's buffers are laid out, as far as these checks need.
#[derive(Clone, Copy)]
enum ArrowFieldLayout {
    /// Validity bitmap, then values of this many bits each.
    Fixed(u64),
    /// Validity bitmap, offsets of this many bytes each, then values.
    Variable(u64),
    /// Validity bitmap and offsets of this many bytes each, then one child.
    List(u64),
    /// Validity bitmap, then children.
    Nested,
    /// No buffers.
    Null,
}

impl ArrowBatchLayout<'_> {
    /// Checks `field` and its children against the next nodes and buffers.
    ///
    /// Returns `false`, checking nothing further, when the field uses a type
    /// these checks don't model, since later fields can't be located after it.
    fn check_field(
        &mut self,
        field: arrow_ipc::Field<'_>,
        dictionary_values: bool,
    ) -> io::Result<bool> {
        let Some(layout) = arrow_field_layout(&field, dictionary_values)? else {
            return Ok(false);
        };

        let (length, null_count) = self.next_node()?;
        if let ArrowFieldLayout::Null = layout {
            return Ok(true);
        }
        let validity = self.next_buffer()?;
        if null_count > 0 && validity < length.div_ceil(8) {
            return Err(invalid_data("Arrow validity buffer is too short"));
        }
        let too_short = || invalid_data("Arrow buffer is too short for its column length");
        match layout {
            ArrowFieldLayout::Fixed(bits) => {
                let needed = length.checked_mul(bits).ok_or_else(too_short)?.div_ceil(8);
                if self.next_buffer()? < needed {
                    return Err(too_short());
                }
            }
            ArrowFieldLayout::Variable(width) | ArrowFieldLayout::List(width) => {
                let offsets = self.next_buffer()?;
                let needed = length
                    .checked_add(1)
                    .and_then(|count| count.checked_mul(width))
                    .ok_or_else(too_short)?;
                if length > 0 && offsets < needed {
                    return Err(too_short());
                }
                if let ArrowFieldLayout::Variable(_) = layout {
                    self.next_buffer()?;
                }
            }
            ArrowFieldLayout::Nested | ArrowFieldLayout::Null => {}
        }
        if let ArrowFieldLayout::List(_) | ArrowFieldLayout::Nested = layout {
            for child in field.children().into_iter().flatten() {
                if !self.check_field(child, false)? {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    fn next_node(&mut self) -> io::Result<(u64, u64)> {
        let node = self
            .nodes
            .get(self.next_node)
            .ok_or_else(|| invalid_data("Arrow batch has fewer field nodes than its schema"))?;
        self.next_node += 1;
        // Both counts were checked to be non-negative before the walk.
        Ok((node.length() as u64, node.null_count() as u64))
    }

    /// Returns the next buffer's length; its extent was checked against the body.
    fn next_buffer(&mut self) -> io::Result<u64> {
        let buffer = self
            .buffers
            .get(self.next_buffer)
            .ok_or_else(|| invalid_data("Arrow batch has fewer buffers than its schema"))?;
        self.next_buffer += 1;
        Ok(buffer.length() as u64)
    }
}

/// Whether checking `field` bounds its length by the bytes present: its
/// layout has a values or offsets buffer that grows with the length.
fn column_bounds_its_length(field: &arrow_ipc::Field<'_>) -> io::Result<bool> {
    Ok(matches!(
        arrow_field_layout(field, false)?,
        Some(
            ArrowFieldLayout::Fixed(1..)
                | ArrowFieldLayout::Variable(_)
                | ArrowFieldLayout::List(_)
        )
    ))
}

/// Returns `field`'s buffer layout, or `None` for a type the checks don't model.
///
/// A dictionary-encoded field is laid out as its index type in record batches
/// and as its value type in its own dictionary batches.
fn arrow_field_layout(
    field: &arrow_ipc::Field<'_>,
    dictionary_values: bool,
) -> io::Result<Option<ArrowFieldLayout>> {
    use arrow_ipc::{DateUnit, IntervalUnit, Precision, Type};

    let bits =
        |width: i32| u64::try_from(width).map_err(|_| invalid_data("Arrow type width is negative"));
    if !dictionary_values && let Some(encoding) = field.dictionary() {
        let index_bits = encoding.indexType().map_or(32, |index| index.bitWidth());
        return Ok(Some(ArrowFieldLayout::Fixed(bits(index_bits)?)));
    }
    let layout = match field.type_type() {
        Type::Null => ArrowFieldLayout::Null,
        Type::Bool => ArrowFieldLayout::Fixed(1),
        Type::Int => {
            ArrowFieldLayout::Fixed(bits(field.type_as_int().map_or(0, |int| int.bitWidth()))?)
        }
        Type::FloatingPoint => ArrowFieldLayout::Fixed(
            match field
                .type_as_floating_point()
                .map(|float| float.precision())
            {
                Some(Precision::HALF) => 16,
                Some(Precision::SINGLE) => 32,
                _ => 64,
            },
        ),
        Type::Decimal => ArrowFieldLayout::Fixed(bits(
            field
                .type_as_decimal()
                .map_or(128, |decimal| decimal.bitWidth()),
        )?),
        Type::Date => ArrowFieldLayout::Fixed(match field.type_as_date().map(|date| date.unit()) {
            Some(DateUnit::DAY) => 32,
            _ => 64,
        }),
        Type::Time => ArrowFieldLayout::Fixed(bits(
            field.type_as_time().map_or(32, |time| time.bitWidth()),
        )?),
        Type::Timestamp | Type::Duration => ArrowFieldLayout::Fixed(64),
        Type::Interval => ArrowFieldLayout::Fixed(
            match field.type_as_interval().map(|interval| interval.unit()) {
                Some(IntervalUnit::YEAR_MONTH) => 32,
                Some(IntervalUnit::DAY_TIME) => 64,
                _ => 128,
            },
        ),
        Type::FixedSizeBinary => ArrowFieldLayout::Fixed(
            bits(
                field
                    .type_as_fixed_size_binary()
                    .map_or(0, |binary| binary.byteWidth()),
            )?
            .checked_mul(8)
            .ok_or_else(|| invalid_data("Arrow type width overflow"))?,
        ),
        Type::Utf8 | Type::Binary => ArrowFieldLayout::Variable(4),
        Type::LargeUtf8 | Type::LargeBinary => ArrowFieldLayout::Variable(8),
        Type::List | Type::Map => ArrowFieldLayout::List(4),
        Type::LargeList => ArrowFieldLayout::List(8),
        Type::FixedSizeList | Type::Struct_ => ArrowFieldLayout::Nested,
        _ => return Ok(None),
    };
    Ok(Some(layout))
}

fn read_arrow_footer<R: Read + Seek>(
    shared: SharedReader<R>,
    offset: u64,
    length: u64,
) -> io::Result<Vec<u8>> {
    if length < 10 {
        return Err(invalid_data("embedded Arrow file is too short"));
    }

    let trailer = read_exact_at(shared.clone(), offset + length - 10, 10)?;
    if trailer[4..] != ARROW_MAGIC[..] {
        return Err(invalid_data(
            "embedded Arrow file is missing trailing magic",
        ));
    }

    let footer_len = i32::from_le_bytes(slice_to_array::<4>(&trailer[..4]));
    if footer_len < 0 {
        return Err(invalid_data("Arrow footer length is negative"));
    }
    let footer_len = footer_len as u64;
    if footer_len > length - 10 {
        return Err(invalid_data("Arrow footer length exceeds section"));
    }

    read_exact_at(
        shared,
        offset + length - 10 - footer_len,
        usize::try_from(footer_len)
            .map_err(|_| invalid_data("Arrow footer length exceeds usize"))?,
    )
}

fn parse_arrow_message(data: &[u8]) -> io::Result<arrow_ipc::Message<'_>> {
    if data.len() < 4 {
        return Err(invalid_data("Arrow message is truncated"));
    }

    let message = if data[..4] == [0xff; 4] {
        if data.len() < 8 {
            return Err(invalid_data("Arrow continuation message is truncated"));
        }
        &data[8..]
    } else {
        &data[4..]
    };

    root_as_message(message)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))
}

fn read_exact_at<R: Read + Seek>(
    shared: SharedReader<R>,
    offset: u64,
    len: usize,
) -> io::Result<Vec<u8>> {
    let mut data = vec![0; len];
    let mut reader = shared.borrow_mut();
    reader.seek(SeekFrom::Start(offset))?;
    reader.read_exact(&mut data)?;
    Ok(data)
}

fn checked_padded_len(length: u64) -> io::Result<u64> {
    length
        .checked_add((8 - length % 8) % 8)
        .ok_or_else(|| invalid_data("POD5 padded section length overflow"))
}

fn fb_root_table(data: &[u8]) -> io::Result<usize> {
    let offset = fb_u32_at(data, 0)? as usize;
    if offset >= data.len() {
        return Err(invalid_data("POD5 footer root table is out of bounds"));
    }
    Ok(offset)
}

fn fb_table_vector_field(data: &[u8], table: usize, slot: u16) -> io::Result<Vec<usize>> {
    let Some(field) = fb_field_position(data, table, slot)? else {
        return Ok(Vec::new());
    };
    let vector = fb_uoffset_target(data, field)?;
    let len = fb_u32_at(data, vector)? as usize;
    let elements = vector
        .checked_add(4)
        .ok_or_else(|| invalid_data("FlatBuffer vector offset overflow"))?;
    if len > data.len().saturating_sub(elements) / 4 {
        return Err(invalid_data("FlatBuffer vector is longer than its data"));
    }
    let mut tables = Vec::with_capacity(len);

    for index in 0..len {
        let element = elements
            .checked_add(
                index
                    .checked_mul(4)
                    .ok_or_else(|| invalid_data("FlatBuffer vector offset overflow"))?,
            )
            .ok_or_else(|| invalid_data("FlatBuffer vector offset overflow"))?;
        tables.push(fb_uoffset_target(data, element)?);
    }

    Ok(tables)
}

fn fb_string_field(data: &[u8], table: usize, slot: u16) -> io::Result<Option<String>> {
    let Some(field) = fb_field_position(data, table, slot)? else {
        return Ok(None);
    };
    let string = fb_uoffset_target(data, field)?;
    let len = fb_u32_at(data, string)? as usize;
    let start = string
        .checked_add(4)
        .ok_or_else(|| invalid_data("FlatBuffer string offset overflow"))?;
    let end = start
        .checked_add(len)
        .ok_or_else(|| invalid_data("FlatBuffer string length overflow"))?;
    let bytes = data
        .get(start..end)
        .ok_or_else(|| invalid_data("FlatBuffer string is out of bounds"))?;

    String::from_utf8(bytes.to_vec())
        .map(Some)
        .map_err(|_| invalid_data("FlatBuffer string is not valid UTF-8"))
}

fn fb_i64_field(data: &[u8], table: usize, slot: u16) -> io::Result<Option<i64>> {
    let Some(position) = fb_field_position(data, table, slot)? else {
        return Ok(None);
    };
    Ok(Some(i64::from_le_bytes(slice_to_array::<8>(
        data.get(position..position + 8)
            .ok_or_else(|| invalid_data("FlatBuffer i64 is out of bounds"))?,
    ))))
}

fn fb_i16_field(data: &[u8], table: usize, slot: u16) -> io::Result<Option<i16>> {
    let Some(position) = fb_field_position(data, table, slot)? else {
        return Ok(None);
    };
    Ok(Some(i16::from_le_bytes(slice_to_array::<2>(
        data.get(position..position + 2)
            .ok_or_else(|| invalid_data("FlatBuffer i16 is out of bounds"))?,
    ))))
}

fn fb_field_position(data: &[u8], table: usize, slot: u16) -> io::Result<Option<usize>> {
    let vtable = fb_vtable_position(data, table)?;
    let vtable_len = fb_u16_at(data, vtable)?;
    if slot.saturating_add(2) > vtable_len {
        return Ok(None);
    }
    let field_offset = fb_u16_at(data, vtable + slot as usize)?;
    if field_offset == 0 {
        return Ok(None);
    }
    let position = table
        .checked_add(field_offset as usize)
        .ok_or_else(|| invalid_data("FlatBuffer field offset overflow"))?;
    if position >= data.len() {
        return Err(invalid_data("FlatBuffer field is out of bounds"));
    }
    Ok(Some(position))
}

fn fb_vtable_position(data: &[u8], table: usize) -> io::Result<usize> {
    let offset = fb_i32_at(data, table)? as isize;
    let vtable = table as isize - offset;
    if vtable < 0 {
        return Err(invalid_data("FlatBuffer vtable is out of bounds"));
    }
    let vtable = vtable as usize;
    if vtable >= data.len() {
        return Err(invalid_data("FlatBuffer vtable is out of bounds"));
    }
    Ok(vtable)
}

fn fb_uoffset_target(data: &[u8], position: usize) -> io::Result<usize> {
    let offset = fb_u32_at(data, position)? as usize;
    let target = position
        .checked_add(offset)
        .ok_or_else(|| invalid_data("FlatBuffer offset overflow"))?;
    if target >= data.len() {
        return Err(invalid_data("FlatBuffer offset is out of bounds"));
    }
    Ok(target)
}

fn fb_u16_at(data: &[u8], position: usize) -> io::Result<u16> {
    Ok(u16::from_le_bytes(slice_to_array::<2>(
        data.get(position..position + 2)
            .ok_or_else(|| invalid_data("FlatBuffer u16 is out of bounds"))?,
    )))
}

fn fb_u32_at(data: &[u8], position: usize) -> io::Result<u32> {
    Ok(u32::from_le_bytes(slice_to_array::<4>(
        data.get(position..position + 4)
            .ok_or_else(|| invalid_data("FlatBuffer u32 is out of bounds"))?,
    )))
}

fn fb_i32_at(data: &[u8], position: usize) -> io::Result<i32> {
    Ok(i32::from_le_bytes(slice_to_array::<4>(
        data.get(position..position + 4)
            .ok_or_else(|| invalid_data("FlatBuffer i32 is out of bounds"))?,
    )))
}

/// Signal table rows per Arrow record batch, as official POD5 writers use.
///
/// Official readers find a row's batch by dividing its number by the first
/// batch's row count, so every batch but the last has this many rows.
const SIGNAL_BATCH_ROWS: usize = 100;

/// A payload that passed every check on its contents, with its Run Info and
/// Reads tables encoded. Writing it still compresses uncompressed signal rows
/// and encodes the Signal table, a batch at a time.
struct PreparedPod5<'a> {
    signals: &'a [Pod5Signal],
    metadata: Pod5WriterMetadata,
    marker: [u8; 16],
    run_info: Vec<u8>,
    reads: Vec<u8>,
}

fn prepare_pod5(pod5: &Pod5) -> io::Result<PreparedPod5<'_>> {
    validate_pod5_for_writing(pod5)?;

    let metadata = pod5_writer_metadata(pod5);
    let run_info =
        write_arrow_section(build_run_info_batch(&pod5.run_infos, &metadata)?, &metadata)?;
    let reads = write_arrow_section(build_reads_batch(&pod5.records, &metadata)?, &metadata)?;
    let marker = pod5_section_marker(pod5, &metadata)?;

    Ok(PreparedPod5 {
        signals: &pod5.signals,
        metadata,
        marker,
        run_info,
        reads,
    })
}

impl PreparedPod5<'_> {
    /// Writes the POD5 file, the Signal table a batch at a time.
    fn write_to<W: Write>(&self, writer: W) -> io::Result<()> {
        let mut output = CountingWriter {
            inner: BufWriter::new(writer),
            written: 0,
        };
        let mut footer_entries = Vec::with_capacity(3);

        output.write_all(POD5_MAGIC)?;
        output.write_all(&self.marker)?;

        for (kind, content_type) in [
            (Pod5SectionKind::Signal, 1),
            (Pod5SectionKind::RunInfo, 4),
            (Pod5SectionKind::Reads, 0),
        ] {
            let offset = output.written;
            match kind {
                Pod5SectionKind::Signal => {
                    write_signal_section(self.signals, &self.metadata, &mut output)?;
                }
                Pod5SectionKind::RunInfo => output.write_all(&self.run_info)?,
                _ => output.write_all(&self.reads)?,
            }
            let length = output.written - offset;
            output.pad_to_8()?;
            output.write_all(&self.marker)?;
            footer_entries.push(Pod5FooterEntry {
                offset: file_offset(offset)?,
                length: file_offset(length)?,
                content_type,
            });
        }

        let footer = build_pod5_footer(&self.metadata, &footer_entries);
        output.write_all(POD5_FOOTER_MAGIC)?;
        let footer_payload_start = output.written;
        output.write_all(&footer)?;
        output.pad_to_8()?;
        let footer_len = output.written - footer_payload_start;
        output.write_all(&footer_len.to_le_bytes())?;
        output.write_all(&self.marker)?;
        output.write_all(POD5_MAGIC)?;

        output
            .inner
            .into_inner()
            .map(drop)
            .map_err(io::IntoInnerError::into_error)
    }
}

fn file_offset(value: u64) -> io::Result<i64> {
    i64::try_from(value).map_err(|_| invalid_data("POD5 output exceeds the footer's offset range"))
}

/// Counts the bytes written through it, for the POD5 footer's offsets.
struct CountingWriter<W> {
    inner: W,
    written: u64,
}

impl<W: Write> CountingWriter<W> {
    fn pad_to_8(&mut self) -> io::Result<()> {
        let padding = (8 - self.written % 8) % 8;
        self.write_all(&[0; 8][..padding as usize])
    }
}

impl<W: Write> Write for CountingWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let written = self.inner.write(buf)?;
        self.written += written as u64;
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

/// Writes the Signal table as an Arrow file of [`SIGNAL_BATCH_ROWS`]-row
/// batches, building one batch at a time.
fn write_signal_section<W: Write>(
    signals: &[Pod5Signal],
    metadata: &Pod5WriterMetadata,
    output: W,
) -> io::Result<()> {
    let schema = Arc::new(signal_schema(metadata));
    let mut writer = FileWriter::try_new(output, schema.as_ref()).map_err(arrow_write_error)?;
    for (key, value) in arrow_metadata(metadata) {
        writer.write_metadata(key, value);
    }
    // An empty table is still written as one empty batch.
    for rows in signals
        .chunks(SIGNAL_BATCH_ROWS)
        .chain(signals.is_empty().then_some(signals))
    {
        writer
            .write(&build_signal_batch(rows, &schema)?)
            .map_err(arrow_write_error)?;
    }
    writer.finish().map_err(arrow_write_error)
}

fn validate_pod5_for_writing(pod5: &Pod5) -> io::Result<()> {
    // Empty tables are valid: official pod5 writes files with no reads, and
    // reads with no samples have no signal rows.
    if let Some(version) = &pod5.header.pod5_version {
        validate_pod5_version(version)?;
    }
    if let Some(identifier) = &pod5.header.file_identifier
        && !is_canonical_uuid(identifier)
    {
        return Err(invalid_data(format!(
            "MINKNOW:file_identifier must be a UUID, got {identifier}"
        )));
    }

    let mut signal_read_ids = Vec::with_capacity(pod5.signals.len());
    for signal in &pod5.signals {
        signal_read_ids.push(uuid_string_to_bytes(&signal.read_id)?);
        validate_signal_payload(signal)?;
    }

    for record in &pod5.records {
        let read_id = uuid_string_to_bytes(&record.read_id)?;
        if !POD5_END_REASONS.contains(&record.end_reason.as_str()) {
            return Err(invalid_data(format!(
                "POD5 end_reason {:?} is not a known end reason",
                record.end_reason
            )));
        }
        let mut total = 0u64;
        for &index in &record.signal_rows {
            let signal = pod5
                .signals
                .get(index as usize)
                .ok_or_else(|| invalid_data("POD5 signal row index is out of bounds"))?;
            if signal_read_ids[index as usize] != read_id {
                return Err(invalid_data("POD5 signal row read_id does not match read"));
            }
            total = total
                .checked_add(u64::from(signal.samples))
                .ok_or_else(|| invalid_data("POD5 signal sample count overflow"))?;
        }
        if total != record.num_samples {
            return Err(invalid_data(
                "POD5 read num_samples does not match referenced signal rows",
            ));
        }
        if !pod5
            .run_infos
            .iter()
            .any(|run_info| run_info.acquisition_id == record.run_info)
        {
            return Err(invalid_data("POD5 read references missing run info row"));
        }
    }

    Ok(())
}

/// End reasons official POD5 readers accept in the Reads `end_reason` column.
const POD5_END_REASONS: &[&str] = &[
    "unknown",
    "mux_change",
    "unblock_mux_change",
    "data_service_unblock_mux_change",
    "signal_positive",
    "signal_negative",
    "api_request",
    "device_data_error",
    "analysis_config_change",
    "paused",
];

/// Checks that a signal row's payload decodes to exactly its declared samples.
fn validate_signal_payload(signal: &Pod5Signal) -> io::Result<()> {
    match &signal.payload {
        Pod5SignalPayload::Uncompressed(samples) if samples.len() != signal.samples as usize => {
            Err(invalid_data(
                "POD5 signal row samples does not match its payload length",
            ))
        }
        Pod5SignalPayload::Uncompressed(_) => Ok(()),
        Pod5SignalPayload::Vbz(data) => {
            decompress_vbz_signal(data, signal.samples as usize).map(drop)
        }
    }
}

/// Whether `value` is a hyphenated `8-4-4-4-12` hex UUID.
fn is_canonical_uuid(value: &str) -> bool {
    value.len() == 36
        && value.bytes().enumerate().all(|(index, byte)| match index {
            8 | 13 | 18 | 23 => byte == b'-',
            _ => byte.is_ascii_hexdigit(),
        })
}

fn pod5_writer_metadata(pod5: &Pod5) -> Pod5WriterMetadata {
    Pod5WriterMetadata {
        file_identifier: pod5
            .header
            .file_identifier
            .clone()
            .unwrap_or_else(|| "00000000-0000-0000-0000-000000000000".to_string()),
        software: pod5
            .header
            .software
            .clone()
            .unwrap_or_else(|| "brust pod5 writer".to_string()),
        pod5_version: pod5
            .header
            .pod5_version
            .clone()
            .unwrap_or_else(|| "0.3.34".to_string()),
    }
}

#[derive(Debug)]
struct Pod5WriterMetadata {
    file_identifier: String,
    software: String,
    pod5_version: String,
}

#[derive(Debug)]
struct Pod5FooterEntry {
    offset: i64,
    length: i64,
    content_type: i16,
}

fn pod5_section_marker(pod5: &Pod5, metadata: &Pod5WriterMetadata) -> io::Result<[u8; 16]> {
    if pod5.header.section_marker != [0; 16] {
        return Ok(pod5.header.section_marker);
    }
    let marker = *b"BRUSTPOD5WRITER!";
    if marker == uuid_string_to_bytes(&metadata.file_identifier)? {
        return Err(invalid_data(
            "POD5 section marker collides with file identifier",
        ));
    }
    Ok(marker)
}

fn write_arrow_section(batch: RecordBatch, metadata: &Pod5WriterMetadata) -> io::Result<Vec<u8>> {
    let mut data = Vec::new();
    let schema = batch.schema();
    let mut writer = FileWriter::try_new(&mut data, schema.as_ref()).map_err(arrow_error)?;
    for (key, value) in arrow_metadata(metadata) {
        writer.write_metadata(key, value);
    }
    writer.write(&batch).map_err(arrow_error)?;
    writer.finish().map_err(arrow_error)?;
    drop(writer);
    Ok(data)
}

fn signal_schema(metadata: &Pod5WriterMetadata) -> Schema {
    Schema::new_with_metadata(
        vec![
            uuid_field("read_id"),
            vbz_field("signal"),
            Field::new("samples", DataType::UInt32, false),
        ],
        arrow_metadata(metadata),
    )
}

fn build_signal_batch(signals: &[Pod5Signal], schema: &SchemaRef) -> io::Result<RecordBatch> {
    let mut read_id = FixedSizeBinaryBuilder::with_capacity(signals.len(), 16);
    let total_signal_bytes = signals
        .iter()
        .filter_map(Pod5Signal::compressed_bytes)
        .map(<[u8]>::len)
        .sum::<usize>();
    let mut signal = LargeBinaryBuilder::with_capacity(signals.len(), total_signal_bytes);
    let mut samples = Vec::with_capacity(signals.len());

    for row in signals {
        read_id
            .append_value(uuid_string_to_bytes(&row.read_id)?)
            .map_err(arrow_error)?;
        let payload = match &row.payload {
            Pod5SignalPayload::Vbz(data) => data.clone(),
            Pod5SignalPayload::Uncompressed(samples) => compress_vbz_signal(samples)?,
        };
        signal.append_value(payload);
        samples.push(row.samples);
    }

    let read_id = Arc::new(read_id.finish()) as ArrayRef;
    let signal = Arc::new(signal.finish()) as ArrayRef;
    let samples = Arc::new(UInt32Array::from(samples)) as ArrayRef;

    RecordBatch::try_new(schema.clone(), vec![read_id, signal, samples]).map_err(arrow_error)
}

fn build_run_info_batch(
    run_infos: &[Pod5RunInfo],
    metadata: &Pod5WriterMetadata,
) -> io::Result<RecordBatch> {
    let strings = |field: fn(&Pod5RunInfo) -> &str| -> ArrayRef {
        Arc::new(StringArray::from(
            run_infos.iter().map(field).collect::<Vec<_>>(),
        ))
    };
    let timestamps = |field: fn(&Pod5RunInfo) -> i64| -> ArrayRef {
        Arc::new(
            TimestampMillisecondArray::from(run_infos.iter().map(field).collect::<Vec<_>>())
                .with_timezone("UTC"),
        )
    };
    let int16s = |field: fn(&Pod5RunInfo) -> i16| -> ArrayRef {
        Arc::new(Int16Array::from(
            run_infos.iter().map(field).collect::<Vec<_>>(),
        ))
    };

    let arrays: Vec<ArrayRef> = vec![
        strings(|run_info| &run_info.acquisition_id),
        timestamps(|run_info| run_info.acquisition_start_time),
        int16s(|run_info| run_info.adc_max),
        int16s(|run_info| run_info.adc_min),
        string_map_array(run_infos.iter().map(|run_info| &run_info.context_tags))?,
        strings(|run_info| &run_info.experiment_name),
        strings(|run_info| &run_info.flow_cell_id),
        strings(|run_info| &run_info.flow_cell_product_code),
        strings(|run_info| &run_info.protocol_name),
        strings(|run_info| &run_info.protocol_run_id),
        timestamps(|run_info| run_info.protocol_start_time),
        strings(|run_info| &run_info.sample_id),
        Arc::new(UInt16Array::from(
            run_infos
                .iter()
                .map(|run_info| run_info.sample_rate)
                .collect::<Vec<_>>(),
        )),
        strings(|run_info| &run_info.sequencing_kit),
        strings(|run_info| &run_info.sequencer_position),
        strings(|run_info| &run_info.sequencer_position_type),
        strings(|run_info| &run_info.software),
        strings(|run_info| &run_info.system_name),
        strings(|run_info| &run_info.system_type),
        string_map_array(run_infos.iter().map(|run_info| &run_info.tracking_id))?,
    ];

    let fields = RUN_INFO_TABLE_FIELDS
        .iter()
        .zip(arrays.iter())
        .map(|(name, array)| Field::new(*name, array.data_type().clone(), false))
        .collect::<Vec<_>>();
    let schema = Schema::new_with_metadata(fields, arrow_metadata(metadata));

    RecordBatch::try_new(Arc::new(schema), arrays).map_err(arrow_error)
}

fn build_reads_batch(
    records: &[Pod5Record],
    metadata: &Pod5WriterMetadata,
) -> io::Result<RecordBatch> {
    let mut read_id = FixedSizeBinaryBuilder::with_capacity(records.len(), 16);
    let mut signal = ListBuilder::new(UInt64Builder::new());
    let mut pore_type = StringDictionaryBuilder::<Int16Type>::new();
    let mut end_reason = StringDictionaryBuilder::<Int16Type>::new();
    let mut run_info = StringDictionaryBuilder::<Int16Type>::new();

    for record in records {
        read_id
            .append_value(uuid_string_to_bytes(&record.read_id)?)
            .map_err(arrow_error)?;
        for &row in &record.signal_rows {
            signal.values().append_value(row);
        }
        signal.append(true);
        // Dictionary keys are 16-bit, as in official POD5 files.
        let key_overflow =
            |_| invalid_data("POD5 column has more distinct values than 16-bit keys can index");
        pore_type.append(&record.pore_type).map_err(key_overflow)?;
        end_reason
            .append(&record.end_reason)
            .map_err(key_overflow)?;
        run_info.append(&record.run_info).map_err(key_overflow)?;
    }

    let mut arrays: Vec<ArrayRef> = vec![
        Arc::new(read_id.finish()),
        Arc::new(signal.finish()),
        Arc::new(UInt32Array::from(
            records
                .iter()
                .map(|record| record.read_number)
                .collect::<Vec<_>>(),
        )),
        Arc::new(UInt64Array::from(
            records
                .iter()
                .map(|record| record.start_sample)
                .collect::<Vec<_>>(),
        )),
        Arc::new(Float32Array::from(
            records
                .iter()
                .map(|record| record.median_before)
                .collect::<Vec<_>>(),
        )),
        Arc::new(UInt64Array::from(
            records
                .iter()
                .map(|record| record.num_minknow_events)
                .collect::<Vec<_>>(),
        )),
        Arc::new(Float32Array::from(
            records
                .iter()
                .map(|record| record.tracked_scaling_scale)
                .collect::<Vec<_>>(),
        )),
        Arc::new(Float32Array::from(
            records
                .iter()
                .map(|record| record.tracked_scaling_shift)
                .collect::<Vec<_>>(),
        )),
        Arc::new(Float32Array::from(
            records
                .iter()
                .map(|record| record.predicted_scaling_scale)
                .collect::<Vec<_>>(),
        )),
        Arc::new(Float32Array::from(
            records
                .iter()
                .map(|record| record.predicted_scaling_shift)
                .collect::<Vec<_>>(),
        )),
        Arc::new(UInt32Array::from(
            records
                .iter()
                .map(|record| record.num_reads_since_mux_change)
                .collect::<Vec<_>>(),
        )),
        Arc::new(Float32Array::from(
            records
                .iter()
                .map(|record| record.time_since_mux_change)
                .collect::<Vec<_>>(),
        )),
        Arc::new(UInt64Array::from(
            records
                .iter()
                .map(|record| record.num_samples)
                .collect::<Vec<_>>(),
        )),
        Arc::new(UInt16Array::from(
            records
                .iter()
                .map(|record| record.channel)
                .collect::<Vec<_>>(),
        )),
        Arc::new(UInt8Array::from(
            records.iter().map(|record| record.well).collect::<Vec<_>>(),
        )),
        Arc::new(pore_type.finish()),
        Arc::new(Float32Array::from(
            records
                .iter()
                .map(|record| record.calibration_offset)
                .collect::<Vec<_>>(),
        )),
        Arc::new(Float32Array::from(
            records
                .iter()
                .map(|record| record.calibration_scale)
                .collect::<Vec<_>>(),
        )),
        Arc::new(end_reason.finish()),
        Arc::new(BooleanArray::from(
            records
                .iter()
                .map(|record| record.end_reason_forced)
                .collect::<Vec<_>>(),
        )),
        Arc::new(run_info.finish()),
    ];
    let mut names = READS_TABLE_FIELDS.to_vec();
    if writes_open_pore_level(&metadata.pod5_version) {
        arrays.push(Arc::new(Float32Array::from(
            records
                .iter()
                .map(|record| record.open_pore_level.unwrap_or(f32::NAN))
                .collect::<Vec<_>>(),
        )));
        names.push(OPEN_PORE_LEVEL_FIELD);
    } else if records
        .iter()
        .any(|record| record.open_pore_level.is_some())
    {
        return Err(invalid_data(format!(
            "POD5 open_pore_level needs MINKNOW:pod5_version {OPEN_PORE_LEVEL_MIN_VERSION} or later, got {}",
            metadata.pod5_version
        )));
    }

    let fields = names
        .iter()
        .zip(arrays.iter())
        .map(|(name, array)| {
            if *name == "read_id" {
                uuid_field("read_id")
            } else {
                Field::new(*name, array.data_type().clone(), false)
            }
        })
        .collect::<Vec<_>>();
    let schema = Schema::new_with_metadata(fields, arrow_metadata(metadata));

    RecordBatch::try_new(Arc::new(schema), arrays).map_err(arrow_error)
}

/// Reads column added in POD5 0.3.30; official readers require it from then on.
const OPEN_PORE_LEVEL_FIELD: &str = "open_pore_level";
/// First POD5 version whose Reads table has [`OPEN_PORE_LEVEL_FIELD`].
const OPEN_PORE_LEVEL_MIN_VERSION: &str = "0.3.30";

/// Whether a Reads table stamped with `version` must carry `open_pore_level`.
///
/// A version that doesn't parse as `major.minor.patch` gets the newest schema.
fn writes_open_pore_level(version: &str) -> bool {
    match (
        pod5_version_number(version),
        pod5_version_number(OPEN_PORE_LEVEL_MIN_VERSION),
    ) {
        (Some(version), Some(minimum)) => version >= minimum,
        _ => true,
    }
}

fn pod5_version_number(version: &str) -> Option<(u64, u64, u64)> {
    let mut parts = version.split('.').map(|part| {
        if part.is_empty() || !part.bytes().all(|byte| byte.is_ascii_digit()) {
            None
        } else {
            part.parse::<u64>().ok()
        }
    });
    let number = (parts.next()??, parts.next()??, parts.next()??);
    parts.next().is_none().then_some(number)
}

fn string_map_array<'a>(
    rows: impl Iterator<Item = &'a Vec<(String, String)>>,
) -> io::Result<ArrayRef> {
    let mut builder = MapBuilder::new(None, StringBuilder::new(), StringBuilder::new());
    for entries in rows {
        for (key, value) in entries {
            builder.keys().append_value(key);
            builder.values().append_value(value);
        }
        builder.append(true).map_err(arrow_error)?;
    }
    Ok(Arc::new(builder.finish()))
}

fn uuid_field(name: &str) -> Field {
    Field::new(name, DataType::FixedSizeBinary(16), false)
        .with_metadata(extension_metadata("minknow.uuid"))
}

fn vbz_field(name: &str) -> Field {
    Field::new(name, DataType::LargeBinary, false).with_metadata(extension_metadata("minknow.vbz"))
}

fn extension_metadata(name: &str) -> HashMap<String, String> {
    HashMap::from([
        ("ARROW:extension:name".to_string(), name.to_string()),
        ("ARROW:extension:metadata".to_string(), String::new()),
    ])
}

fn arrow_metadata(metadata: &Pod5WriterMetadata) -> HashMap<String, String> {
    HashMap::from([
        (
            "MINKNOW:file_identifier".to_string(),
            metadata.file_identifier.clone(),
        ),
        ("MINKNOW:software".to_string(), metadata.software.clone()),
        (
            "MINKNOW:pod5_version".to_string(),
            metadata.pod5_version.clone(),
        ),
    ])
}

fn build_pod5_footer(metadata: &Pod5WriterMetadata, entries: &[Pod5FooterEntry]) -> Vec<u8> {
    let mut builder = flatbuffers::FlatBufferBuilder::new();
    let file_identifier = builder.create_string(&metadata.file_identifier);
    let software = builder.create_string(&metadata.software);
    let pod5_version = builder.create_string(&metadata.pod5_version);
    let mut embedded = Vec::with_capacity(entries.len());

    for entry in entries {
        let start = builder.start_table();
        builder.push_slot_always::<i64>(4, entry.offset);
        builder.push_slot_always::<i64>(6, entry.length);
        builder.push_slot_always::<i16>(8, 0);
        builder.push_slot_always::<i16>(10, entry.content_type);
        embedded.push(builder.end_table(start));
    }

    let contents = builder.create_vector(&embedded);
    let start = builder.start_table();
    builder.push_slot_always(4, file_identifier);
    builder.push_slot_always(6, software);
    builder.push_slot_always(8, pod5_version);
    builder.push_slot_always(10, contents);
    let footer = builder.end_table(start);
    builder.finish(footer, None);
    builder.finished_data().to_vec()
}

fn infer_section_kind(schema: &SchemaRef) -> Pod5SectionKind {
    if has_fields(schema, READS_TABLE_FIELDS) {
        Pod5SectionKind::Reads
    } else if has_fields(schema, SIGNAL_TABLE_FIELDS)
        && !schema
            .fields()
            .iter()
            .any(|field| field.name() == "channel")
    {
        Pod5SectionKind::Signal
    } else if has_fields(schema, RUN_INFO_TABLE_FIELDS) {
        Pod5SectionKind::RunInfo
    } else {
        Pod5SectionKind::Unknown
    }
}

fn has_fields(schema: &SchemaRef, expected: &[&str]) -> bool {
    expected
        .iter()
        .all(|name| schema.fields().iter().any(|field| field.name() == *name))
}

fn parse_signal_batch(batch: &RecordBatch) -> io::Result<Vec<Pod5Signal>> {
    let read_id = fixed_binary_column(batch, "read_id")?;
    let samples = uint32_column(batch, "samples")?;
    let signal = column(batch, "signal")?;

    if let Some(signal) = signal.as_any().downcast_ref::<LargeBinaryArray>() {
        let mut rows = Vec::with_capacity(batch.num_rows());
        for row in 0..batch.num_rows() {
            rows.push(Pod5Signal {
                read_id: uuid_bytes_to_string(read_id.value(row))?,
                samples: samples.value(row),
                payload: Pod5SignalPayload::Vbz(signal.value(row).to_vec()),
            });
        }
        return Ok(rows);
    }

    if let Some(signal) = signal.as_any().downcast_ref::<LargeListArray>() {
        let mut rows = Vec::with_capacity(batch.num_rows());
        for row in 0..batch.num_rows() {
            let value = signal.value(row);
            let value = value
                .as_any()
                .downcast_ref::<Int16Array>()
                .ok_or_else(|| invalid_data("POD5 uncompressed signal has unexpected type"))?;
            let decoded = (0..value.len())
                .map(|index| value.value(index))
                .collect::<Vec<_>>();
            if decoded.len() != samples.value(row) as usize {
                return Err(invalid_data(
                    "POD5 uncompressed signal length does not match samples column",
                ));
            }

            rows.push(Pod5Signal {
                read_id: uuid_bytes_to_string(read_id.value(row))?,
                samples: samples.value(row),
                payload: Pod5SignalPayload::Uncompressed(decoded),
            });
        }
        return Ok(rows);
    }

    Err(invalid_data(
        "POD5 Signal table has unsupported signal column type",
    ))
}

fn parse_reads_batch(batch: &RecordBatch) -> io::Result<Vec<Pod5Record>> {
    let read_id = fixed_binary_column(batch, "read_id")?;
    let signal = list_column(batch, "signal")?;
    let read_number = uint32_column(batch, "read_number")?;
    let start = uint64_column(batch, "start")?;
    let median_before = float32_column(batch, "median_before")?;
    let num_minknow_events = uint64_column(batch, "num_minknow_events")?;
    let tracked_scaling_scale = float32_column(batch, "tracked_scaling_scale")?;
    let tracked_scaling_shift = float32_column(batch, "tracked_scaling_shift")?;
    let predicted_scaling_scale = float32_column(batch, "predicted_scaling_scale")?;
    let predicted_scaling_shift = float32_column(batch, "predicted_scaling_shift")?;
    let num_reads_since_mux_change = uint32_column(batch, "num_reads_since_mux_change")?;
    let time_since_mux_change = float32_column(batch, "time_since_mux_change")?;
    let num_samples = uint64_column(batch, "num_samples")?;
    let channel = uint16_column(batch, "channel")?;
    let well = uint8_column(batch, "well")?;
    let pore_type = dict_string_column(batch, "pore_type")?;
    let calibration_offset = float32_column(batch, "calibration_offset")?;
    let calibration_scale = float32_column(batch, "calibration_scale")?;
    let end_reason = dict_string_column(batch, "end_reason")?;
    let end_reason_forced = bool_column(batch, "end_reason_forced")?;
    let run_info = dict_string_column(batch, "run_info")?;
    let open_pore_level = match batch.column_by_name(OPEN_PORE_LEVEL_FIELD) {
        Some(column) => Some(
            column
                .as_any()
                .downcast_ref::<Float32Array>()
                .ok_or_else(|| invalid_data("POD5 column has unexpected type"))?,
        ),
        None => None,
    };

    let mut records = Vec::with_capacity(batch.num_rows());
    for row in 0..batch.num_rows() {
        let id = uuid_bytes_to_string(read_id.value(row))?;
        let signal_rows = uint64_list_value(signal, row)?;

        records.push(Pod5Record {
            read_id: id,
            signal_rows,
            read_number: read_number.value(row),
            start_sample: start.value(row),
            median_before: median_before.value(row),
            num_minknow_events: num_minknow_events.value(row),
            tracked_scaling_scale: tracked_scaling_scale.value(row),
            tracked_scaling_shift: tracked_scaling_shift.value(row),
            predicted_scaling_scale: predicted_scaling_scale.value(row),
            predicted_scaling_shift: predicted_scaling_shift.value(row),
            num_reads_since_mux_change: num_reads_since_mux_change.value(row),
            time_since_mux_change: time_since_mux_change.value(row),
            num_samples: num_samples.value(row),
            channel: channel.value(row),
            well: well.value(row),
            pore_type: dictionary_string_value(pore_type, row)?,
            calibration_offset: calibration_offset.value(row),
            calibration_scale: calibration_scale.value(row),
            end_reason: dictionary_string_value(end_reason, row)?,
            end_reason_forced: end_reason_forced.value(row),
            run_info: dictionary_string_value(run_info, row)?,
            open_pore_level: open_pore_level
                .filter(|column| column.is_valid(row))
                .map(|column| column.value(row))
                .filter(|level| !level.is_nan()),
        });
    }

    Ok(records)
}

fn parse_run_info_batch(batch: &RecordBatch) -> io::Result<Vec<Pod5RunInfo>> {
    let acquisition_id = string_column(batch, "acquisition_id")?;
    let sample_id = string_column(batch, "sample_id")?;
    let experiment_name = string_column(batch, "experiment_name")?;
    let flow_cell_id = string_column(batch, "flow_cell_id")?;
    let sequencing_kit = string_column(batch, "sequencing_kit")?;
    let sample_rate = uint16_column(batch, "sample_rate")?;
    let software = string_column(batch, "software")?;
    let acquisition_start_time = timestamp_ms_column(batch, "acquisition_start_time")?;
    let adc_max = int16_column(batch, "adc_max")?;
    let adc_min = int16_column(batch, "adc_min")?;
    let context_tags = map_column(batch, "context_tags")?;
    let flow_cell_product_code = string_column(batch, "flow_cell_product_code")?;
    let protocol_name = string_column(batch, "protocol_name")?;
    let protocol_run_id = string_column(batch, "protocol_run_id")?;
    let protocol_start_time = timestamp_ms_column(batch, "protocol_start_time")?;
    let sequencer_position = string_column(batch, "sequencer_position")?;
    let sequencer_position_type = string_column(batch, "sequencer_position_type")?;
    let system_name = string_column(batch, "system_name")?;
    let system_type = string_column(batch, "system_type")?;
    let tracking_id = map_column(batch, "tracking_id")?;

    let mut run_infos = Vec::with_capacity(batch.num_rows());
    for row in 0..batch.num_rows() {
        run_infos.push(Pod5RunInfo {
            acquisition_id: acquisition_id.value(row).to_string(),
            sample_id: sample_id.value(row).to_string(),
            experiment_name: experiment_name.value(row).to_string(),
            flow_cell_id: flow_cell_id.value(row).to_string(),
            sequencing_kit: sequencing_kit.value(row).to_string(),
            sample_rate: sample_rate.value(row),
            software: software.value(row).to_string(),
            acquisition_start_time: acquisition_start_time.value(row),
            adc_max: adc_max.value(row),
            adc_min: adc_min.value(row),
            context_tags: string_map_value(context_tags, row)?,
            flow_cell_product_code: flow_cell_product_code.value(row).to_string(),
            protocol_name: protocol_name.value(row).to_string(),
            protocol_run_id: protocol_run_id.value(row).to_string(),
            protocol_start_time: protocol_start_time.value(row),
            sequencer_position: sequencer_position.value(row).to_string(),
            sequencer_position_type: sequencer_position_type.value(row).to_string(),
            system_name: system_name.value(row).to_string(),
            system_type: system_type.value(row).to_string(),
            tracking_id: string_map_value(tracking_id, row)?,
        });
    }

    Ok(run_infos)
}

#[derive(Default)]
struct Pod5Metadata {
    file_identifier: Option<String>,
    software: Option<String>,
    pod5_version: Option<String>,
}

fn update_metadata(
    arrow_metadata: &std::collections::HashMap<String, String>,
    metadata: &mut Pod5Metadata,
) -> io::Result<()> {
    if let Some(version) = arrow_metadata.get("MINKNOW:pod5_version") {
        validate_pod5_version(version)?;
    }

    merge_metadata(
        &mut metadata.file_identifier,
        arrow_metadata.get("MINKNOW:file_identifier"),
        "MINKNOW:file_identifier",
    )?;
    merge_metadata(
        &mut metadata.software,
        arrow_metadata.get("MINKNOW:software"),
        "MINKNOW:software",
    )?;
    merge_metadata(
        &mut metadata.pod5_version,
        arrow_metadata.get("MINKNOW:pod5_version"),
        "MINKNOW:pod5_version",
    )?;

    Ok(())
}

fn validate_pod5_version(version: &str) -> io::Result<()> {
    let parts = version.split('.').collect::<Vec<_>>();
    if parts.len() != 3
        || parts
            .iter()
            .any(|part| part.is_empty() || !part.bytes().all(|byte| byte.is_ascii_digit()))
    {
        return Err(invalid_data(format!(
            "MINKNOW:pod5_version must be semantic major.minor.patch, got {version}"
        )));
    }

    Ok(())
}

fn merge_metadata(
    existing: &mut Option<String>,
    value: Option<&String>,
    key: &'static str,
) -> io::Result<()> {
    let Some(value) = value else {
        return Ok(());
    };

    match existing {
        Some(existing) if existing != value => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("inconsistent POD5 Arrow metadata for {key}"),
        )),
        Some(_) => Ok(()),
        None => {
            *existing = Some(value.clone());
            Ok(())
        }
    }
}

fn fixed_binary_column<'a>(
    batch: &'a RecordBatch,
    name: &str,
) -> io::Result<&'a FixedSizeBinaryArray> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<FixedSizeBinaryArray>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

fn list_column<'a>(batch: &'a RecordBatch, name: &str) -> io::Result<&'a ListArray> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<ListArray>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

fn string_column<'a>(batch: &'a RecordBatch, name: &str) -> io::Result<&'a StringArray> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

fn int16_column<'a>(batch: &'a RecordBatch, name: &str) -> io::Result<&'a Int16Array> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<Int16Array>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

fn timestamp_ms_column<'a>(
    batch: &'a RecordBatch,
    name: &str,
) -> io::Result<&'a TimestampMillisecondArray> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<TimestampMillisecondArray>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

fn map_column<'a>(batch: &'a RecordBatch, name: &str) -> io::Result<&'a MapArray> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<MapArray>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

/// Returns one row of a `map<string, string>` column as key/value pairs.
fn string_map_value(map: &MapArray, row: usize) -> io::Result<Vec<(String, String)>> {
    if map.is_null(row) {
        return Ok(Vec::new());
    }
    let keys = map
        .keys()
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| invalid_data("POD5 map keys have unexpected type"))?;
    let values = map
        .values()
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| invalid_data("POD5 map values have unexpected type"))?;
    let offsets = map.value_offsets();
    let (start, end) = (offsets[row] as usize, offsets[row + 1] as usize);
    Ok((start..end)
        .map(|entry| {
            (
                keys.value(entry).to_string(),
                values.value(entry).to_string(),
            )
        })
        .collect())
}

fn uint8_column<'a>(batch: &'a RecordBatch, name: &str) -> io::Result<&'a UInt8Array> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<UInt8Array>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

fn uint16_column<'a>(batch: &'a RecordBatch, name: &str) -> io::Result<&'a UInt16Array> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<UInt16Array>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

fn uint32_column<'a>(batch: &'a RecordBatch, name: &str) -> io::Result<&'a UInt32Array> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<UInt32Array>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

fn uint64_column<'a>(batch: &'a RecordBatch, name: &str) -> io::Result<&'a UInt64Array> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<UInt64Array>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

fn float32_column<'a>(batch: &'a RecordBatch, name: &str) -> io::Result<&'a Float32Array> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<Float32Array>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

fn bool_column<'a>(batch: &'a RecordBatch, name: &str) -> io::Result<&'a BooleanArray> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<BooleanArray>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

fn dict_string_column<'a>(
    batch: &'a RecordBatch,
    name: &str,
) -> io::Result<&'a DictionaryArray<Int16Type>> {
    column(batch, name)?
        .as_any()
        .downcast_ref::<DictionaryArray<Int16Type>>()
        .ok_or_else(|| invalid_data("POD5 column has unexpected type"))
}

fn column<'a>(batch: &'a RecordBatch, name: &str) -> io::Result<&'a dyn Array> {
    let index = batch
        .schema()
        .index_of(name)
        .map_err(|_| invalid_data("POD5 column is missing"))?;
    Ok(batch.column(index).as_ref())
}

fn uint64_list_value(array: &ListArray, row: usize) -> io::Result<Vec<u64>> {
    if array.is_null(row) {
        return Ok(Vec::new());
    }

    let values = array.value(row);
    let values = values
        .as_any()
        .downcast_ref::<UInt64Array>()
        .ok_or_else(|| invalid_data("POD5 signal row list has unexpected type"))?;

    Ok((0..values.len()).map(|index| values.value(index)).collect())
}

fn dictionary_string_value(array: &DictionaryArray<Int16Type>, row: usize) -> io::Result<String> {
    if array.is_null(row) {
        return Ok(String::new());
    }

    let key = array.keys().value(row);
    if key < 0 {
        return Err(invalid_data("POD5 dictionary key is negative"));
    }

    let values = array
        .values()
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| invalid_data("POD5 dictionary values have unexpected type"))?;
    let index = key as usize;
    if index >= values.len() {
        return Err(invalid_data("POD5 dictionary key is out of bounds"));
    }

    Ok(values.value(index).to_string())
}

fn uuid_bytes_to_string(bytes: &[u8]) -> io::Result<String> {
    if bytes.len() != 16 {
        return Err(invalid_data("POD5 read_id is not 16 bytes"));
    }

    Ok(format!(
        "{:02x}{:02x}{:02x}{:02x}-{:02x}{:02x}-{:02x}{:02x}-{:02x}{:02x}-{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}",
        bytes[0],
        bytes[1],
        bytes[2],
        bytes[3],
        bytes[4],
        bytes[5],
        bytes[6],
        bytes[7],
        bytes[8],
        bytes[9],
        bytes[10],
        bytes[11],
        bytes[12],
        bytes[13],
        bytes[14],
        bytes[15]
    ))
}

fn uuid_string_to_bytes(value: &str) -> io::Result<[u8; 16]> {
    let mut hex = String::with_capacity(32);
    for byte in value.bytes() {
        if byte != b'-' {
            hex.push(byte as char);
        }
    }
    if hex.len() != 32 || !hex.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err(invalid_data("POD5 UUID string is malformed"));
    }

    let mut bytes = [0u8; 16];
    for index in 0..16 {
        bytes[index] = u8::from_str_radix(&hex[index * 2..index * 2 + 2], 16)
            .map_err(|_| invalid_data("POD5 UUID string is malformed"))?;
    }

    Ok(bytes)
}

/// Compresses raw int16 ADC samples into a POD5 VBZ signal blob.
///
/// The transform implemented here is:
///
/// 1. first-order delta encode with wrapping `i16` subtraction,
/// 2. zigzag encode signed deltas into `u16`,
/// 3. SVB16 encode with one control bit per value, LSB first,
/// 4. zstd-compress the SVB16 byte stream at level 1.
pub fn compress_vbz_signal(samples: &[i16]) -> io::Result<Vec<u8>> {
    if samples.is_empty() {
        return Ok(Vec::new());
    }

    let inner = encode_vbz_inner(samples);
    zstd::bulk::compress(&inner, 1)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))
}

/// Decompresses a POD5 VBZ signal blob into raw int16 ADC samples.
///
/// `num_samples` is required because the SVB16 control stream length is derived
/// from the expected number of samples.
pub fn decompress_vbz_signal(data: &[u8], num_samples: usize) -> io::Result<Vec<i16>> {
    if num_samples == 0 {
        if data.is_empty() {
            return Ok(Vec::new());
        }
        return Err(invalid_data("empty VBZ signal expected for zero samples"));
    }

    let max_inner_len = num_samples
        .div_ceil(8)
        .checked_add(
            num_samples
                .checked_mul(2)
                .ok_or_else(|| invalid_data("VBZ decompressed length overflow"))?,
        )
        .ok_or_else(|| invalid_data("VBZ decompressed length overflow"))?;
    let inner = decompress_vbz_frame(data, max_inner_len)?;

    decode_vbz_inner(&inner, num_samples)
}

/// Highest compression ratio for which a frame's declared size is used to
/// size the output up front; above it the output grows as data decodes.
const VBZ_TRUSTED_RATIO: u64 = 64;

/// Decompresses a VBZ zstd frame that may hold at most `max_len` bytes.
///
/// `max_len` and the frame header come from the file, so neither is trusted to
/// size a large allocation before any data has decoded.
fn decompress_vbz_frame(data: &[u8], max_len: usize) -> io::Result<Vec<u8>> {
    let zstd_error =
        |error: io::Error| io::Error::new(io::ErrorKind::InvalidData, error.to_string());
    let declared = zstd::zstd_safe::get_frame_content_size(data).ok().flatten();
    if let Some(size) = declared
        && size <= max_len as u64
        && size <= (data.len() as u64).saturating_mul(VBZ_TRUSTED_RATIO)
    {
        return zstd::bulk::decompress(data, size as usize).map_err(zstd_error);
    }

    let mut inner = Vec::new();
    zstd::stream::read::Decoder::with_buffer(data)
        .map_err(zstd_error)?
        .take((max_len as u64).saturating_add(1))
        .read_to_end(&mut inner)
        .map_err(zstd_error)?;
    if inner.len() > max_len {
        return Err(invalid_data(
            "VBZ data decompresses past its declared samples",
        ));
    }
    Ok(inner)
}

fn encode_vbz_inner(samples: &[i16]) -> Vec<u8> {
    let control_len = samples.len().div_ceil(8);
    let mut output = vec![0u8; control_len];
    let mut previous = 0i16;

    for (index, &sample) in samples.iter().enumerate() {
        let delta = sample.wrapping_sub(previous);
        previous = sample;
        let code = zigzag_encode_i16(delta);

        if code <= u8::MAX as u16 {
            output.push(code as u8);
        } else {
            output[index / 8] |= 1 << (index % 8);
            output.extend_from_slice(&code.to_le_bytes());
        }
    }

    output
}

fn decode_vbz_inner(data: &[u8], num_samples: usize) -> io::Result<Vec<i16>> {
    let control_len = num_samples.div_ceil(8);
    if data.len() < control_len {
        return Err(invalid_data("VBZ control stream is truncated"));
    }

    let control = &data[..control_len];
    let values = &data[control_len..];
    // Every sample takes at least one value byte, so this bounds the output.
    if values.len() < num_samples {
        return Err(invalid_data("VBZ data stream is truncated"));
    }
    let mut value_offset = 0usize;
    let mut previous = 0i16;
    let mut output = Vec::with_capacity(num_samples);

    for index in 0..num_samples {
        let bit = (control[index / 8] >> (index % 8)) & 1;
        let code = if bit == 0 {
            let Some(&value) = values.get(value_offset) else {
                return Err(invalid_data("VBZ data stream is truncated"));
            };
            value_offset += 1;
            u16::from(value)
        } else {
            let bytes = values
                .get(value_offset..value_offset + 2)
                .ok_or_else(|| invalid_data("VBZ data stream is truncated"))?;
            value_offset += 2;
            u16::from_le_bytes([bytes[0], bytes[1]])
        };

        let delta = zigzag_decode_i16(code);
        previous = previous.wrapping_add(delta);
        output.push(previous);
    }

    if value_offset != values.len() {
        return Err(invalid_data("VBZ data stream has trailing bytes"));
    }

    Ok(output)
}

fn zigzag_encode_i16(value: i16) -> u16 {
    ((value as u16) << 1) ^ ((value >> 15) as u16)
}

fn zigzag_decode_i16(value: u16) -> i16 {
    ((value >> 1) ^ 0u16.wrapping_sub(value & 1)) as i16
}

fn slice_to_array<const N: usize>(slice: &[u8]) -> [u8; N] {
    let mut output = [0u8; N];
    output.copy_from_slice(slice);
    output
}

fn arrow_error(error: ArrowError) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, error.to_string())
}

/// Like [`arrow_error`], but a failure of the output stream is returned as
/// it is, so callers can tell a full disk or a closed pipe from bad data.
fn arrow_write_error(error: ArrowError) -> io::Error {
    match error {
        ArrowError::IoError(_, error) => error,
        error => arrow_error(error),
    }
}

fn invalid_data(message: impl Into<String>) -> io::Error {
    Error::invalid(Format::Pod5, message).into()
}

#[cfg(test)]
mod pod5_tests {
    use super::*;
    use std::cell::Cell;

    const A_100_POD5: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/A_100.pod5");
    const RECORD_COUNT: usize = 100;
    const FIRST_READ_ID: &str = "1cadb1e9-592f-4e22-9285-4626f2b7da9f";
    const LAST_READ_ID: &str = "6ae6c2e9-fe1a-4b1d-befb-ac98a5a16c9a";
    const RUN_ID: &str = "1f03c20c2da347cc99b58b232181a7464126f4cb";

    struct CountingReadSeek<R> {
        inner: R,
        bytes_read: Rc<Cell<u64>>,
    }

    impl<R> CountingReadSeek<R> {
        fn new(inner: R) -> (Self, Rc<Cell<u64>>) {
            let bytes_read = Rc::new(Cell::new(0));
            (
                Self {
                    inner,
                    bytes_read: bytes_read.clone(),
                },
                bytes_read,
            )
        }
    }

    impl<R: Read> Read for CountingReadSeek<R> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            let bytes = self.inner.read(buf)?;
            self.bytes_read.set(self.bytes_read.get() + bytes as u64);
            Ok(bytes)
        }
    }

    impl<R: Seek> Seek for CountingReadSeek<R> {
        fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
            self.inner.seek(pos)
        }
    }

    fn assert_records_equivalent(actual: &[Pod5Record], expected: &[Pod5Record]) {
        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.iter().zip(expected) {
            assert_eq!(actual.read_id, expected.read_id);
            assert_eq!(actual.signal_rows, expected.signal_rows);
            assert_eq!(actual.read_number, expected.read_number);
            assert_eq!(actual.start_sample, expected.start_sample);
            assert_eq!(
                actual.median_before.to_bits(),
                expected.median_before.to_bits()
            );
            assert_eq!(actual.num_minknow_events, expected.num_minknow_events);
            assert_eq!(
                actual.tracked_scaling_scale.to_bits(),
                expected.tracked_scaling_scale.to_bits()
            );
            assert_eq!(
                actual.tracked_scaling_shift.to_bits(),
                expected.tracked_scaling_shift.to_bits()
            );
            assert_eq!(
                actual.predicted_scaling_scale.to_bits(),
                expected.predicted_scaling_scale.to_bits()
            );
            assert_eq!(
                actual.predicted_scaling_shift.to_bits(),
                expected.predicted_scaling_shift.to_bits()
            );
            assert_eq!(
                actual.num_reads_since_mux_change,
                expected.num_reads_since_mux_change
            );
            assert_eq!(
                actual.time_since_mux_change.to_bits(),
                expected.time_since_mux_change.to_bits()
            );
            assert_eq!(actual.num_samples, expected.num_samples);
            assert_eq!(actual.channel, expected.channel);
            assert_eq!(actual.well, expected.well);
            assert_eq!(actual.pore_type, expected.pore_type);
            assert_eq!(
                actual.open_pore_level.map(f32::to_bits),
                expected.open_pore_level.map(f32::to_bits)
            );
            assert_eq!(
                actual.calibration_offset.to_bits(),
                expected.calibration_offset.to_bits()
            );
            assert_eq!(
                actual.calibration_scale.to_bits(),
                expected.calibration_scale.to_bits()
            );
            assert_eq!(actual.end_reason, expected.end_reason);
            assert_eq!(actual.end_reason_forced, expected.end_reason_forced);
            assert_eq!(actual.run_info, expected.run_info);
        }
    }

    #[test]
    fn reads_one_record_at_a_time_from_path() {
        let mut reader = Pod5Reader::from_path(A_100_POD5).unwrap();

        assert_eq!(reader.header.magic, *POD5_MAGIC);
        assert_eq!(reader.header.read_count(), RECORD_COUNT);
        assert_eq!(reader.header.signal_count(), RECORD_COUNT);
        assert_eq!(reader.header.run_info_count(), 1);
        assert_eq!(
            reader.header.software.as_deref(),
            Some("brust test fixture")
        );
        assert_eq!(reader.header.pod5_version.as_deref(), Some("0.3.28"));
        assert_eq!(reader.run_infos[0].acquisition_id, RUN_ID);
        assert_eq!(reader.run_infos[0].sample_rate, 5000);

        let first = reader.read_record().unwrap().unwrap();
        assert_eq!(first.read_id, FIRST_READ_ID);
        assert_eq!(first.channel, 269);
        assert_eq!(first.well, 3);
        assert_eq!(first.read_number, 5243);
        assert_eq!(first.start_sample, 18311270);
        assert_eq!(first.num_samples, 10162);
        assert_eq!(first.num_minknow_events, 1252);
        assert_eq!(first.pore_type, "not_set");
        assert_eq!(first.end_reason, "signal_positive");
        assert!(!first.end_reason_forced);
        assert_eq!(first.run_info, RUN_ID);
        assert_eq!(first.signal_rows, vec![0]);

        let signal = reader.signal_for_record(&first).unwrap();
        assert_eq!(signal.len(), first.num_samples as usize);
        assert_eq!(
            &signal[..20],
            &[
                617, 450, 473, 454, 454, 460, 463, 462, 467, 476, 464, 464, 464, 473, 490, 471,
                465, 442, 454, 470,
            ]
        );

        let mut count = 1;
        let mut last = first;
        while let Some(record) = reader.read_record().unwrap() {
            count += 1;
            last = record;
        }

        assert_eq!(count, RECORD_COUNT);
        assert_eq!(last.read_id, LAST_READ_ID);
        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn reader_from_seekable_source_does_not_buffer_entire_file() {
        let file = File::open(A_100_POD5).expect("fixture should open");
        let file_len = file.metadata().unwrap().len();
        let (file, bytes_read) = CountingReadSeek::new(file);
        let mut reader = Pod5Reader::from_reader(file).unwrap();

        assert!(
            bytes_read.get() < file_len,
            "Pod5Reader construction read {} bytes from a {file_len}-byte file",
            bytes_read.get()
        );

        let first = reader.read_record().unwrap().unwrap();
        assert_eq!(first.read_id, FIRST_READ_ID);
    }

    #[test]
    fn reader_caches_decompressed_signal_rows() {
        let file = File::open(A_100_POD5).expect("fixture should open");
        let (file, bytes_read) = CountingReadSeek::new(file);
        let mut reader = Pod5Reader::from_reader(file).unwrap();
        let first = reader.read_record().unwrap().unwrap();

        let signal = reader.signal_for_record(&first).unwrap();
        let bytes_after_first_signal = bytes_read.get();
        let cached_signal = reader.signal_for_record(&first).unwrap();

        assert_eq!(cached_signal, signal);
        assert_eq!(bytes_read.get(), bytes_after_first_signal);
    }

    #[test]
    fn read_all_materializes_sample_file() {
        let pod5 = Pod5::from_path(A_100_POD5).unwrap();

        assert_eq!(pod5.records.len(), RECORD_COUNT);
        assert_eq!(pod5.records[0].read_id, FIRST_READ_ID);
        assert_eq!(pod5.records[RECORD_COUNT - 1].read_id, LAST_READ_ID);
        assert_eq!(pod5.run_infos.len(), 1);
        assert_eq!(pod5.signals.len(), RECORD_COUNT);
        assert_eq!(pod5.header.sections.len(), 4);
        assert!(
            pod5.header
                .sections
                .iter()
                .any(|section| section.kind == Pod5SectionKind::Footer)
        );
        assert!(pod5.records.iter().all(|record| record.num_samples > 0));
        assert_eq!(
            pod5.signal_by_read_id(FIRST_READ_ID)
                .unwrap()
                .unwrap()
                .len(),
            pod5.records[0].num_samples as usize
        );
    }

    #[test]
    fn records_iterator_streams_until_eof() {
        let mut reader = Pod5Reader::from_path(A_100_POD5).unwrap();
        let records = reader.records().collect::<io::Result<Vec<_>>>().unwrap();

        assert_eq!(records.len(), RECORD_COUNT);
        assert_eq!(records[0].read_id, FIRST_READ_ID);
        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn pod5_from_reader_materializes_stream() {
        let file = File::open(A_100_POD5).expect("fixture should open");
        let pod5 = Pod5::from_reader(file).expect("POD5 should materialize from reader");

        assert_eq!(pod5.records.len(), RECORD_COUNT);
        assert_eq!(pod5.records[0].read_id, FIRST_READ_ID);
    }

    #[test]
    fn materialized_pod5_can_be_deep_cloned() {
        let original = Pod5::from_path(A_100_POD5).expect("POD5 should materialize");
        let mut cloned = original.clone();
        assert_eq!(cloned.records.len(), original.records.len());
        assert_eq!(cloned.records[0].read_id, original.records[0].read_id);
        assert_eq!(
            cloned.run_infos[0].acquisition_id,
            original.run_infos[0].acquisition_id
        );

        cloned.records[0].read_id.push_str("_clone");
        cloned.run_infos[0].sample_id.push_str("_clone");

        assert_ne!(cloned.records[0].read_id, original.records[0].read_id);
        assert_ne!(
            cloned.run_infos[0].sample_id,
            original.run_infos[0].sample_id
        );
        assert_eq!(original.records[0].read_id, FIRST_READ_ID);
        assert_eq!(original.run_infos[0].sample_id, "A_NB01");
    }

    #[test]
    fn vbz_transform_round_trips_samples() {
        let samples = [0i16, 1, -1, 255, 256, -256, i16::MAX, i16::MIN, -3, 4, 4, 5];
        let compressed = compress_vbz_signal(&samples).unwrap();
        let decoded = decompress_vbz_signal(&compressed, samples.len()).unwrap();

        assert_eq!(decoded, samples);
    }

    #[test]
    fn fixture_vbz_blobs_decompress_and_recompress_identically() {
        let pod5 = Pod5::from_path(A_100_POD5).unwrap();

        for signal in &pod5.signals {
            let samples = signal.decompress().unwrap();
            assert_eq!(samples.len(), signal.samples as usize);

            if let Some(original) = signal.compressed_bytes() {
                let recompressed = signal.compress().unwrap();
                assert_eq!(recompressed, original);
            }
        }
    }

    #[test]
    fn sample_summaries_saturate_instead_of_overflowing() {
        // num_samples comes from the file; summing crafted values must not panic.
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        for record in &mut pod5.records {
            record.num_samples = u64::MAX / 2 + 1;
        }

        assert_eq!(pod5.total_samples(), u64::MAX);
        assert!(
            pod5.channel_summaries()
                .iter()
                .all(|summary| summary.sample_count >= u64::MAX / 2)
        );
        assert_eq!(pod5.run_info_summaries()[0].sample_count, u64::MAX);
    }

    #[test]
    fn writer_round_trips_an_empty_payload() {
        // Official pod5 writes files with no reads; brust reads them.
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.records.clear();
        pod5.signals.clear();
        pod5.run_infos.clear();
        let mut output = Vec::new();
        pod5.to_writer(&mut output).unwrap();

        let round_tripped = Pod5::from_reader(&output[..]).unwrap();
        assert!(round_tripped.records.is_empty());
        assert!(round_tripped.signals.is_empty());
        assert!(round_tripped.run_infos.is_empty());
    }

    #[test]
    fn writer_round_trips_a_zero_sample_read() {
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.records.truncate(1);
        pod5.records[0].signal_rows.clear();
        pod5.records[0].num_samples = 0;
        pod5.signals.clear();
        let mut output = Vec::new();
        pod5.to_writer(&mut output).unwrap();

        let round_tripped = Pod5::from_reader(&output[..]).unwrap();
        assert_eq!(round_tripped.records[0].num_samples, 0);
        assert_eq!(
            round_tripped
                .signal_for_record(&round_tripped.records[0])
                .unwrap(),
            Vec::<i16>::new()
        );
    }

    #[test]
    fn too_many_distinct_dictionary_values_are_an_error_not_a_panic() {
        // Dictionary keys are 16-bit, as in official POD5 files.
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        let template = pod5.records[0].clone();
        pod5.records = (0..32_769)
            .map(|index| Pod5Record {
                pore_type: format!("pore-{index}"),
                ..template.clone()
            })
            .collect();
        let mut output = Vec::new();

        let error = pod5.to_writer(&mut output).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(output.is_empty());
    }

    #[test]
    fn writer_round_trips_materialized_pod5() {
        let pod5 = Pod5::from_path(A_100_POD5).unwrap();
        let mut output = Vec::new();
        pod5.to_writer(&mut output).unwrap();
        let round_tripped = Pod5::from_reader(&output[..]).unwrap();

        assert_records_equivalent(&round_tripped.records, &pod5.records);
        assert_eq!(round_tripped.run_infos, pod5.run_infos);
        assert_eq!(round_tripped.signals, pod5.signals);
        assert_eq!(
            round_tripped.header.file_identifier,
            pod5.header.file_identifier
        );
        assert_eq!(round_tripped.header.software, pod5.header.software);
        assert_eq!(round_tripped.header.pod5_version, pod5.header.pod5_version);
        assert_eq!(round_tripped.header.read_count(), RECORD_COUNT);
        assert_eq!(round_tripped.header.signal_count(), RECORD_COUNT);
        assert_eq!(round_tripped.header.run_info_count(), 1);
    }

    fn written_reads_field_names(pod5: &Pod5) -> Vec<String> {
        let mut output = Vec::new();
        pod5.to_writer(&mut output).unwrap();
        let reader = Pod5Reader::from_reader(Cursor::new(output.clone())).unwrap();
        let section = reader
            .header
            .sections
            .iter()
            .find(|section| section.kind == Pod5SectionKind::Reads)
            .unwrap();
        let start = section.offset as usize;
        let arrow = output[start..start + section.length as usize].to_vec();
        let arrow_reader = FileReader::try_new(Cursor::new(arrow), None).unwrap();
        arrow_reader
            .schema()
            .fields()
            .iter()
            .map(|field| field.name().clone())
            .collect()
    }

    #[test]
    fn writer_emits_open_pore_level_for_versions_that_require_it() {
        // Official pod5 >= 0.3.30 refuses Reads tables without this column.
        for version in [None, Some("0.3.30"), Some("0.3.34"), Some("0.3.39")] {
            let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
            pod5.header.pod5_version = version.map(str::to_string);
            let fields = written_reads_field_names(&pod5);
            assert!(
                fields.iter().any(|name| name == "open_pore_level"),
                "version {version:?} wrote {fields:?}"
            );
        }
    }

    fn assert_fixture_run_info_columns(run_info: &Pod5RunInfo) {
        let pair = |key: &str, value: &str| (key.to_string(), value.to_string());
        assert_eq!(run_info.acquisition_start_time, 1_688_015_479_792);
        assert_eq!(run_info.adc_max, 4095);
        assert_eq!(run_info.adc_min, -4096);
        assert_eq!(run_info.context_tags.len(), 8);
        assert_eq!(run_info.context_tags[0], pair("barcoding_enabled", "0"));
        assert_eq!(
            run_info.context_tags[7],
            pair("sequencing_kit", "sqk-nbd114-24")
        );
        assert_eq!(run_info.flow_cell_product_code, "FLO-MIN114");
        assert_eq!(
            run_info.protocol_name,
            "sequencing/sequencing_MIN114_DNA_e8_2_400K:FLO-MIN114:SQK-NBD114-24:400"
        );
        assert_eq!(
            run_info.protocol_run_id,
            "f51f2183-4844-4a44-9be9-d7455407b836"
        );
        assert_eq!(run_info.protocol_start_time, 1_688_015_166_811);
        assert_eq!(run_info.sequencer_position, "MN40918");
        assert_eq!(run_info.sequencer_position_type, "MinION Mk1B");
        assert_eq!(run_info.system_name, "grhl-c214-07");
        assert_eq!(run_info.system_type, "Darwin 21.6.0");
        assert_eq!(run_info.tracking_id.len(), 29);
        assert_eq!(run_info.tracking_id[0], pair("asic_id", "616835618"));
        assert_eq!(run_info.tracking_id[28], pair("version", "5.5.3"));
    }

    #[test]
    fn fixture_run_info_reads_every_column() {
        let pod5 = Pod5::from_path(A_100_POD5).unwrap();

        assert_fixture_run_info_columns(&pod5.run_infos[0]);
    }

    #[test]
    fn writer_preserves_every_run_info_column() {
        let pod5 = Pod5::from_path(A_100_POD5).unwrap();
        let mut output = Vec::new();
        pod5.to_writer(&mut output).unwrap();
        let round_tripped = Pod5::from_reader(&output[..]).unwrap();

        assert_fixture_run_info_columns(&round_tripped.run_infos[0]);
    }

    #[test]
    fn open_pore_level_version_threshold_is_numeric() {
        for (version, expected) in [
            ("0.3.29", false),
            ("0.2.99", false),
            ("0.3.30", true),
            ("0.3.100", true),
            ("0.10.0", true),
            ("1.0.0", true),
            // Unparseable versions get the newest schema.
            ("0.3", true),
            ("0.3.29.1", true),
            ("+0.3.29", true),
            ("", true),
        ] {
            assert_eq!(writes_open_pore_level(version), expected, "{version}");
        }
    }

    #[test]
    fn writer_omits_open_pore_level_for_older_versions() {
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.header.pod5_version = Some("0.3.29".to_string());
        let fields = written_reads_field_names(&pod5);

        assert!(!fields.iter().any(|name| name == "open_pore_level"));
        assert_eq!(fields.len(), READS_TABLE_FIELDS.len());
    }

    #[test]
    fn fixture_without_open_pore_level_reads_as_none() {
        let pod5 = Pod5::from_path(A_100_POD5).unwrap();

        assert!(
            pod5.records
                .iter()
                .all(|record| record.open_pore_level.is_none())
        );
    }

    #[test]
    fn open_pore_level_round_trips() {
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.header.pod5_version = Some("0.3.39".to_string());
        pod5.records[0].open_pore_level = Some(222.5);
        let mut output = Vec::new();
        pod5.to_writer(&mut output).unwrap();
        let round_tripped = Pod5::from_reader(&output[..]).unwrap();

        assert_eq!(round_tripped.records[0].open_pore_level, Some(222.5));
        assert_eq!(round_tripped.records[1].open_pore_level, None);
    }

    #[test]
    fn writer_rejects_open_pore_level_for_versions_without_the_column() {
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.header.pod5_version = Some("0.3.28".to_string());
        pod5.records[0].open_pore_level = Some(222.5);
        let mut output = Vec::new();

        let error = pod5.to_writer(&mut output).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(output.is_empty());
    }

    #[test]
    fn writer_compresses_uncompressed_signal_rows() {
        let read_id = "00000000-0000-0000-0000-000000000001".to_string();
        let run_id = "run-1".to_string();
        let samples = vec![10, 12, -3, 400, 399, i16::MIN, i16::MAX];
        let pod5 = Pod5 {
            header: Pod5Header {
                magic: *POD5_MAGIC,
                section_marker: [0; 16],
                sections: Vec::new(),
                file_identifier: Some("00000000-0000-0000-0000-000000000002".to_string()),
                software: Some("brust pod5 writer test".to_string()),
                pod5_version: Some("0.3.34".to_string()),
            },
            run_infos: vec![Pod5RunInfo {
                acquisition_id: run_id.clone(),
                sample_id: "sample".to_string(),
                experiment_name: "experiment".to_string(),
                flow_cell_id: "flow-cell".to_string(),
                sequencing_kit: "kit".to_string(),
                sample_rate: 5000,
                software: "software".to_string(),
                acquisition_start_time: 1_700_000_000_123,
                adc_max: 2047,
                adc_min: -2048,
                context_tags: vec![("experiment_type".to_string(), "rna".to_string())],
                flow_cell_product_code: "FLO-PRO114M".to_string(),
                protocol_name: "protocol".to_string(),
                protocol_run_id: "protocol-run".to_string(),
                protocol_start_time: 1_699_999_999_000,
                sequencer_position: "1A".to_string(),
                sequencer_position_type: "PromethION".to_string(),
                system_name: "host".to_string(),
                system_type: "Linux".to_string(),
                tracking_id: vec![
                    ("asic_id".to_string(), "1".to_string()),
                    ("device_id".to_string(), "PC24B".to_string()),
                ],
            }],
            signals: vec![Pod5Signal {
                read_id: read_id.clone(),
                samples: samples.len() as u32,
                payload: Pod5SignalPayload::Uncompressed(samples.clone()),
            }],
            records: vec![Pod5Record {
                read_id: read_id.clone(),
                signal_rows: vec![0],
                read_number: 1,
                start_sample: 42,
                median_before: 0.5,
                num_minknow_events: 7,
                tracked_scaling_scale: 1.0,
                tracked_scaling_shift: 2.0,
                predicted_scaling_scale: 3.0,
                predicted_scaling_shift: 4.0,
                num_reads_since_mux_change: 5,
                time_since_mux_change: 6.0,
                num_samples: samples.len() as u64,
                channel: 10,
                well: 2,
                pore_type: "not_set".to_string(),
                calibration_offset: 8.0,
                calibration_scale: 9.0,
                end_reason: "signal_positive".to_string(),
                end_reason_forced: false,
                run_info: run_id,
                open_pore_level: None,
            }],
        };

        let mut output = Vec::new();
        pod5.to_writer(&mut output).unwrap();
        let round_tripped = Pod5::from_reader(&output[..]).unwrap();

        assert!(round_tripped.signals[0].is_vbz_compressed());
        assert_eq!(
            round_tripped.signal_by_read_id(&read_id).unwrap().unwrap(),
            samples
        );
        assert_records_equivalent(&round_tripped.records, &pod5.records);
        assert_eq!(round_tripped.run_infos, pod5.run_infos);
    }

    /// Where an Arrow footer `Block` of the first `kind` section sits in `data`.
    fn arrow_block_position(
        data: &[u8],
        kind: Pod5SectionKind,
        dictionary: bool,
        index: usize,
    ) -> usize {
        let reader = Pod5Reader::from_reader(Cursor::new(data.to_vec())).unwrap();
        let section = reader
            .header
            .sections
            .iter()
            .find(|section| section.kind == kind)
            .unwrap();
        let (start, end) = (
            section.offset as usize,
            (section.offset + section.length) as usize,
        );
        let footer_len = i32::from_le_bytes(data[end - 10..end - 6].try_into().unwrap()) as usize;
        let footer_start = end - 10 - footer_len;
        let footer = &data[footer_start..end - 10];
        let parsed = root_as_footer(footer).unwrap();
        let blocks = if dictionary {
            parsed.dictionaries().unwrap()
        } else {
            parsed.recordBatches().unwrap()
        };
        let block: &arrow_ipc::Block = blocks.get(index);
        assert!(start < footer_start);
        footer_start + (block as *const arrow_ipc::Block as usize - footer.as_ptr() as usize)
    }

    /// The fixture with one Arrow footer block's `bodyLength` replaced.
    fn fixture_with_block_body_length(
        kind: Pod5SectionKind,
        dictionary: bool,
        body: i64,
    ) -> Vec<u8> {
        let mut data = std::fs::read(A_100_POD5).unwrap();
        // Block layout: offset i64, metaDataLength i32, padding, bodyLength i64.
        let position = arrow_block_position(&data, kind, dictionary, 0) + 16;
        data[position..position + 8].copy_from_slice(&body.to_le_bytes());
        data
    }

    fn assert_invalid_pod5(data: Vec<u8>) {
        let error = Pod5Reader::from_reader(Cursor::new(data.clone()))
            .and_then(Pod5Reader::read_all)
            .expect_err("malformed POD5 should be rejected");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData, "{error}");
        assert!(Pod5::from_reader(&data[..]).is_err());
    }

    #[test]
    fn negative_dictionary_body_length_is_rejected() {
        assert_invalid_pod5(fixture_with_block_body_length(
            Pod5SectionKind::Reads,
            true,
            -1,
        ));
    }

    #[test]
    fn negative_record_batch_body_length_is_rejected() {
        assert_invalid_pod5(fixture_with_block_body_length(
            Pod5SectionKind::Reads,
            false,
            -8,
        ));
    }

    #[test]
    fn record_batch_body_past_its_section_is_rejected() {
        assert_invalid_pod5(fixture_with_block_body_length(
            Pod5SectionKind::Signal,
            false,
            1 << 44,
        ));
    }

    #[test]
    fn record_batch_length_that_disagrees_with_its_columns_is_rejected() {
        let mut data = std::fs::read(A_100_POD5).unwrap();
        let position = arrow_batch_length_position(&data, Pod5SectionKind::Signal);
        data[position..position + 8].copy_from_slice(&(1i64 << 40).to_le_bytes());

        assert_invalid_pod5(data);
    }

    /// Where the Message of an Arrow footer block starts in `data`.
    fn arrow_message_start(data: &[u8], kind: Pod5SectionKind, dictionary: bool) -> usize {
        let block = arrow_block_position(data, kind, dictionary, 0);
        let section = Pod5Reader::from_reader(Cursor::new(data.to_vec()))
            .unwrap()
            .header
            .sections
            .into_iter()
            .find(|section| section.kind == kind)
            .unwrap();
        let block_offset = i64::from_le_bytes(data[block..block + 8].try_into().unwrap());
        // Block metadata is a 0xFFFFFFFF continuation, a length, then the Message.
        section.offset as usize + block_offset as usize + 8
    }

    /// Where a batch's first field node (length i64, then null_count i64) sits.
    fn arrow_first_node_position(data: &[u8], kind: Pod5SectionKind, dictionary: bool) -> usize {
        let start = arrow_message_start(data, kind, dictionary);
        let message = root_as_message(&data[start..]).unwrap();
        let batch = if dictionary {
            message
                .header_as_dictionary_batch()
                .unwrap()
                .data()
                .unwrap()
        } else {
            message.header_as_record_batch().unwrap()
        };
        let node: &arrow_ipc::FieldNode = batch.nodes().unwrap().get(0);
        start + (node as *const arrow_ipc::FieldNode as usize - data[start..].as_ptr() as usize)
    }

    /// Where the first record batch's `length` field sits in `data`.
    fn arrow_batch_length_position(data: &[u8], kind: Pod5SectionKind) -> usize {
        let start = arrow_message_start(data, kind, false);
        let message = &data[start..];
        let table = fb_root_table(message).unwrap();
        let header_field = fb_field_position(message, table, 8).unwrap().unwrap();
        let header = fb_uoffset_target(message, header_field).unwrap();
        start + fb_field_position(message, header, 4).unwrap().unwrap()
    }

    #[test]
    fn unknown_arrow_column_type_is_rejected() {
        // arrow-ipc panics ("Type NONE not supported") converting this schema.
        let mut data = std::fs::read(A_100_POD5).unwrap();
        let section = Pod5Reader::from_reader(Cursor::new(data.clone()))
            .unwrap()
            .header
            .sections
            .into_iter()
            .find(|section| section.kind == Pod5SectionKind::Reads)
            .unwrap();
        let end = (section.offset + section.length) as usize;
        let footer_len = i32::from_le_bytes(data[end - 10..end - 6].try_into().unwrap()) as usize;
        let footer_start = end - 10 - footer_len;
        let footer = &data[footer_start..end - 10];
        let footer_table = fb_root_table(footer).unwrap();
        let schema_field = fb_field_position(footer, footer_table, 6).unwrap().unwrap();
        let schema = fb_uoffset_target(footer, schema_field).unwrap();
        let signal_field = fb_table_vector_field(footer, schema, 6).unwrap()[1];
        let type_type = fb_field_position(footer, signal_field, 8).unwrap().unwrap();
        data[footer_start + type_type] = 0;

        assert_invalid_pod5(data);
    }

    #[test]
    fn signal_rows_with_an_unchecked_first_column_are_rejected() {
        // A first column of a type the layout checks don't model leaves the
        // row count unbounded, so the section is rejected when opened.
        let mut data = std::fs::read(A_100_POD5).unwrap();
        let batch_length = arrow_batch_length_position(&data, Pod5SectionKind::Signal);
        let node_length = arrow_first_node_position(&data, Pod5SectionKind::Signal, false);
        let section = Pod5Reader::from_reader(Cursor::new(data.clone()))
            .unwrap()
            .header
            .sections
            .into_iter()
            .find(|section| section.kind == Pod5SectionKind::Signal)
            .unwrap();
        let end = (section.offset + section.length) as usize;
        let footer_len = i32::from_le_bytes(data[end - 10..end - 6].try_into().unwrap()) as usize;
        let footer_start = end - 10 - footer_len;
        let footer = &data[footer_start..end - 10];
        let footer_table = fb_root_table(footer).unwrap();
        let schema_field = fb_field_position(footer, footer_table, 6).unwrap().unwrap();
        let schema = fb_uoffset_target(footer, schema_field).unwrap();
        let read_id_field = fb_table_vector_field(footer, schema, 6).unwrap()[0];
        let type_type = footer_start
            + fb_field_position(footer, read_id_field, 8)
                .unwrap()
                .unwrap();
        data[type_type] = arrow_ipc::Type::BinaryView.0;
        data[batch_length..batch_length + 8].copy_from_slice(&i64::MAX.to_le_bytes());
        data[node_length..node_length + 8].copy_from_slice(&i64::MAX.to_le_bytes());

        let error = Pod5Reader::from_reader(Cursor::new(data))
            .map(|reader| reader.header.signal_count())
            .expect_err("unbounded Signal row count should be rejected");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn zero_width_first_column_cannot_vouch_for_a_row_count() {
        // FixedSizeBinary(0) needs no bytes however many rows it claims.
        let mut data = std::fs::read(A_100_POD5).unwrap();
        let batch_length = arrow_batch_length_position(&data, Pod5SectionKind::Signal);
        let node_length = arrow_first_node_position(&data, Pod5SectionKind::Signal, false);
        let section = Pod5Reader::from_reader(Cursor::new(data.clone()))
            .unwrap()
            .header
            .sections
            .into_iter()
            .find(|section| section.kind == Pod5SectionKind::Signal)
            .unwrap();
        let end = (section.offset + section.length) as usize;
        let footer_len = i32::from_le_bytes(data[end - 10..end - 6].try_into().unwrap()) as usize;
        let footer_start = end - 10 - footer_len;
        let footer = &data[footer_start..end - 10];
        let footer_table = fb_root_table(footer).unwrap();
        let schema_field = fb_field_position(footer, footer_table, 6).unwrap().unwrap();
        let schema = fb_uoffset_target(footer, schema_field).unwrap();
        let read_id_field = fb_table_vector_field(footer, schema, 6).unwrap()[0];
        let type_field = fb_field_position(footer, read_id_field, 10)
            .unwrap()
            .unwrap();
        let binary_type = fb_uoffset_target(footer, type_field).unwrap();
        let byte_width = footer_start + fb_field_position(footer, binary_type, 4).unwrap().unwrap();
        assert_eq!(data[byte_width..byte_width + 4], 16i32.to_le_bytes());
        data[byte_width..byte_width + 4].copy_from_slice(&0i32.to_le_bytes());
        data[batch_length..batch_length + 8].copy_from_slice(&i64::MAX.to_le_bytes());
        data[node_length..node_length + 8].copy_from_slice(&i64::MAX.to_le_bytes());

        let error = Pod5Reader::from_reader(Cursor::new(data))
            .map(|reader| reader.header.signal_count())
            .expect_err("unbounded Signal row count should be rejected");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn compressed_arrow_bodies_are_rejected() {
        use arrow_ipc::CompressionType;
        use arrow_ipc::writer::IpcWriteOptions;

        let pod5 = Pod5::from_path(A_100_POD5).unwrap();
        let metadata = pod5_writer_metadata(&pod5);
        let schema = Arc::new(signal_schema(&metadata));
        let batch = build_signal_batch(&pod5.signals, &schema).unwrap();
        let options = IpcWriteOptions::default()
            .try_with_compression(Some(CompressionType::ZSTD))
            .unwrap();
        let mut data = Vec::new();
        let mut writer =
            FileWriter::try_new_with_options(&mut data, batch.schema().as_ref(), options).unwrap();
        writer.write(&batch).unwrap();
        writer.finish().unwrap();
        drop(writer);

        let length = data.len() as u64;
        let shared = Rc::new(RefCell::new(Cursor::new(data)));
        let error = validate_arrow_section(shared, 0, length).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn null_count_without_a_validity_buffer_is_rejected() {
        for (kind, dictionary) in [
            (Pod5SectionKind::Reads, false),
            (Pod5SectionKind::Signal, false),
            (Pod5SectionKind::Reads, true),
        ] {
            let mut data = std::fs::read(A_100_POD5).unwrap();
            let null_count = arrow_first_node_position(&data, kind, dictionary) + 8;
            data[null_count..null_count + 8].copy_from_slice(&1i64.to_le_bytes());

            assert_invalid_pod5(data);
        }
    }

    #[test]
    fn batch_and_column_lengths_beyond_their_buffers_are_rejected() {
        let mut data = std::fs::read(A_100_POD5).unwrap();
        let batch_length = arrow_batch_length_position(&data, Pod5SectionKind::Signal);
        let node_length = arrow_first_node_position(&data, Pod5SectionKind::Signal, false);
        data[batch_length..batch_length + 8].copy_from_slice(&i64::MAX.to_le_bytes());
        data[node_length..node_length + 8].copy_from_slice(&i64::MAX.to_le_bytes());

        // Opening must fail: header counts are reported from the open alone.
        let error = Pod5Reader::from_reader(Cursor::new(data.clone()))
            .map(|reader| reader.header.signal_count())
            .expect_err("bogus Signal row count should be rejected");
        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert_invalid_pod5(data);
    }

    #[test]
    fn largest_representable_vbz_limit_does_not_overflow() {
        // div_ceil(n, 8) + 2n == usize::MAX exactly for this n on 64-bit.
        assert!(decompress_vbz_signal(&[], usize::MAX / 17 * 8).is_err());
    }

    #[test]
    fn huge_footer_vector_length_is_rejected_without_allocating_it() {
        let metadata = Pod5WriterMetadata {
            file_identifier: "00000000-0000-0000-0000-000000000002".to_string(),
            software: "test".to_string(),
            pod5_version: "0.3.34".to_string(),
        };
        let entry = Pod5FooterEntry {
            offset: 24,
            length: 8,
            content_type: 0,
        };
        let mut footer = build_pod5_footer(&metadata, &[entry]);
        let table = fb_root_table(&footer).unwrap();
        let field = fb_field_position(&footer, table, 10).unwrap().unwrap();
        let vector = fb_uoffset_target(&footer, field).unwrap();
        footer[vector..vector + 4].copy_from_slice(&u32::MAX.to_le_bytes());

        assert!(parse_pod5_footer(&footer).is_err());
    }

    #[test]
    fn huge_declared_sample_counts_are_errors_not_allocations() {
        let samples = [1i16, 2, 3];
        let compressed = compress_vbz_signal(&samples).unwrap();
        assert!(decompress_vbz_signal(&compressed, 1 << 40).is_err());

        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.records[0].num_samples = 1 << 62;
        assert!(pod5.signal_for_record(&pod5.records[0]).is_err());
        let mut reader = Pod5Reader::from_path(A_100_POD5).unwrap();
        let mut record = reader.read_record().unwrap().unwrap();
        record.num_samples = 1 << 62;
        assert!(reader.signal_for_record(&record).is_err());
    }

    /// An Arrow Signal table with one record batch per entry of `batches`,
    /// and each batch's row count. Rows are `(last read ID byte, samples,
    /// declared sample count)`.
    #[allow(clippy::type_complexity)]
    fn signal_table(batches: &[&[(u8, &[i16], u32)]]) -> (Vec<u8>, Vec<usize>) {
        use arrow_array::builder::{Int16Builder, LargeListBuilder};

        let schema = Arc::new(Schema::new(vec![
            uuid_field("read_id"),
            Field::new(
                "signal",
                DataType::LargeList(Arc::new(Field::new_list_field(DataType::Int16, true))),
                false,
            ),
            Field::new("samples", DataType::UInt32, false),
        ]));
        let batch = |rows: &[(u8, &[i16], u32)]| {
            let mut read_id = FixedSizeBinaryBuilder::new(16);
            let mut signal = LargeListBuilder::new(Int16Builder::new());
            for (id, values, _) in rows {
                let mut bytes = [0u8; 16];
                bytes[15] = *id;
                read_id.append_value(bytes).unwrap();
                signal.values().append_slice(values);
                signal.append(true);
            }
            let samples = UInt32Array::from(rows.iter().map(|row| row.2).collect::<Vec<_>>());
            let signal = signal
                .finish()
                .into_data()
                .into_builder()
                .data_type(schema.field(1).data_type().clone())
                .build()
                .unwrap();
            RecordBatch::try_new(
                schema.clone(),
                vec![
                    Arc::new(read_id.finish()),
                    Arc::new(arrow_array::LargeListArray::from(signal)),
                    Arc::new(samples),
                ],
            )
            .unwrap()
        };

        let mut data = Vec::new();
        let mut writer = FileWriter::try_new(&mut data, &schema).unwrap();
        for rows in batches {
            writer.write(&batch(rows)).unwrap();
        }
        writer.finish().unwrap();
        drop(writer);
        (data, batches.iter().map(|rows| rows.len()).collect())
    }

    fn signal_cursor_over<R: Read + Seek>(
        reader: R,
        length: u64,
        batch_rows: Vec<usize>,
    ) -> Pod5SignalCursor<R> {
        let section = Pod5Section {
            kind: Pod5SectionKind::Signal,
            offset: 0,
            length,
            padded_length: length,
            row_count: batch_rows.iter().sum(),
        };
        Pod5SignalCursor::new(Rc::new(RefCell::new(reader)), vec![(section, batch_rows)])
    }

    /// A Signal table of three batches: `[row 0]`, `[row 1]`, `[rows 2, 3]`.
    /// Row 1 declares three samples but stores two, so its batch fails to parse.
    fn signal_cursor_with_bad_middle_batch() -> Pod5SignalCursor<Cursor<Vec<u8>>> {
        let (data, batch_rows) = signal_table(&[
            &[(1, &[10, 11], 2)],
            &[(2, &[20, 21], 3)],
            &[(3, &[30, 31], 2), (4, &[40, 41], 2)],
        ]);
        let length = data.len() as u64;
        signal_cursor_over(Cursor::new(data), length, batch_rows)
    }

    #[test]
    fn signal_rows_are_read_from_their_own_batch() {
        // Out-of-order access used to restart at row 0 and read every batch
        // before the one wanted, which is quadratic over a whole file.
        let samples: Vec<i16> = (0..2_000).collect();
        let rows: Vec<(u8, &[i16], u32)> = (0..50).map(|id| (id, &samples[..], 2_000)).collect();
        let batches: Vec<&[(u8, &[i16], u32)]> = rows.chunks(1).collect();
        let (data, batch_rows) = signal_table(&batches);
        let length = data.len() as u64;
        let batch_bytes = length / 50;
        let (reader, bytes_read) = CountingReadSeek::new(Cursor::new(data));
        let mut cursor = signal_cursor_over(reader, length, batch_rows);

        for row in [49, 0, 25, 24, 48] {
            let before = bytes_read.get();
            let signal = cursor.signal_row_at(row).unwrap();
            assert_eq!(signal.read_id[34..], format!("{row:02x}"));
            assert!(bytes_read.get() - before < 3 * batch_bytes, "row {row}");
        }
    }

    #[test]
    fn rows_of_the_batch_just_read_are_served_without_reading_again() {
        let samples = [1i16, 2, 3];
        let rows: Vec<(u8, &[i16], u32)> = (0..10).map(|id| (id, &samples[..], 3)).collect();
        let (data, batch_rows) = signal_table(&[&rows]);
        let length = data.len() as u64;
        let (reader, bytes_read) = CountingReadSeek::new(Cursor::new(data));
        let mut cursor = signal_cursor_over(reader, length, batch_rows);

        cursor.signal_row_at(9).unwrap();
        let before = bytes_read.get();
        for row in [0, 5, 3] {
            let signal = cursor.signal_row_at(row).unwrap();
            assert_eq!(signal.read_id[34..], format!("{row:02x}"));
        }
        assert_eq!(bytes_read.get(), before);

        // A row already handed out is read again with its batch.
        assert_eq!(cursor.signal_row_at(9).unwrap().read_id[34..], *"09");
        assert!(bytes_read.get() > before);
        assert_eq!(cursor.signal_row_at(8).unwrap().read_id[34..], *"08");
    }

    fn row_read_id(cursor: &mut Pod5SignalCursor<Cursor<Vec<u8>>>, row: u64) -> io::Result<String> {
        cursor.signal_row_at(row).map(|signal| signal.read_id)
    }

    #[test]
    fn signal_rows_after_a_bad_batch_keep_their_numbers() {
        let mut cursor = signal_cursor_with_bad_middle_batch();

        assert_eq!(
            row_read_id(&mut cursor, 0).unwrap(),
            "00000000-0000-0000-0000-000000000001"
        );
        assert!(row_read_id(&mut cursor, 1).is_err());
        assert_eq!(
            row_read_id(&mut cursor, 2).unwrap(),
            "00000000-0000-0000-0000-000000000003"
        );
        assert!(row_read_id(&mut cursor, 1).is_err());
        assert_eq!(
            row_read_id(&mut cursor, 3).unwrap(),
            "00000000-0000-0000-0000-000000000004"
        );
    }

    #[test]
    fn signal_rows_past_a_bad_batch_are_reachable_on_a_fresh_cursor() {
        let mut cursor = signal_cursor_with_bad_middle_batch();
        assert_eq!(
            row_read_id(&mut cursor, 3).unwrap(),
            "00000000-0000-0000-0000-000000000004"
        );

        let mut cursor = signal_cursor_with_bad_middle_batch();
        assert_eq!(
            row_read_id(&mut cursor, 2).unwrap(),
            "00000000-0000-0000-0000-000000000003"
        );
        assert!(row_read_id(&mut cursor, 4).is_err());
    }

    fn assert_writer_rejects(pod5: &Pod5) {
        let mut output = Vec::new();
        let error = pod5.to_writer(&mut output).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData, "{error}");
        assert!(output.is_empty());
    }

    #[test]
    fn writer_rejects_uncompressed_payload_shorter_than_declared_samples() {
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        let mut samples = pod5.signals[0].decompress().unwrap();
        samples.pop();
        pod5.signals[0].payload = Pod5SignalPayload::Uncompressed(samples);

        assert_writer_rejects(&pod5);
    }

    #[test]
    fn writer_rejects_vbz_payload_that_disagrees_with_declared_samples() {
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.signals[0].samples += 1;
        let record = pod5
            .records
            .iter_mut()
            .find(|record| record.signal_rows == [0])
            .unwrap();
        record.num_samples += 1;

        assert_writer_rejects(&pod5);
    }

    #[test]
    fn writer_rejects_reads_that_reference_another_reads_signal() {
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        let other_row = pod5
            .records
            .iter()
            .find(|record| record.read_id != pod5.records[0].read_id)
            .unwrap()
            .signal_rows[0];
        pod5.records[0].signal_rows = vec![other_row];
        pod5.records[0].num_samples = u64::from(pod5.signals[other_row as usize].samples);

        assert_writer_rejects(&pod5);
    }

    #[test]
    fn writer_rejects_malformed_header_metadata() {
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.header.pod5_version = Some("0.3".to_string());
        assert_writer_rejects(&pod5);

        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.header.file_identifier = Some("garbage".to_string());
        assert_writer_rejects(&pod5);

        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.header.file_identifier = Some("1cadb1e9592f4e2292854626f2b7da9f".to_string());
        assert_writer_rejects(&pod5);
    }

    #[test]
    fn writer_rejects_unknown_end_reason() {
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.records[0].end_reason = "garbage".to_string();

        assert_writer_rejects(&pod5);
    }

    #[test]
    fn writer_accepts_every_official_end_reason() {
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        let reasons = [
            "unknown",
            "mux_change",
            "unblock_mux_change",
            "data_service_unblock_mux_change",
            "signal_positive",
            "signal_negative",
            "api_request",
            "device_data_error",
            "analysis_config_change",
            "paused",
        ];
        for (record, reason) in pod5.records.iter_mut().zip(reasons.iter().cycle()) {
            record.end_reason = reason.to_string();
        }
        let mut output = Vec::new();
        pod5.to_writer(&mut output).unwrap();
        let round_tripped = Pod5::from_reader(&output[..]).unwrap();

        assert_eq!(round_tripped.records[9].end_reason, "paused");
    }

    #[test]
    fn writer_rejects_invalid_signal_row_references() {
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.records[0].signal_rows = vec![pod5.signals.len() as u64];
        let mut output = Vec::new();

        assert!(pod5.to_writer(&mut output).is_err());
    }

    #[test]
    fn from_reader_accepts_section_marker_bytes_inside_data() {
        // Payloads without a marker get brust's fixed one; data that happens
        // to contain it must not split a section.
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.header.section_marker = [0; 16];
        pod5.run_infos[0].experiment_name = "BRUSTPOD5WRITER!".to_string();
        let mut output = Vec::new();
        pod5.to_writer(&mut output).unwrap();

        let buffered = Pod5::from_reader(&output[..]).unwrap();

        assert_eq!(buffered.records.len(), RECORD_COUNT);
        assert_eq!(buffered.run_infos[0].experiment_name, "BRUSTPOD5WRITER!");
    }

    #[test]
    fn from_reader_rejects_a_zeroed_footer() {
        let mut data = std::fs::read(A_100_POD5).unwrap();
        let footer = Pod5::from_path(A_100_POD5)
            .unwrap()
            .header
            .sections
            .into_iter()
            .find(|section| section.kind == Pod5SectionKind::Footer)
            .unwrap();
        let payload = footer.offset as usize + POD5_FOOTER_MAGIC.len();
        data[payload..payload + footer.length as usize].fill(0);

        assert!(
            Pod5Reader::from_reader(Cursor::new(data.clone()))
                .and_then(Pod5Reader::read_all)
                .is_err()
        );
        assert!(Pod5::from_reader(&data[..]).is_err());
    }

    #[test]
    fn malformed_wrappers_return_errors() {
        assert!(Pod5::from_reader(&b""[..]).is_err());

        let mut data = std::fs::read(A_100_POD5).unwrap();
        data[0] = 0;
        assert!(Pod5::from_reader(&data[..]).is_err());

        let mut data = std::fs::read(A_100_POD5).unwrap();
        let len = data.len();
        data[len - 1] = 0;
        assert!(Pod5::from_reader(&data[..]).is_err());
    }
}
