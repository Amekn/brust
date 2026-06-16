//! Minimal POD5 parsing primitives.
//!
//! The crate exposes two ownership models:
//!
//! - [`Pod5Reader`] is a cursor over POD5 reads. It owns the parsed POD5
//!   wrapper metadata and records, and returns one [`Pod5Record`] at a time.
//! - [`Pod5`] is a fully materialized POD5 payload containing the header,
//!   run-info rows, and read records. It is [`Clone`] for deep-copy workflows.
//!
//! POD5 files are a POD5 wrapper around Apache Arrow IPC/Feather V2 tables.
//! This crate validates the POD5 container markers and reads the embedded
//! Arrow Reads and Run Info tables through Arrow's IPC reader.

use arrow_array::builder::{
    FixedSizeBinaryBuilder, LargeBinaryBuilder, ListBuilder, MapBuilder, StringBuilder,
    StringDictionaryBuilder, TimestampMillisecondBuilder, UInt64Builder,
};
use arrow_array::types::Int16Type;
use arrow_array::{
    Array, ArrayRef, BooleanArray, DictionaryArray, FixedSizeBinaryArray, Float32Array, Int16Array,
    LargeBinaryArray, LargeListArray, ListArray, RecordBatch, StringArray,
    TimestampMillisecondArray, UInt8Array, UInt16Array, UInt32Array, UInt64Array,
};
use arrow_ipc::reader::FileReader;
use arrow_ipc::writer::FileWriter;
use arrow_schema::{ArrowError, DataType, Field, Schema, SchemaRef};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Cursor, Read, Write};
use std::path::Path;
use std::sync::Arc;

/// POD5 file signature at the start and end of every POD5 file.
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

type ParsedPod5 = (
    Pod5Header,
    Vec<Pod5RunInfo>,
    Vec<Pod5Signal>,
    Vec<Pod5Record>,
);

/// A fully materialized POD5 file.
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
/// POD5 stores its data as embedded Arrow tables with a footer, so this reader
/// materializes the table metadata during construction and then exposes the
/// same `read_record`/`records` surface as the text-format crates.
pub struct Pod5Reader<R: Read = File> {
    /// POD5 wrapper and table metadata.
    pub header: Pod5Header,
    /// Run Info table rows parsed during construction.
    pub run_infos: Vec<Pod5RunInfo>,
    /// Signal table rows parsed during construction.
    pub signals: Vec<Pod5Signal>,
    records: Vec<Pod5Record>,
    position: usize,
    _reader: std::marker::PhantomData<R>,
}

/// POD5 file writer over any writable byte stream.
///
/// POD5 is a container of Arrow IPC files plus a footer, so this writer emits a
/// complete materialized [`Pod5`] value at a time rather than streaming
/// individual reads.
pub struct Pod5Writer<W: Write = File> {
    writer: W,
}

/// POD5 wrapper and table metadata.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5Header {
    /// Leading POD5 magic bytes.
    pub magic: [u8; 8],
    /// Per-file section marker used between embedded files.
    pub section_marker: [u8; 16],
    /// Embedded Arrow/footer sections in file order.
    pub sections: Vec<Pod5Section>,
    /// `MINKNOW:file_identifier` Arrow schema metadata, when present.
    pub file_identifier: Option<String>,
    /// `MINKNOW:software` Arrow schema metadata, when present.
    pub software: Option<String>,
    /// `MINKNOW:pod5_version` Arrow schema metadata, when present.
    pub pod5_version: Option<String>,
}

/// One POD5 wrapper section.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5Section {
    /// Section kind inferred from the embedded Arrow schema or footer magic.
    pub kind: Pod5SectionKind,
    /// Absolute byte offset of the section payload, just after the section marker.
    pub offset: u64,
    /// Unpadded payload length in bytes.
    pub length: u64,
    /// Padded payload length up to the next section marker.
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

/// A Run Info table row.
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
}

/// A POD5 Reads table row.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5Record {
    /// Read UUID.
    pub read_id: String,
    /// Signal table row indices referenced by this read.
    pub signal_rows: Vec<u64>,
    /// Read number.
    pub read_number: u32,
    /// Samples recorded on this channel since run start to the first sample.
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
    /// Run Info acquisition identifier.
    pub run_info: String,
}

/// A POD5 Signal table row.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5Signal {
    /// Read UUID associated with this signal row.
    pub read_id: String,
    /// Number of decoded ADC samples in this row.
    pub samples: u32,
    /// Signal payload, compressed or uncompressed depending on the table schema.
    pub payload: Pod5SignalPayload,
}

/// Signal payload representation.
#[derive(Debug, Clone, PartialEq)]
pub enum Pod5SignalPayload {
    /// VBZ-compressed signal bytes from the POD5 Signal table.
    Vbz(Vec<u8>),
    /// Uncompressed int16 ADC samples.
    Uncompressed(Vec<i16>),
}

impl Pod5 {
    /// Opens a POD5 file and materializes all records into memory.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Pod5Reader::from_path(path)?.read_all()
    }

    /// Materializes all records from a POD5 byte stream.
    pub fn from_reader<R: Read>(reader: R) -> io::Result<Self> {
        Pod5Reader::from_reader(reader)?.read_all()
    }

    /// Writes this POD5 payload to a filesystem path.
    pub fn to_path<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut writer = Pod5Writer::from_path(path)?;
        writer.write_all(self)?;
        writer.flush()
    }

    /// Writes this POD5 payload to a writable byte stream.
    pub fn to_writer<W: Write>(&self, writer: W) -> io::Result<()> {
        let mut writer = Pod5Writer::from_writer(writer);
        writer.write_all(self)?;
        writer.flush()
    }

    /// Decompresses and concatenates signal rows referenced by `record`.
    pub fn signal_for_record(&self, record: &Pod5Record) -> io::Result<Vec<i16>> {
        record.signal(&self.signals)
    }

    /// Finds a read by ID and returns its decompressed signal, if present.
    pub fn signal_by_read_id(&self, read_id: &str) -> io::Result<Option<Vec<i16>>> {
        self.records
            .iter()
            .find(|record| record.read_id == read_id)
            .map(|record| self.signal_for_record(record))
            .transpose()
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

impl<R: Read> Pod5Reader<R> {
    /// Creates a POD5 reader from a byte stream.
    ///
    /// POD5 requires footer and table access, so the stream is read into memory
    /// during construction.
    pub fn from_reader(mut reader: R) -> io::Result<Self> {
        let mut data = Vec::new();
        reader.read_to_end(&mut data)?;
        let (header, run_infos, signals, records) = parse_pod5(&data)?;

        Ok(Self {
            header,
            run_infos,
            signals,
            records,
            position: 0,
            _reader: std::marker::PhantomData,
        })
    }

    /// Reads the next POD5 record.
    ///
    /// Returns `Ok(None)` when all parsed reads have been returned.
    pub fn read_record(&mut self) -> io::Result<Option<Pod5Record>> {
        let Some(record) = self.records.get(self.position).cloned() else {
            return Ok(None);
        };

        self.position += 1;
        Ok(Some(record))
    }

    /// Reads the next POD5 record.
    ///
    /// This is a compatibility alias for [`Pod5Reader::read_record`].
    pub fn read(&mut self) -> io::Result<Option<Pod5Record>> {
        self.read_record()
    }

    /// Returns an iterator over records in the POD5 file.
    pub fn records(&mut self) -> Pod5Records<'_, R> {
        Pod5Records { reader: self }
    }

    /// Consumes this reader and materializes the remaining POD5 records.
    pub fn read_all(mut self) -> io::Result<Pod5> {
        let mut records = Vec::new();
        while let Some(record) = self.read_record()? {
            records.push(record);
        }

        Ok(Pod5 {
            header: self.header,
            run_infos: self.run_infos,
            signals: self.signals,
            records,
        })
    }
}

impl Pod5Writer<File> {
    /// Creates or truncates a POD5 file at a filesystem path.
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

impl<W: Write> Pod5Writer<W> {
    /// Creates a POD5 writer from a writable byte stream.
    pub fn from_writer(writer: W) -> Self {
        Self { writer }
    }

    /// Writes a complete materialized POD5 payload.
    pub fn write_all(&mut self, pod5: &Pod5) -> io::Result<()> {
        let data = encode_pod5(pod5)?;
        self.writer.write_all(&data)
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
    pub fn signal(&self, signals: &[Pod5Signal]) -> io::Result<Vec<i16>> {
        let total = self
            .signal_rows
            .iter()
            .filter_map(|&index| signals.get(index as usize))
            .map(|signal| signal.samples as usize)
            .sum();
        let mut samples = Vec::with_capacity(total);

        for &index in &self.signal_rows {
            let signal = signals
                .get(index as usize)
                .ok_or_else(|| invalid_data("POD5 signal row index is out of bounds"))?;
            samples.extend(signal.decompress()?);
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
    pub fn compress(&self) -> io::Result<Vec<u8>> {
        compress_vbz_signal(&self.decompress()?)
    }
}

/// Iterator over records from a [`Pod5Reader`].
pub struct Pod5Records<'a, R: Read> {
    reader: &'a mut Pod5Reader<R>,
}

impl<R: Read> Iterator for Pod5Records<'_, R> {
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

fn encode_pod5(pod5: &Pod5) -> io::Result<Vec<u8>> {
    validate_pod5_for_writing(pod5)?;

    let metadata = pod5_writer_metadata(pod5);
    let signal = write_arrow_section(build_signal_batch(&pod5.signals, &metadata)?, &metadata)?;
    let run_info =
        write_arrow_section(build_run_info_batch(&pod5.run_infos, &metadata)?, &metadata)?;
    let reads = write_arrow_section(build_reads_batch(&pod5.records, &metadata)?, &metadata)?;
    let sections = [
        (Pod5SectionKind::Signal, signal),
        (Pod5SectionKind::RunInfo, run_info),
        (Pod5SectionKind::Reads, reads),
    ];

    let marker = pod5_section_marker(pod5, &metadata)?;
    let mut output = Vec::new();
    let mut footer_entries = Vec::with_capacity(sections.len());

    output.extend_from_slice(POD5_MAGIC);
    output.extend_from_slice(&marker);

    for (kind, section) in sections {
        let offset = output.len() as i64;
        let length = section.len() as i64;
        output.extend_from_slice(&section);
        pad_to_8(&mut output);
        output.extend_from_slice(&marker);
        footer_entries.push(Pod5FooterEntry {
            offset,
            length,
            content_type: match kind {
                Pod5SectionKind::Reads => 0,
                Pod5SectionKind::Signal => 1,
                Pod5SectionKind::RunInfo => 4,
                Pod5SectionKind::Unknown | Pod5SectionKind::Footer => {
                    return Err(invalid_data("POD5 writer cannot encode this section kind"));
                }
            },
        });
    }

    let footer = build_pod5_footer(&metadata, &footer_entries);
    output.extend_from_slice(POD5_FOOTER_MAGIC);
    let footer_payload_start = output.len();
    output.extend_from_slice(&footer);
    pad_to_8(&mut output);
    let footer_len = (output.len() - footer_payload_start) as u64;
    output.extend_from_slice(&footer_len.to_le_bytes());
    output.extend_from_slice(&marker);
    output.extend_from_slice(POD5_MAGIC);

    Ok(output)
}

fn validate_pod5_for_writing(pod5: &Pod5) -> io::Result<()> {
    if pod5.signals.is_empty() {
        return Err(invalid_data("POD5 writer requires at least one signal row"));
    }
    if pod5.run_infos.is_empty() {
        return Err(invalid_data(
            "POD5 writer requires at least one run info row",
        ));
    }

    for signal in &pod5.signals {
        uuid_string_to_bytes(&signal.read_id)?;
    }

    for record in &pod5.records {
        uuid_string_to_bytes(&record.read_id)?;
        let mut total = 0u64;
        for &index in &record.signal_rows {
            let signal = pod5
                .signals
                .get(index as usize)
                .ok_or_else(|| invalid_data("POD5 signal row index is out of bounds"))?;
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

fn build_signal_batch(
    signals: &[Pod5Signal],
    metadata: &Pod5WriterMetadata,
) -> io::Result<RecordBatch> {
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
    let schema = Schema::new_with_metadata(
        vec![
            uuid_field("read_id"),
            vbz_field("signal"),
            Field::new("samples", DataType::UInt32, false),
        ],
        arrow_metadata(metadata),
    );

    RecordBatch::try_new(Arc::new(schema), vec![read_id, signal, samples]).map_err(arrow_error)
}

fn build_run_info_batch(
    run_infos: &[Pod5RunInfo],
    metadata: &Pod5WriterMetadata,
) -> io::Result<RecordBatch> {
    let rows = run_infos.len();
    let context_tags = Arc::new(empty_string_map_array(rows)?) as ArrayRef;
    let tracking_id = Arc::new(empty_string_map_array(rows)?) as ArrayRef;
    let acquisition_start_time = Arc::new(timestamp_utc_zeros(rows)) as ArrayRef;
    let protocol_start_time = Arc::new(timestamp_utc_zeros(rows)) as ArrayRef;

    let arrays: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from(
            run_infos
                .iter()
                .map(|run_info| run_info.acquisition_id.clone())
                .collect::<Vec<_>>(),
        )),
        acquisition_start_time.clone(),
        Arc::new(Int16Array::from(vec![0i16; rows])),
        Arc::new(Int16Array::from(vec![0i16; rows])),
        context_tags.clone(),
        Arc::new(StringArray::from(
            run_infos
                .iter()
                .map(|run_info| run_info.experiment_name.clone())
                .collect::<Vec<_>>(),
        )),
        Arc::new(StringArray::from(
            run_infos
                .iter()
                .map(|run_info| run_info.flow_cell_id.clone())
                .collect::<Vec<_>>(),
        )),
        Arc::new(StringArray::from(vec![String::new(); rows])),
        Arc::new(StringArray::from(vec![String::new(); rows])),
        Arc::new(StringArray::from(vec![String::new(); rows])),
        protocol_start_time.clone(),
        Arc::new(StringArray::from(
            run_infos
                .iter()
                .map(|run_info| run_info.sample_id.clone())
                .collect::<Vec<_>>(),
        )),
        Arc::new(UInt16Array::from(
            run_infos
                .iter()
                .map(|run_info| run_info.sample_rate)
                .collect::<Vec<_>>(),
        )),
        Arc::new(StringArray::from(
            run_infos
                .iter()
                .map(|run_info| run_info.sequencing_kit.clone())
                .collect::<Vec<_>>(),
        )),
        Arc::new(StringArray::from(vec![String::new(); rows])),
        Arc::new(StringArray::from(vec![String::new(); rows])),
        Arc::new(StringArray::from(
            run_infos
                .iter()
                .map(|run_info| run_info.software.clone())
                .collect::<Vec<_>>(),
        )),
        Arc::new(StringArray::from(vec![String::new(); rows])),
        Arc::new(StringArray::from(vec![String::new(); rows])),
        tracking_id.clone(),
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
        pore_type.append_value(&record.pore_type);
        end_reason.append_value(&record.end_reason);
        run_info.append_value(&record.run_info);
    }

    let arrays: Vec<ArrayRef> = vec![
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

    let fields = READS_TABLE_FIELDS
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

fn empty_string_map_array(rows: usize) -> io::Result<arrow_array::MapArray> {
    let mut builder = MapBuilder::new(None, StringBuilder::new(), StringBuilder::new());
    for _ in 0..rows {
        builder.append(true).map_err(arrow_error)?;
    }
    Ok(builder.finish())
}

fn timestamp_utc_zeros(rows: usize) -> TimestampMillisecondArray {
    let mut builder = TimestampMillisecondBuilder::new().with_timezone("UTC");
    for _ in 0..rows {
        builder.append_value(0);
    }
    builder.finish()
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

fn pad_to_8(data: &mut Vec<u8>) {
    let padding = (8 - data.len() % 8) % 8;
    data.extend(std::iter::repeat_n(0, padding));
}

fn parse_pod5(data: &[u8]) -> io::Result<ParsedPod5> {
    if data.len() < 8 + 16 + 16 + 8 {
        return Err(invalid_data("POD5 file is too short"));
    }

    let magic = slice_to_array::<8>(
        data.get(..8)
            .ok_or_else(|| invalid_data("POD5 magic is missing"))?,
    );

    if &magic != POD5_MAGIC {
        return Err(invalid_data("invalid POD5 leading magic"));
    }

    if data.get(data.len() - 8..).unwrap() != POD5_MAGIC {
        return Err(invalid_data("invalid POD5 trailing magic"));
    }

    let section_marker = slice_to_array::<16>(
        data.get(8..24)
            .ok_or_else(|| invalid_data("POD5 section marker is missing"))?,
    );
    let final_marker_offset = data.len() - 8 - 16;
    if data.get(final_marker_offset..final_marker_offset + 16) != Some(&section_marker) {
        return Err(invalid_data("POD5 final section marker is missing"));
    }

    let marker_positions = find_section_markers(data, &section_marker);
    if marker_positions.len() < 3 {
        return Err(invalid_data("POD5 file has too few section markers"));
    }

    if marker_positions[0] != 8 || *marker_positions.last().unwrap() != final_marker_offset {
        return Err(invalid_data("POD5 section markers are malformed"));
    }

    let mut sections = Vec::new();
    let mut run_infos = Vec::new();
    let mut signals = Vec::new();
    let mut records = Vec::new();
    let mut metadata = Pod5Metadata::default();
    let mut seen_footer = false;

    for window in marker_positions.windows(2) {
        let section_start = window[0] + 16;
        let section_end = window[1];

        if section_start > section_end {
            return Err(invalid_data("POD5 section offsets are malformed"));
        }

        let section = &data[section_start..section_end];
        if section.starts_with(POD5_FOOTER_MAGIC) {
            let footer_len = parse_footer_section(section)?;
            sections.push(Pod5Section {
                kind: Pod5SectionKind::Footer,
                offset: section_start as u64,
                length: footer_len as u64,
                padded_length: section.len() as u64,
                row_count: 0,
            });
            seen_footer = true;
            continue;
        }

        let arrow_len = arrow_payload_len(section)?;
        let payload = &section[..arrow_len];
        let mut reader = FileReader::try_new(Cursor::new(payload), None).map_err(arrow_error)?;
        update_metadata(reader.custom_metadata(), &mut metadata)?;
        let schema = reader.schema();
        let kind = infer_section_kind(&schema);
        let mut row_count = 0usize;

        for batch in &mut reader {
            let batch = batch.map_err(arrow_error)?;
            row_count += batch.num_rows();

            match kind {
                Pod5SectionKind::Reads => records.extend(parse_reads_batch(&batch)?),
                Pod5SectionKind::Signal => signals.extend(parse_signal_batch(&batch)?),
                Pod5SectionKind::RunInfo => run_infos.extend(parse_run_info_batch(&batch)?),
                Pod5SectionKind::Unknown | Pod5SectionKind::Footer => {}
            }
        }

        sections.push(Pod5Section {
            kind,
            offset: section_start as u64,
            length: arrow_len as u64,
            padded_length: section.len() as u64,
            row_count,
        });
    }

    if !seen_footer {
        return Err(invalid_data("POD5 footer section is missing"));
    }

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

    Ok((header, run_infos, signals, records))
}

fn find_section_markers(data: &[u8], section_marker: &[u8; 16]) -> Vec<usize> {
    let mut positions = Vec::new();
    let mut offset = 0usize;

    while let Some(position) = data[offset..]
        .windows(section_marker.len())
        .position(|window| window == section_marker)
    {
        let absolute = offset + position;
        positions.push(absolute);
        offset = absolute + section_marker.len();
    }

    positions
}

fn parse_footer_section(section: &[u8]) -> io::Result<usize> {
    if section.len() < POD5_FOOTER_MAGIC.len() + 8 {
        return Err(invalid_data("POD5 footer section is too short"));
    }

    let footer_len_offset = section.len() - 8;
    let footer_len = u64::from_le_bytes(
        section[footer_len_offset..]
            .try_into()
            .map_err(|_| invalid_data("POD5 footer length is malformed"))?,
    ) as usize;
    let footer_start = POD5_FOOTER_MAGIC.len();
    let footer_end = footer_start
        .checked_add(footer_len)
        .ok_or_else(|| invalid_data("POD5 footer length overflow"))?;

    if footer_end > footer_len_offset {
        return Err(invalid_data("POD5 footer length exceeds section"));
    }

    if !section[footer_end..footer_len_offset]
        .iter()
        .all(|byte| *byte == 0)
    {
        return Err(invalid_data("POD5 footer padding is not zeroed"));
    }

    Ok(footer_len)
}

fn arrow_payload_len(section: &[u8]) -> io::Result<usize> {
    if !section.starts_with(ARROW_MAGIC) {
        return Err(invalid_data("POD5 section is neither Arrow nor footer"));
    }

    let end_magic_offset = section
        .windows(ARROW_MAGIC.len())
        .rposition(|window| window == ARROW_MAGIC)
        .ok_or_else(|| invalid_data("embedded Arrow file is missing trailing magic"))?;
    let arrow_len = end_magic_offset + ARROW_MAGIC.len();

    if arrow_len == ARROW_MAGIC.len() {
        return Err(invalid_data("embedded Arrow file has no payload"));
    }

    if !section[arrow_len..].iter().all(|byte| *byte == 0) {
        return Err(invalid_data("POD5 Arrow section padding is not zeroed"));
    }

    Ok(arrow_len)
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
/// The local transform is:
///
/// 1. first-order delta encode with wrapping `i16` subtraction,
/// 2. zigzag encode signed deltas into `u16`,
/// 3. SVB16 encode with one control bit per value, LSB first,
/// 4. zstd-compress the SVB16 byte stream at the POD5-compatible level.
pub fn compress_vbz_signal(samples: &[i16]) -> io::Result<Vec<u8>> {
    if samples.is_empty() {
        return Ok(Vec::new());
    }

    let inner = encode_vbz_inner(samples);
    zstd::bulk::compress(&inner, 1)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))
}

/// Decompresses a POD5 VBZ signal blob into raw int16 ADC samples.
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
    let inner = zstd::bulk::decompress(data, max_inner_len)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))?;

    decode_vbz_inner(&inner, num_samples)
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

fn invalid_data(message: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

#[cfg(test)]
mod pod5_tests {
    use super::*;

    const A_100_POD5: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/A_100.pod5");
    const RECORD_COUNT: usize = 100;
    const FIRST_READ_ID: &str = "1cadb1e9-592f-4e22-9285-4626f2b7da9f";
    const LAST_READ_ID: &str = "6ae6c2e9-fe1a-4b1d-befb-ac98a5a16c9a";
    const RUN_ID: &str = "1f03c20c2da347cc99b58b232181a7464126f4cb";

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
        assert_eq!(reader.signals.len(), RECORD_COUNT);

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

        let signal = first.signal(&reader.signals).unwrap();
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

    #[test]
    fn writer_rejects_invalid_signal_row_references() {
        let mut pod5 = Pod5::from_path(A_100_POD5).unwrap();
        pod5.records[0].signal_rows = vec![pod5.signals.len() as u64];
        let mut output = Vec::new();

        assert!(pod5.to_writer(&mut output).is_err());
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
