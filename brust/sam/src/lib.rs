//! Minimal SAM parsing primitives.
//!
//! The crate exposes two ownership models:
//!
//! - [`SamReader`] is a streaming parser over a SAM byte stream. It owns
//!   reader state, parses the optional header up front, and returns one
//!   [`SamRecord`] at a time.
//! - [`Sam`] is a fully materialized SAM payload containing the header and all
//!   alignment records. It is [`Clone`] for deep-copy workflows.
//!
//! SAM alignment lines are tab-delimited and contain the 11 mandatory fields
//! defined by SAM v1.6, followed by zero or more optional `TAG:TYPE:VALUE`
//! fields.

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

/// A fully materialized SAM file.
///
/// `Sam` owns its header and all alignment records, so cloning this type
/// performs a deep copy of the parsed SAM data.
#[derive(Debug, Clone, PartialEq)]
pub struct Sam {
    /// Optional SAM header section.
    pub header: SamHeader,
    /// All alignment records loaded from the stream.
    pub records: Vec<SamRecord>,
}

/// Streaming SAM parser over any readable byte stream.
///
/// `SamReader` owns the buffered reader and parser cursor, so it is not
/// cloneable. Use [`SamReader::read_all`] or [`Sam::from_path`] when a
/// deep-copyable, materialized representation is needed.
pub struct SamReader<R: Read = File> {
    /// Header records parsed before the first alignment line.
    pub header: SamHeader,
    reader: BufReader<R>,
    // Alignment line already read while parsing the header.
    pending_record: Option<String>,
    // One-based line counter maintained for diagnostics.
    line_number: usize,
}

/// Streaming SAM writer over any writable byte stream.
///
/// `SamWriter` emits header records followed by alignment records. Use
/// [`SamWriter::write_all`] for a materialized [`Sam`] value, or call
/// [`SamWriter::write_header`] once before streaming records manually.
pub struct SamWriter<W: Write = File> {
    writer: W,
}

/// SAM header section.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct SamHeader {
    /// All header records in file order.
    pub records: Vec<SamHeaderRecord>,
}

/// A single SAM header record.
#[derive(Debug, Clone, PartialEq)]
pub struct SamHeaderRecord {
    /// Two-character record type without the leading `@`, such as `HD` or `SQ`.
    pub record_type: String,
    /// `TAG:VALUE` fields for all non-comment header records.
    pub fields: Vec<SamHeaderField>,
    /// Comment text for `@CO` header records.
    pub comment: Option<String>,
}

/// A `TAG:VALUE` field from a SAM header record.
#[derive(Debug, Clone, PartialEq)]
pub struct SamHeaderField {
    /// Two-character header field tag.
    pub tag: String,
    /// Header field value.
    pub value: String,
}

/// A SAM alignment record.
#[derive(Debug, Clone, PartialEq)]
pub struct SamRecord {
    /// Query template name.
    pub qname: String,
    /// Bitwise SAM flag.
    pub flag: u16,
    /// Reference sequence name.
    pub rname: String,
    /// One-based leftmost mapping position, or `0` when unavailable.
    pub pos: u32,
    /// Mapping quality, with `255` representing unavailable.
    pub mapq: u8,
    /// CIGAR string, or `*` when unavailable.
    pub cigar: String,
    /// Reference name of the mate/next read.
    pub rnext: String,
    /// One-based position of the mate/next read, or `0` when unavailable.
    pub pnext: u32,
    /// Observed template length.
    pub tlen: i32,
    /// Segment sequence, or `*` when not stored.
    pub seq: String,
    /// Base quality string, or `*` when not stored.
    pub qual: String,
    /// Optional `TAG:TYPE:VALUE` fields.
    pub optional: Vec<SamOptionalField>,
}

/// Optional `TAG:TYPE:VALUE` field from a SAM alignment record.
#[derive(Debug, Clone, PartialEq)]
pub struct SamOptionalField {
    /// Two-character optional field tag.
    pub tag: String,
    /// Parsed optional value.
    pub value: SamOptionalValue,
}

/// Parsed SAM optional field value.
#[derive(Debug, Clone, PartialEq)]
pub enum SamOptionalValue {
    /// Printable character (`A`).
    Character(char),
    /// Signed integer (`i`).
    Integer(i64),
    /// Single-precision floating point (`f`).
    Float(f32),
    /// Printable string (`Z`).
    String(String),
    /// Hex-encoded byte array (`H`).
    Hex(Vec<u8>),
    /// Integer or floating point array (`B`).
    Array(SamOptionalArray),
}

/// Parsed SAM optional array value.
#[derive(Debug, Clone, PartialEq)]
pub enum SamOptionalArray {
    /// Signed 8-bit integer array (`B:c`).
    Int8(Vec<i8>),
    /// Unsigned 8-bit integer array (`B:C`).
    UInt8(Vec<u8>),
    /// Signed 16-bit integer array (`B:s`).
    Int16(Vec<i16>),
    /// Unsigned 16-bit integer array (`B:S`).
    UInt16(Vec<u16>),
    /// Signed 32-bit integer array (`B:i`).
    Int32(Vec<i32>),
    /// Unsigned 32-bit integer array (`B:I`).
    UInt32(Vec<u32>),
    /// Single-precision floating point array (`B:f`).
    Float(Vec<f32>),
}

/// One CIGAR operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SamCigarOp {
    /// Operation length.
    pub len: u32,
    /// Operation code, one of `MIDNSHP=X`.
    pub op: char,
}

impl SamRecord {
    /// Creates a SAM alignment record from parsed fields.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        qname: String,
        flag: u16,
        rname: String,
        pos: u32,
        mapq: u8,
        cigar: String,
        rnext: String,
        pnext: u32,
        tlen: i32,
        seq: String,
        qual: String,
        optional: Vec<SamOptionalField>,
    ) -> Self {
        Self {
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext,
            pnext,
            tlen,
            seq,
            qual,
            optional,
        }
    }

    /// Returns the query/read name.
    pub fn read_name(&self) -> &str {
        &self.qname
    }

    /// Returns `true` when the segment-unmapped flag (`0x4`) is set.
    pub fn is_unmapped(&self) -> bool {
        self.flag & flags::UNMAPPED != 0
    }

    /// Returns `true` when the reverse-complemented flag (`0x10`) is set.
    pub fn is_reverse_complemented(&self) -> bool {
        self.flag & flags::REVERSE_COMPLEMENTED != 0
    }

    /// Returns the optional field matching `tag`, if present.
    pub fn aux(&self, tag: &str) -> Option<&SamOptionalValue> {
        self.optional
            .iter()
            .find(|field| field.tag == tag)
            .map(|field| &field.value)
    }

    /// Parses this record's CIGAR string into operations.
    ///
    /// Returns an empty vector when the CIGAR field is `*`.
    pub fn cigar_ops(&self) -> io::Result<Vec<SamCigarOp>> {
        parse_cigar(&self.cigar)
    }

    /// Returns the query-consuming CIGAR length.
    pub fn query_len_from_cigar(&self) -> io::Result<u64> {
        Ok(self
            .cigar_ops()?
            .iter()
            .filter(|op| matches!(op.op, 'M' | 'I' | 'S' | '=' | 'X'))
            .map(|op| u64::from(op.len))
            .sum())
    }

    /// Returns the reference-consuming CIGAR length.
    pub fn reference_len_from_cigar(&self) -> io::Result<u64> {
        Ok(self
            .cigar_ops()?
            .iter()
            .filter(|op| matches!(op.op, 'M' | 'D' | 'N' | '=' | 'X'))
            .map(|op| u64::from(op.len))
            .sum())
    }
}

impl SamOptionalValue {
    /// Returns the SAM type code for this value.
    pub fn value_type(&self) -> char {
        match self {
            Self::Character(_) => 'A',
            Self::Integer(_) => 'i',
            Self::Float(_) => 'f',
            Self::String(_) => 'Z',
            Self::Hex(_) => 'H',
            Self::Array(_) => 'B',
        }
    }
}

impl Sam {
    /// Opens a SAM file and materializes all records into memory.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        SamReader::from_path(path)?.read_all()
    }

    /// Materializes all records from a SAM byte stream.
    pub fn from_reader<R: Read>(reader: R) -> io::Result<Self> {
        SamReader::from_reader(reader)?.read_all()
    }

    /// Writes this SAM payload to a filesystem path.
    pub fn to_path<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut writer = SamWriter::from_path(path)?;
        writer.write_all(self)?;
        writer.flush()
    }

    /// Writes this SAM payload to a writable byte stream.
    pub fn to_writer<W: Write>(&self, writer: W) -> io::Result<()> {
        let mut writer = SamWriter::from_writer(writer);
        writer.write_all(self)?;
        writer.flush()
    }
}

impl SamReader<File> {
    /// Opens a SAM file from a filesystem path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        Self::from_reader(file)
    }

    /// Opens a SAM file from a filesystem path.
    ///
    /// This is a convenience alias for [`SamReader::from_path`].
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_path(path)
    }
}

impl<R: Read> SamReader<R> {
    /// Creates a streaming parser from a SAM byte stream.
    ///
    /// Header records are parsed immediately and are available in
    /// [`SamReader::header`]. The first alignment line is held for the first
    /// call to [`SamReader::read_record`].
    pub fn from_reader(reader: R) -> io::Result<Self> {
        let mut reader = BufReader::new(reader);
        let mut records = Vec::new();
        let mut pending_record = None;
        let mut line_number = 0;
        let mut line = String::new();

        loop {
            line.clear();
            let bytes = reader.read_line(&mut line)?;

            if bytes == 0 {
                break;
            }

            line_number += 1;

            if line.starts_with('@') {
                records.push(SamHeaderRecord::parse(line.trim_end_matches(['\r', '\n']))?);
            } else {
                pending_record = Some(line.clone());
                break;
            }
        }

        let header = SamHeader { records };
        header.validate()?;

        Ok(Self {
            header,
            reader,
            pending_record,
            line_number,
        })
    }

    /// Reads the next SAM alignment record.
    ///
    /// Returns `Ok(None)` when the stream is at EOF before another alignment
    /// line begins. Header records are consumed during construction and are not
    /// returned by this method.
    pub fn read_record(&mut self) -> io::Result<Option<SamRecord>> {
        let line = match self.pending_record.take() {
            Some(line) => line,
            None => {
                let mut line = String::new();
                let bytes = self.reader.read_line(&mut line)?;

                if bytes == 0 {
                    return Ok(None);
                }

                self.line_number += 1;
                line
            }
        };

        if line.starts_with('@') {
            return Err(invalid_data("header record found after alignment section"));
        }

        SamRecord::parse(line.trim_end_matches(['\r', '\n'])).map(Some)
    }

    /// Reads the next SAM alignment record.
    ///
    /// This is a compatibility alias for [`SamReader::read_record`].
    pub fn read(&mut self) -> io::Result<Option<SamRecord>> {
        self.read_record()
    }

    /// Returns an iterator over alignment records in the stream.
    pub fn records(&mut self) -> SamRecords<'_, R> {
        SamRecords { reader: self }
    }

    /// Consumes this reader and materializes the entire SAM stream.
    pub fn read_all(mut self) -> io::Result<Sam> {
        let mut records = Vec::new();
        while let Some(record) = self.read_record()? {
            records.push(record);
        }

        Ok(Sam {
            header: self.header,
            records,
        })
    }
}

impl SamWriter<File> {
    /// Creates or truncates a SAM file at a filesystem path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self::from_writer(file))
    }

    /// Creates or truncates a SAM file at a filesystem path.
    ///
    /// This is a convenience alias for [`SamWriter::from_path`].
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_path(path)
    }
}

impl<W: Write> SamWriter<W> {
    /// Creates a SAM writer from a writable byte stream.
    pub fn from_writer(writer: W) -> Self {
        Self { writer }
    }

    /// Writes a SAM header section.
    pub fn write_header(&mut self, header: &SamHeader) -> io::Result<()> {
        header.validate()?;
        for record in &header.records {
            let line = format_sam_header_record(record)?;
            self.writer.write_all(line.as_bytes())?;
            self.writer.write_all(b"\n")?;
        }
        Ok(())
    }

    /// Writes one SAM alignment record.
    pub fn write_record(&mut self, record: &SamRecord) -> io::Result<()> {
        let line = format_sam_record(record)?;
        SamRecord::parse(&line)?;
        self.writer.write_all(line.as_bytes())?;
        self.writer.write_all(b"\n")?;
        Ok(())
    }

    /// Writes one SAM alignment record.
    ///
    /// This is a compatibility alias for [`SamWriter::write_record`].
    pub fn write(&mut self, record: &SamRecord) -> io::Result<()> {
        self.write_record(record)
    }

    /// Writes a materialized SAM payload.
    pub fn write_all(&mut self, sam: &Sam) -> io::Result<()> {
        self.write_header(&sam.header)?;
        for record in &sam.records {
            self.write_record(record)?;
        }
        Ok(())
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

/// Iterator over records from a [`SamReader`].
pub struct SamRecords<'a, R: Read> {
    reader: &'a mut SamReader<R>,
}

impl<R: Read> Iterator for SamRecords<'_, R> {
    type Item = io::Result<SamRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(error) => Some(Err(error)),
        }
    }
}

impl SamHeader {
    /// Returns the first header record of `record_type`, if present.
    pub fn first(&self, record_type: &str) -> Option<&SamHeaderRecord> {
        self.records
            .iter()
            .find(|record| record.record_type == record_type)
    }

    /// Returns all header records of `record_type`.
    pub fn records_of_type<'a>(
        &'a self,
        record_type: &'a str,
    ) -> impl Iterator<Item = &'a SamHeaderRecord> {
        self.records
            .iter()
            .filter(move |record| record.record_type == record_type)
    }

    /// Returns `@SQ` sequence names in header order.
    pub fn sequence_names(&self) -> Vec<&str> {
        self.records_of_type("SQ")
            .filter_map(|record| record.value("SN"))
            .collect()
    }

    fn validate(&self) -> io::Result<()> {
        let mut hd_count = 0usize;
        let mut sq_names = HashSet::new();
        let mut rg_ids = HashSet::new();
        let mut pg_ids = HashSet::new();

        for (index, record) in self.records.iter().enumerate() {
            match record.record_type.as_str() {
                "HD" => {
                    hd_count += 1;
                    if index != 0 {
                        return Err(invalid_data("@HD must be the first header record"));
                    }
                    if hd_count > 1 {
                        return Err(invalid_data("only one @HD header record is allowed"));
                    }
                }
                "SQ" => {
                    let name = record
                        .value("SN")
                        .ok_or_else(|| invalid_data("@SQ missing SN"))?;
                    let len = record
                        .value("LN")
                        .ok_or_else(|| invalid_data("@SQ missing LN"))?;
                    let len = len
                        .parse::<u32>()
                        .map_err(|_| invalid_data("invalid @SQ LN"))?;
                    if len == 0 || len > i32::MAX as u32 {
                        return Err(invalid_data("@SQ LN out of range"));
                    }
                    if !is_reference_name(name) {
                        return Err(invalid_data("invalid @SQ SN reference name"));
                    }
                    if !sq_names.insert(name.to_string()) {
                        return Err(invalid_data("duplicate @SQ SN"));
                    }
                }
                "RG" => {
                    let id = record
                        .value("ID")
                        .ok_or_else(|| invalid_data("@RG missing ID"))?;
                    if !rg_ids.insert(id.to_string()) {
                        return Err(invalid_data("duplicate @RG ID"));
                    }
                }
                "PG" => {
                    let id = record
                        .value("ID")
                        .ok_or_else(|| invalid_data("@PG missing ID"))?;
                    if !pg_ids.insert(id.to_string()) {
                        return Err(invalid_data("duplicate @PG ID"));
                    }
                }
                "CO" => {}
                _ => return Err(invalid_data("unknown SAM header record type")),
            }
        }

        Ok(())
    }
}

impl SamHeaderRecord {
    /// Creates a non-comment header record from parsed fields.
    pub fn new(record_type: String, fields: Vec<SamHeaderField>) -> Self {
        Self {
            record_type,
            fields,
            comment: None,
        }
    }

    /// Creates an `@CO` comment header record.
    pub fn comment(comment: String) -> Self {
        Self {
            record_type: "CO".to_string(),
            fields: Vec::new(),
            comment: Some(comment),
        }
    }

    /// Returns the value of the first field matching `tag`, if present.
    pub fn value(&self, tag: &str) -> Option<&str> {
        self.fields
            .iter()
            .find(|field| field.tag == tag)
            .map(|field| field.value.as_str())
    }

    fn parse(line: &str) -> io::Result<Self> {
        let body = line
            .strip_prefix('@')
            .ok_or_else(|| invalid_data("SAM header line must start with @"))?;
        let mut parts = body.split('\t');
        let record_type = parts
            .next()
            .ok_or_else(|| invalid_data("SAM header record type is missing"))?;

        if !matches!(record_type, "HD" | "SQ" | "RG" | "PG" | "CO") {
            return Err(invalid_data("unknown SAM header record type"));
        }

        if record_type == "CO" {
            let comment = body
                .strip_prefix("CO\t")
                .map(str::to_string)
                .unwrap_or_default();
            return Ok(Self::comment(comment));
        }

        let mut seen = HashSet::new();
        let mut fields = Vec::new();

        for field in parts {
            let (tag, value) = field
                .split_once(':')
                .ok_or_else(|| invalid_data("header field must be TAG:VALUE"))?;

            if !is_header_tag(tag) {
                return Err(invalid_data("invalid SAM header field tag"));
            }

            if value.is_empty() || !is_printable_or_space_ascii(value) {
                return Err(invalid_data("invalid SAM header field value"));
            }

            if !seen.insert(tag.to_string()) {
                return Err(invalid_data("duplicate SAM header field tag"));
            }

            fields.push(SamHeaderField {
                tag: tag.to_string(),
                value: value.to_string(),
            });
        }

        if fields.is_empty() {
            return Err(invalid_data("SAM header record has no fields"));
        }

        Ok(Self::new(record_type.to_string(), fields))
    }
}

impl SamHeaderField {
    /// Creates a SAM header field from parsed components.
    pub fn new(tag: String, value: String) -> Self {
        Self { tag, value }
    }
}

impl SamRecord {
    fn parse(line: &str) -> io::Result<Self> {
        if line.is_empty() {
            return Err(invalid_data("empty SAM alignment line"));
        }

        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 11 {
            return Err(invalid_data("SAM alignment line has fewer than 11 fields"));
        }

        let qname = fields[0];
        if !is_query_name(qname) {
            return Err(invalid_data("invalid QNAME"));
        }

        let flag = parse_u16(fields[1], "invalid FLAG")?;
        let rname = fields[2];
        if rname != "*" && !is_reference_name(rname) {
            return Err(invalid_data("invalid RNAME"));
        }

        let pos = parse_u32_max_i32(fields[3], "invalid POS")?;
        let mapq = parse_u8(fields[4], "invalid MAPQ")?;
        let cigar = fields[5];
        let cigar_ops = parse_cigar(cigar)?;
        let rnext = fields[6];
        if !(rnext == "*" || rnext == "=" || is_reference_name(rnext)) {
            return Err(invalid_data("invalid RNEXT"));
        }

        let pnext = parse_u32_max_i32(fields[7], "invalid PNEXT")?;
        let tlen = parse_i32(fields[8], "invalid TLEN")?;
        let seq = fields[9];
        let qual = fields[10];

        if !is_sequence(seq) {
            return Err(invalid_data("invalid SEQ"));
        }

        if !is_quality(qual) {
            return Err(invalid_data("invalid QUAL"));
        }

        if cigar != "*" && seq != "*" {
            let query_len = cigar_ops
                .iter()
                .filter(|op| matches!(op.op, 'M' | 'I' | 'S' | '=' | 'X'))
                .map(|op| u64::from(op.len))
                .sum::<u64>();

            if query_len != seq.len() as u64 {
                return Err(invalid_data("SEQ length does not match CIGAR query length"));
            }
        }

        if qual != "*" && seq == "*" {
            return Err(invalid_data("QUAL is present while SEQ is unavailable"));
        }

        if qual != "*" && qual.len() != seq.len() {
            return Err(invalid_data("QUAL length does not match SEQ length"));
        }

        let mut optional = Vec::new();
        let mut seen_optional = HashSet::new();

        for field in fields.iter().skip(11) {
            let optional_field = SamOptionalField::parse(field)?;
            if !seen_optional.insert(optional_field.tag.clone()) {
                return Err(invalid_data("duplicate SAM optional field tag"));
            }
            optional.push(optional_field);
        }

        Ok(Self::new(
            qname.to_string(),
            flag,
            rname.to_string(),
            pos,
            mapq,
            cigar.to_string(),
            rnext.to_string(),
            pnext,
            tlen,
            seq.to_string(),
            qual.to_string(),
            optional,
        ))
    }
}

impl SamOptionalField {
    /// Creates a SAM optional field from parsed components.
    pub fn new(tag: String, value: SamOptionalValue) -> Self {
        Self { tag, value }
    }

    fn parse(field: &str) -> io::Result<Self> {
        let mut parts = field.splitn(3, ':');
        let tag = parts
            .next()
            .ok_or_else(|| invalid_data("optional field tag missing"))?;
        let value_type = parts
            .next()
            .ok_or_else(|| invalid_data("optional field type missing"))?;
        let value = parts
            .next()
            .ok_or_else(|| invalid_data("optional field value missing"))?;

        if !is_optional_tag(tag) {
            return Err(invalid_data("invalid SAM optional field tag"));
        }

        let value_type = value_type
            .chars()
            .next()
            .filter(|_| value_type.len() == 1)
            .ok_or_else(|| invalid_data("invalid SAM optional field type"))?;

        let value = match value_type {
            'A' => {
                let mut chars = value.chars();
                let character = chars
                    .next()
                    .ok_or_else(|| invalid_data("empty A optional value"))?;
                if chars.next().is_some() || !is_printable_ascii_char(character) {
                    return Err(invalid_data("invalid A optional value"));
                }
                SamOptionalValue::Character(character)
            }
            'i' => SamOptionalValue::Integer(parse_i64(value, "invalid i optional value")?),
            'f' => SamOptionalValue::Float(
                value
                    .parse::<f32>()
                    .map_err(|_| invalid_data("invalid f optional value"))?,
            ),
            'Z' => {
                if !is_printable_or_space_ascii(value) {
                    return Err(invalid_data("invalid Z optional value"));
                }
                SamOptionalValue::String(value.to_string())
            }
            'H' => SamOptionalValue::Hex(parse_hex(value)?),
            'B' => SamOptionalValue::Array(parse_optional_array(value)?),
            _ => return Err(invalid_data("invalid SAM optional field type")),
        };

        Ok(Self::new(tag.to_string(), value))
    }
}

/// SAM bitwise flag constants.
pub mod flags {
    /// Template has multiple segments in sequencing.
    pub const MULTIPLE_SEGMENTS: u16 = 0x1;
    /// Each segment is properly aligned according to the aligner.
    pub const PROPERLY_ALIGNED: u16 = 0x2;
    /// Segment is unmapped.
    pub const UNMAPPED: u16 = 0x4;
    /// Next segment in template is unmapped.
    pub const NEXT_UNMAPPED: u16 = 0x8;
    /// Sequence is reverse complemented.
    pub const REVERSE_COMPLEMENTED: u16 = 0x10;
    /// Sequence of next segment is reverse complemented.
    pub const NEXT_REVERSE_COMPLEMENTED: u16 = 0x20;
    /// First segment in template.
    pub const FIRST_SEGMENT: u16 = 0x40;
    /// Last segment in template.
    pub const LAST_SEGMENT: u16 = 0x80;
    /// Secondary alignment.
    pub const SECONDARY: u16 = 0x100;
    /// Not passing filters.
    pub const FILTERED: u16 = 0x200;
    /// PCR or optical duplicate.
    pub const DUPLICATE: u16 = 0x400;
    /// Supplementary alignment.
    pub const SUPPLEMENTARY: u16 = 0x800;
}

fn parse_cigar(cigar: &str) -> io::Result<Vec<SamCigarOp>> {
    if cigar == "*" {
        return Ok(Vec::new());
    }

    if cigar.is_empty() {
        return Err(invalid_data("empty CIGAR"));
    }

    let mut ops = Vec::new();
    let mut len = 0u32;
    let mut saw_digit = false;

    for byte in cigar.bytes() {
        match byte {
            b'0'..=b'9' => {
                saw_digit = true;
                len = len
                    .checked_mul(10)
                    .and_then(|value| value.checked_add(u32::from(byte - b'0')))
                    .ok_or_else(|| invalid_data("CIGAR operation length overflow"))?;
            }
            b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'=' | b'X' => {
                if !saw_digit || len == 0 {
                    return Err(invalid_data("invalid CIGAR operation length"));
                }
                ops.push(SamCigarOp {
                    len,
                    op: char::from(byte),
                });
                len = 0;
                saw_digit = false;
            }
            _ => return Err(invalid_data("invalid CIGAR operation")),
        }
    }

    if saw_digit {
        return Err(invalid_data("CIGAR ends with length but no operation"));
    }

    if ops.is_empty() {
        return Err(invalid_data("CIGAR contains no operations"));
    }

    validate_cigar_clipping(&ops)?;
    Ok(ops)
}

fn validate_cigar_clipping(ops: &[SamCigarOp]) -> io::Result<()> {
    for (index, op) in ops.iter().enumerate() {
        if op.op == 'H' && index != 0 && index != ops.len() - 1 {
            return Err(invalid_data("H CIGAR operation must be first or last"));
        }

        if op.op == 'S' {
            let before = &ops[..index];
            let after = &ops[index + 1..];
            let near_start = before.iter().all(|op| op.op == 'H');
            let near_end = after.iter().all(|op| op.op == 'H');
            if !near_start && !near_end {
                return Err(invalid_data(
                    "S CIGAR operation may only have H operations between it and an end",
                ));
            }
        }
    }

    Ok(())
}

fn parse_optional_array(value: &str) -> io::Result<SamOptionalArray> {
    let (array_type, values) = value
        .split_once(',')
        .map_or((value, ""), |(array_type, values)| (array_type, values));

    if array_type.len() != 1 {
        return Err(invalid_data("invalid B optional array type"));
    }

    let parts = values.split(',').filter(|part| !part.is_empty());

    match array_type.as_bytes()[0] {
        b'c' => parts
            .map(|part| {
                part.parse::<i8>()
                    .map_err(|_| invalid_data("invalid B:c value"))
            })
            .collect::<io::Result<Vec<_>>>()
            .map(SamOptionalArray::Int8),
        b'C' => parts
            .map(|part| {
                part.parse::<u8>()
                    .map_err(|_| invalid_data("invalid B:C value"))
            })
            .collect::<io::Result<Vec<_>>>()
            .map(SamOptionalArray::UInt8),
        b's' => parts
            .map(|part| {
                part.parse::<i16>()
                    .map_err(|_| invalid_data("invalid B:s value"))
            })
            .collect::<io::Result<Vec<_>>>()
            .map(SamOptionalArray::Int16),
        b'S' => parts
            .map(|part| {
                part.parse::<u16>()
                    .map_err(|_| invalid_data("invalid B:S value"))
            })
            .collect::<io::Result<Vec<_>>>()
            .map(SamOptionalArray::UInt16),
        b'i' => parts
            .map(|part| {
                part.parse::<i32>()
                    .map_err(|_| invalid_data("invalid B:i value"))
            })
            .collect::<io::Result<Vec<_>>>()
            .map(SamOptionalArray::Int32),
        b'I' => parts
            .map(|part| {
                part.parse::<u32>()
                    .map_err(|_| invalid_data("invalid B:I value"))
            })
            .collect::<io::Result<Vec<_>>>()
            .map(SamOptionalArray::UInt32),
        b'f' => parts
            .map(|part| {
                part.parse::<f32>()
                    .map_err(|_| invalid_data("invalid B:f value"))
            })
            .collect::<io::Result<Vec<_>>>()
            .map(SamOptionalArray::Float),
        _ => Err(invalid_data("invalid B optional array type")),
    }
}

fn parse_hex(value: &str) -> io::Result<Vec<u8>> {
    if !value.len().is_multiple_of(2) {
        return Err(invalid_data("hex optional value has odd length"));
    }

    let mut bytes = Vec::with_capacity(value.len() / 2);
    for index in (0..value.len()).step_by(2) {
        let byte = u8::from_str_radix(&value[index..index + 2], 16)
            .map_err(|_| invalid_data("invalid hex optional value"))?;
        bytes.push(byte);
    }

    Ok(bytes)
}

fn format_sam_header_record(record: &SamHeaderRecord) -> io::Result<String> {
    if record.record_type == "CO" {
        let comment = record.comment.as_deref().unwrap_or_default();
        if !is_printable_or_space_ascii(comment) {
            return Err(invalid_data("invalid @CO comment"));
        }
        return Ok(format!("@CO\t{comment}"));
    }

    if !matches!(record.record_type.as_str(), "HD" | "SQ" | "RG" | "PG") {
        return Err(invalid_data("unknown SAM header record type"));
    }
    if record.fields.is_empty() {
        return Err(invalid_data("SAM header record has no fields"));
    }

    let mut line = format!("@{}", record.record_type);
    let mut seen = HashSet::new();
    for field in &record.fields {
        if !is_header_tag(&field.tag) {
            return Err(invalid_data("invalid SAM header field tag"));
        }
        if field.value.is_empty() || !is_printable_or_space_ascii(&field.value) {
            return Err(invalid_data("invalid SAM header field value"));
        }
        if !seen.insert(field.tag.as_str()) {
            return Err(invalid_data("duplicate SAM header field tag"));
        }
        line.push('\t');
        line.push_str(&field.tag);
        line.push(':');
        line.push_str(&field.value);
    }

    Ok(line)
}

fn format_sam_record(record: &SamRecord) -> io::Result<String> {
    let mut fields = vec![
        record.qname.clone(),
        record.flag.to_string(),
        record.rname.clone(),
        record.pos.to_string(),
        record.mapq.to_string(),
        record.cigar.clone(),
        record.rnext.clone(),
        record.pnext.to_string(),
        record.tlen.to_string(),
        record.seq.clone(),
        record.qual.clone(),
    ];

    for optional in &record.optional {
        fields.push(format_sam_optional_field(optional)?);
    }

    Ok(fields.join("\t"))
}

fn format_sam_optional_field(field: &SamOptionalField) -> io::Result<String> {
    if !is_optional_tag(&field.tag) {
        return Err(invalid_data("invalid SAM optional field tag"));
    }

    Ok(format!(
        "{}:{}:{}",
        field.tag,
        field.value.value_type(),
        format_sam_optional_value(&field.value)?
    ))
}

fn format_sam_optional_value(value: &SamOptionalValue) -> io::Result<String> {
    match value {
        SamOptionalValue::Character(value) => {
            if !is_printable_ascii_char(*value) {
                return Err(invalid_data("invalid A optional value"));
            }
            Ok(value.to_string())
        }
        SamOptionalValue::Integer(value) => Ok(value.to_string()),
        SamOptionalValue::Float(value) => Ok(value.to_string()),
        SamOptionalValue::String(value) => {
            if !is_printable_or_space_ascii(value) {
                return Err(invalid_data("invalid Z optional value"));
            }
            Ok(value.clone())
        }
        SamOptionalValue::Hex(value) => Ok(bytes_to_hex(value)),
        SamOptionalValue::Array(value) => Ok(format_sam_optional_array(value)),
    }
}

fn format_sam_optional_array(value: &SamOptionalArray) -> String {
    match value {
        SamOptionalArray::Int8(values) => format_array_values('c', values),
        SamOptionalArray::UInt8(values) => format_array_values('C', values),
        SamOptionalArray::Int16(values) => format_array_values('s', values),
        SamOptionalArray::UInt16(values) => format_array_values('S', values),
        SamOptionalArray::Int32(values) => format_array_values('i', values),
        SamOptionalArray::UInt32(values) => format_array_values('I', values),
        SamOptionalArray::Float(values) => format_array_values('f', values),
    }
}

fn format_array_values<T: ToString>(array_type: char, values: &[T]) -> String {
    let mut output = array_type.to_string();
    for value in values {
        output.push(',');
        output.push_str(&value.to_string());
    }
    output
}

fn bytes_to_hex(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789ABCDEF";
    let mut output = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        output.push(HEX[(byte >> 4) as usize] as char);
        output.push(HEX[(byte & 0x0f) as usize] as char);
    }
    output
}

fn parse_u8(value: &str, message: &'static str) -> io::Result<u8> {
    value.parse::<u8>().map_err(|_| invalid_data(message))
}

fn parse_u16(value: &str, message: &'static str) -> io::Result<u16> {
    value.parse::<u16>().map_err(|_| invalid_data(message))
}

fn parse_u32_max_i32(value: &str, message: &'static str) -> io::Result<u32> {
    let value = value.parse::<u32>().map_err(|_| invalid_data(message))?;
    if value > i32::MAX as u32 {
        return Err(invalid_data(message));
    }

    Ok(value)
}

fn parse_i32(value: &str, message: &'static str) -> io::Result<i32> {
    value.parse::<i32>().map_err(|_| invalid_data(message))
}

fn parse_i64(value: &str, message: &'static str) -> io::Result<i64> {
    value.parse::<i64>().map_err(|_| invalid_data(message))
}

fn is_query_name(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 254
        && value
            .bytes()
            .all(|byte| (33..=126).contains(&byte) && byte != b'@')
}

fn is_reference_name(value: &str) -> bool {
    let mut bytes = value.bytes();
    let Some(first) = bytes.next() else {
        return false;
    };

    if first == b'*' || first == b'=' || !is_reference_name_byte(first) {
        return false;
    }

    bytes.all(is_reference_name_byte)
}

fn is_reference_name_byte(byte: u8) -> bool {
    matches!(
        byte,
        b'0'..=b'9'
            | b'A'..=b'Z'
            | b'a'..=b'z'
            | b'!'
            | b'#'
            | b'$'
            | b'%'
            | b'&'
            | b'*'
            | b'+'
            | b'.'
            | b'/'
            | b':'
            | b';'
            | b'='
            | b'?'
            | b'@'
            | b'^'
            | b'_'
            | b'|'
            | b'~'
            | b'-'
    )
}

fn is_sequence(value: &str) -> bool {
    value == "*"
        || (!value.is_empty()
            && value
                .bytes()
                .all(|byte| byte.is_ascii_alphabetic() || byte == b'=' || byte == b'.'))
}

fn is_quality(value: &str) -> bool {
    value == "*" || (!value.is_empty() && value.bytes().all(|byte| (33..=126).contains(&byte)))
}

fn is_header_tag(value: &str) -> bool {
    let bytes = value.as_bytes();
    bytes.len() == 2 && bytes[0].is_ascii_alphabetic() && bytes[1].is_ascii_alphanumeric()
}

fn is_optional_tag(value: &str) -> bool {
    is_header_tag(value)
}

fn is_printable_ascii_char(value: char) -> bool {
    value.is_ascii() && (33..=126).contains(&(value as u8))
}

fn is_printable_or_space_ascii(value: &str) -> bool {
    value.bytes().all(|byte| (32..=126).contains(&byte))
}

fn invalid_data(message: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

#[cfg(test)]
mod sam_tests {
    use super::*;

    const ALIGNED_SAM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/aligned.sam");
    const UNALIGNED_SAM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/unaligned.sam");
    const RECORD_COUNT: usize = 100;
    const FIRST_ALIGNED_ID: &str = "2f890ab0-63b9-45f6-ba8d-3c367ca26a63";
    const FIRST_UNALIGNED_ID: &str = "d3843747-8d64-429e-bc47-44763b006ad1";

    #[test]
    fn reads_one_record_at_a_time_from_path() {
        let mut reader = SamReader::from_path(ALIGNED_SAM).unwrap();
        assert_eq!(reader.header.first("HD").unwrap().value("VN"), Some("1.6"));
        assert_eq!(reader.header.sequence_names(), vec!["fc_reference"]);

        let first = reader.read_record().unwrap().unwrap();
        assert_eq!(first.qname, FIRST_ALIGNED_ID);
        assert_eq!(first.flag, flags::REVERSE_COMPLEMENTED);
        assert!(!first.is_unmapped());
        assert!(first.is_reverse_complemented());
        assert_eq!(first.rname, "fc_reference");
        assert_eq!(first.pos, 2);
        assert_eq!(first.mapq, 60);
        assert_eq!(first.cigar, "7S333M1I334M37S");
        assert_eq!(
            first.query_len_from_cigar().unwrap(),
            first.seq.len() as u64
        );
        assert_eq!(first.aux("NM"), Some(&SamOptionalValue::Integer(8)));
        assert_eq!(first.aux("tp"), Some(&SamOptionalValue::Character('P')));

        let mut count = 1;
        while reader.read_record().unwrap().is_some() {
            count += 1;
        }
        assert_eq!(count, RECORD_COUNT);
        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn read_all_materializes_sample_file() {
        let sam = Sam::from_path(UNALIGNED_SAM).unwrap();

        assert_eq!(sam.records.len(), RECORD_COUNT);
        assert_eq!(sam.records[0].qname, FIRST_UNALIGNED_ID);
        assert_eq!(sam.header.first("HD").unwrap().value("SO"), Some("unknown"));
        assert!(sam.header.first("RG").is_some());
        assert!(sam.records.iter().all(SamRecord::is_unmapped));
        assert!(sam.records.iter().all(|record| record.cigar == "*"));
        assert!(sam.records.iter().all(|record| record.rname == "*"));
    }

    #[test]
    fn records_iterator_streams_until_eof() {
        let mut reader = SamReader::from_path(ALIGNED_SAM).unwrap();
        let records = reader.records().collect::<io::Result<Vec<_>>>().unwrap();

        assert_eq!(records.len(), RECORD_COUNT);
        assert_eq!(records[0].read_name(), FIRST_ALIGNED_ID);
        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn sam_from_reader_materializes_stream() {
        let file = File::open(ALIGNED_SAM).expect("fixture should open");
        let sam = Sam::from_reader(file).expect("SAM should materialize from reader");

        assert_eq!(sam.records.len(), RECORD_COUNT);
        assert_eq!(sam.records[0].qname, FIRST_ALIGNED_ID);
    }

    #[test]
    fn materialized_sam_can_be_deep_cloned() {
        let original = Sam::from_path(ALIGNED_SAM).expect("SAM should materialize");
        let mut cloned = original.clone();
        assert_eq!(cloned, original);

        cloned.records[0].qname.push_str("_clone");
        cloned.header.records[0].fields[0].value = "9.9".to_string();

        assert_ne!(cloned, original);
        assert_eq!(original.records[0].qname, FIRST_ALIGNED_ID);
        assert_eq!(
            original.header.first("HD").unwrap().value("VN"),
            Some("1.6")
        );
    }

    #[test]
    fn parses_header_optional_fields_and_cigar_from_reader() {
        let input = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ref\tLN:45\nr001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\tIIIIIIIIIIIIIIIII\tNM:i:1\tSA:Z:ref,29,-,6H5M,17,0;\n";
        let sam = Sam::from_reader(&input[..]).unwrap();

        assert_eq!(sam.header.sequence_names(), vec!["ref"]);
        assert_eq!(sam.records.len(), 1);
        assert_eq!(sam.records[0].flag, 99);
        assert_eq!(
            sam.records[0].cigar_ops().unwrap(),
            vec![
                SamCigarOp { len: 8, op: 'M' },
                SamCigarOp { len: 2, op: 'I' },
                SamCigarOp { len: 4, op: 'M' },
                SamCigarOp { len: 1, op: 'D' },
                SamCigarOp { len: 3, op: 'M' },
            ]
        );
        assert_eq!(sam.records[0].reference_len_from_cigar().unwrap(), 16);
        assert_eq!(
            sam.records[0].aux("NM"),
            Some(&SamOptionalValue::Integer(1))
        );
    }

    #[test]
    fn writer_round_trips_materialized_records() {
        let sam = Sam::from_path(ALIGNED_SAM).unwrap();
        let mut output = Vec::new();
        sam.to_writer(&mut output).unwrap();

        assert_eq!(Sam::from_reader(&output[..]).unwrap(), sam);
    }

    #[test]
    fn writer_round_trips_optional_value_types() {
        let sam = Sam {
            header: SamHeader {
                records: vec![
                    SamHeaderRecord::new(
                        "HD".to_string(),
                        vec![SamHeaderField::new("VN".to_string(), "1.6".to_string())],
                    ),
                    SamHeaderRecord::new(
                        "SQ".to_string(),
                        vec![
                            SamHeaderField::new("SN".to_string(), "ref".to_string()),
                            SamHeaderField::new("LN".to_string(), "45".to_string()),
                        ],
                    ),
                    SamHeaderRecord::comment("writer comment".to_string()),
                ],
            },
            records: vec![SamRecord::new(
                "r001".to_string(),
                0,
                "ref".to_string(),
                1,
                60,
                "4M".to_string(),
                "*".to_string(),
                0,
                0,
                "ACGT".to_string(),
                "IIII".to_string(),
                vec![
                    SamOptionalField::new("AA".to_string(), SamOptionalValue::Character('x')),
                    SamOptionalField::new("NM".to_string(), SamOptionalValue::Integer(1)),
                    SamOptionalField::new("AS".to_string(), SamOptionalValue::Float(2.5)),
                    SamOptionalField::new(
                        "ZZ".to_string(),
                        SamOptionalValue::String("value".to_string()),
                    ),
                    SamOptionalField::new(
                        "HH".to_string(),
                        SamOptionalValue::Hex(vec![0, 15, 255]),
                    ),
                    SamOptionalField::new(
                        "BI".to_string(),
                        SamOptionalValue::Array(SamOptionalArray::Int16(vec![-2, 3])),
                    ),
                ],
            )],
        };
        let mut output = Vec::new();
        sam.to_writer(&mut output).unwrap();

        assert_eq!(Sam::from_reader(&output[..]).unwrap(), sam);
        assert!(String::from_utf8(output).unwrap().contains("HH:H:000FFF"));
    }

    #[test]
    fn writer_rejects_invalid_record() {
        let record = SamRecord::new(
            "bad name".to_string(),
            0,
            "*".to_string(),
            0,
            0,
            "*".to_string(),
            "*".to_string(),
            0,
            0,
            "*".to_string(),
            "*".to_string(),
            Vec::new(),
        );
        let mut output = Vec::new();
        let mut writer = SamWriter::from_writer(&mut output);

        assert!(writer.write_record(&record).is_err());
    }

    #[test]
    fn malformed_records_return_errors() {
        assert!(Sam::from_reader(&b"bad\tline\n"[..]).is_err());
        assert!(Sam::from_reader(&b"r1\t0\tref\t1\t60\t2M\t*\t0\t0\tA\tI\n"[..]).is_err());
        assert!(Sam::from_reader(&b"r1\t0\tref\t1\t60\t1M1S1M\t*\t0\t0\tAAA\tIII\n"[..]).is_err());
        assert!(
            Sam::from_reader(&b"@SQ\tSN:ref\tLN:0\nr1\t0\tref\t1\t60\t1M\t*\t0\t0\tA\tI\n"[..])
                .is_err()
        );
    }
}
