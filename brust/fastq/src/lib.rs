//! FASTQ reader and writer primitives.
//!
//! The crate exposes both streaming and materialized APIs:
//!
//! - [`FastqReader`] reads one [`FastqRecord`] at a time from any [`Read`]
//!   source.
//! - [`Fastq`] owns every parsed record and can be cloned or written back out.
//!
//! FASTQ headers are split at the first whitespace character: the first token
//! after `@` is the record ID, and the remaining text is stored as the optional
//! description. Sequence lines are concatenated until the `+` separator, then
//! quality lines are concatenated until their length exactly matches the
//! sequence length. The writer emits a canonical four-line record with a bare
//! `+` separator and validates that sequence and quality lengths match.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

/// Streaming FASTQ parser over any readable byte stream.
///
/// `FastqReader` keeps the underlying reader and parser cursor in one place, so
/// it is intentionally not cloneable. Use [`FastqReader::read_all`] or
/// [`Fastq::from_path`] when a cloneable, in-memory representation is needed.
pub struct FastqReader<R: Read = File> {
    reader: BufReader<R>,
    // Number of lines consumed from the underlying reader.
    line_number: usize,
}

/// Streaming FASTQ writer over any writable byte stream.
///
/// `FastqWriter` emits one [`FastqRecord`] at a time. It validates record shape
/// before writing and always writes sequence and quality as single lines.
pub struct FastqWriter<W: Write = File> {
    writer: W,
}

impl FastqReader<File> {
    /// Opens a FASTQ file from a filesystem path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        Self::from_reader(file)
    }

    /// Opens a FASTQ file from a filesystem path.
    ///
    /// This is a convenience alias for [`FastqReader::from_path`].
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_path(path)
    }
}

impl<R: Read> FastqReader<R> {
    /// Creates a streaming parser from a FASTQ byte stream.
    pub fn from_reader(reader: R) -> io::Result<Self> {
        Ok(Self {
            reader: BufReader::new(reader),
            line_number: 0,
        })
    }

    /// Reads the next FASTQ record from the stream.
    ///
    /// Returns `Ok(None)` when the stream is at EOF before another header. A
    /// malformed or truncated record returns `InvalidData`. Sequence line
    /// endings are removed until a `+` separator line is reached, and quality
    /// line endings are removed until the accumulated quality length equals the
    /// parsed sequence length.
    pub fn read_record(&mut self) -> io::Result<Option<FastqRecord>> {
        let mut line = String::new();
        let byte = self.reader.read_line(&mut line)?;
        if byte == 0 {
            return Ok(None);
        }
        self.line_number += 1;

        let header = line
            .strip_prefix('@')
            .ok_or_else(|| invalid_data("FASTQ header line must start with @"))?
            .trim();
        let mut parts = header.splitn(2, char::is_whitespace);

        let id = parts.next().unwrap_or("").to_string();
        let description = parts
            .next()
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(String::from);

        let mut sequence = String::new();

        loop {
            line.clear();
            let byte = self.reader.read_line(&mut line)?;
            if byte == 0 {
                return Err(invalid_data("FASTQ sequence ended before + separator"));
            }
            self.line_number += 1;
            if line.starts_with('+') {
                break;
            }
            sequence.push_str(line.trim_end());
        }

        let mut quality = String::new();
        loop {
            line.clear();
            let byte = self.reader.read_line(&mut line)?;
            if byte == 0 {
                return Err(invalid_data("FASTQ quality ended before sequence length"));
            }
            self.line_number += 1;
            quality.push_str(line.trim_end());
            if quality.len() > sequence.len() {
                return Err(invalid_data("FASTQ quality length exceeds sequence length"));
            }
            if quality.len() == sequence.len() {
                break;
            }
        }

        Ok(Some(FastqRecord {
            id,
            description,
            sequence,
            quality,
        }))
    }

    /// Reads the next FASTQ record.
    ///
    /// This is a compatibility alias for [`FastqReader::read_record`].
    pub fn read(&mut self) -> io::Result<Option<FastqRecord>> {
        self.read_record()
    }

    /// Returns an iterator over records in the stream.
    pub fn records(&mut self) -> FastqRecords<'_, R> {
        FastqRecords { reader: self }
    }

    /// Consumes this reader and materializes the entire FASTQ stream.
    pub fn read_all(mut self) -> io::Result<Fastq> {
        let mut records = Vec::new();
        while let Some(record) = self.read()? {
            records.push(record);
        }
        Ok(Fastq { records })
    }
}

impl FastqWriter<File> {
    /// Creates or truncates a FASTQ file at a filesystem path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self::from_writer(file))
    }

    /// Creates or truncates a FASTQ file at a filesystem path.
    ///
    /// This is a convenience alias for [`FastqWriter::from_path`].
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_path(path)
    }
}

impl<W: Write> FastqWriter<W> {
    /// Creates a FASTQ writer from a writable byte stream.
    pub fn from_writer(writer: W) -> Self {
        Self { writer }
    }

    /// Writes one FASTQ record.
    ///
    /// The record ID must be non-empty and contain no whitespace. IDs,
    /// descriptions, sequences, and qualities must not contain line endings,
    /// and sequence and quality lengths must match.
    pub fn write_record(&mut self, record: &FastqRecord) -> io::Result<()> {
        validate_fastq_record(record)?;

        write!(self.writer, "@{}", record.id)?;
        if let Some(description) = &record.description
            && !description.is_empty()
        {
            write!(self.writer, " {description}")?;
        }
        self.writer.write_all(b"\n")?;
        self.writer.write_all(record.sequence.as_bytes())?;
        self.writer.write_all(b"\n+\n")?;
        self.writer.write_all(record.quality.as_bytes())?;
        self.writer.write_all(b"\n")?;
        Ok(())
    }

    /// Writes one FASTQ record.
    ///
    /// This is a compatibility alias for [`FastqWriter::write_record`].
    pub fn write(&mut self, record: &FastqRecord) -> io::Result<()> {
        self.write_record(record)
    }

    /// Writes all records from a materialized FASTQ value.
    ///
    /// Each record is validated by [`FastqWriter::write_record`].
    pub fn write_all(&mut self, fastq: &Fastq) -> io::Result<()> {
        for record in &fastq.records {
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

/// A fully materialized FASTQ payload.
///
/// `Fastq` owns all records from the stream, so cloning this type performs a
/// deep copy of the parsed FASTQ data.
#[derive(Debug, Clone, PartialEq)]
pub struct Fastq {
    /// All sequence records loaded from the stream.
    pub records: Vec<FastqRecord>,
}

impl Fastq {
    /// Opens a FASTQ file and materializes all records into memory.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        FastqReader::from_path(path)?.read_all()
    }

    /// Materializes all records from a FASTQ byte stream.
    pub fn from_reader<R: Read>(reader: R) -> io::Result<Self> {
        FastqReader::from_reader(reader)?.read_all()
    }

    /// Writes this FASTQ payload to a filesystem path.
    pub fn to_path<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut writer = FastqWriter::from_path(path)?;
        writer.write_all(self)?;
        writer.flush()
    }

    /// Writes this FASTQ payload to a writable byte stream.
    pub fn to_writer<W: Write>(&self, writer: W) -> io::Result<()> {
        let mut writer = FastqWriter::from_writer(writer);
        writer.write_all(self)?;
        writer.flush()
    }
}

/// A FASTQ sequence record.
///
/// The ID is the first token after `@`, the description is the remaining header
/// text after the first whitespace, the sequence is stored without line
/// endings, and the quality string is stored without line endings.
#[derive(Debug, Clone, PartialEq)]
pub struct FastqRecord {
    /// Record identifier from the FASTQ header.
    pub id: String,
    /// Optional description text from the FASTQ header.
    pub description: Option<String>,
    /// Sequence bases with record line endings removed.
    pub sequence: String,
    /// Base quality characters with record line endings removed.
    pub quality: String,
}

impl FastqRecord {
    /// Creates a FASTQ sequence record from parsed fields.
    pub fn new(id: String, description: Option<String>, sequence: String, quality: String) -> Self {
        Self {
            id,
            description,
            sequence,
            quality,
        }
    }
}

/// Iterator over records from a [`FastqReader`].
pub struct FastqRecords<'a, R: Read> {
    reader: &'a mut FastqReader<R>,
}

impl<R: Read> Iterator for FastqRecords<'_, R> {
    type Item = io::Result<FastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(error) => Some(Err(error)),
        }
    }
}

fn validate_fastq_record(record: &FastqRecord) -> io::Result<()> {
    if record.id.is_empty() || record.id.chars().any(char::is_whitespace) {
        return Err(invalid_data(
            "FASTQ record ID must be non-empty and contain no whitespace",
        ));
    }
    if record.id.contains(['\r', '\n']) {
        return Err(invalid_data(
            "FASTQ record ID must not contain line endings",
        ));
    }
    if record
        .description
        .as_deref()
        .is_some_and(|description| description.contains(['\r', '\n']))
    {
        return Err(invalid_data(
            "FASTQ record description must not contain line endings",
        ));
    }
    if record.sequence.contains(['\r', '\n']) || record.quality.contains(['\r', '\n']) {
        return Err(invalid_data(
            "FASTQ sequence and quality must not contain line endings",
        ));
    }
    if record.sequence.len() != record.quality.len() {
        return Err(invalid_data(
            "FASTQ sequence and quality lengths must match",
        ));
    }

    Ok(())
}

fn invalid_data(message: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

#[cfg(test)]
mod fastq_tests {
    use super::*;

    const UDP0057_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/UDP0057_sub100.fastq");
    const UDP0057_RECORD_COUNT: usize = 100;
    const UDP0057_READ_LENGTH: usize = 687;
    const FIRST_ID: &str = "FT100112191L1C001R00100001474:UMI_GATGTTAT_AGTAGTTC";
    const MIDDLE_ID: &str = "FT100112191L1C001R00100008713:UMI_CCAGGTCT_GTTTCTTG";
    const LAST_ID: &str = "FT100112191L1C001R00100016533:UMI_GATGCTTA_AGGTTGTC";

    fn assert_sample_record_shape(record: &FastqRecord) {
        assert_eq!(record.description.as_deref(), None);
        assert_eq!(record.sequence.len(), UDP0057_READ_LENGTH);
        assert_eq!(record.quality.len(), UDP0057_READ_LENGTH);
        assert!(!record.sequence.contains('\n'));
        assert!(!record.sequence.contains('\r'));
        assert!(!record.quality.contains('\n'));
        assert!(!record.quality.contains('\r'));
    }

    #[test]
    fn reads_one_record_at_a_time_from_path() {
        let mut reader = FastqReader::from_path(UDP0057_PATH).unwrap();
        let mut records = Vec::new();

        while let Some(record) = reader.read_record().unwrap() {
            assert_sample_record_shape(&record);
            records.push(record);
        }

        assert_eq!(records.len(), UDP0057_RECORD_COUNT);
        assert_eq!(records[0].id, FIRST_ID);
        assert_eq!(records[UDP0057_RECORD_COUNT / 2 - 1].id, MIDDLE_ID);
        assert_eq!(records[UDP0057_RECORD_COUNT - 1].id, LAST_ID);
        assert!(records[0].sequence.starts_with("TCGCCATCCTGCTGTCATCCGCGT"));
        assert!(records[0].sequence.ends_with("GGTTCTGGTTCCGGTGATTTTGAT"));
        assert!(records[0].quality.starts_with("IIIIIFGIIIIIIIIIIIIIIIII"));
        assert!(records[0].quality.ends_with("96IIII'B/)IICGHIIC/+/III"));
        assert!(
            reader
                .read_record()
                .expect("second EOF read should succeed")
                .is_none()
        );
    }

    #[test]
    fn read_all_materializes_sample_file() {
        let fastq = Fastq::from_path(UDP0057_PATH).unwrap();
        let ids = fastq
            .records
            .iter()
            .map(|record| record.id.as_str())
            .collect::<Vec<_>>();

        assert_eq!(fastq.records.len(), UDP0057_RECORD_COUNT);
        assert_eq!(ids[0], FIRST_ID);
        assert_eq!(ids[UDP0057_RECORD_COUNT / 2 - 1], MIDDLE_ID);
        assert_eq!(ids[UDP0057_RECORD_COUNT - 1], LAST_ID);
        assert!(
            fastq
                .records
                .iter()
                .all(|record| record.sequence.len() == record.quality.len())
        );
    }

    #[test]
    fn records_iterates_over_sample_file() {
        let mut reader = FastqReader::from_path(UDP0057_PATH).unwrap();
        let records = reader.records().collect::<io::Result<Vec<_>>>().unwrap();

        assert_eq!(records.len(), UDP0057_RECORD_COUNT);
        assert_eq!(records[0].id, FIRST_ID);
        assert_eq!(records[UDP0057_RECORD_COUNT - 1].id, LAST_ID);
        assert!(
            reader
                .read_record()
                .expect("post-iterator EOF should succeed")
                .is_none()
        );
    }

    #[test]
    fn fastq_from_reader_materializes_stream() {
        let file = File::open(UDP0057_PATH).expect("fixture should open");
        let fastq = Fastq::from_reader(file).expect("Fastq should materialize from reader");

        assert_eq!(fastq.records.len(), UDP0057_RECORD_COUNT);
        assert_eq!(fastq.records[0].id, FIRST_ID);
        assert_eq!(fastq.records[UDP0057_RECORD_COUNT - 1].id, LAST_ID);
    }

    #[test]
    fn materialized_fastq_can_be_deep_cloned() {
        let original = Fastq::from_path(UDP0057_PATH).expect("Fastq should materialize");
        let mut cloned = original.clone();
        assert_eq!(cloned, original);

        cloned.records[0].id.push_str("_clone");
        cloned.records[0].sequence.push('A');
        cloned.records[0].quality.push('I');

        assert_ne!(cloned, original);
        assert_eq!(original.records[0].id, FIRST_ID);
        assert_eq!(original.records[0].sequence.len(), UDP0057_READ_LENGTH);
        assert_eq!(original.records[0].quality.len(), UDP0057_READ_LENGTH);
    }

    #[test]
    fn parses_description_and_multiline_payload() {
        let input = b"@seq1 sample description\nAC\nGT\n+\nIIII\n@seq2\nTT\nAA\n+\nHH\nHH\n";
        let fastq = Fastq::from_reader(&input[..]).unwrap();

        assert_eq!(
            fastq.records,
            vec![
                FastqRecord::new(
                    "seq1".to_string(),
                    Some("sample description".to_string()),
                    "ACGT".to_string(),
                    "IIII".to_string(),
                ),
                FastqRecord::new(
                    "seq2".to_string(),
                    None,
                    "TTAA".to_string(),
                    "HHHH".to_string()
                ),
            ]
        );
    }

    #[test]
    fn writer_round_trips_materialized_records() {
        let fastq =
            Fastq::from_reader(&b"@seq1 description\nACGT\n+\nIIII\n@seq2\nTTAA\n+\nHHHH\n"[..])
                .unwrap();
        let mut output = Vec::new();
        fastq.to_writer(&mut output).unwrap();

        assert_eq!(Fastq::from_reader(&output[..]).unwrap(), fastq);
    }

    #[test]
    fn writer_emits_expected_record_shape() {
        let record = FastqRecord::new(
            "seq1".to_string(),
            Some("description".to_string()),
            "ACGT".to_string(),
            "IIII".to_string(),
        );
        let mut output = Vec::new();
        let mut writer = FastqWriter::from_writer(&mut output);
        writer.write_record(&record).unwrap();

        assert_eq!(
            String::from_utf8(output).unwrap(),
            "@seq1 description\nACGT\n+\nIIII\n"
        );
    }

    #[test]
    fn writer_rejects_invalid_records() {
        let record = FastqRecord::new(
            "seq1".to_string(),
            None,
            "ACGT".to_string(),
            "III".to_string(),
        );
        let mut output = Vec::new();
        let mut writer = FastqWriter::from_writer(&mut output);

        assert!(writer.write_record(&record).is_err());
    }

    #[test]
    fn malformed_truncated_record_returns_error() {
        assert!(Fastq::from_reader(&b"@seq1\nACGT\n"[..]).is_err());
        assert!(Fastq::from_reader(&b"@seq1\nACGT\n+\nII\n"[..]).is_err());
    }

    #[test]
    fn malformed_overlong_quality_returns_error() {
        let err = Fastq::from_reader(&b"@seq1\nACGT\n+\nIIIII\n"[..]).unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert_eq!(
            err.to_string(),
            "FASTQ quality length exceeds sequence length"
        );
    }
}
