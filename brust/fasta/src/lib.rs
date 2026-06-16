//! FASTA reader and writer primitives.
//!
//! The crate exposes both streaming and materialized APIs:
//!
//! - [`FastaReader`] reads one [`FastaRecord`] at a time from any [`Read`]
//!   source.
//! - [`Fasta`] owns every parsed record and can be cloned or written back out.
//!
//! Header lines are split at the first whitespace character. The first token
//! after `>` becomes the record ID; the remaining header text becomes the
//! optional description. Sequence lines are concatenated with line endings
//! removed. The parser does not validate sequence alphabets; the writer only
//! rejects records whose IDs, descriptions, or sequences cannot be represented
//! as valid FASTA lines.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

/// Streaming FASTA parser over any readable byte stream.
///
/// `FastaReader` keeps the underlying reader and parser cursor in one place, so
/// it is intentionally not cloneable. Use [`FastaReader::read_all`] or
/// [`Fasta::from_path`] when a cloneable, in-memory representation is needed.
pub struct FastaReader<R: Read = File> {
    reader: BufReader<R>,
    // Header line consumed while finishing the previous record.
    pending_header: Option<String>,
    // Number of lines consumed from the underlying reader.
    line_number: usize,
}

/// Streaming FASTA writer over any writable byte stream.
///
/// `FastaWriter` emits one [`FastaRecord`] at a time. By default it writes each
/// sequence as a single line; call [`FastaWriter::set_line_width`] to wrap
/// sequence output at a fixed width.
pub struct FastaWriter<W: Write = File> {
    writer: W,
    line_width: Option<usize>,
}

impl FastaReader<File> {
    /// Opens a FASTA file from a filesystem path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        Self::from_reader(file)
    }

    /// Opens a FASTA file from a filesystem path.
    ///
    /// This is a convenience alias for [`FastaReader::from_path`].
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_path(path)
    }
}

impl<R: Read> FastaReader<R> {
    /// Creates a streaming parser from a FASTA byte stream.
    pub fn from_reader(reader: R) -> io::Result<Self> {
        Ok(Self {
            reader: BufReader::new(reader),
            pending_header: None,
            line_number: 0,
        })
    }

    /// Reads the next FASTA record from the stream.
    ///
    /// Lines before the next header are skipped. Returns `Ok(None)` when EOF is
    /// reached before another header line begins. Sequence lines are
    /// concatenated until the following header or EOF, with record line endings
    /// removed.
    pub fn read_record(&mut self) -> io::Result<Option<FastaRecord>> {
        let header_line = match self.pending_header.take() {
            Some(header) => header,
            None => {
                let mut line = String::new();
                loop {
                    line.clear();
                    let bytes = self.reader.read_line(&mut line)?;

                    if bytes == 0 {
                        return Ok(None);
                    }

                    self.line_number += 1;

                    if line.starts_with('>') {
                        break line;
                    }
                }
            }
        };

        let header = header_line.strip_prefix('>').unwrap().trim();
        let mut parts = header.splitn(2, char::is_whitespace);

        let id = parts.next().unwrap_or("").to_string();
        let description = parts
            .next()
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(String::from);

        let mut sequence = String::new();
        let mut line = String::new();

        loop {
            line.clear();
            let bytes = self.reader.read_line(&mut line)?;

            if bytes == 0 {
                break;
            }

            self.line_number += 1;

            if line.starts_with('>') {
                self.pending_header = Some(line);
                break;
            }

            sequence.push_str(line.trim_end());
        }

        Ok(Some(FastaRecord::new(id, description, sequence)))
    }

    /// Reads the next FASTA record.
    ///
    /// This is a compatibility alias for [`FastaReader::read_record`].
    pub fn read(&mut self) -> io::Result<Option<FastaRecord>> {
        self.read_record()
    }

    /// Returns an iterator over records in the stream.
    pub fn records(&mut self) -> FastaRecords<'_, R> {
        FastaRecords { reader: self }
    }

    /// Consumes this reader and materializes the entire FASTA stream.
    pub fn read_all(mut self) -> io::Result<Fasta> {
        let mut records = Vec::new();
        while let Some(record) = self.read_record()? {
            records.push(record);
        }

        Ok(Fasta { records })
    }
}

impl FastaWriter<File> {
    /// Creates or truncates a FASTA file at a filesystem path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self::from_writer(file))
    }

    /// Creates or truncates a FASTA file at a filesystem path.
    ///
    /// This is a convenience alias for [`FastaWriter::from_path`].
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_path(path)
    }
}

impl<W: Write> FastaWriter<W> {
    /// Creates a FASTA writer from a writable byte stream.
    pub fn from_writer(writer: W) -> Self {
        Self {
            writer,
            line_width: None,
        }
    }

    /// Sets the sequence line wrapping width for future records.
    ///
    /// A width of `0` disables wrapping and writes each sequence on one line.
    pub fn set_line_width(&mut self, width: usize) {
        self.line_width = (width > 0).then_some(width);
    }

    /// Writes one FASTA record.
    ///
    /// The record ID must be non-empty and contain no whitespace. IDs,
    /// descriptions, and sequences must not contain line endings.
    pub fn write_record(&mut self, record: &FastaRecord) -> io::Result<()> {
        validate_fasta_record(record)?;

        write!(self.writer, ">{}", record.id)?;
        if let Some(description) = &record.description
            && !description.is_empty()
        {
            write!(self.writer, " {description}")?;
        }
        self.writer.write_all(b"\n")?;

        match self.line_width {
            Some(width) => {
                for chunk in record.sequence.as_bytes().chunks(width) {
                    self.writer.write_all(chunk)?;
                    self.writer.write_all(b"\n")?;
                }
            }
            None => {
                self.writer.write_all(record.sequence.as_bytes())?;
                self.writer.write_all(b"\n")?;
            }
        }

        Ok(())
    }

    /// Writes one FASTA record.
    ///
    /// This is a compatibility alias for [`FastaWriter::write_record`].
    pub fn write(&mut self, record: &FastaRecord) -> io::Result<()> {
        self.write_record(record)
    }

    /// Writes all records from a materialized FASTA file.
    ///
    /// Each record is validated by [`FastaWriter::write_record`].
    pub fn write_all(&mut self, fasta: &Fasta) -> io::Result<()> {
        for record in &fasta.records {
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

/// A fully materialized FASTA payload.
///
/// `Fasta` owns all records from the stream, so cloning this type performs a
/// deep copy of the parsed FASTA data.
#[derive(Debug, Clone, PartialEq)]
pub struct Fasta {
    /// All sequence records loaded from the stream.
    pub records: Vec<FastaRecord>,
}

impl Fasta {
    /// Opens a FASTA file and materializes all records into memory.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        FastaReader::from_path(path)?.read_all()
    }

    /// Materializes all records from a FASTA byte stream.
    pub fn from_reader<R: Read>(reader: R) -> io::Result<Self> {
        FastaReader::from_reader(reader)?.read_all()
    }

    /// Writes this FASTA payload to a filesystem path.
    pub fn to_path<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut writer = FastaWriter::from_path(path)?;
        writer.write_all(self)?;
        writer.flush()
    }

    /// Writes this FASTA payload to a writable byte stream.
    pub fn to_writer<W: Write>(&self, writer: W) -> io::Result<()> {
        let mut writer = FastaWriter::from_writer(writer);
        writer.write_all(self)?;
        writer.flush()
    }
}

/// A FASTA sequence record.
///
/// The ID is the first token after `>`, the description is the remaining header
/// text after the first whitespace, and the sequence is stored as concatenated
/// sequence lines without line endings.
#[derive(Debug, Clone, PartialEq)]
pub struct FastaRecord {
    /// Record identifier from the FASTA header.
    pub id: String,
    /// Optional description text from the FASTA header.
    pub description: Option<String>,
    /// Sequence bases or residues with record line endings removed.
    pub sequence: String,
}

impl FastaRecord {
    /// Creates a FASTA sequence record from parsed fields.
    pub fn new(id: String, description: Option<String>, sequence: String) -> Self {
        Self {
            id,
            description,
            sequence,
        }
    }
}

/// Iterator over records from a [`FastaReader`].
pub struct FastaRecords<'a, R: Read> {
    reader: &'a mut FastaReader<R>,
}

impl<R: Read> Iterator for FastaRecords<'_, R> {
    type Item = io::Result<FastaRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(error) => Some(Err(error)),
        }
    }
}

fn validate_fasta_record(record: &FastaRecord) -> io::Result<()> {
    if record.id.is_empty() || record.id.chars().any(char::is_whitespace) {
        return Err(invalid_data(
            "FASTA record ID must be non-empty and contain no whitespace",
        ));
    }
    if record.id.contains(['\r', '\n']) {
        return Err(invalid_data(
            "FASTA record ID must not contain line endings",
        ));
    }
    if record
        .description
        .as_deref()
        .is_some_and(|description| description.contains(['\r', '\n']))
    {
        return Err(invalid_data(
            "FASTA record description must not contain line endings",
        ));
    }
    if record.sequence.contains(['\r', '\n']) {
        return Err(invalid_data("FASTA sequence must not contain line endings"));
    }

    Ok(())
}

fn invalid_data(message: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

#[cfg(test)]
mod fasta_tests {
    use super::*;

    const ACE2_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/ace2_fragments.fasta");
    const ACE2_DESCRIPTION: &str = "pTWIST-FcaRI_FLAG_GST_His";
    const ACE2_IDS: [&str; 6] = [
        "ACE2_pangolin_fragment",
        "ACE2_mouse_fragment",
        "ACE2_mink_fragment",
        "ACE2_hu_fragment",
        "ACE2_civet_fragment",
        "ACE2_bat_fragment",
    ];

    #[test]
    fn reads_one_record_at_a_time_from_path() {
        let mut reader = FastaReader::from_path(ACE2_PATH).unwrap();

        for expected_id in ACE2_IDS {
            let record = reader.read_record().unwrap().unwrap();
            assert_eq!(record.id, expected_id);
            assert_eq!(record.description.as_deref(), Some(ACE2_DESCRIPTION));
            assert!(record.sequence.starts_with("AAAGGCGGCCGC"));
            assert!(record.sequence.ends_with("CTGGTT"));
            assert!(!record.sequence.contains('\n'));
            assert!(!record.sequence.contains('\r'));
        }

        assert_eq!(reader.read_record().unwrap(), None);
    }

    #[test]
    fn read_all_materializes_sample_file() {
        let fasta = Fasta::from_path(ACE2_PATH).unwrap();
        let ids = fasta
            .records
            .iter()
            .map(|record| record.id.as_str())
            .collect::<Vec<_>>();

        assert_eq!(fasta.records.len(), ACE2_IDS.len());
        assert_eq!(ids, ACE2_IDS.to_vec());
        assert_eq!(fasta.records[0].id, "ACE2_pangolin_fragment");
        assert_eq!(fasta.records[5].id, "ACE2_bat_fragment");
    }

    #[test]
    fn records_iterates_over_sample_file() {
        let mut reader = FastaReader::from_path(ACE2_PATH).unwrap();
        let records = reader.records().collect::<io::Result<Vec<_>>>().unwrap();
        let ids = records
            .iter()
            .map(|record| record.id.as_str())
            .collect::<Vec<_>>();

        assert_eq!(records.len(), ACE2_IDS.len());
        assert_eq!(ids, ACE2_IDS.to_vec());
    }

    #[test]
    fn parses_header_without_description() {
        let input = b">seq1\nACGT\n>seq2\nTTAA\n";
        let fasta = Fasta::from_reader(&input[..]).unwrap();

        assert_eq!(
            fasta.records,
            vec![
                FastaRecord::new("seq1".to_string(), None, "ACGT".to_string()),
                FastaRecord::new("seq2".to_string(), None, "TTAA".to_string()),
            ]
        );
    }

    #[test]
    fn reads_final_record_without_trailing_header() {
        let input = b">seq1\nACGT\n>seq2 desc\nTTAA";
        let fasta = Fasta::from_reader(&input[..]).unwrap();

        assert_eq!(
            fasta.records,
            vec![
                FastaRecord::new("seq1".to_string(), None, "ACGT".to_string()),
                FastaRecord::new(
                    "seq2".to_string(),
                    Some("desc".to_string()),
                    "TTAA".to_string()
                ),
            ]
        );
    }

    #[test]
    fn writer_round_trips_materialized_records() {
        let fasta = Fasta::from_reader(&b">seq1 description\nACGT\n>seq2\nTTAA\n"[..]).unwrap();
        let mut output = Vec::new();
        fasta.to_writer(&mut output).unwrap();

        assert_eq!(Fasta::from_reader(&output[..]).unwrap(), fasta);
    }

    #[test]
    fn writer_wraps_sequence_lines_when_configured() {
        let record = FastaRecord::new("seq1".to_string(), None, "ACGTACGT".to_string());
        let mut output = Vec::new();
        let mut writer = FastaWriter::from_writer(&mut output);
        writer.set_line_width(3);
        writer.write_record(&record).unwrap();

        assert_eq!(String::from_utf8(output).unwrap(), ">seq1\nACG\nTAC\nGT\n");
    }

    #[test]
    fn writer_rejects_invalid_records() {
        let record = FastaRecord::new("bad id".to_string(), None, "ACGT".to_string());
        let mut output = Vec::new();
        let mut writer = FastaWriter::from_writer(&mut output);

        assert!(writer.write_record(&record).is_err());
    }
}
