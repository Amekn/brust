use brust_core::Compression;
use brust_fastq::{Fastq, FastqRecord, FastqWriter};
use std::io::{self, BufWriter, Write};

fn record(sequence: &str, quality: &str) -> FastqRecord {
    FastqRecord::new("r1".into(), None, sequence.into(), quality.into())
}

fn write(record: &FastqRecord) -> io::Result<Vec<u8>> {
    let mut writer = FastqWriter::from_writer(Vec::new());
    writer.write_record(record)?;
    writer.finish()
}

#[test]
fn zero_length_read_round_trips() {
    let empty = record("", "");
    let output = write(&empty).unwrap();

    assert_eq!(Fastq::from_reader(&output[..]).unwrap().records, [empty]);
}

#[test]
fn writer_rejects_records_its_reader_would_misread() {
    // A sequence line starting with '+' reads as the separator, and trailing
    // whitespace is trimmed on reading.
    for (sequence, quality) in [("+ACG", "IIII"), ("ACG ", "IIII"), ("ACGT", "III ")] {
        let error = write(&record(sequence, quality)).unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidData, "{sequence:?}");
    }
}

#[test]
fn writer_rejects_quality_outside_printable_ascii() {
    let error = write(&record("ACGT", "II\u{1}I")).unwrap_err();

    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}

struct FullDisk;

impl Write for FullDisk {
    fn write(&mut self, _: &[u8]) -> io::Result<usize> {
        Err(io::Error::other("disk full"))
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

#[test]
fn gzip_finish_reports_errors_from_a_buffered_sink() {
    for compression in [Compression::Uncompressed, Compression::Gzip] {
        let mut writer =
            FastqWriter::from_writer_with_compression(BufWriter::new(FullDisk), compression);
        writer.write_record(&record("ACGT", "IIII")).unwrap();

        assert!(writer.finish().is_err(), "{compression:?}");
    }
}
