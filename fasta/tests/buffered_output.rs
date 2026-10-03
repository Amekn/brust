//! `from_path` output is buffered: `finish` writes the rest and reports errors.

use brust_fasta::{FastaRecord, FastaWriter};
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "brust-fasta-buffered-{}-{name}",
        std::process::id()
    ))
}

fn record() -> FastaRecord {
    FastaRecord::new("r1".into(), None, "ACGT".into())
}

#[test]
fn from_path_writes_through_a_buffer_that_finish_flushes() {
    let path = temp_path("finish.fa");
    let mut writer = FastaWriter::from_path(&path).unwrap();
    writer.write_record(&record()).unwrap();
    let file: BufWriter<File> = writer.finish().unwrap();

    assert!(file.buffer().is_empty());
    assert_eq!(fs::read(&path).unwrap(), b">r1\nACGT\n");
    fs::remove_file(path).unwrap();
}

#[test]
fn dropping_a_from_path_writer_still_writes_its_records() {
    let path = temp_path("drop.fa");
    let mut writer = FastaWriter::from_path(&path).unwrap();
    writer.write_record(&record()).unwrap();
    drop(writer);

    assert_eq!(fs::read(&path).unwrap(), b">r1\nACGT\n");
    fs::remove_file(path).unwrap();
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
fn finish_reports_errors_from_a_buffered_sink() {
    let mut writer = FastaWriter::from_writer(BufWriter::new(FullDisk));
    writer.write_record(&record()).unwrap();

    assert!(writer.finish().is_err());
}
