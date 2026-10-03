//! `from_path` output is buffered: `finish` writes the rest and reports errors.

use brust_sam::{Sam, SamWriter};
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

const SAM: &[u8] = b"@HD\tVN:1.6\nr1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";

fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("brust-sam-buffered-{}-{name}", std::process::id()))
}

#[test]
fn from_path_writes_through_a_buffer_that_finish_flushes() {
    let path = temp_path("finish.sam");
    let mut writer = SamWriter::from_path(&path).unwrap();
    writer.write_all(&Sam::from_reader(SAM).unwrap()).unwrap();
    let file: BufWriter<File> = writer.finish().unwrap();

    assert!(file.buffer().is_empty());
    assert_eq!(fs::read(&path).unwrap(), SAM);
    fs::remove_file(path).unwrap();
}

#[test]
fn dropping_a_from_path_writer_still_writes_its_records() {
    let path = temp_path("drop.sam");
    let mut writer = SamWriter::from_path(&path).unwrap();
    writer.write_all(&Sam::from_reader(SAM).unwrap()).unwrap();
    drop(writer);

    assert_eq!(fs::read(&path).unwrap(), SAM);
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
    let mut writer = SamWriter::from_writer(BufWriter::new(FullDisk));
    writer.write_all(&Sam::from_reader(SAM).unwrap()).unwrap();

    assert!(writer.finish().is_err());
}
