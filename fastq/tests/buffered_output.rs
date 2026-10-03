//! `from_path` output is buffered: `finish` writes the rest.

use brust_core::Compression;
use brust_fastq::{FastqRecord, FastqWriter};
use std::fs::{self, File};
use std::io::BufWriter;
use std::path::PathBuf;

fn temp_path(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "brust-fastq-buffered-{}-{name}",
        std::process::id()
    ))
}

fn record() -> FastqRecord {
    FastqRecord::new("r1".into(), None, "ACGT".into(), "IIII".into())
}

#[test]
fn from_path_writes_through_a_buffer_that_finish_flushes() {
    for (name, compression) in [
        ("plain.fastq", Compression::Uncompressed),
        ("packed.fastq.gz", Compression::Gzip),
    ] {
        let path = temp_path(name);
        let mut writer = FastqWriter::from_path(&path).unwrap();
        writer.write_record(&record()).unwrap();
        let file: BufWriter<File> = writer.finish().unwrap();

        assert!(file.buffer().is_empty(), "{compression:?}");
        let fastq = brust_fastq::Fastq::from_path(&path).unwrap();
        assert_eq!(fastq.records, [record()], "{compression:?}");
        fs::remove_file(path).unwrap();
    }
}

#[test]
fn from_path_with_compression_is_buffered_too() {
    let path = temp_path("explicit.tmp");
    let mut writer = FastqWriter::from_path_with_compression(&path, Compression::Gzip).unwrap();
    writer.write_record(&record()).unwrap();
    let _: BufWriter<File> = writer.finish().unwrap();

    let reader = brust_fastq::FastqReader::from_reader_with_compression(
        File::open(&path).unwrap(),
        Compression::Gzip,
    );
    assert_eq!(reader.unwrap().read_all().unwrap().records, [record()]);
    fs::remove_file(path).unwrap();
}

#[test]
fn dropping_a_from_path_writer_still_writes_its_records() {
    for name in ["drop.fastq", "drop.fastq.gz"] {
        let path = temp_path(name);
        let mut writer = FastqWriter::from_path(&path).unwrap();
        writer.write_record(&record()).unwrap();
        drop(writer);

        let fastq = brust_fastq::Fastq::from_path(&path).unwrap();
        assert_eq!(fastq.records, [record()], "{name}");
        fs::remove_file(path).unwrap();
    }
}
