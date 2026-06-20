use brust_core::Compression;
use brust_fastq::{Fastq, FastqReader, FastqWriter};
use std::fs;
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};

const RECORD_COUNT: usize = 100;

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

struct TempFile(PathBuf);

impl TempFile {
    fn new(extension: &str) -> Self {
        let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        Self(std::env::temp_dir().join(format!(
            "brust-fastq-compression-{}-{counter}.{extension}",
            process::id()
        )))
    }

    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TempFile {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(name)
}

#[test]
fn compressed_fixture_streams_through_path_reader() {
    // The repository fixture exercises real multi-record gzip decompression.
    let mut reader = FastqReader::from_path(fixture("UDP0057_sub100.fastq.gz")).unwrap();
    let records = reader
        .records()
        .collect::<std::io::Result<Vec<_>>>()
        .unwrap();

    assert_eq!(records.len(), RECORD_COUNT);
    assert_eq!(records[0].sequence.len(), 687);
    assert_eq!(records[0].sequence.len(), records[0].quality.len());
}

#[test]
fn path_writer_compresses_fq_gz_and_round_trips() {
    // A `.fq.gz` destination selects gzip while preserving the streaming API.
    let expected = Fastq::from_path(fixture("UDP0057_sub100.fastq")).unwrap();
    let output = TempFile::new("fq.gz");
    let mut writer = FastqWriter::from_path(output.path()).unwrap();
    for record in &expected.records {
        writer.write_record(record).unwrap();
    }
    writer.finish().unwrap();

    assert_eq!(&fs::read(output.path()).unwrap()[..2], &[0x1f, 0x8b]);
    assert_eq!(Fastq::from_path(output.path()).unwrap(), expected);
    assert_eq!(Compression::from_path(output.path()), Compression::Gzip);
}
