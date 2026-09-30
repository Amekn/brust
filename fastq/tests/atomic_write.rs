use brust_fastq::{Fastq, FastqWriter};
use std::fs;
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};

const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/UDP0057_sub100.fastq");

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn new(name: &str) -> Self {
        let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("brust-fastq-{name}-{}-{counter}", process::id()));
        fs::create_dir_all(&path).expect("temporary test directory should be created");
        Self { path }
    }

    fn join(&self, name: &str) -> PathBuf {
        self.path.join(name)
    }

    fn path(&self) -> &Path {
        &self.path
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

#[test]
fn atomic_writer_matches_from_path() {
    let dir = TempDir::new("writer-matches");
    let fixture = Fastq::from_path(FIXTURE).unwrap();

    let mut plain = FastqWriter::from_path(dir.join("plain.fastq")).unwrap();
    plain.write_all(&fixture).unwrap();
    plain.finish().unwrap();
    let mut atomic = FastqWriter::from_path_atomic(dir.join("atomic.fastq")).unwrap();
    atomic.write_all(&fixture).unwrap();
    atomic.commit().unwrap();

    assert_eq!(
        fs::read(dir.join("plain.fastq")).unwrap(),
        fs::read(dir.join("atomic.fastq")).unwrap()
    );
    assert_eq!(names(dir.path()), ["atomic.fastq", "plain.fastq"]);
}

#[test]
fn dropped_atomic_writer_leaves_nothing() {
    let dir = TempDir::new("writer-dropped");
    let fixture = Fastq::from_path(FIXTURE).unwrap();

    let mut atomic = FastqWriter::from_path_atomic(dir.join("out.fastq")).unwrap();
    atomic.write_all(&fixture).unwrap();
    drop(atomic);

    assert!(names(dir.path()).is_empty());
}

#[test]
fn to_path_atomic_matches_to_path() {
    let dir = TempDir::new("to-path");
    let fixture = Fastq::from_path(FIXTURE).unwrap();

    for (plain, atomic) in [
        ("plain.fastq", "atomic.fastq"),
        ("plain.fastq.gz", "atomic.fastq.gz"),
    ] {
        fixture.to_path(dir.join(plain)).unwrap();
        fixture.to_path_atomic(dir.join(atomic)).unwrap();
        assert_eq!(
            fs::read(dir.join(plain)).unwrap(),
            fs::read(dir.join(atomic)).unwrap()
        );
    }
    assert_eq!(
        names(dir.path()),
        [
            "atomic.fastq",
            "atomic.fastq.gz",
            "plain.fastq",
            "plain.fastq.gz"
        ]
    );
}

#[test]
fn atomic_gzip_target_writes_gzip() {
    let dir = TempDir::new("gzip-target");
    let fixture = Fastq::from_path(FIXTURE).unwrap();

    let mut writer = FastqWriter::from_path_atomic(dir.join("reads.fq.gz")).unwrap();
    writer.write_all(&fixture).unwrap();
    writer.commit().unwrap();

    assert_eq!(
        &fs::read(dir.join("reads.fq.gz")).unwrap()[..2],
        &[0x1f, 0x8b]
    );
    assert_eq!(
        Fastq::from_path(dir.join("reads.fq.gz")).unwrap().records,
        fixture.records
    );
    assert_eq!(names(dir.path()), ["reads.fq.gz"]);
}

#[test]
fn atomic_gzip_detection_uses_target_name() {
    let dir = TempDir::new("gzip-target-name");
    let fixture = Fastq::from_path(FIXTURE).unwrap();
    let name = "réads 1.FQ.GZ";

    let mut writer = FastqWriter::from_path_atomic(dir.join(name)).unwrap();
    writer.write_all(&fixture).unwrap();
    writer.commit().unwrap();

    assert_eq!(&fs::read(dir.join(name)).unwrap()[..2], &[0x1f, 0x8b]);
    assert_eq!(names(dir.path()), [name]);
}
