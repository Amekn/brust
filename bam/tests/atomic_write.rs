use brust_bam::{Bam, BamWriter};
use std::fs;
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};

const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/aligned.bam");

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn new(name: &str) -> Self {
        let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("brust-bam-{name}-{}-{counter}", process::id()));
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
    let fixture = Bam::from_path(FIXTURE).unwrap();

    let mut plain = BamWriter::from_path(dir.join("plain.bam")).unwrap();
    plain.write_all(&fixture).unwrap();
    plain.finish().unwrap();
    let mut atomic = BamWriter::from_path_atomic(dir.join("atomic.bam")).unwrap();
    atomic.write_all(&fixture).unwrap();
    atomic.commit().unwrap();

    assert_eq!(
        fs::read(dir.join("plain.bam")).unwrap(),
        fs::read(dir.join("atomic.bam")).unwrap()
    );
    assert_eq!(names(dir.path()), ["atomic.bam", "plain.bam"]);
}

#[test]
fn atomic_threaded_writer_matches_single_thread() {
    let dir = TempDir::new("threaded-matches");
    let fixture = Bam::from_path(FIXTURE).unwrap();

    let mut single = BamWriter::from_path_atomic(dir.join("single.bam")).unwrap();
    single.write_all(&fixture).unwrap();
    single.commit().unwrap();
    let mut threaded =
        BamWriter::from_path_atomic_with_threads(dir.join("threaded.bam"), 4).unwrap();
    threaded.write_all(&fixture).unwrap();
    threaded.commit().unwrap();

    assert_eq!(
        fs::read(dir.join("single.bam")).unwrap(),
        fs::read(dir.join("threaded.bam")).unwrap()
    );
    assert_eq!(names(dir.path()), ["single.bam", "threaded.bam"]);
}

#[test]
fn dropped_atomic_writer_leaves_nothing() {
    let dir = TempDir::new("writer-dropped");
    let fixture = Bam::from_path(FIXTURE).unwrap();

    let mut atomic = BamWriter::from_path_atomic(dir.join("out.bam")).unwrap();
    atomic.write_all(&fixture).unwrap();
    drop(atomic);

    assert!(names(dir.path()).is_empty());
}

#[test]
fn dropped_threaded_writer_leaves_nothing() {
    let dir = TempDir::new("threaded-dropped");
    let fixture = Bam::from_path(FIXTURE).unwrap();

    let mut writer = BamWriter::from_path_atomic_with_threads(dir.join("out.bam"), 4).unwrap();
    writer.write_all(&fixture).unwrap();
    drop(writer);

    assert!(names(dir.path()).is_empty());
}

#[test]
fn to_path_atomic_matches_to_path() {
    let dir = TempDir::new("to-path");
    let fixture = Bam::from_path(FIXTURE).unwrap();

    fixture.to_path(dir.join("plain.bam")).unwrap();
    fixture.to_path_atomic(dir.join("atomic.bam")).unwrap();

    assert_eq!(
        fs::read(dir.join("plain.bam")).unwrap(),
        fs::read(dir.join("atomic.bam")).unwrap()
    );
    assert_eq!(names(dir.path()), ["atomic.bam", "plain.bam"]);
}
