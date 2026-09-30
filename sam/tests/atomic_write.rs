use brust_sam::{Sam, SamWriter};
use std::fs;
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};

const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/aligned.sam");

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn new(name: &str) -> Self {
        let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("brust-sam-{name}-{}-{counter}", process::id()));
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
    let payload = Sam::from_path(FIXTURE).unwrap();

    let mut plain = SamWriter::from_path(dir.join("plain.sam")).unwrap();
    plain.write_all(&payload).unwrap();
    plain.flush().unwrap();
    let mut atomic = SamWriter::from_path_atomic(dir.join("atomic.sam")).unwrap();
    atomic.write_all(&payload).unwrap();
    atomic.commit().unwrap();

    assert_eq!(
        fs::read(dir.join("plain.sam")).unwrap(),
        fs::read(dir.join("atomic.sam")).unwrap()
    );
    assert_eq!(names(dir.path()), ["atomic.sam", "plain.sam"]);
}

#[test]
fn dropped_atomic_writer_leaves_nothing() {
    let dir = TempDir::new("writer-dropped");
    let payload = Sam::from_path(FIXTURE).unwrap();

    let mut atomic = SamWriter::from_path_atomic(dir.join("out.sam")).unwrap();
    atomic.write_all(&payload).unwrap();
    drop(atomic);

    assert!(names(dir.path()).is_empty());
}

#[test]
fn to_path_atomic_matches_to_path() {
    let dir = TempDir::new("to-path");
    let payload = Sam::from_path(FIXTURE).unwrap();

    payload.to_path(dir.join("plain.sam")).unwrap();
    payload.to_path_atomic(dir.join("atomic.sam")).unwrap();

    assert_eq!(
        fs::read(dir.join("plain.sam")).unwrap(),
        fs::read(dir.join("atomic.sam")).unwrap()
    );
    assert_eq!(names(dir.path()), ["atomic.sam", "plain.sam"]);
}
