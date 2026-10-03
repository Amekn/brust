use brust_pod5::{Pod5, Pod5Writer};
use std::fs;
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};

const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/A_100.pod5");

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn new(name: &str) -> Self {
        let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("brust-pod5-{name}-{}-{counter}", process::id()));
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
    let payload = Pod5::from_path(FIXTURE).unwrap();

    let mut plain = Pod5Writer::from_path(dir.join("plain.pod5")).unwrap();
    plain.write_all(&payload).unwrap();
    plain.flush().unwrap();
    let mut atomic = Pod5Writer::from_path_atomic(dir.join("atomic.pod5")).unwrap();
    atomic.write_all(&payload).unwrap();
    atomic.commit().unwrap();

    assert_eq!(
        fs::read(dir.join("plain.pod5")).unwrap(),
        fs::read(dir.join("atomic.pod5")).unwrap()
    );
    assert_eq!(names(dir.path()), ["atomic.pod5", "plain.pod5"]);
}

#[test]
fn dropped_atomic_writer_leaves_nothing() {
    let dir = TempDir::new("writer-dropped");
    let payload = Pod5::from_path(FIXTURE).unwrap();

    let mut atomic = Pod5Writer::from_path_atomic(dir.join("out.pod5")).unwrap();
    atomic.write_all(&payload).unwrap();
    drop(atomic);

    assert!(names(dir.path()).is_empty());
}

#[test]
fn to_path_atomic_matches_to_path() {
    let dir = TempDir::new("to-path");
    let payload = Pod5::from_path(FIXTURE).unwrap();

    payload.to_path(dir.join("plain.pod5")).unwrap();
    payload.to_path_atomic(dir.join("atomic.pod5")).unwrap();

    assert_eq!(
        fs::read(dir.join("plain.pod5")).unwrap(),
        fs::read(dir.join("atomic.pod5")).unwrap()
    );
    assert_eq!(names(dir.path()), ["atomic.pod5", "plain.pod5"]);
}

#[test]
fn writer_rejects_a_second_payload() {
    // Each payload is a whole POD5 container; appending a second one used to
    // corrupt the file while both calls returned Ok.
    let payload = Pod5::from_path(FIXTURE).unwrap();
    let mut single = Vec::new();
    payload.to_writer(&mut single).unwrap();

    let mut writer = Pod5Writer::from_writer(Vec::new());
    writer.write_all(&payload).unwrap();
    let error = writer.write_all(&payload).unwrap_err();

    assert_eq!(error.kind(), std::io::ErrorKind::InvalidInput);
    assert_eq!(writer.into_inner(), single);
}

#[test]
fn atomic_commit_without_a_payload_publishes_nothing() {
    let dir = TempDir::new("writer-empty");
    let writer = Pod5Writer::from_path_atomic(dir.join("out.pod5")).unwrap();

    let error = writer.commit().unwrap_err();

    assert_eq!(error.kind(), std::io::ErrorKind::InvalidInput);
    assert!(names(dir.path()).is_empty());
}

#[test]
fn to_path_keeps_existing_file_when_payload_is_invalid() {
    let dir = TempDir::new("to-path-invalid");
    let target = dir.join("out.pod5");
    fs::write(&target, b"old contents").unwrap();
    let mut payload = Pod5::from_path(FIXTURE).unwrap();
    payload.records[0].signal_rows = vec![payload.signals.len() as u64];

    assert!(payload.to_path(&target).is_err());
    assert_eq!(fs::read(&target).unwrap(), b"old contents");
}
