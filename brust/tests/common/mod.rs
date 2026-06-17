use std::fs;
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

pub fn fixture(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join(relative)
}

pub struct TempDir {
    path: PathBuf,
}

impl TempDir {
    pub fn new(name: &str) -> Self {
        let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path =
            std::env::temp_dir().join(format!("brust-test-{name}-{}-{counter}", process::id()));
        fs::create_dir_all(&path).expect("temporary test directory should be created");
        Self { path }
    }

    pub fn join(&self, path: &str) -> PathBuf {
        self.path.join(path)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}
