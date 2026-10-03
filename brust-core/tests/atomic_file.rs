use brust_core::AtomicFile;
use std::fs;
use std::io::{ErrorKind, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn new(name: &str) -> Self {
        let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "brust-core-atomic-{name}-{}-{counter}",
            process::id()
        ));
        fs::create_dir_all(&path).expect("temporary test directory should be created");
        Self { path }
    }

    fn join(&self, path: &str) -> PathBuf {
        self.path.join(path)
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
fn target_is_untouched_until_commit() {
    let dir = TempDir::new("untouched");
    fs::write(dir.join("old.txt"), b"old").unwrap();
    let mut a = AtomicFile::create(dir.join("new.txt")).unwrap();
    let mut b = AtomicFile::create(dir.join("old.txt")).unwrap();
    a.write_all(b"new").unwrap();
    a.flush().unwrap();
    b.write_all(b"new").unwrap();
    b.flush().unwrap();
    assert!(!dir.join("new.txt").exists());
    assert_eq!(fs::read(dir.join("old.txt")).unwrap(), b"old");
}

#[test]
fn commit_publishes_bytes_and_leaves_no_temp() {
    for existing in [false, true] {
        let dir = TempDir::new("commit");
        if existing {
            fs::write(dir.join("out.txt"), b"old").unwrap();
        }
        let mut file = AtomicFile::create(dir.join("out.txt")).unwrap();
        file.write_all(b"new\n").unwrap();
        file.commit().unwrap();
        assert_eq!(fs::read(dir.join("out.txt")).unwrap(), b"new\n");
        assert_eq!(names(dir.path()), ["out.txt"]);
    }
}

#[test]
fn drop_without_commit_removes_temp() {
    let dir = TempDir::new("drop-missing");
    {
        let mut file = AtomicFile::create(dir.join("out.txt")).unwrap();
        file.write_all(b"partial").unwrap();
    }
    assert!(names(dir.path()).is_empty());

    let dir = TempDir::new("drop-existing");
    fs::write(dir.join("out.txt"), b"old").unwrap();
    {
        let mut file = AtomicFile::create(dir.join("out.txt")).unwrap();
        file.write_all(b"partial").unwrap();
        file.flush().unwrap();
    }
    assert_eq!(names(dir.path()), ["out.txt"]);
    assert_eq!(fs::read(dir.join("out.txt")).unwrap(), b"old");
}

#[test]
fn panic_before_commit_removes_temp() {
    let dir = TempDir::new("panic");
    let target = dir.join("out.txt");
    let result = std::panic::catch_unwind(|| {
        let mut file = AtomicFile::create(&target).unwrap();
        file.write_all(b"partial").unwrap();
        panic!("writer failed");
    });
    assert!(result.is_err());
    assert!(names(dir.path()).is_empty());
}

#[test]
fn temp_name_is_hidden_and_keeps_target_name() {
    let dir = TempDir::new("temp-name");
    let _file = AtomicFile::create(dir.join("reads.fq.gz")).unwrap();
    let names = names(dir.path());
    assert_eq!(names.len(), 1);
    assert!(
        names[0].starts_with('.')
            && names[0].contains(".tmp.")
            && names[0].ends_with("reads.fq.gz")
    );
}

#[test]
fn path_without_file_name_is_invalid_input() {
    let dir = TempDir::new("no-file-name");
    let error = AtomicFile::create(dir.path().join("..")).unwrap_err();
    assert_eq!(error.kind(), ErrorKind::InvalidInput);
    assert!(names(dir.path()).is_empty());
}

#[test]
fn missing_parent_folder_is_not_found() {
    let dir = TempDir::new("missing-parent");
    let error = AtomicFile::create(dir.join("missing/out.txt")).unwrap_err();
    assert_eq!(error.kind(), ErrorKind::NotFound);
    assert!(names(dir.path()).is_empty());
}

#[test]
fn failed_rename_removes_temp_and_keeps_target_folder() {
    let dir = TempDir::new("failed-rename");
    let mut file = AtomicFile::create(dir.join("out")).unwrap();
    // A folder appearing at the target after `create` makes the rename fail.
    fs::create_dir(dir.join("out")).unwrap();
    fs::write(dir.join("out/keep.txt"), b"keep").unwrap();
    file.write_all(b"new").unwrap();
    assert!(file.commit().is_err());
    assert_eq!(names(dir.path()), ["out"]);
    assert_eq!(fs::read(dir.join("out/keep.txt")).unwrap(), b"keep");
}

#[test]
fn two_files_for_one_target_never_mix() {
    let dir = TempDir::new("two-files");
    let mut a = AtomicFile::create(dir.join("out.bin")).unwrap();
    let mut b = AtomicFile::create(dir.join("out.bin")).unwrap();
    a.write_all(&vec![b'a'; 65_536]).unwrap();
    b.write_all(&vec![b'b'; 65_536]).unwrap();
    assert_eq!(names(dir.path()).len(), 2);
    a.commit().unwrap();
    b.commit().unwrap();
    assert_eq!(fs::read(dir.join("out.bin")).unwrap(), vec![b'b'; 65_536]);
    assert_eq!(names(dir.path()), ["out.bin"]);
}

#[cfg(unix)]
#[test]
fn commit_replaces_read_only_target() {
    use std::os::unix::fs::PermissionsExt;

    let dir = TempDir::new("read-only");
    fs::write(dir.join("out.txt"), b"old").unwrap();
    fs::set_permissions(dir.join("out.txt"), fs::Permissions::from_mode(0o444)).unwrap();
    let mut file = AtomicFile::create(dir.join("out.txt")).unwrap();
    file.write_all(b"new").unwrap();
    file.commit().unwrap();
    assert_eq!(fs::read(dir.join("out.txt")).unwrap(), b"new");
    assert_eq!(names(dir.path()), ["out.txt"]);
}

#[cfg(unix)]
#[test]
fn folder_sync_failure_reports_published_file() {
    use std::os::unix::fs::PermissionsExt;

    let dir = TempDir::new("folder-sync");
    let mut file = AtomicFile::create(dir.join("out.txt")).unwrap();
    file.write_all(b"new").unwrap();
    // Write + search without read: rename works, opening the folder does not.
    fs::set_permissions(dir.path(), fs::Permissions::from_mode(0o300)).unwrap();
    if fs::File::open(dir.path()).is_ok() {
        // Running as root: permissions are not enforced, so the failure can't be provoked.
        fs::set_permissions(dir.path(), fs::Permissions::from_mode(0o700)).unwrap();
        return;
    }
    let result = file.commit();
    fs::set_permissions(dir.path(), fs::Permissions::from_mode(0o700)).unwrap();
    let error = result.unwrap_err();
    assert_eq!(error.kind(), ErrorKind::PermissionDenied);
    assert!(error.to_string().contains("into place"));
    assert_eq!(fs::read(dir.join("out.txt")).unwrap(), b"new");
    assert_eq!(names(dir.path()), ["out.txt"]);
}

#[test]
fn atomic_file_is_send_and_sync() {
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<AtomicFile>();
}

#[test]
fn longest_legal_file_name_can_be_written() {
    // The temporary name used to prefix the whole target name, so a legal
    // 255-byte name failed with "File name too long".
    let dir = TempDir::new("long-name");
    for name in [format!("{}.fastq.gz", "a".repeat(246)), "é".repeat(127)] {
        assert!(name.len() <= 255);
        let mut file = AtomicFile::create(dir.join(&name)).unwrap();
        file.write_all(b"data").unwrap();
        file.commit().unwrap();
        assert_eq!(fs::read(dir.join(&name)).unwrap(), b"data");
    }
}

#[test]
fn existing_directory_target_is_refused_up_front() {
    let dir = TempDir::new("directory-target");
    fs::create_dir(dir.join("out")).unwrap();

    let error = AtomicFile::create(dir.join("out")).unwrap_err();

    assert_eq!(error.kind(), ErrorKind::InvalidInput);
    assert_eq!(names(dir.path()), ["out"]);
}

#[cfg(unix)]
#[test]
fn special_file_target_is_refused_instead_of_replaced() {
    // Renaming over a FIFO, socket or device would silently replace it with a
    // regular file.
    let dir = TempDir::new("special-target");
    let socket = dir.join("out.sock");
    let _listener = std::os::unix::net::UnixListener::bind(&socket).unwrap();

    let error = AtomicFile::create(&socket).unwrap_err();

    assert_eq!(error.kind(), ErrorKind::InvalidInput);
    use std::os::unix::fs::FileTypeExt;
    assert!(fs::metadata(&socket).unwrap().file_type().is_socket());
    assert_eq!(names(dir.path()), ["out.sock"]);
}
