use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

/// A buffered file sink that publishes its contents at the target path only
/// when [`AtomicFile::commit`] succeeds.
///
/// [`AtomicFile::create`] opens a hidden temporary file beside the target,
/// named `.{pid}.{counter}.tmp.{file_name}`. The target is not opened or
/// touched. Writes go through a `BufWriter`, and `flush` only pushes the buffer
/// to the temporary file; it does not sync it to disk.
///
/// [`AtomicFile::commit`] runs in this order:
///
/// 1. flush the buffer;
/// 2. `sync_all` the temporary file;
/// 3. close it;
/// 4. rename it over the target;
/// 5. on Unix, open the parent folder and `sync_all` it.
///
/// The new contents appear at the target only when `commit` reaches the
/// rename. If flush, sync or rename fails, the temporary file is removed, the
/// target is left as it was, and the error is returned. `commit` can also
/// return an error after the rename: if opening or syncing the parent folder
/// fails, the new file is already in place. That error keeps the original
/// [`io::ErrorKind`] and its message says the file was renamed but the folder
/// sync failed.
///
/// Dropping an `AtomicFile` without committing it (an early return, `?`, or a
/// panic that unwinds) removes the temporary file, so nothing appears at the
/// target and an existing target is unchanged. A kill, Ctrl-C or power loss
/// before `commit` leaves the temporary file behind; the target is still
/// untouched, and the `.tmp.` in the name makes the leftover easy to spot.
///
/// Known limits:
///
/// - The rename replaces a symlink at the target with a plain file.
/// - The new file gets default permissions (umask), not the old file's.
/// - Hard links to the old file keep the old contents.
/// - A read-only existing target is replaced, because `rename` does not need
///   write access to the old file. This is deliberate and matches `mv`.
/// - Paths are kept as given, so a relative path is resolved against the
///   current folder at each step. Don't change the current folder between
///   `create` and `commit`.
/// - The rename is atomic only within one filesystem. Putting the temporary
///   file beside the target guarantees that.
/// - Only Unix syncs the folder after the rename. On other platforms a power
///   cut soon after `commit` can undo the rename. The target then holds the
///   old file, or is missing if it was new, but never a partial one.
#[derive(Debug)]
pub struct AtomicFile {
    writer: Option<BufWriter<File>>,
    temp_path: PathBuf,
    target: PathBuf,
    /// Set only once the rename has succeeded.
    published: bool,
}

impl AtomicFile {
    /// Creates a temporary file beside `path`, ready for writing.
    ///
    /// Returns `InvalidInput` when `path` has no file name, and the open error
    /// (for example `NotFound` for a missing parent folder) otherwise. Nothing
    /// is created on failure.
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<AtomicFile> {
        let target = path.as_ref();
        let file_name = target.file_name().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("output path {} must include a file name", target.display()),
            )
        })?;
        let parent = parent_of(target);

        for _ in 0..100 {
            let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
            // Keep the complete target name at the end so suffix checks such
            // as `.gz` and `.bam` still work on the temporary path.
            let mut temp_name = OsString::from(format!(".{}.{}.tmp.", process::id(), counter));
            temp_name.push(file_name);
            let temp_path = parent.join(temp_name);

            match OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temp_path)
            {
                Ok(file) => {
                    return Ok(AtomicFile {
                        writer: Some(BufWriter::new(file)),
                        temp_path,
                        target: target.to_path_buf(),
                        published: false,
                    });
                }
                Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(error),
            }
        }

        Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!(
                "could not create temporary output beside {}",
                target.display()
            ),
        ))
    }

    /// Flushes, syncs and renames the temporary file over the target.
    ///
    /// On an error from the flush, sync or rename, the temporary file is
    /// removed and the target is unchanged. An error about syncing the parent
    /// folder (Unix only) comes after the rename: the new file is already in
    /// place. See the type documentation for details.
    pub fn commit(mut self) -> io::Result<()> {
        // Any early return drops `self`, and `Drop` removes the temp file.
        let writer = self.writer.take().expect("writer is present until commit");
        let file = match writer.into_inner() {
            Ok(file) => file,
            Err(error) => {
                // Take the buffer apart so it is not flushed again on drop.
                let (error, writer) = error.into_parts();
                let _ = writer.into_parts();
                return Err(error);
            }
        };
        file.sync_all()?;
        // Closing is unchecked: std ignores close errors, and the data is synced.
        drop(file);

        fs::rename(&self.temp_path, &self.target)?;
        self.published = true;

        #[cfg(unix)]
        {
            let parent = parent_of(&self.target);
            File::open(parent)
                .and_then(|folder| folder.sync_all())
                .map_err(|error| {
                    io::Error::new(
                        error.kind(),
                        format!(
                            "renamed {} into place but syncing folder {} failed: {error}",
                            self.target.display(),
                            parent.display()
                        ),
                    )
                })?;
        }

        Ok(())
    }

    fn writer(&mut self) -> &mut BufWriter<File> {
        self.writer
            .as_mut()
            .expect("writer is present until commit")
    }
}

/// The folder holding `path`, or `.` when the path has no parent part.
fn parent_of(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

impl Write for AtomicFile {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.writer().write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.writer().flush()
    }
}

impl Drop for AtomicFile {
    fn drop(&mut self) {
        if let Some(writer) = self.writer.take() {
            // Discard the buffer without flushing; the file closes here.
            let _ = writer.into_parts();
        }
        if !self.published {
            let _ = fs::remove_file(&self.temp_path);
        }
    }
}
