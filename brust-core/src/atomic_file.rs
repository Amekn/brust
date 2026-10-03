use std::ffi::{OsStr, OsString};
use std::fs::{self, File, OpenOptions};
use std::hash::BuildHasher;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

/// A buffered file sink that publishes its contents at the target path only
/// when [`AtomicFile::commit`] succeeds.
///
/// [`AtomicFile::create`] opens a hidden temporary file beside the target,
/// named `.{pid}.{counter}.tmp.{file_name}`, or `.{pid}.{counter}.tmp.{hash}`
/// when that would be longer than 255 bytes. The target is not opened
/// or touched, but an existing target that is not a regular file or symlink
/// (a folder, FIFO, socket or device) is refused, since the rename would
/// replace it. Writes go through a `BufWriter`, and `flush` only pushes the buffer
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
/// rename. If flush, sync or rename fails, the temporary file is removed (best
/// effort), the target is left as it was, and the error is returned. `commit` can also
/// return an error after the rename: if opening or syncing the parent folder
/// fails, the new file is already in place. That error keeps the original
/// [`io::ErrorKind`] and its message says the file was renamed but the folder
/// sync failed.
///
/// Dropping an `AtomicFile` without committing it (an early return, `?`, or a
/// panic that unwinds) removes the temporary file, so nothing appears at the
/// target and an existing target is unchanged. A kill, Ctrl-C or power loss
/// before `commit` leaves the temporary file behind; the target is still
/// untouched, and the `.{pid}.{counter}.tmp` prefix makes the leftover easy to
/// spot.
///
/// Known limits:
///
/// - The rename replaces a symlink at the target with a plain file.
/// - The new file gets default permissions (umask), not the old file's.
/// - Hard links to the old file keep the old contents.
/// - On Unix, a read-only existing target is replaced, because `rename` does
///   not need write access to the old file. This is deliberate and matches
///   `mv`. Other platforms may refuse the rename and keep the old file.
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
    /// Returns `InvalidInput` when `path` has no file name or names an existing
    /// folder, FIFO, socket or device, and the open error (for example
    /// `NotFound` for a missing parent folder) otherwise. Nothing is created on
    /// failure.
    pub fn create<P: AsRef<Path>>(path: P) -> io::Result<AtomicFile> {
        let target = path.as_ref();
        let file_name = target.file_name().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("output path {} must include a file name", target.display()),
            )
        })?;
        // The rename would replace a folder, FIFO, socket or device with a
        // regular file (or fail only after all the output was written).
        if let Ok(metadata) = fs::symlink_metadata(target) {
            let file_type = metadata.file_type();
            if !file_type.is_file() && !file_type.is_symlink() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "output path {} exists and is not a regular file",
                        target.display()
                    ),
                ));
            }
        }
        let parent = parent_of(target);

        for _ in 0..100 {
            let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
            let temp_path = parent.join(temp_file_name(file_name, counter));

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
    /// removed (best effort) and the target is unchanged. An error about syncing the parent
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

/// Longest file name, in bytes, that common filesystems accept.
const MAX_FILE_NAME_BYTES: usize = 255;

/// Returns `.{pid}.{counter}.tmp.{file_name}`, or `.{pid}.{counter}.tmp.{hash}`
/// when that would be longer than [`MAX_FILE_NAME_BYTES`].
///
/// A name that fits can't spell the target's own name on any filesystem: it is
/// the target's name after a visible ASCII prefix. A shortened name could
/// (case, Unicode normalisation, ignorable characters and Windows trailing
/// dots all make different spellings one name), so a long name is replaced by
/// 16 hex digits of a hash with a random per-process key, which no target
/// name can be made to match.
fn temp_file_name(file_name: &OsStr, counter: u64) -> OsString {
    let mut name = OsString::from(format!(".{}.{counter}.tmp.", process::id()));
    if name.len() + file_name.len() <= MAX_FILE_NAME_BYTES {
        name.push(file_name);
    } else {
        name.push(format!("{:016x}", RANDOM_STATE.hash_one(file_name)));
    }
    name
}

/// Randomly keyed hasher for shortened temporary names.
static RANDOM_STATE: std::sync::LazyLock<std::hash::RandomState> =
    std::sync::LazyLock::new(std::hash::RandomState::new);

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

#[cfg(test)]
mod tests {
    use super::*;

    /// A 255-byte target name that starts with the temporary prefix for
    /// `counter`, with the marker in the given case.
    fn aliasing_name(counter: u64, marker: &str) -> String {
        let prefix = format!(".{}.{counter}.{marker}.", process::id());
        format!("{prefix}{}", "a".repeat(MAX_FILE_NAME_BYTES - prefix.len()))
    }

    #[test]
    fn long_target_names_become_a_keyed_hash_in_temporary_names() {
        // A shortened copy of a long name could spell the target's own name.
        let prefix = format!(".{}.7.tmp.", process::id());
        for marker in ["tmp", "TMP"] {
            let name = aliasing_name(7, marker);
            let temp = temp_file_name(OsStr::new(&name), 7).into_string().unwrap();
            let hash = temp.strip_prefix(&prefix).unwrap();
            assert_eq!(hash.len(), 16, "{temp}");
            assert!(hash.bytes().all(|byte| byte.is_ascii_hexdigit()), "{temp}");
        }
        let short = temp_file_name(OsStr::new("out.fastq.gz"), 7);
        assert_eq!(short, OsString::from(format!("{prefix}out.fastq.gz")));
    }

    #[test]
    fn output_never_appears_at_a_long_target_before_commit() {
        // This is the only test in this binary that calls AtomicFile::create,
        // so the counter it reads is the one create will use next.
        let dir = std::env::temp_dir().join(format!("brust-core-alias-{}", process::id()));
        fs::create_dir_all(&dir).unwrap();
        let target = dir.join(aliasing_name(TEMP_COUNTER.load(Ordering::Relaxed), "tmp"));

        let mut file = AtomicFile::create(&target).unwrap();
        file.write_all(b"data").unwrap();
        file.flush().unwrap();
        let visible_early = target.exists();
        file.commit().unwrap();
        let contents = fs::read(&target).unwrap();
        fs::remove_dir_all(&dir).unwrap();

        assert!(
            !visible_early,
            "output appeared at the target before commit"
        );
        assert_eq!(contents, b"data");
    }
}
