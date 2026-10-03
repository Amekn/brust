//! Write failures keep their own error, so callers can tell a full disk or a
//! closed pipe from bad data.

mod common;

use brust_pod5::Pod5Writer;
use std::io::{self, Write};

/// An OS error code to inject. What it means differs by platform (EPIPE on
/// Linux and macOS), so the tests compare against its own kind.
const OS_ERROR: i32 = 32;

struct FailingWrite {
    /// Bytes to accept before failing.
    room: usize,
}

impl Write for FailingWrite {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        if self.room == 0 {
            return Err(io::Error::from_raw_os_error(OS_ERROR));
        }
        let accepted = buf.len().min(self.room);
        self.room -= accepted;
        Ok(accepted)
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

struct FailingFlush;

impl Write for FailingFlush {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Err(io::Error::from_raw_os_error(OS_ERROR))
    }
}

fn assert_injected_error(error: io::Error) {
    let injected = io::Error::from_raw_os_error(OS_ERROR);
    assert_eq!(error.kind(), injected.kind(), "{error}");
    assert_eq!(error.raw_os_error(), Some(OS_ERROR), "{error}");
}

#[test]
fn write_errors_inside_the_signal_table_keep_their_kind() {
    let pod5 = common::pod5_with_reads(300);
    // Fail at the start, inside the Signal table, and near the end.
    for room in [0, 100_000, 2_000_000, 2_800_000] {
        let mut writer = Pod5Writer::from_writer(FailingWrite { room });

        assert_injected_error(writer.write_all(&pod5).unwrap_err());
    }
}

#[test]
fn flush_errors_keep_their_kind() {
    let mut writer = Pod5Writer::from_writer(FailingFlush);

    assert_injected_error(writer.write_all(&common::pod5_with_reads(10)).unwrap_err());
}
