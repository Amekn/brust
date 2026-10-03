//! Shared POD5 test payloads.

use brust_pod5::Pod5;

pub const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/A_100.pod5");

/// A payload of `reads` reads made by cycling through the fixture's reads,
/// each copy with its own read ID and its own Signal rows.
pub fn pod5_with_reads(reads: usize) -> Pod5 {
    let base = Pod5::from_path(FIXTURE).unwrap();
    let mut pod5 = Pod5 {
        signals: Vec::new(),
        records: Vec::new(),
        ..base.clone()
    };
    for (number, record) in base.records.iter().cycle().take(reads).enumerate() {
        let read_id = format!("{number:08x}{}", &record.read_id[8..]);
        let mut record = record.clone();
        for row in &mut record.signal_rows {
            let mut signal = base.signals[*row as usize].clone();
            signal.read_id = read_id.clone();
            pod5.signals.push(signal);
            *row = pod5.signals.len() as u64 - 1;
        }
        record.read_id = read_id;
        pod5.records.push(record);
    }
    pod5
}
