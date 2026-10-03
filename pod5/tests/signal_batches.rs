//! The Signal table is written in batches of 100 rows, as official POD5
//! writers do, and every row reads back in any order.

mod common;

use arrow_ipc::reader::FileReader;
use brust_pod5::{Pod5, Pod5Reader, Pod5SectionKind};
use std::io::Cursor;

fn signal_batch_rows(bytes: &[u8]) -> Vec<usize> {
    let reader = Pod5Reader::from_reader(Cursor::new(bytes)).unwrap();
    let section = reader
        .header
        .sections
        .iter()
        .find(|section| section.kind == Pod5SectionKind::Signal)
        .unwrap();
    let start = section.offset as usize;
    let arrow = &bytes[start..start + section.length as usize];
    FileReader::try_new(Cursor::new(arrow), None)
        .unwrap()
        .map(|batch| batch.unwrap().num_rows())
        .collect()
}

#[test]
fn signal_rows_are_written_in_batches_of_100() {
    // Official readers find a row's batch by dividing its number by the first
    // batch's row count, so every batch but the last must be that size.
    for (reads, batches) in [
        (0, vec![0]),
        (1, vec![1]),
        (100, vec![100]),
        (250, vec![100, 100, 50]),
    ] {
        let mut bytes = Vec::new();
        common::pod5_with_reads(reads)
            .to_writer(&mut bytes)
            .unwrap();

        assert_eq!(signal_batch_rows(&bytes), batches, "{reads} reads");
    }
}

#[test]
fn every_read_reads_back_from_a_batched_signal_table() {
    let pod5 = common::pod5_with_reads(250);
    let mut bytes = Vec::new();
    pod5.to_writer(&mut bytes).unwrap();

    let back = Pod5::from_reader(&bytes[..]).unwrap();
    assert_eq!(back.signals, pod5.signals);
    let mut reader = Pod5Reader::from_reader(Cursor::new(&bytes)).unwrap();
    for record in pod5.records.iter().rev() {
        assert_eq!(
            reader.signal_for_record(record).unwrap(),
            pod5.signal_for_record(record).unwrap(),
            "{}",
            record.read_id
        );
    }
}
