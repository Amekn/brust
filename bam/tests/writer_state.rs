//! A BamWriter writes one header and records that fit its reference
//! dictionary, and an atomic writer publishes only a complete BAM.

use brust_bam::{Bam, BamWriter, SamToBamConverter};
use sam::Sam;
use std::fs;
use std::io;

/// A one-record BAM whose only reference is `name`.
fn bam_on(name: &str) -> Bam {
    let text = format!("@SQ\tSN:{name}\tLN:100\nr1\t0\t{name}\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n");
    let sam = Sam::from_reader(text.as_bytes()).unwrap();
    let converter = SamToBamConverter::new(&sam.header).unwrap();
    Bam {
        header: converter.header().clone(),
        refs: converter.refs().to_vec(),
        records: vec![converter.convert_record(&sam.records[0]).unwrap()],
    }
}

#[test]
fn write_all_refuses_a_payload_with_another_reference_dictionary() {
    // Its records' reference IDs would be read against the first header.
    let mut writer = BamWriter::from_writer(Vec::new());
    writer.write_all(&bam_on("chrA")).unwrap();

    let error = writer.write_all(&bam_on("chrB")).unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidInput);

    // The same dictionary may be written again, appending its records.
    writer.write_all(&bam_on("chrA")).unwrap();
    let bytes = writer.finish().unwrap();
    let bam = Bam::from_reader(&bytes[..]).unwrap();
    assert_eq!(bam.records.len(), 2);
    assert_eq!(bam.refs[0].name, "chrA");
}

#[test]
fn write_record_refuses_reference_ids_outside_the_written_dictionary() {
    let bam = bam_on("chrA");
    let mut writer = BamWriter::from_writer(Vec::new());
    writer.write_header(&bam.header, &bam.refs).unwrap();
    for (ref_id, next_ref_id) in [(1, -1), (0, 1), (-2, -1)] {
        let mut record = bam.records[0].clone();
        record.fixed.ref_id = ref_id;
        record.fixed.next_ref_id = next_ref_id;

        let error = writer.write_record(&record).unwrap_err();
        assert_eq!(
            error.kind(),
            io::ErrorKind::InvalidInput,
            "{ref_id} {next_ref_id}"
        );
    }
}

#[test]
fn atomic_commit_without_a_header_publishes_nothing() {
    // An EOF block alone is not a BAM.
    let dir = std::env::temp_dir().join(format!("brust-bam-empty-commit-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let target = dir.join("out.bam");
    let writer = BamWriter::from_path_atomic(&target).unwrap();

    let error = writer.commit().unwrap_err();
    let published = target.exists();
    fs::remove_dir_all(&dir).unwrap();

    assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
    assert!(!published);
}
