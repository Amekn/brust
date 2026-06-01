use bam::*;

#[test]
fn reader_opens_bam_header() {
    let reader = BamReader::new("aligned.bam").expect("failed to create BamReader");
    assert_eq!(reader.header.magic, [0x42, 0x41, 0x4D, 0x01]);
    assert_eq!(reader.header.n_ref, 1);
}

#[test]
fn bam_materializes_and_clones_records() {
    let bam = Bam::from_path("aligned.bam").expect("failed to materialize Bam");
    let cloned = bam.clone();

    assert_eq!(bam.records.len(), 100);
    assert_eq!(cloned, bam);
}
