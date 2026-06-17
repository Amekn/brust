use brust_bam::{Bam, BamReader, BgzfVirtualOffset};
use sam::{Sam, SamHeader, SamHeaderField, SamHeaderRecord, SamRecord};

const ALIGNED_BAM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/aligned.bam");

#[test]
fn reader_reports_virtual_offsets_for_records() {
    let mut reader = BamReader::from_path(ALIGNED_BAM).unwrap();
    let first = reader.read_record_with_virtual_offset().unwrap().unwrap();
    let second = reader.read_record_with_virtual_offset().unwrap().unwrap();

    assert_eq!(
        first.virtual_offset,
        BgzfVirtualOffset::from_raw(first.virtual_offset.raw())
    );
    assert!(second.virtual_offset > first.virtual_offset);
    assert!(first.virtual_offset.compressed_offset() > 0);
}

#[test]
fn bam_record_formats_as_sam_line() {
    let mut reader = BamReader::from_path(ALIGNED_BAM).unwrap();
    let record = reader.read_record().unwrap().unwrap();
    let line = record.to_sam_line(&reader.refs).unwrap();

    assert!(line.starts_with(
        "2f890ab0-63b9-45f6-ba8d-3c367ca26a63\t16\tfc_reference\t2\t60\t7S333M1I334M37S"
    ));
    assert!(line.contains("\tNM:i:8\t"));
}

#[test]
fn supported_sam_payload_converts_to_bam_and_back_to_sam() {
    let sam = Sam {
        header: SamHeader {
            records: vec![
                SamHeaderRecord::new(
                    "HD".to_string(),
                    vec![SamHeaderField::new("VN".to_string(), "1.6".to_string())],
                ),
                SamHeaderRecord::new(
                    "SQ".to_string(),
                    vec![
                        SamHeaderField::new("SN".to_string(), "ref".to_string()),
                        SamHeaderField::new("LN".to_string(), "10".to_string()),
                    ],
                ),
            ],
        },
        records: vec![SamRecord::new(
            "r1".to_string(),
            0,
            "ref".to_string(),
            1,
            60,
            "4M".to_string(),
            "*".to_string(),
            0,
            0,
            "ACGT".to_string(),
            "IIII".to_string(),
            Vec::new(),
        )],
    };

    let bam = Bam::from_sam(&sam).unwrap();
    assert_eq!(bam.refs[0].name, "ref");
    assert_eq!(
        bam.records[0].to_sam_record(&bam.refs).unwrap(),
        sam.records[0]
    );
}
