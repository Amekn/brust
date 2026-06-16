use pod5::Pod5;

const A_100_POD5: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/A_100.pod5");
const FIRST_READ_ID: &str = "1cadb1e9-592f-4e22-9285-4626f2b7da9f";

#[test]
fn read_lookup_and_summary_report_fixture_shape() {
    let pod5 = Pod5::from_path(A_100_POD5).unwrap();
    let lookup = pod5.read_lookup();
    let summary = pod5.summary();

    assert_eq!(lookup.len(), 100);
    assert!(lookup.contains_key(FIRST_READ_ID));
    assert_eq!(summary.read_count, 100);
    assert_eq!(summary.signal_count, 100);
    assert_eq!(summary.run_info_count, 1);
    assert_eq!(summary.total_samples, pod5.total_samples());
    assert_eq!(
        summary
            .channels
            .iter()
            .map(|channel| channel.read_count)
            .sum::<usize>(),
        100
    );
    assert_eq!(summary.run_infos[0].read_count, 100);
}

#[test]
fn materialized_signal_cache_reuses_decompressed_rows() {
    let pod5 = Pod5::from_path(A_100_POD5).unwrap();
    let record = pod5.read_by_id(FIRST_READ_ID).unwrap();
    let cache = pod5.signal_cache();

    let first = cache.signal_for_record(record).unwrap();
    let second = cache.signal_for_record(record).unwrap();

    assert_eq!(first, second);
    assert_eq!(cache.cached_row_count(), record.signal_rows.len());
}
