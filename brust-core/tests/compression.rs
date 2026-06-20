use brust_core::Compression;
use std::path::Path;

#[test]
fn public_compression_api_handles_fastq_paths() {
    // This integration test keeps the path contract visible to downstream crates.
    assert_eq!(
        Compression::from_path(Path::new("run/reads.fq.gz")),
        Compression::Gzip
    );
    assert_eq!(
        Compression::from_path(Path::new("run/reads.fastq")),
        Compression::Uncompressed
    );
}
