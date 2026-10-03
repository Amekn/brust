//! Validation helpers for the `brust` facade.
//!
//! Validation is intentionally parser-backed: each function streams through the
//! target format reader and relies on the lower-level crate to perform the
//! format checks it already owns. This avoids duplicating validation logic while
//! still surfacing structured `brust::Error` diagnostics.

use crate::{Error, Format, Result};
use std::collections::HashSet;
use std::path::Path;

/// Validates a file for a format selected at runtime.
pub fn validate<P: AsRef<Path>>(format: Format, input: P) -> Result<()> {
    match format {
        Format::Fasta => validate_fasta(input),
        Format::Fastq => validate_fastq(input),
        Format::Sam => validate_sam(input),
        Format::Bam => validate_bam(input),
        Format::Pod5 => validate_pod5(input),
    }
}

/// Validates plain or gzip FASTQ by streaming every record through the parser.
pub fn validate_fastq<P: AsRef<Path>>(input: P) -> Result<()> {
    let mut reader = fastq::FastqReader::from_path(input.as_ref())?;
    while let Some(_record) = reader.read_record()? {
        // Reading is validation; parser errors carry line/context where available.
    }
    Ok(())
}

/// Validates FASTA by streaming every record through the parser.
pub fn validate_fasta<P: AsRef<Path>>(input: P) -> Result<()> {
    let mut reader = fasta::FastaReader::from_path(input.as_ref())?;
    while let Some(_record) = reader.read_record()? {
        // Reading is validation; parser errors carry line/context where available.
    }
    Ok(())
}

/// Validates SAM headers and records by streaming through the parser.
///
/// When the header has `@SQ` lines, each record's RNAME and RNEXT (other than
/// `*` and `=`) must name one of them, as SAMv1 requires.
pub fn validate_sam<P: AsRef<Path>>(input: P) -> Result<()> {
    let mut reader = sam::SamReader::from_path(input.as_ref())?;
    let references = reader
        .header
        .records
        .iter()
        .filter(|record| record.record_type == "SQ")
        .filter_map(|record| record.value("SN").map(str::to_string))
        .collect::<HashSet<_>>();
    while let Some(record) = reader.read_record()? {
        // Reading is validation; parser errors carry line/context where available.
        if references.is_empty() {
            continue;
        }
        for name in [&record.rname, &record.rnext] {
            if name != "*" && name != "=" && !references.contains(name.as_str()) {
                return Err(Error::invalid(
                    Format::Sam,
                    format!(
                        "SAM record {} names reference {name}, which is missing from @SQ",
                        record.qname
                    ),
                ));
            }
        }
    }
    Ok(())
}

/// Validates BAM headers, references, BGZF blocks, and records by streaming,
/// and requires the BGZF EOF block.
///
/// The header text must make the SAM header `bam-to-sam` writes, each
/// record's read name must end in exactly the NUL its stored length counts,
/// and each record must make a valid SAM line (QUAL within 0-93, tag names,
/// read name characters and so on), so a BAM that validates converts to SAM.
pub fn validate_bam<P: AsRef<Path>>(input: P) -> Result<()> {
    let mut reader = bam::BamReader::from_path(input.as_ref())?;
    reader.set_require_eof_block(true);
    let refs = reader.refs.clone();
    // The header bam-to-sam would write, with the writer's checks.
    let header = crate::convert::sam_header_for_bam(&reader.header.text, &refs)?;
    sam::SamWriter::from_writer(std::io::sink()).write_header(&header)?;
    while let Some(record) = reader.read_record()? {
        // Reading is validation; parser errors carry context where available.
        let name = record.read_name();
        if usize::from(record.fixed.l_read_name) != name.len() + 1 {
            return Err(Error::invalid(
                Format::Bam,
                format!("BAM read name {name} does not end in the one NUL its length counts"),
            ));
        }
        record
            .to_sam_record(&refs)
            .and_then(|sam| sam.validate())
            .map_err(|error| {
                Error::invalid(
                    Format::Bam,
                    format!("BAM record {name} does not make a valid SAM line: {error}"),
                )
            })?;
    }
    Ok(())
}

/// Validates POD5 wrapper metadata and read rows by streaming through the parser.
pub fn validate_pod5<P: AsRef<Path>>(input: P) -> Result<()> {
    let mut reader = pod5::Pod5Reader::from_path(input.as_ref())?;
    while let Some(_record) = reader.read_record()? {
        // Reading is validation; parser errors carry context where available.
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEST_DIR_COUNTER: AtomicU64 = AtomicU64::new(0);

    struct TestDir {
        path: PathBuf,
    }

    impl TestDir {
        fn new() -> Self {
            let counter = TEST_DIR_COUNTER.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir().join(format!(
                "brust-validate-test-{}-{counter}",
                std::process::id()
            ));
            fs::create_dir(&path).unwrap();
            Self { path }
        }

        fn path(&self, name: &str) -> PathBuf {
            self.path.join(name)
        }
    }

    impl Drop for TestDir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.path);
        }
    }

    #[test]
    fn validate_dispatch_accepts_valid_text_formats() {
        // Dispatch should route to the same parser-backed validation as direct calls.
        let dir = TestDir::new();
        let fasta = dir.path("seqs.fasta");
        let fastq = dir.path("reads.fastq");
        let sam = dir.path("reads.sam");
        fs::write(&fasta, ">seq1 description\nACGT\n").unwrap();
        fs::write(&fastq, "@read1\nACGT\n+\nIIII\n").unwrap();
        fs::write(&sam, "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n").unwrap();

        validate(Format::Fasta, &fasta).unwrap();
        validate(Format::Fastq, &fastq).unwrap();
        validate(Format::Sam, &sam).unwrap();
    }

    #[test]
    fn validate_sam_requires_declared_references() {
        // SAMv1: when @SQ lines are present, RNAME and RNEXT must be among them.
        let dir = TestDir::new();
        let header = "@SQ\tSN:chr1\tLN:1000\n";
        for record in [
            "r1\t0\tchrX\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
            "r1\t0\tchr1\t1\t60\t4M\tchrY\t5\t0\tACGT\tIIII\n",
        ] {
            let input = dir.path("undeclared.sam");
            fs::write(&input, format!("{header}{record}")).unwrap();

            let error = validate_sam(&input).unwrap_err();
            assert_eq!(error.format(), Some(Format::Sam), "{record:?}");
        }

        let input = dir.path("declared.sam");
        fs::write(
            &input,
            format!("{header}r1\t0\tchr1\t1\t60\t4M\t=\t5\t0\tACGT\tIIII\n"),
        )
        .unwrap();
        validate_sam(&input).unwrap();
    }

    #[test]
    fn validate_sam_allows_any_reference_without_sq_lines() {
        let dir = TestDir::new();
        let input = dir.path("headerless.sam");
        fs::write(&input, "r1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n").unwrap();

        validate_sam(&input).unwrap();
    }

    /// Writes a one-record BAM after letting `edit` change the record.
    fn bam_with(dir: &TestDir, name: &str, edit: impl FnOnce(&mut bam::BamRecord)) -> PathBuf {
        let text = b"@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n";
        let sam = sam::Sam::from_reader(&text[..]).unwrap();
        let converter = bam::SamToBamConverter::new(&sam.header).unwrap();
        let mut record = converter.convert_record(&sam.records[0]).unwrap();
        edit(&mut record);
        record.fixed.block_size = 0;
        let path = dir.path(name);
        let mut writer = bam::BamWriter::from_path(&path).unwrap();
        writer
            .write_header(converter.header(), converter.refs())
            .unwrap();
        writer.write_record(&record).unwrap();
        writer.finish().unwrap();
        path
    }

    #[test]
    fn validate_bam_rejects_records_bam_to_sam_cannot_convert() {
        // samtools reads these, but they break SAM's rules for QUAL and tags.
        let dir = TestDir::new();
        let quality = bam_with(&dir, "quality.bam", |record| record.variable.qual[0] = 94);
        let tag = bam_with(&dir, "tag.bam", |record| {
            record.auxiliary.push(bam::BamRecordAuxiliary {
                tag: "1A".to_string(),
                value: bam::BamAuxValue::i(1),
            })
        });
        for path in [quality, tag] {
            let error = validate_bam(&path).unwrap_err();
            assert_eq!(error.format(), Some(Format::Bam), "{}", path.display());
        }
        validate_bam(bam_with(&dir, "ok.bam", |_| {})).unwrap();
    }

    #[test]
    fn validate_bam_rejects_header_text_bam_to_sam_cannot_convert() {
        let dir = TestDir::new();
        let text = "@SQ\tSN:ref\tLN:100\n@SQ\tSN:ref\tLN:100\n";
        let input = dir.path("duplicate-sq.bam");
        let sam = sam::Sam::from_reader(&b"@SQ\tSN:ref\tLN:100\n"[..]).unwrap();
        let converter = bam::SamToBamConverter::new(&sam.header).unwrap();
        let mut header = converter.header().clone();
        header.text = text.as_bytes().to_vec();
        header.l_text = 0;
        let mut writer = bam::BamWriter::from_path(&input).unwrap();
        writer.write_header(&header, converter.refs()).unwrap();
        writer.finish().unwrap();

        assert!(validate_bam(&input).is_err());
    }

    #[test]
    fn validate_bam_rejects_a_read_name_without_its_nul() {
        use std::io::{Read, Write};
        let dir = TestDir::new();
        let good = bam_with(&dir, "good.bam", |_| {});
        let mut raw = Vec::new();
        bam::BgzfReader::new(fs::File::open(&good).unwrap())
            .read_to_end(&mut raw)
            .unwrap();
        let name = raw.windows(3).position(|window| window == b"r1\0").unwrap();
        raw[name + 2] = b'X';
        let mut writer = bam::BgzfWriter::new(Vec::new());
        writer.write_all(&raw).unwrap();
        let input = dir.path("no-nul.bam");
        fs::write(&input, writer.finish().unwrap()).unwrap();

        let error = validate_bam(&input).unwrap_err();
        assert_eq!(error.format(), Some(Format::Bam));
    }

    #[test]
    fn validate_fastq_reports_parser_errors() {
        // The lower-level FASTQ parser should surface malformed records as FASTQ errors.
        let dir = TestDir::new();
        let input = dir.path("broken.fastq");
        fs::write(&input, "@read1\nACGT\n+\nIII\n").unwrap();

        let error = validate_fastq(&input).unwrap_err();

        assert_eq!(error.format(), Some(Format::Fastq));
    }

    #[test]
    fn validate_fastq_accepts_gzip_streams() {
        // Path-based validation inherits transparent decompression from FastqReader.
        let dir = TestDir::new();
        let input = dir.path("reads.fastq.gz");
        let mut writer = fastq::FastqWriter::from_path(&input).unwrap();
        writer
            .write_record(&fastq::FastqRecord::new(
                "read1".to_string(),
                None,
                "ACGT".to_string(),
                "IIII".to_string(),
            ))
            .unwrap();
        writer.finish().unwrap();

        validate_fastq(input).unwrap();
    }

    #[test]
    fn validate_sam_reports_parser_errors() {
        // Invalid mandatory fields are rejected while streaming alignment records.
        let dir = TestDir::new();
        let input = dir.path("broken.sam");
        fs::write(
            &input,
            "read1\tnot-a-flag\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
        )
        .unwrap();

        let error = validate_sam(&input).unwrap_err();

        assert_eq!(error.format(), Some(Format::Sam));
    }
}
