//! Format conversion helpers for the `brust` facade.
//!
//! All public conversions stream records through the relevant readers and
//! writers, avoiding whole-file materialization in the facade layer. Outputs are
//! written to a temporary file beside the target. On success the file is synced
//! to disk and renamed over the target, so an existing output file is never
//! replaced by a partial conversion when parsing or writing fails. A `.fq.gz` or
//! `.fastq.gz` target selects streaming gzip compression; FASTA and SAM
//! targets can't be compressed, so a `.gz` name for them is refused. BAM input
//! must end with the BGZF EOF block; a BAM without it is rejected as possibly
//! truncated.

use crate::{AtomicFile, Compression, Error, Format, Result};
use std::path::Path;

/// Supported file conversion paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Conversion {
    /// Convert FASTQ reads to FASTA records by dropping qualities.
    FastqToFasta,
    /// Convert FASTQ reads to unmapped SAM records.
    FastqToSam,
    /// Convert FASTQ reads to unmapped BAM records.
    FastqToBam,
    /// Convert SAM records to BAM.
    SamToBam,
    /// Convert BAM records to SAM.
    BamToSam,
    /// Convert SAM records with stored sequence and qualities to FASTQ.
    ///
    /// Reverse-strand reads are written in their original orientation,
    /// secondary and supplementary records are skipped, and read names are
    /// unchanged. Unmapped, QC-fail and duplicate records are still written.
    /// See [`sam_to_fastq`].
    SamToFastq,
    /// Convert BAM records with stored sequence and qualities to FASTQ.
    ///
    /// Reverse-strand reads are written in their original orientation,
    /// secondary and supplementary records are skipped, and read names are
    /// unchanged. Unmapped, QC-fail and duplicate records are still written.
    /// See [`bam_to_fastq`].
    BamToFastq,
}

impl Conversion {
    /// Stable kebab-case name used by the CLI.
    pub fn name(self) -> &'static str {
        match self {
            Self::FastqToFasta => "fastq-to-fasta",
            Self::FastqToSam => "fastq-to-sam",
            Self::FastqToBam => "fastq-to-bam",
            Self::SamToBam => "sam-to-bam",
            Self::BamToSam => "bam-to-sam",
            Self::SamToFastq => "sam-to-fastq",
            Self::BamToFastq => "bam-to-fastq",
        }
    }

    /// Input format consumed by this conversion.
    pub fn input_format(self) -> Format {
        match self {
            Self::FastqToFasta | Self::FastqToSam | Self::FastqToBam => Format::Fastq,
            Self::SamToBam | Self::SamToFastq => Format::Sam,
            Self::BamToSam | Self::BamToFastq => Format::Bam,
        }
    }

    /// Output format produced by this conversion.
    pub fn output_format(self) -> Format {
        match self {
            Self::FastqToFasta => Format::Fasta,
            Self::FastqToSam | Self::BamToSam => Format::Sam,
            Self::FastqToBam | Self::SamToBam => Format::Bam,
            Self::SamToFastq | Self::BamToFastq => Format::Fastq,
        }
    }
}

/// Options for [`convert_with`].
///
/// Build one with [`ConvertOptions::default`] and the builder methods; the
/// struct is `#[non_exhaustive]`, so it cannot be built with a struct literal
/// outside this crate and new options can be added without breaking callers.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct ConvertOptions {
    /// BGZF compression threads used when the output is BAM.
    ///
    /// `0` and `1` compress on the calling thread. This has no effect on
    /// conversions that do not write BAM. Output bytes are identical for any
    /// thread count.
    pub threads: usize,
}

impl Default for ConvertOptions {
    /// One thread, so compression runs inline on the calling thread.
    fn default() -> Self {
        Self { threads: 1 }
    }
}

impl ConvertOptions {
    /// Sets the number of BGZF compression threads.
    ///
    /// `threads` affects only conversions that write BAM
    /// ([`Conversion::FastqToBam`] and [`Conversion::SamToBam`]); every other
    /// conversion ignores it. `0` and `1` compress on the calling thread. The
    /// output is byte-identical for any thread count.
    #[must_use]
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }
}

/// Converts between supported formats selected at runtime.
///
/// This is [`convert_with`] using [`ConvertOptions::default`].
pub fn convert<I: AsRef<Path>, O: AsRef<Path>>(
    conversion: Conversion,
    input: I,
    output: O,
) -> Result<()> {
    convert_with(conversion, input, output, &ConvertOptions::default())
}

/// Converts between supported formats selected at runtime, with options.
///
/// `options.threads` affects only conversions that write BAM
/// ([`Conversion::FastqToBam`] and [`Conversion::SamToBam`]); the other
/// conversions ignore it. The output is byte-identical for any thread count.
/// As with [`convert`], the output is written to a temporary file, synced to
/// disk and renamed on success, so a conversion that fails before the rename
/// leaves an existing output file untouched. After the rename, only syncing
/// the output's folder can still fail (see [`AtomicFile::commit`]); the new
/// file is then already in place. BAM input must end with the BGZF EOF block.
pub fn convert_with<I: AsRef<Path>, O: AsRef<Path>>(
    conversion: Conversion,
    input: I,
    output: O,
    options: &ConvertOptions,
) -> Result<()> {
    match conversion {
        Conversion::FastqToFasta => fastq_to_fasta(input, output),
        Conversion::FastqToSam => fastq_to_sam(input, output),
        Conversion::FastqToBam => fastq_to_bam_with(input, output, options.threads),
        Conversion::SamToBam => sam_to_bam_with(input, output, options.threads),
        Conversion::BamToSam => bam_to_sam(input, output),
        Conversion::SamToFastq => sam_to_fastq(input, output),
        Conversion::BamToFastq => bam_to_fastq(input, output),
    }
}

/// Converts FASTQ records to FASTA records, dropping quality scores.
pub fn fastq_to_fasta<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    require_uncompressed_output(output.as_ref(), Format::Fasta)?;
    let input = input.as_ref();
    let mut reader = fastq::FastqReader::from_path(input)?;
    let mut writer = fasta::FastaWriter::from_path_atomic(output.as_ref())?;
    while let Some(record) = reader.read_record()? {
        writer.write_record(&record.to_fasta_record())?;
    }
    writer.commit()?;
    Ok(())
}

/// Converts FASTQ reads to unmapped SAM records.
pub fn fastq_to_sam<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    require_uncompressed_output(output.as_ref(), Format::Sam)?;
    let input = input.as_ref();
    let mut reader = fastq::FastqReader::from_path(input)?;
    let mut writer = sam::SamWriter::from_path_atomic(output.as_ref())?;
    while let Some(record) = reader.read_record()? {
        writer.write_record(&fastq_record_to_unmapped_sam(&record))?;
    }
    writer.commit()?;
    Ok(())
}

/// Converts FASTQ reads to unmapped BAM records.
pub fn fastq_to_bam<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    fastq_to_bam_with(input, output, 1)
}

fn fastq_to_bam_with<I: AsRef<Path>, O: AsRef<Path>>(
    input: I,
    output: O,
    threads: usize,
) -> Result<()> {
    let input = input.as_ref();
    let header = sam::SamHeader::default();
    let converter = bam::SamToBamConverter::new(&header)?;
    let mut reader = fastq::FastqReader::from_path(input)?;
    let mut writer = bam::BamWriter::from_path_atomic_with_threads(output.as_ref(), threads)?;

    writer.write_header(converter.header(), converter.refs())?;
    while let Some(record) = reader.read_record()? {
        let mut bam_record = converter.convert_record(&fastq_record_to_unmapped_sam(&record))?;
        // In SAM text a one-base quality of `*` (Q9) means "no quality", so
        // take the scores from the FASTQ record itself.
        if !record.quality.is_empty() {
            bam_record.variable.qual = record.quality.bytes().map(|score| score - 33).collect();
        }
        writer.write_record(&bam_record)?;
    }
    writer.commit()?;
    Ok(())
}

/// Converts a supported SAM payload to BAM.
pub fn sam_to_bam<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    sam_to_bam_with(input, output, 1)
}

fn sam_to_bam_with<I: AsRef<Path>, O: AsRef<Path>>(
    input: I,
    output: O,
    threads: usize,
) -> Result<()> {
    let input = input.as_ref();
    let mut reader = sam::SamReader::from_path(input)?;
    let converter = bam::SamToBamConverter::new(&reader.header)?;
    let mut writer = bam::BamWriter::from_path_atomic_with_threads(output.as_ref(), threads)?;

    writer.write_header(converter.header(), converter.refs())?;
    while let Some(record) = reader.read_record()? {
        let record = converter.convert_record(&record)?;
        writer.write_record(&record)?;
    }
    writer.commit()?;
    Ok(())
}

/// Converts BAM to SAM.
///
/// The input must end with the BGZF EOF block.
pub fn bam_to_sam<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    require_uncompressed_output(output.as_ref(), Format::Sam)?;
    let input = input.as_ref();
    let mut reader = bam::BamReader::from_path(input)?;
    reader.set_require_eof_block(true);
    let refs = reader.refs.clone();
    let header = sam_header_for_bam(&reader.header.text, &refs)?;
    let mut writer = sam::SamWriter::from_writer(AtomicFile::create(output.as_ref())?);
    writer.write_header(&header)?;
    while let Some(record) = reader.read_record()? {
        writer.write_record(&record.to_sam_record(&refs)?)?;
    }
    writer.into_inner().commit()?;
    Ok(())
}

/// Converts SAM records with stored sequence and qualities to FASTQ.
///
/// Reverse-strand reads (flag `0x10`) are written in their original
/// orientation: SEQ is reverse-complemented and QUAL reversed, keeping case.
/// Secondary (`0x100`) and supplementary (`0x800`) records are skipped, even
/// when their SEQ or QUAL is `*`. Unmapped (`0x4`), QC-fail (`0x200`) and
/// duplicate (`0x400`) records are still written. Read names are unchanged,
/// with no `/1` or `/2` suffix. Any other record without SEQ or QUAL is an
/// error.
///
/// A `.fq.gz` or `.fastq.gz` output path is compressed while records stream.
pub fn sam_to_fastq<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    let input = input.as_ref();
    let mut reader = sam::SamReader::from_path(input)?;
    let mut writer = fastq::FastqWriter::from_path_atomic(output.as_ref())?;
    while let Some(record) = reader.read_record()? {
        if let Some(read) = sam_record_to_fastq(&record, record.qual == "*", Format::Sam)? {
            writer.write_record(&read)?;
        }
    }
    writer.commit()?;
    Ok(())
}

/// Converts BAM records with stored sequence and qualities to FASTQ.
///
/// Reverse-strand reads (flag `0x10`) are written in their original
/// orientation: SEQ is reverse-complemented and QUAL reversed, keeping case.
/// Secondary (`0x100`) and supplementary (`0x800`) records are skipped, even
/// when their SEQ or QUAL is `*`. Unmapped (`0x4`), QC-fail (`0x200`) and
/// duplicate (`0x400`) records are still written. Read names are unchanged,
/// with no `/1` or `/2` suffix. Any other record without SEQ or QUAL is an
/// error.
///
/// A `.fq.gz` or `.fastq.gz` output path is compressed while records stream.
pub fn bam_to_fastq<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    let input = input.as_ref();
    let mut reader = bam::BamReader::from_path(input)?;
    reader.set_require_eof_block(true);
    let refs = reader.refs.clone();
    let mut writer = fastq::FastqWriter::from_path_atomic(output.as_ref())?;
    while let Some(record) = reader.read_record()? {
        // Missing quality is all 0xFF in BAM; the SAM text `*` would also
        // match a one-base read of quality 9.
        let quality_missing = record.variable.qual.iter().all(|&score| score == 0xff);
        let record = record.to_sam_record(&refs)?;
        if let Some(read) = sam_record_to_fastq(&record, quality_missing, Format::Bam)? {
            writer.write_record(&read)?;
        }
    }
    writer.commit()?;
    Ok(())
}

/// Builds the SAM header for a BAM: its header text, parsed and checked like
/// any SAM header, plus `@SQ` lines from the binary reference list when the
/// text has none (the list is authoritative, and SAM records need them).
pub(crate) fn sam_header_for_bam(text: &[u8], refs: &[bam::BamRef]) -> Result<sam::SamHeader> {
    let invalid = |message: String| Error::invalid(Format::Bam, message);
    // SAM text must be UTF-8; the BAM header keeps its bytes as stored.
    let text = std::str::from_utf8(text)
        .map_err(|error| invalid(format!("BAM header text is not UTF-8: {error}")))?;
    let parsed = sam::Sam::from_reader(text.as_bytes())
        .map_err(|error| invalid(format!("BAM header text is not a SAM header: {error}")))?;
    if !parsed.records.is_empty() {
        return Err(invalid(
            "BAM header text holds a line that is not a header line".to_string(),
        ));
    }
    let mut header = parsed.header;
    if !header
        .records
        .iter()
        .any(|record| record.record_type == "SQ")
    {
        let at = usize::from(
            header
                .records
                .first()
                .is_some_and(|record| record.record_type == "HD"),
        );
        let sequences = refs.iter().map(|reference| {
            sam::SamHeaderRecord::new(
                "SQ".to_string(),
                vec![
                    sam::SamHeaderField::new("SN".to_string(), reference.name.clone()),
                    sam::SamHeaderField::new("LN".to_string(), reference.l_seq.to_string()),
                ],
            )
        });
        header.records.splice(at..at, sequences);
    }
    Ok(header)
}

/// Refuses a `.gz` output name for a format whose writer can't compress, so
/// a file named `.gz` never holds plain text.
fn require_uncompressed_output(output: &Path, format: Format) -> Result<()> {
    if Compression::from_path(output) == Compression::Gzip {
        return Err(Error::Io(format!(
            "{} output can't be gzip-compressed; remove the .gz suffix from {}",
            format.name(),
            output.display()
        )));
    }
    Ok(())
}

fn fastq_record_to_unmapped_sam(record: &fastq::FastqRecord) -> sam::SamRecord {
    // SAM has no empty SEQ or QUAL field; a zero-length read stores `*`.
    let (sequence, quality) = if record.sequence.is_empty() {
        ("*".to_string(), "*".to_string())
    } else {
        (record.sequence.clone(), record.quality.clone())
    };
    sam::SamRecord::new(
        record.id.clone(),
        sam::flags::UNMAPPED,
        "*".to_string(),
        0,
        255,
        "*".to_string(),
        "*".to_string(),
        0,
        0,
        sequence,
        quality,
        Vec::new(),
    )
}

/// Builds the FASTQ read for a SAM record, or `None` for a secondary or
/// supplementary record, which duplicates the read written from its primary
/// record.
fn sam_record_to_fastq(
    record: &sam::SamRecord,
    quality_missing: bool,
    error_format: Format,
) -> Result<Option<fastq::FastqRecord>> {
    if record.is_secondary() || record.is_supplementary() {
        return Ok(None);
    }
    if record.seq == "*" {
        return Err(Error::invalid(
            error_format,
            "cannot convert record without SEQ to FASTQ",
        ));
    }
    if quality_missing {
        return Err(Error::invalid(
            error_format,
            "cannot convert record without QUAL to FASTQ",
        ));
    }

    Ok(Some(fastq::FastqRecord::new(
        record.qname.clone(),
        None,
        record.original_seq(),
        record.original_qual(),
    )))
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
                "brust-convert-test-{}-{counter}",
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
    fn conversion_metadata_matches_supported_paths() {
        // The CLI relies on these stable names and format pairs for dispatch.
        assert_eq!(Conversion::FastqToFasta.name(), "fastq-to-fasta");
        assert_eq!(Conversion::FastqToFasta.input_format(), Format::Fastq);
        assert_eq!(Conversion::FastqToFasta.output_format(), Format::Fasta);
        assert_eq!(Conversion::SamToBam.name(), "sam-to-bam");
        assert_eq!(Conversion::SamToBam.input_format(), Format::Sam);
        assert_eq!(Conversion::SamToBam.output_format(), Format::Bam);
        assert_eq!(Conversion::BamToFastq.name(), "bam-to-fastq");
        assert_eq!(Conversion::BamToFastq.input_format(), Format::Bam);
        assert_eq!(Conversion::BamToFastq.output_format(), Format::Fastq);
    }

    #[test]
    fn runtime_convert_fastq_to_fasta_drops_quality_scores() {
        // A tiny fixture keeps the conversion contract readable in the assertion.
        let dir = TestDir::new();
        let input = dir.path("reads.fastq");
        let output = dir.path("reads.fasta");
        fs::write(&input, "@read1 sample\nACGT\n+\nIIII\n").unwrap();

        convert(Conversion::FastqToFasta, &input, &output).unwrap();

        assert_eq!(fs::read_to_string(output).unwrap(), ">read1 sample\nACGT\n");
    }

    #[test]
    fn zero_length_fastq_reads_convert_to_sam_and_bam_without_seq() {
        // SAM stores an empty read as SEQ and QUAL `*`, as samtools import does.
        let dir = TestDir::new();
        let input = dir.path("reads.fastq");
        fs::write(&input, "@empty\n\n+\n\n@read1\nACGT\n+\nIIII\n").unwrap();

        convert(Conversion::FastqToSam, &input, dir.path("reads.sam")).unwrap();
        convert(Conversion::FastqToBam, &input, dir.path("reads.bam")).unwrap();

        let sam = sam::Sam::from_path(dir.path("reads.sam")).unwrap();
        assert_eq!(
            (sam.records[0].seq.as_str(), sam.records[0].qual.as_str()),
            ("*", "*")
        );
        assert_eq!(sam.records[1].seq, "ACGT");
        let bam = bam::Bam::from_path(dir.path("reads.bam")).unwrap();
        assert_eq!(bam.records[0].fixed.l_seq, 0);
        assert_eq!(bam.records[1].fixed.l_seq, 4);
    }

    #[test]
    fn gzip_named_fasta_and_sam_outputs_are_refused() {
        // These writers can't compress, so a `.gz` name would hold plain text.
        let dir = TestDir::new();
        let input = dir.path("reads.fastq");
        fs::write(&input, "@read1\nACGT\n+\nIIII\n").unwrap();

        for (conversion, name) in [
            (Conversion::FastqToFasta, "reads.fa.gz"),
            (Conversion::FastqToSam, "reads.sam.gz"),
        ] {
            let output = dir.path(name);
            let error = convert(conversion, &input, &output).unwrap_err();
            assert!(matches!(error, Error::Io(_)), "{error:?}");
            assert!(!output.exists(), "{name}");
        }
        convert(Conversion::FastqToFasta, &input, dir.path("reads.fa")).unwrap();
    }

    /// Writes a BAM with the given header text, one `ref` reference of length
    /// 100, and one record placed on it.
    fn bam_with_header_text(dir: &TestDir, text: &[u8]) -> PathBuf {
        let header = sam::Sam::from_reader(&b"@SQ\tSN:ref\tLN:100\n"[..])
            .unwrap()
            .header;
        let converter = bam::SamToBamConverter::new(&header).unwrap();
        let record = sam::Sam::from_reader(
            &b"@SQ\tSN:ref\tLN:100\nr1\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"[..],
        )
        .unwrap()
        .records
        .remove(0);
        let path = dir.path("in.bam");
        let mut writer = bam::BamWriter::from_path(&path).unwrap();
        let bam_header = bam::BamHeader {
            magic: *b"BAM\x01",
            l_text: 0,
            text: text.to_vec(),
            n_ref: 1,
        };
        writer.write_header(&bam_header, converter.refs()).unwrap();
        writer
            .write_record(&converter.convert_record(&record).unwrap())
            .unwrap();
        writer.finish().unwrap();
        path
    }

    #[test]
    fn bam_to_sam_declares_references_missing_from_the_header_text() {
        // The binary reference list is authoritative; SAM needs @SQ lines.
        let dir = TestDir::new();
        for (text, expected) in [
            ("", "@SQ\tSN:ref\tLN:100\n"),
            ("@HD\tVN:1.6\n", "@HD\tVN:1.6\n@SQ\tSN:ref\tLN:100\n"),
        ] {
            let input = bam_with_header_text(&dir, text.as_bytes());
            let output = dir.path("out.sam");

            bam_to_sam(&input, &output).unwrap();

            let sam = fs::read_to_string(&output).unwrap();
            assert!(sam.starts_with(expected), "{sam:?}");
            assert!(crate::validate::validate_sam(&output).is_ok());
        }
    }

    #[test]
    fn bam_to_sam_refuses_header_text_that_is_not_utf8() {
        // SAM text must be UTF-8; the bytes used to be replaced with U+FFFD.
        let dir = TestDir::new();
        let input = bam_with_header_text(&dir, b"@HD\tVN:1.6\n@CO\tcaf\xe9\n");
        let output = dir.path("latin1.sam");

        let error = bam_to_sam(&input, &output).unwrap_err();
        assert_eq!(error.format(), Some(Format::Bam), "{error}");
        assert!(!output.exists());
    }

    #[test]
    fn bam_to_sam_rejects_header_text_that_is_not_a_sam_header() {
        let dir = TestDir::new();
        for text in [
            "@HD\tVN:1.6\n\n@SQ\tSN:ref\tLN:100\n",
            "@SQ\tSN:ref\tLN:100\n@SQ\tSN:ref\tLN:100\n",
            "@SQ\tSN:ref\tLN:100\nnot a header line\n",
        ] {
            let input = bam_with_header_text(&dir, text.as_bytes());
            let output = dir.path("bad.sam");

            assert!(bam_to_sam(&input, &output).is_err(), "{text:?}");
            assert!(!output.exists());
        }
    }

    #[test]
    fn single_base_q9_reads_keep_their_quality_through_bam() {
        // In SAM text a lone quality `*` (Q9) reads as "no quality".
        let dir = TestDir::new();
        let input = dir.path("q9.fastq");
        fs::write(&input, "@q9\nA\n+\n*\n").unwrap();
        let bam_path = dir.path("q9.bam");
        let output = dir.path("back.fastq");

        convert(Conversion::FastqToBam, &input, &bam_path).unwrap();
        let bam = bam::Bam::from_path(&bam_path).unwrap();
        assert_eq!(bam.records[0].variable.qual, [9]);
        convert(Conversion::BamToFastq, &bam_path, &output).unwrap();

        assert_eq!(fs::read_to_string(output).unwrap(), "@q9\nA\n+\n*\n");
    }

    #[test]
    fn sam_to_fastq_compresses_gzip_destination() {
        // Atomic temporary output retains `.gz`, so the FASTQ writer compresses it.
        let dir = TestDir::new();
        let input = dir.path("reads.sam");
        let output = dir.path("reads.fq.gz");
        fs::write(&input, "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n").unwrap();

        sam_to_fastq(&input, &output).unwrap();

        assert_eq!(&fs::read(&output).unwrap()[..2], &[0x1f, 0x8b]);
        let fastq = fastq::Fastq::from_path(output).unwrap();
        assert_eq!(fastq.records.len(), 1);
        assert_eq!(fastq.records[0].sequence, "ACGT");
    }

    #[test]
    fn sam_to_fastq_rejects_records_without_quality_scores() {
        // FASTQ output requires both SEQ and QUAL, even when SAM parsing succeeds.
        let dir = TestDir::new();
        let input = dir.path("reads.sam");
        let output = dir.path("reads.fastq");
        fs::write(&input, "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n").unwrap();

        let error = sam_to_fastq(&input, &output).unwrap_err();

        assert_eq!(error.format(), Some(Format::Sam));
        assert!(!output.exists());
    }

    #[test]
    fn failed_conversion_preserves_existing_output_file() {
        // Atomic writes must leave the original output untouched on parse failure.
        let dir = TestDir::new();
        let input = dir.path("broken.fastq");
        let output = dir.path("reads.fasta");
        fs::write(&input, "@read1\nACGT\n+\nIII\n").unwrap();
        fs::write(&output, "keep me\n").unwrap();

        let error = fastq_to_fasta(&input, &output).unwrap_err();

        assert_eq!(error.format(), Some(Format::Fastq));
        assert_eq!(fs::read_to_string(output).unwrap(), "keep me\n");
    }
}
