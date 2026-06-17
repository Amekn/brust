//! Format conversion helpers for the `brust` facade.
//!
//! All public conversions stream records through the relevant readers and
//! writers, avoiding whole-file materialization in the facade layer. Outputs are
//! written through a temporary file that is renamed only after successful
//! completion, protecting existing output files from being replaced by a partial
//! conversion when parsing or writing fails.

use crate::{Error, Format, Result};
use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

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
    SamToFastq,
    /// Convert BAM records with stored sequence and qualities to FASTQ.
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

/// Converts between supported formats selected at runtime.
pub fn convert<I: AsRef<Path>, O: AsRef<Path>>(
    conversion: Conversion,
    input: I,
    output: O,
) -> Result<()> {
    match conversion {
        Conversion::FastqToFasta => fastq_to_fasta(input, output),
        Conversion::FastqToSam => fastq_to_sam(input, output),
        Conversion::FastqToBam => fastq_to_bam(input, output),
        Conversion::SamToBam => sam_to_bam(input, output),
        Conversion::BamToSam => bam_to_sam(input, output),
        Conversion::SamToFastq => sam_to_fastq(input, output),
        Conversion::BamToFastq => bam_to_fastq(input, output),
    }
}

/// Converts FASTQ records to FASTA records, dropping quality scores.
pub fn fastq_to_fasta<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    let input = input.as_ref();
    write_atomic(output.as_ref(), |temp_output| {
        let mut reader = fastq::FastqReader::from_path(input)?;
        let mut writer = fasta::FastaWriter::from_path(temp_output)?;
        while let Some(record) = reader.read_record()? {
            writer.write_record(&record.to_fasta_record())?;
        }
        writer.flush()?;
        Ok(())
    })
}

/// Converts FASTQ reads to unmapped SAM records.
pub fn fastq_to_sam<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    let input = input.as_ref();
    write_atomic(output.as_ref(), |temp_output| {
        let mut reader = fastq::FastqReader::from_path(input)?;
        let mut writer = sam::SamWriter::from_path(temp_output)?;
        while let Some(record) = reader.read_record()? {
            writer.write_record(&fastq_record_to_unmapped_sam(&record))?;
        }
        writer.flush()?;
        Ok(())
    })
}

/// Converts FASTQ reads to unmapped BAM records.
pub fn fastq_to_bam<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    let input = input.as_ref();
    write_atomic(output.as_ref(), |temp_output| {
        let header = sam::SamHeader::default();
        let converter = bam::SamToBamConverter::new(&header)?;
        let mut reader = fastq::FastqReader::from_path(input)?;
        let mut writer = bam::BamWriter::from_path(temp_output)?;

        writer.write_header(converter.header(), converter.refs())?;
        while let Some(record) = reader.read_record()? {
            let record = fastq_record_to_unmapped_sam(&record);
            let record = converter.convert_record(&record)?;
            writer.write_record(&record)?;
        }
        writer.finish()?;
        Ok(())
    })
}

/// Converts a supported SAM payload to BAM.
pub fn sam_to_bam<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    let input = input.as_ref();
    write_atomic(output.as_ref(), |temp_output| {
        let mut reader = sam::SamReader::from_path(input)?;
        let converter = bam::SamToBamConverter::new(&reader.header)?;
        let mut writer = bam::BamWriter::from_path(temp_output)?;

        writer.write_header(converter.header(), converter.refs())?;
        while let Some(record) = reader.read_record()? {
            let record = converter.convert_record(&record)?;
            writer.write_record(&record)?;
        }
        writer.finish()?;
        Ok(())
    })
}

/// Converts BAM to SAM.
pub fn bam_to_sam<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    let input = input.as_ref();
    write_atomic(output.as_ref(), |temp_output| {
        let mut reader = bam::BamReader::from_path(input)?;
        let mut file = File::create(temp_output)?;
        if !reader.header.text.is_empty() {
            file.write_all(reader.header.text.as_bytes())?;
            if !reader.header.text.ends_with('\n') {
                file.write_all(b"\n")?;
            }
        }
        let refs = reader.refs.clone();
        let mut writer = sam::SamWriter::from_writer(file);
        while let Some(record) = reader.read_record()? {
            writer.write_record(&record.to_sam_record(&refs)?)?;
        }
        writer.flush()?;
        Ok(())
    })
}

/// Converts SAM records with stored sequence and qualities to FASTQ.
pub fn sam_to_fastq<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    let input = input.as_ref();
    write_atomic(output.as_ref(), |temp_output| {
        let mut reader = sam::SamReader::from_path(input)?;
        let mut writer = fastq::FastqWriter::from_path(temp_output)?;
        while let Some(record) = reader.read_record()? {
            writer.write_record(&sam_record_to_fastq(&record, Format::Sam)?)?;
        }
        writer.flush()?;
        Ok(())
    })
}

/// Converts BAM records with stored sequence and qualities to FASTQ.
pub fn bam_to_fastq<I: AsRef<Path>, O: AsRef<Path>>(input: I, output: O) -> Result<()> {
    let input = input.as_ref();
    write_atomic(output.as_ref(), |temp_output| {
        let mut reader = bam::BamReader::from_path(input)?;
        let refs = reader.refs.clone();
        let mut writer = fastq::FastqWriter::from_path(temp_output)?;
        while let Some(record) = reader.read_record()? {
            let record = record.to_sam_record(&refs)?;
            writer.write_record(&sam_record_to_fastq(&record, Format::Bam)?)?;
        }
        writer.flush()?;
        Ok(())
    })
}

fn fastq_record_to_unmapped_sam(record: &fastq::FastqRecord) -> sam::SamRecord {
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
        record.sequence.clone(),
        record.quality.clone(),
        Vec::new(),
    )
}

fn sam_record_to_fastq(
    record: &sam::SamRecord,
    error_format: Format,
) -> Result<fastq::FastqRecord> {
    if record.seq == "*" {
        return Err(Error::invalid(
            error_format,
            "cannot convert record without SEQ to FASTQ",
        ));
    }
    if record.qual == "*" {
        return Err(Error::invalid(
            error_format,
            "cannot convert record without QUAL to FASTQ",
        ));
    }

    Ok(fastq::FastqRecord::new(
        record.qname.clone(),
        None,
        record.seq.clone(),
        record.qual.clone(),
    ))
}

fn write_atomic(output: &Path, write: impl FnOnce(&Path) -> Result<()>) -> Result<()> {
    let temp_output = create_temp_output_path(output)?;

    match write(&temp_output) {
        Ok(()) => {
            if let Err(error) = fs::rename(&temp_output, output) {
                let _ = fs::remove_file(&temp_output);
                return Err(error.into());
            }
            Ok(())
        }
        Err(error) => {
            let _ = fs::remove_file(&temp_output);
            Err(error)
        }
    }
}

fn create_temp_output_path(output: &Path) -> Result<PathBuf> {
    let parent = output
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let file_name = output.file_name().ok_or_else(|| {
        Error::Io(format!(
            "output path {} must include a file name",
            output.display()
        ))
    })?;

    for _ in 0..100 {
        let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let mut temp_name = OsString::from(".");
        temp_name.push(file_name);
        temp_name.push(format!(".{}.{}.tmp", process::id(), counter));
        let temp_path = parent.join(temp_name);

        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temp_path)
        {
            Ok(_) => return Ok(temp_path),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(error.into()),
        }
    }

    Err(Error::Io(format!(
        "could not create temporary output beside {}",
        output.display()
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
