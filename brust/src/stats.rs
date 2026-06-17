//! Structured statistics for the formats exposed by the `brust` facade.
//!
//! The functions in this module stream records through the format readers and
//! return typed summary structs. The lower-level parsers remain responsible for
//! format validation; stats collection avoids repeating expensive validation
//! work that the readers already perform.

use crate::{Format, Result};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;
use std::fs;
use std::path::Path;

/// Statistics for any supported Brust format.
#[derive(Debug, Clone, PartialEq)]
// Keep variants unboxed so callers can pattern-match directly on typed stats.
#[allow(clippy::large_enum_variant)]
pub enum Stats {
    /// FASTA sequence statistics.
    Fasta(FastaStats),
    /// FASTQ sequence and quality statistics.
    Fastq(FastqStats),
    /// SAM header and alignment statistics.
    Sam(SamStats),
    /// BAM header, alignment, and record-size statistics.
    Bam(BamStats),
    /// POD5 read, run, channel, and signal-row metadata statistics.
    Pod5(Pod5Stats),
}

impl Stats {
    /// Returns the file format represented by this stats value.
    pub fn format(&self) -> Format {
        match self {
            Self::Fasta(_) => Format::Fasta,
            Self::Fastq(_) => Format::Fastq,
            Self::Sam(_) => Format::Sam,
            Self::Bam(_) => Format::Bam,
            Self::Pod5(_) => Format::Pod5,
        }
    }

    /// Returns a display adapter for human-readable CLI output.
    pub fn display(&self) -> StatsDisplay<'_> {
        StatsDisplay(self)
    }
}

/// Display adapter returned by [`Stats::display`].
pub struct StatsDisplay<'a>(&'a Stats);

impl fmt::Display for StatsDisplay<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(self.0, formatter)
    }
}

impl fmt::Display for Stats {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Fasta(stats) => fmt::Display::fmt(stats, formatter),
            Self::Fastq(stats) => fmt::Display::fmt(stats, formatter),
            Self::Sam(stats) => fmt::Display::fmt(stats, formatter),
            Self::Bam(stats) => fmt::Display::fmt(stats, formatter),
            Self::Pod5(stats) => fmt::Display::fmt(stats, formatter),
        }
    }
}

/// Computes stats for a format selected at runtime.
pub fn stats<P: AsRef<Path>>(format: Format, input: P) -> Result<Stats> {
    match format {
        Format::Fasta => fasta_stats(input).map(Stats::Fasta),
        Format::Fastq => fastq_stats(input).map(Stats::Fastq),
        Format::Sam => sam_stats(input).map(Stats::Sam),
        Format::Bam => bam_stats(input).map(Stats::Bam),
        Format::Pod5 => pod5_stats(input).map(Stats::Pod5),
    }
}

/// Summary of unsigned lengths/count-like values.
#[derive(Debug, Clone, PartialEq)]
pub struct LengthStats {
    /// Number of observed values.
    pub count: u64,
    /// Sum of all observed values.
    pub total: u64,
    /// Smallest observed value.
    pub min: Option<u64>,
    /// Largest observed value.
    pub max: Option<u64>,
    /// Arithmetic mean of observed values.
    pub mean: Option<f64>,
    /// N50 value, computed from values sorted descending by length.
    pub n50: Option<u64>,
    /// N90 value, computed from values sorted descending by length.
    pub n90: Option<u64>,
}

/// Summary of signed numeric values.
#[derive(Debug, Clone, PartialEq)]
pub struct NumberStats {
    /// Number of observed values.
    pub count: u64,
    /// Smallest observed value.
    pub min: Option<i64>,
    /// Largest observed value.
    pub max: Option<i64>,
    /// Arithmetic mean of observed values.
    pub mean: Option<f64>,
}

/// Summary of floating-point values.
#[derive(Debug, Clone, PartialEq)]
pub struct FloatStats {
    /// Number of finite observed values included in min/max/mean.
    pub count: u64,
    /// Number of NaN or infinite values skipped from min/max/mean.
    pub non_finite_count: u64,
    /// Smallest finite observed value.
    pub min: Option<f64>,
    /// Largest finite observed value.
    pub max: Option<f64>,
    /// Arithmetic mean of finite observed values.
    pub mean: Option<f64>,
}

/// Base/residue composition for sequence-bearing formats.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BaseComposition {
    /// Total number of observed sequence symbols.
    pub total: u64,
    /// Count of `A` symbols.
    pub a: u64,
    /// Count of `C` symbols.
    pub c: u64,
    /// Count of `G` symbols.
    pub g: u64,
    /// Count of `T` symbols.
    pub t: u64,
    /// Count of `U` symbols.
    pub u: u64,
    /// Count of `N` symbols.
    pub n: u64,
    /// Count of IUPAC ambiguous nucleotide symbols other than `N`.
    pub ambiguous: u64,
    /// Count of gap symbols (`-` and `.`).
    pub gap: u64,
    /// Count of symbols outside the common nucleotide/IUPAC buckets.
    pub other: u64,
}

impl BaseComposition {
    /// Returns the number of `G` plus `C` symbols.
    pub fn gc_count(&self) -> u64 {
        self.g + self.c
    }

    /// Returns the GC fraction, or `None` when no bases were observed.
    pub fn gc_fraction(&self) -> Option<f64> {
        (self.total > 0).then(|| self.gc_count() as f64 / self.total as f64)
    }

    fn add_sequence(&mut self, sequence: &str) {
        for byte in sequence.bytes() {
            self.total += 1;
            match byte.to_ascii_uppercase() {
                b'A' => self.a += 1,
                b'C' => self.c += 1,
                b'G' => self.g += 1,
                b'T' => self.t += 1,
                b'U' => self.u += 1,
                b'N' => self.n += 1,
                b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' => {
                    self.ambiguous += 1;
                }
                b'-' | b'.' => self.gap += 1,
                _ => self.other += 1,
            }
        }
    }
}

/// Quality-score summary for FASTQ-like Phred+33 qualities.
#[derive(Debug, Clone, PartialEq)]
pub struct QualityStats {
    /// Number of reads contributing quality data.
    pub read_count: u64,
    /// Number of quality bases observed.
    pub base_count: u64,
    /// Lowest Phred score observed after Phred+33 decoding.
    pub min_phred: Option<u8>,
    /// Highest Phred score observed after Phred+33 decoding.
    pub max_phred: Option<u8>,
    /// Mean Phred score across all bases.
    pub mean_phred: Option<f64>,
    /// Number of bases with Phred score at least 20.
    pub q20_bases: u64,
    /// Number of bases with Phred score at least 30.
    pub q30_bases: u64,
    /// Fraction of bases with Phred score at least 20.
    pub q20_fraction: Option<f64>,
    /// Fraction of bases with Phred score at least 30.
    pub q30_fraction: Option<f64>,
    /// Summary of per-read mean Phred scores.
    pub per_read_mean_phred: FloatStats,
}

/// FASTA-specific statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct FastaStats {
    /// Input file size from filesystem metadata.
    pub file_size_bytes: u64,
    /// Number of FASTA records.
    pub records: u64,
    /// Distribution of sequence lengths.
    pub sequence_lengths: LengthStats,
    /// Base/residue composition across all records.
    pub bases: BaseComposition,
    /// Number of records with non-empty description text.
    pub records_with_description: u64,
    /// Number of records whose sequence is empty.
    pub empty_records: u64,
    /// Number of distinct record IDs.
    pub unique_ids: u64,
    /// Number of records whose ID was seen previously.
    pub duplicate_id_records: u64,
}

/// FASTQ-specific statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct FastqStats {
    /// Input file size from filesystem metadata.
    pub file_size_bytes: u64,
    /// Number of FASTQ reads.
    pub reads: u64,
    /// Distribution of read sequence lengths.
    pub read_lengths: LengthStats,
    /// Base composition across all read sequences.
    pub bases: BaseComposition,
    /// Quality-score summary across all reads.
    pub qualities: QualityStats,
    /// Number of reads with non-empty description text.
    pub reads_with_description: u64,
    /// Number of reads containing at least one `N` base.
    pub reads_with_n: u64,
    /// Number of distinct read IDs.
    pub unique_ids: u64,
    /// Number of reads whose ID was seen previously.
    pub duplicate_id_records: u64,
}

/// SAM header statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct SamHeaderStats {
    /// Total number of SAM header records.
    pub header_records: u64,
    /// Number of `@HD` records.
    pub hd_records: u64,
    /// Number of `@SQ` records.
    pub sq_records: u64,
    /// Number of `@RG` records.
    pub rg_records: u64,
    /// Number of `@PG` records.
    pub pg_records: u64,
    /// Number of `@CO` records.
    pub co_records: u64,
    /// `@HD VN` value, when present.
    pub version: Option<String>,
    /// `@HD SO` value, when present.
    pub sort_order: Option<String>,
    /// `@HD GO` value, when present.
    pub group_order: Option<String>,
    /// Reference sequences declared by `@SQ` records.
    pub references: Vec<ReferenceSequenceStats>,
    /// Sum of declared reference lengths from `@SQ LN` fields.
    pub declared_reference_bases: u64,
}

/// BAM header/reference-dictionary statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct BamHeaderStats {
    /// Length in bytes of the parsed SAM header text.
    pub header_text_bytes: u64,
    /// Original BAM `l_text` value, including any NUL padding.
    pub binary_header_text_bytes: u64,
    /// Number of lines in the parsed SAM header text.
    pub header_lines: u64,
    /// Number of entries in the BAM reference dictionary.
    pub reference_count: u64,
    /// Reference sequence dictionary entries.
    pub references: Vec<ReferenceSequenceStats>,
    /// Sum of declared reference lengths from the BAM reference dictionary.
    pub declared_reference_bases: u64,
}

/// Reference sequence metadata.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceSequenceStats {
    /// Reference sequence name.
    pub name: String,
    /// Declared reference sequence length.
    pub length: u64,
}

/// Shared SAM/BAM alignment statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentStats {
    /// Total number of alignment records.
    pub records: u64,
    /// Number of primary alignments.
    pub primary: u64,
    /// Number of secondary alignments.
    pub secondary: u64,
    /// Number of supplementary alignments.
    pub supplementary: u64,
    /// Number of records not marked unmapped.
    pub mapped: u64,
    /// Number of records marked unmapped.
    pub unmapped: u64,
    /// Number of records marked as part of a multi-segment template.
    pub paired: u64,
    /// Number of records marked as properly aligned.
    pub properly_paired: u64,
    /// Number of records whose mate/next segment is marked unmapped.
    pub mate_unmapped: u64,
    /// Number of records on the reverse strand.
    pub reverse_complemented: u64,
    /// Number of records whose mate/next segment is on the reverse strand.
    pub mate_reverse_complemented: u64,
    /// Number of records marked as the first segment.
    pub first_segment: u64,
    /// Number of records marked as the last segment.
    pub last_segment: u64,
    /// Number of records marked duplicate.
    pub duplicate: u64,
    /// Number of records failing platform/vendor quality checks.
    pub qc_fail: u64,
    /// Number of records with mapping quality `0`.
    pub mapq_zero: u64,
    /// Number of records with mapping quality `255`.
    pub mapq_unavailable: u64,
    /// Mapping quality summary.
    pub mapq: NumberStats,
    /// Template length (`TLEN`) summary.
    pub template_length: NumberStats,
    /// Query/read length summary.
    pub query_lengths: LengthStats,
    /// Reference-consuming alignment length summary.
    pub reference_lengths: LengthStats,
    /// Total CIGAR operation lengths keyed by operation code.
    pub cigar_ops: BTreeMap<char, u64>,
    /// Count of optional/auxiliary tags observed across records.
    pub optional_tag_counts: BTreeMap<String, u64>,
    /// Count of records per reference name, with `*` for unmapped/unavailable.
    pub records_by_reference: BTreeMap<String, u64>,
}

/// SAM-specific statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct SamStats {
    /// Input file size from filesystem metadata.
    pub file_size_bytes: u64,
    /// SAM header summary.
    pub header: SamHeaderStats,
    /// Alignment record summary.
    pub alignments: AlignmentStats,
}

/// BAM-specific statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct BamStats {
    /// Input file size from filesystem metadata.
    pub file_size_bytes: u64,
    /// BAM header and reference dictionary summary.
    pub header: BamHeaderStats,
    /// Alignment record summary.
    pub alignments: AlignmentStats,
    /// Distribution of BAM record block sizes.
    pub record_block_sizes: LengthStats,
    /// Number of records whose qualities are all unavailable (`0xff`).
    pub unavailable_quality_records: u64,
}

/// POD5 per-channel statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5ChannelStats {
    /// One-indexed channel number.
    pub channel: u16,
    /// Number of reads observed on this channel.
    pub read_count: u64,
    /// Sum of `num_samples` for reads on this channel.
    pub sample_count: u64,
    /// Estimated read duration for this channel, when sample rate is known.
    pub duration_seconds: Option<f64>,
}

/// POD5 per-run statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5RunStats {
    /// Acquisition/run identifier.
    pub acquisition_id: String,
    /// User-supplied sample identifier.
    pub sample_id: String,
    /// User-supplied experiment name.
    pub experiment_name: String,
    /// Flow-cell identifier.
    pub flow_cell_id: String,
    /// Sequencing kit name.
    pub sequencing_kit: String,
    /// Samples per second for this run.
    pub sample_rate: u16,
    /// MinKNOW/software string from the run-info row.
    pub software: String,
    /// Number of reads associated with this run.
    pub read_count: u64,
    /// Sum of `num_samples` for reads associated with this run.
    pub sample_count: u64,
    /// Estimated read duration for this run, when sample rate is known.
    pub duration_seconds: Option<f64>,
}

/// POD5-specific statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct Pod5Stats {
    /// Input file size from filesystem metadata.
    pub file_size_bytes: u64,
    /// POD5 file identifier metadata, when present.
    pub file_identifier: Option<String>,
    /// File-level software metadata, when present.
    pub software: Option<String>,
    /// POD5 version metadata, when present.
    pub pod5_version: Option<String>,
    /// Number of Reads table rows.
    pub read_count: u64,
    /// Number of Signal table rows from wrapper metadata.
    pub signal_count: u64,
    /// Number of Run Info table rows.
    pub run_info_count: u64,
    /// Sum of `num_samples` across reads.
    pub total_samples: u64,
    /// Estimated total read duration, when sample rates are known.
    pub duration_seconds: Option<f64>,
    /// Distribution of per-read `num_samples` values.
    pub read_sample_lengths: LengthStats,
    /// Distribution of signal-row reference counts per read.
    pub signal_rows_per_read: LengthStats,
    /// Distribution of per-read MinKNOW event counts.
    pub minknow_events: LengthStats,
    /// Per-channel read and sample summaries.
    pub channels: Vec<Pod5ChannelStats>,
    /// Per-run read and sample summaries.
    pub runs: Vec<Pod5RunStats>,
    /// Count of reads per one-indexed well/mux.
    pub wells: BTreeMap<u8, u64>,
    /// Count of reads per pore type string.
    pub pore_types: BTreeMap<String, u64>,
    /// Count of reads per end reason string.
    pub end_reasons: BTreeMap<String, u64>,
    /// Number of reads whose end reason was forced.
    pub forced_end_reason_reads: u64,
    /// Number of reads whose run-info reference could not be resolved.
    pub reads_with_missing_run_info: u64,
    /// Summary of `median_before` values.
    pub median_before: FloatStats,
    /// Summary of tracked scaling scale values.
    pub tracked_scaling_scale: FloatStats,
    /// Summary of tracked scaling shift values.
    pub tracked_scaling_shift: FloatStats,
    /// Summary of predicted scaling scale values.
    pub predicted_scaling_scale: FloatStats,
    /// Summary of predicted scaling shift values.
    pub predicted_scaling_shift: FloatStats,
    /// Summary of calibration offset values.
    pub calibration_offset: FloatStats,
    /// Summary of calibration scale values.
    pub calibration_scale: FloatStats,
    /// Summary of seconds since the last mux change.
    pub time_since_mux_change: FloatStats,
}

/// Computes typed FASTA statistics.
pub fn fasta_stats<P: AsRef<Path>>(input: P) -> Result<FastaStats> {
    let path = input.as_ref();
    let file_size_bytes = fs::metadata(path)?.len();
    let mut reader = fasta::FastaReader::from_path(path)?;
    let mut sequence_lengths = LengthAccumulator::default();
    let mut bases = BaseComposition::default();
    let mut ids = HashSet::new();
    let mut records = 0;
    let mut records_with_description = 0;
    let mut empty_records = 0;
    let mut duplicate_id_records = 0;

    while let Some(record) = reader.read_record()? {
        records += 1;
        let len = record.sequence.len() as u64;
        sequence_lengths.push(len);
        bases.add_sequence(&record.sequence);

        if record
            .description
            .as_deref()
            .is_some_and(|value| !value.is_empty())
        {
            records_with_description += 1;
        }
        if record.sequence.is_empty() {
            empty_records += 1;
        }
        if !ids.insert(record.id) {
            duplicate_id_records += 1;
        }
    }

    Ok(FastaStats {
        file_size_bytes,
        records,
        sequence_lengths: sequence_lengths.finish(),
        bases,
        records_with_description,
        empty_records,
        unique_ids: ids.len() as u64,
        duplicate_id_records,
    })
}

/// Computes typed FASTQ statistics.
pub fn fastq_stats<P: AsRef<Path>>(input: P) -> Result<FastqStats> {
    let path = input.as_ref();
    let file_size_bytes = fs::metadata(path)?.len();
    let mut reader = fastq::FastqReader::from_path(path)?;
    let mut read_lengths = LengthAccumulator::default();
    let mut bases = BaseComposition::default();
    let mut qualities = QualityAccumulator::default();
    let mut ids = HashSet::new();
    let mut reads = 0;
    let mut reads_with_description = 0;
    let mut reads_with_n = 0;
    let mut duplicate_id_records = 0;

    while let Some(record) = reader.read_record()? {
        reads += 1;
        read_lengths.push(record.sequence.len() as u64);
        bases.add_sequence(&record.sequence);
        qualities.push(&record.quality);

        if record
            .description
            .as_deref()
            .is_some_and(|value| !value.is_empty())
        {
            reads_with_description += 1;
        }
        if record
            .sequence
            .bytes()
            .any(|byte| byte.eq_ignore_ascii_case(&b'n'))
        {
            reads_with_n += 1;
        }
        if !ids.insert(record.id) {
            duplicate_id_records += 1;
        }
    }

    Ok(FastqStats {
        file_size_bytes,
        reads,
        read_lengths: read_lengths.finish(),
        bases,
        qualities: qualities.finish(),
        reads_with_description,
        reads_with_n,
        unique_ids: ids.len() as u64,
        duplicate_id_records,
    })
}

/// Computes typed SAM statistics.
pub fn sam_stats<P: AsRef<Path>>(input: P) -> Result<SamStats> {
    let path = input.as_ref();
    let file_size_bytes = fs::metadata(path)?.len();
    let mut reader = sam::SamReader::from_path(path)?;
    let header = summarize_sam_header(&reader.header);
    let mut alignments = AlignmentAccumulator::default();

    while let Some(record) = reader.read_record()? {
        alignments.push_sam(&record)?;
    }

    Ok(SamStats {
        file_size_bytes,
        header,
        alignments: alignments.finish(),
    })
}

/// Computes typed BAM statistics.
pub fn bam_stats<P: AsRef<Path>>(input: P) -> Result<BamStats> {
    let path = input.as_ref();
    let file_size_bytes = fs::metadata(path)?.len();
    let mut reader = bam::BamReader::from_path(path)?;
    let header = summarize_bam_header(&reader.header, &reader.refs);
    let mut alignments = AlignmentAccumulator::default();
    let mut record_block_sizes = LengthAccumulator::default();
    let mut unavailable_quality_records = 0;

    while let Some(record) = reader.read_record()? {
        record_block_sizes.push(u64::from(record.fixed.block_size));
        if record
            .variable
            .qual
            .iter()
            .all(|quality| *quality == u8::MAX)
        {
            unavailable_quality_records += 1;
        }
        alignments.push_bam(&record, &reader.refs);
    }

    Ok(BamStats {
        file_size_bytes,
        header,
        alignments: alignments.finish(),
        record_block_sizes: record_block_sizes.finish(),
        unavailable_quality_records,
    })
}

/// Computes typed POD5 statistics without decompressing signal payloads.
pub fn pod5_stats<P: AsRef<Path>>(input: P) -> Result<Pod5Stats> {
    let path = input.as_ref();
    let file_size_bytes = fs::metadata(path)?.len();
    let mut reader = pod5::Pod5Reader::from_path(path)?;
    let header = reader.header.clone();
    let run_infos = reader.run_infos.clone();
    let run_info_count = header.run_info_count() as u64;
    let signal_count = header.signal_count() as u64;

    let sample_rates = run_infos
        .iter()
        .map(|run_info| (run_info.acquisition_id.clone(), run_info.sample_rate))
        .collect::<HashMap<_, _>>();
    let mut run_counts = BTreeMap::<String, (u64, u64, f64)>::new();
    let mut channel_counts = BTreeMap::<u16, (u64, u64, f64)>::new();
    let mut wells = BTreeMap::<u8, u64>::new();
    let mut pore_types = BTreeMap::<String, u64>::new();
    let mut end_reasons = BTreeMap::<String, u64>::new();
    let mut read_sample_lengths = LengthAccumulator::default();
    let mut signal_rows_per_read = LengthAccumulator::default();
    let mut minknow_events = LengthAccumulator::default();
    let mut median_before = FloatAccumulator::default();
    let mut tracked_scaling_scale = FloatAccumulator::default();
    let mut tracked_scaling_shift = FloatAccumulator::default();
    let mut predicted_scaling_scale = FloatAccumulator::default();
    let mut predicted_scaling_shift = FloatAccumulator::default();
    let mut calibration_offset = FloatAccumulator::default();
    let mut calibration_scale = FloatAccumulator::default();
    let mut time_since_mux_change = FloatAccumulator::default();
    let mut read_count = 0;
    let mut total_samples = 0;
    let mut duration_seconds = 0.0;
    let mut duration_count = 0;
    let mut forced_end_reason_reads = 0;
    let mut reads_with_missing_run_info = 0;

    while let Some(record) = reader.read_record()? {
        read_count += 1;
        total_samples += record.num_samples;
        read_sample_lengths.push(record.num_samples);
        signal_rows_per_read.push(record.signal_rows.len() as u64);
        minknow_events.push(record.num_minknow_events);
        median_before.push(f64::from(record.median_before));
        tracked_scaling_scale.push(f64::from(record.tracked_scaling_scale));
        tracked_scaling_shift.push(f64::from(record.tracked_scaling_shift));
        predicted_scaling_scale.push(f64::from(record.predicted_scaling_scale));
        predicted_scaling_shift.push(f64::from(record.predicted_scaling_shift));
        calibration_offset.push(f64::from(record.calibration_offset));
        calibration_scale.push(f64::from(record.calibration_scale));
        time_since_mux_change.push(f64::from(record.time_since_mux_change));

        *wells.entry(record.well).or_default() += 1;
        *pore_types.entry(record.pore_type.clone()).or_default() += 1;
        *end_reasons.entry(record.end_reason.clone()).or_default() += 1;
        if record.end_reason_forced {
            forced_end_reason_reads += 1;
        }

        let read_duration = sample_rates
            .get(&record.run_info)
            .copied()
            .filter(|sample_rate| *sample_rate > 0)
            .map(|sample_rate| record.num_samples as f64 / f64::from(sample_rate));

        if let Some(read_duration) = read_duration {
            duration_seconds += read_duration;
            duration_count += 1;
        } else {
            reads_with_missing_run_info += 1;
        }

        let run_entry = run_counts.entry(record.run_info.clone()).or_default();
        run_entry.0 += 1;
        run_entry.1 += record.num_samples;
        if let Some(read_duration) = read_duration {
            run_entry.2 += read_duration;
        }

        let channel_entry = channel_counts.entry(record.channel).or_default();
        channel_entry.0 += 1;
        channel_entry.1 += record.num_samples;
        if let Some(read_duration) = read_duration {
            channel_entry.2 += read_duration;
        }
    }

    let channels = channel_counts
        .into_iter()
        .map(
            |(channel, (read_count, sample_count, duration))| Pod5ChannelStats {
                channel,
                read_count,
                sample_count,
                duration_seconds: (duration > 0.0).then_some(duration),
            },
        )
        .collect();

    let runs = run_infos
        .into_iter()
        .map(|run_info| {
            let (read_count, sample_count, duration) = run_counts
                .get(&run_info.acquisition_id)
                .copied()
                .unwrap_or_default();
            Pod5RunStats {
                acquisition_id: run_info.acquisition_id,
                sample_id: run_info.sample_id,
                experiment_name: run_info.experiment_name,
                flow_cell_id: run_info.flow_cell_id,
                sequencing_kit: run_info.sequencing_kit,
                sample_rate: run_info.sample_rate,
                software: run_info.software,
                read_count,
                sample_count,
                duration_seconds: (duration > 0.0).then_some(duration),
            }
        })
        .collect();

    Ok(Pod5Stats {
        file_size_bytes,
        file_identifier: header.file_identifier,
        software: header.software,
        pod5_version: header.pod5_version,
        read_count,
        signal_count,
        run_info_count,
        total_samples,
        duration_seconds: (duration_count > 0).then_some(duration_seconds),
        read_sample_lengths: read_sample_lengths.finish(),
        signal_rows_per_read: signal_rows_per_read.finish(),
        minknow_events: minknow_events.finish(),
        channels,
        runs,
        wells,
        pore_types,
        end_reasons,
        forced_end_reason_reads,
        reads_with_missing_run_info,
        median_before: median_before.finish(),
        tracked_scaling_scale: tracked_scaling_scale.finish(),
        tracked_scaling_shift: tracked_scaling_shift.finish(),
        predicted_scaling_scale: predicted_scaling_scale.finish(),
        predicted_scaling_shift: predicted_scaling_shift.finish(),
        calibration_offset: calibration_offset.finish(),
        calibration_scale: calibration_scale.finish(),
        time_since_mux_change: time_since_mux_change.finish(),
    })
}

#[derive(Default)]
struct LengthAccumulator {
    lengths: Vec<u64>,
    total: u64,
}

impl LengthAccumulator {
    fn push(&mut self, length: u64) {
        self.total += length;
        self.lengths.push(length);
    }

    fn finish(mut self) -> LengthStats {
        self.lengths.sort_unstable_by(|left, right| right.cmp(left));
        let count = self.lengths.len() as u64;
        let min = self.lengths.last().copied();
        let max = self.lengths.first().copied();
        let mean = (count > 0).then(|| self.total as f64 / count as f64);
        let n50 = nx(&self.lengths, self.total, 50, 100);
        let n90 = nx(&self.lengths, self.total, 90, 100);

        LengthStats {
            count,
            total: self.total,
            min,
            max,
            mean,
            n50,
            n90,
        }
    }
}

#[derive(Default)]
struct NumberAccumulator {
    count: u64,
    sum: i128,
    min: Option<i64>,
    max: Option<i64>,
}

impl NumberAccumulator {
    fn push(&mut self, value: i64) {
        self.count += 1;
        self.sum += i128::from(value);
        self.min = Some(self.min.map_or(value, |current| current.min(value)));
        self.max = Some(self.max.map_or(value, |current| current.max(value)));
    }

    fn finish(self) -> NumberStats {
        NumberStats {
            count: self.count,
            min: self.min,
            max: self.max,
            mean: (self.count > 0).then(|| self.sum as f64 / self.count as f64),
        }
    }
}

#[derive(Default)]
struct FloatAccumulator {
    count: u64,
    non_finite_count: u64,
    sum: f64,
    min: Option<f64>,
    max: Option<f64>,
}

impl FloatAccumulator {
    fn push(&mut self, value: f64) {
        if !value.is_finite() {
            self.non_finite_count += 1;
            return;
        }

        self.count += 1;
        self.sum += value;
        self.min = Some(self.min.map_or(value, |current| current.min(value)));
        self.max = Some(self.max.map_or(value, |current| current.max(value)));
    }

    fn finish(self) -> FloatStats {
        FloatStats {
            count: self.count,
            non_finite_count: self.non_finite_count,
            min: self.min,
            max: self.max,
            mean: (self.count > 0).then(|| self.sum / self.count as f64),
        }
    }
}

#[derive(Default)]
struct QualityAccumulator {
    read_count: u64,
    base_count: u64,
    sum_phred: u128,
    min_phred: Option<u8>,
    max_phred: Option<u8>,
    q20_bases: u64,
    q30_bases: u64,
    per_read_mean_phred: FloatAccumulator,
}

impl QualityAccumulator {
    fn push(&mut self, quality: &str) {
        self.read_count += 1;
        let mut read_sum = 0u64;
        let mut read_bases = 0u64;

        for byte in quality.bytes() {
            let phred = byte.saturating_sub(33);
            self.base_count += 1;
            self.sum_phred += u128::from(phred);
            self.min_phred = Some(self.min_phred.map_or(phred, |current| current.min(phred)));
            self.max_phred = Some(self.max_phred.map_or(phred, |current| current.max(phred)));
            if phred >= 20 {
                self.q20_bases += 1;
            }
            if phred >= 30 {
                self.q30_bases += 1;
            }
            read_sum += u64::from(phred);
            read_bases += 1;
        }

        if read_bases > 0 {
            self.per_read_mean_phred
                .push(read_sum as f64 / read_bases as f64);
        }
    }

    fn finish(self) -> QualityStats {
        QualityStats {
            read_count: self.read_count,
            base_count: self.base_count,
            min_phred: self.min_phred,
            max_phred: self.max_phred,
            mean_phred: (self.base_count > 0)
                .then(|| self.sum_phred as f64 / self.base_count as f64),
            q20_bases: self.q20_bases,
            q30_bases: self.q30_bases,
            q20_fraction: (self.base_count > 0)
                .then(|| self.q20_bases as f64 / self.base_count as f64),
            q30_fraction: (self.base_count > 0)
                .then(|| self.q30_bases as f64 / self.base_count as f64),
            per_read_mean_phred: self.per_read_mean_phred.finish(),
        }
    }
}

#[derive(Default)]
struct AlignmentAccumulator {
    records: u64,
    primary: u64,
    secondary: u64,
    supplementary: u64,
    mapped: u64,
    unmapped: u64,
    paired: u64,
    properly_paired: u64,
    mate_unmapped: u64,
    reverse_complemented: u64,
    mate_reverse_complemented: u64,
    first_segment: u64,
    last_segment: u64,
    duplicate: u64,
    qc_fail: u64,
    mapq_zero: u64,
    mapq_unavailable: u64,
    mapq: NumberAccumulator,
    template_length: NumberAccumulator,
    query_lengths: LengthAccumulator,
    reference_lengths: LengthAccumulator,
    cigar_ops: BTreeMap<char, u64>,
    optional_tag_counts: BTreeMap<String, u64>,
    records_by_reference: BTreeMap<String, u64>,
}

impl AlignmentAccumulator {
    fn push_sam(&mut self, record: &sam::SamRecord) -> Result<()> {
        self.push_flags(
            record.flag,
            record.rname.as_str(),
            record.mapq,
            i64::from(record.tlen),
        );

        for field in &record.optional {
            *self
                .optional_tag_counts
                .entry(field.tag.clone())
                .or_default() += 1;
        }

        if record.seq != "*" {
            self.query_lengths.push(record.seq.len() as u64);
        }

        if record.cigar != "*" {
            let ops = record.cigar_ops()?;
            let mut query_len = 0u64;
            let mut reference_len = 0u64;
            for op in ops {
                let len = u64::from(op.len);
                *self.cigar_ops.entry(op.op).or_default() += len;
                if matches!(op.op, 'M' | 'I' | 'S' | '=' | 'X') {
                    query_len += len;
                }
                if matches!(op.op, 'M' | 'D' | 'N' | '=' | 'X') {
                    reference_len += len;
                }
            }
            if record.seq == "*" {
                self.query_lengths.push(query_len);
            }
            self.reference_lengths.push(reference_len);
        }

        Ok(())
    }

    fn push_bam(&mut self, record: &bam::BamRecord, refs: &[bam::BamRef]) {
        let reference_name = bam_reference_name(record, refs);
        self.push_flags(
            record.fixed.flag,
            &reference_name,
            record.fixed.mapq,
            i64::from(record.fixed.tlen),
        );
        self.query_lengths.push(u64::from(record.fixed.l_seq));

        for field in &record.auxiliary {
            *self
                .optional_tag_counts
                .entry(field.tag.clone())
                .or_default() += 1;
        }

        if !record.variable.cigar.is_empty() {
            let mut reference_len = 0u64;
            for packed in &record.variable.cigar {
                let len = u64::from(packed >> 4);
                let Some(op) = bam_cigar_op(*packed) else {
                    continue;
                };
                *self.cigar_ops.entry(op).or_default() += len;
                if matches!(op, 'M' | 'D' | 'N' | '=' | 'X') {
                    reference_len += len;
                }
            }
            self.reference_lengths.push(reference_len);
        }
    }

    fn push_flags(&mut self, flag: u16, reference_name: &str, mapq: u8, template_length: i64) {
        self.records += 1;
        self.mapq.push(i64::from(mapq));
        self.template_length.push(template_length);
        *self
            .records_by_reference
            .entry(reference_name.to_string())
            .or_default() += 1;

        if mapq == 0 {
            self.mapq_zero += 1;
        }
        if mapq == 255 {
            self.mapq_unavailable += 1;
        }
        if flag & sam::flags::SECONDARY != 0 {
            self.secondary += 1;
        }
        if flag & sam::flags::SUPPLEMENTARY != 0 {
            self.supplementary += 1;
        }
        if flag & sam::flags::SECONDARY == 0 && flag & sam::flags::SUPPLEMENTARY == 0 {
            self.primary += 1;
        }
        if flag & sam::flags::UNMAPPED != 0 {
            self.unmapped += 1;
        } else {
            self.mapped += 1;
        }
        if flag & sam::flags::MULTIPLE_SEGMENTS != 0 {
            self.paired += 1;
        }
        if flag & sam::flags::PROPERLY_ALIGNED != 0 {
            self.properly_paired += 1;
        }
        if flag & sam::flags::NEXT_UNMAPPED != 0 {
            self.mate_unmapped += 1;
        }
        if flag & sam::flags::REVERSE_COMPLEMENTED != 0 {
            self.reverse_complemented += 1;
        }
        if flag & sam::flags::NEXT_REVERSE_COMPLEMENTED != 0 {
            self.mate_reverse_complemented += 1;
        }
        if flag & sam::flags::FIRST_SEGMENT != 0 {
            self.first_segment += 1;
        }
        if flag & sam::flags::LAST_SEGMENT != 0 {
            self.last_segment += 1;
        }
        if flag & sam::flags::DUPLICATE != 0 {
            self.duplicate += 1;
        }
        if flag & sam::flags::FILTERED != 0 {
            self.qc_fail += 1;
        }
    }

    fn finish(self) -> AlignmentStats {
        AlignmentStats {
            records: self.records,
            primary: self.primary,
            secondary: self.secondary,
            supplementary: self.supplementary,
            mapped: self.mapped,
            unmapped: self.unmapped,
            paired: self.paired,
            properly_paired: self.properly_paired,
            mate_unmapped: self.mate_unmapped,
            reverse_complemented: self.reverse_complemented,
            mate_reverse_complemented: self.mate_reverse_complemented,
            first_segment: self.first_segment,
            last_segment: self.last_segment,
            duplicate: self.duplicate,
            qc_fail: self.qc_fail,
            mapq_zero: self.mapq_zero,
            mapq_unavailable: self.mapq_unavailable,
            mapq: self.mapq.finish(),
            template_length: self.template_length.finish(),
            query_lengths: self.query_lengths.finish(),
            reference_lengths: self.reference_lengths.finish(),
            cigar_ops: self.cigar_ops,
            optional_tag_counts: self.optional_tag_counts,
            records_by_reference: self.records_by_reference,
        }
    }
}

fn summarize_sam_header(header: &sam::SamHeader) -> SamHeaderStats {
    let mut hd_records = 0;
    let mut sq_records = 0;
    let mut rg_records = 0;
    let mut pg_records = 0;
    let mut co_records = 0;
    let mut references = Vec::new();

    for record in &header.records {
        match record.record_type.as_str() {
            "HD" => hd_records += 1,
            "SQ" => {
                sq_records += 1;
                if let (Some(name), Some(length)) = (record.value("SN"), record.value("LN"))
                    && let Ok(length) = length.parse::<u64>()
                {
                    references.push(ReferenceSequenceStats {
                        name: name.to_string(),
                        length,
                    });
                }
            }
            "RG" => rg_records += 1,
            "PG" => pg_records += 1,
            "CO" => co_records += 1,
            _ => {}
        }
    }

    let declared_reference_bases = references.iter().map(|reference| reference.length).sum();
    let hd = header.first("HD");

    SamHeaderStats {
        header_records: header.records.len() as u64,
        hd_records,
        sq_records,
        rg_records,
        pg_records,
        co_records,
        version: hd.and_then(|record| record.value("VN")).map(str::to_string),
        sort_order: hd.and_then(|record| record.value("SO")).map(str::to_string),
        group_order: hd.and_then(|record| record.value("GO")).map(str::to_string),
        references,
        declared_reference_bases,
    }
}

fn summarize_bam_header(header: &bam::BamHeader, refs: &[bam::BamRef]) -> BamHeaderStats {
    let references = refs
        .iter()
        .map(|reference| ReferenceSequenceStats {
            name: reference.name.clone(),
            length: u64::from(reference.l_seq),
        })
        .collect::<Vec<_>>();
    let declared_reference_bases = references.iter().map(|reference| reference.length).sum();

    BamHeaderStats {
        header_text_bytes: header.text.len() as u64,
        binary_header_text_bytes: u64::from(header.l_text),
        header_lines: header.text.lines().count() as u64,
        reference_count: refs.len() as u64,
        references,
        declared_reference_bases,
    }
}

fn bam_reference_name(record: &bam::BamRecord, refs: &[bam::BamRef]) -> String {
    usize::try_from(record.fixed.ref_id)
        .ok()
        .and_then(|index| refs.get(index))
        .map(|reference| reference.name.clone())
        .unwrap_or_else(|| "*".to_string())
}

fn bam_cigar_op(packed: u32) -> Option<char> {
    const OPS: &[u8; 9] = b"MIDNSHP=X";
    OPS.get((packed & 0x0f) as usize).copied().map(char::from)
}

fn nx(lengths_descending: &[u64], total: u64, numerator: u128, denominator: u128) -> Option<u64> {
    if lengths_descending.is_empty() {
        return None;
    }
    let threshold = (u128::from(total) * numerator).div_ceil(denominator);
    let mut cumulative = 0u128;

    for length in lengths_descending {
        cumulative += u128::from(*length);
        if cumulative >= threshold {
            return Some(*length);
        }
    }

    lengths_descending.last().copied()
}

impl fmt::Display for FastaStats {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(formatter, "FASTA statistics")?;
        writeln!(formatter, "  file_size_bytes: {}", self.file_size_bytes)?;
        writeln!(formatter, "  records: {}", self.records)?;
        writeln!(
            formatter,
            "  records_with_description: {}",
            self.records_with_description
        )?;
        writeln!(formatter, "  empty_records: {}", self.empty_records)?;
        writeln!(formatter, "  unique_ids: {}", self.unique_ids)?;
        writeln!(
            formatter,
            "  duplicate_id_records: {}",
            self.duplicate_id_records
        )?;
        write_length_stats(formatter, "sequence_lengths", &self.sequence_lengths)?;
        write_base_composition(formatter, &self.bases)
    }
}

impl fmt::Display for FastqStats {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(formatter, "FASTQ statistics")?;
        writeln!(formatter, "  file_size_bytes: {}", self.file_size_bytes)?;
        writeln!(formatter, "  reads: {}", self.reads)?;
        writeln!(
            formatter,
            "  reads_with_description: {}",
            self.reads_with_description
        )?;
        writeln!(formatter, "  reads_with_n: {}", self.reads_with_n)?;
        writeln!(formatter, "  unique_ids: {}", self.unique_ids)?;
        writeln!(
            formatter,
            "  duplicate_id_records: {}",
            self.duplicate_id_records
        )?;
        write_length_stats(formatter, "read_lengths", &self.read_lengths)?;
        write_base_composition(formatter, &self.bases)?;
        write_quality_stats(formatter, &self.qualities)
    }
}

impl fmt::Display for SamStats {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(formatter, "SAM statistics")?;
        writeln!(formatter, "  file_size_bytes: {}", self.file_size_bytes)?;
        write_sam_header_stats(formatter, &self.header)?;
        write_alignment_stats(formatter, &self.alignments)
    }
}

impl fmt::Display for BamStats {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(formatter, "BAM statistics")?;
        writeln!(formatter, "  file_size_bytes: {}", self.file_size_bytes)?;
        write_bam_header_stats(formatter, &self.header)?;
        writeln!(
            formatter,
            "  unavailable_quality_records: {}",
            self.unavailable_quality_records
        )?;
        write_length_stats(formatter, "record_block_sizes", &self.record_block_sizes)?;
        write_alignment_stats(formatter, &self.alignments)
    }
}

impl fmt::Display for Pod5Stats {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(formatter, "POD5 statistics")?;
        writeln!(formatter, "  file_size_bytes: {}", self.file_size_bytes)?;
        writeln!(
            formatter,
            "  file_identifier: {}",
            self.file_identifier.as_deref().unwrap_or("-")
        )?;
        writeln!(
            formatter,
            "  software: {}",
            self.software.as_deref().unwrap_or("-")
        )?;
        writeln!(
            formatter,
            "  pod5_version: {}",
            self.pod5_version.as_deref().unwrap_or("-")
        )?;
        writeln!(formatter, "  read_count: {}", self.read_count)?;
        writeln!(formatter, "  signal_count: {}", self.signal_count)?;
        writeln!(formatter, "  run_info_count: {}", self.run_info_count)?;
        writeln!(formatter, "  total_samples: {}", self.total_samples)?;
        writeln!(
            formatter,
            "  duration_seconds: {}",
            display_f64(self.duration_seconds)
        )?;
        writeln!(
            formatter,
            "  forced_end_reason_reads: {}",
            self.forced_end_reason_reads
        )?;
        writeln!(
            formatter,
            "  reads_with_missing_run_info: {}",
            self.reads_with_missing_run_info
        )?;
        write_length_stats(formatter, "read_sample_lengths", &self.read_sample_lengths)?;
        write_length_stats(
            formatter,
            "signal_rows_per_read",
            &self.signal_rows_per_read,
        )?;
        write_length_stats(formatter, "minknow_events", &self.minknow_events)?;
        write_float_stats(formatter, "median_before", &self.median_before)?;
        write_float_stats(
            formatter,
            "tracked_scaling_scale",
            &self.tracked_scaling_scale,
        )?;
        write_float_stats(
            formatter,
            "tracked_scaling_shift",
            &self.tracked_scaling_shift,
        )?;
        write_float_stats(
            formatter,
            "predicted_scaling_scale",
            &self.predicted_scaling_scale,
        )?;
        write_float_stats(
            formatter,
            "predicted_scaling_shift",
            &self.predicted_scaling_shift,
        )?;
        write_float_stats(formatter, "calibration_offset", &self.calibration_offset)?;
        write_float_stats(formatter, "calibration_scale", &self.calibration_scale)?;
        write_float_stats(
            formatter,
            "time_since_mux_change",
            &self.time_since_mux_change,
        )?;
        write_map(formatter, "wells", &self.wells)?;
        write_map(formatter, "pore_types", &self.pore_types)?;
        write_map(formatter, "end_reasons", &self.end_reasons)?;
        write_pod5_channels(formatter, &self.channels)?;
        write_pod5_runs(formatter, &self.runs)
    }
}

fn write_length_stats(
    formatter: &mut fmt::Formatter<'_>,
    name: &str,
    stats: &LengthStats,
) -> fmt::Result {
    writeln!(formatter, "  {name}:")?;
    writeln!(formatter, "    count: {}", stats.count)?;
    writeln!(formatter, "    total: {}", stats.total)?;
    writeln!(formatter, "    min: {}", display_u64(stats.min))?;
    writeln!(formatter, "    max: {}", display_u64(stats.max))?;
    writeln!(formatter, "    mean: {}", display_f64(stats.mean))?;
    writeln!(formatter, "    n50: {}", display_u64(stats.n50))?;
    writeln!(formatter, "    n90: {}", display_u64(stats.n90))
}

fn write_number_stats(
    formatter: &mut fmt::Formatter<'_>,
    name: &str,
    stats: &NumberStats,
) -> fmt::Result {
    writeln!(formatter, "  {name}:")?;
    writeln!(formatter, "    count: {}", stats.count)?;
    writeln!(formatter, "    min: {}", display_i64(stats.min))?;
    writeln!(formatter, "    max: {}", display_i64(stats.max))?;
    writeln!(formatter, "    mean: {}", display_f64(stats.mean))
}

fn write_float_stats(
    formatter: &mut fmt::Formatter<'_>,
    name: &str,
    stats: &FloatStats,
) -> fmt::Result {
    writeln!(formatter, "  {name}:")?;
    writeln!(formatter, "    count: {}", stats.count)?;
    writeln!(
        formatter,
        "    non_finite_count: {}",
        stats.non_finite_count
    )?;
    writeln!(formatter, "    min: {}", display_f64(stats.min))?;
    writeln!(formatter, "    max: {}", display_f64(stats.max))?;
    writeln!(formatter, "    mean: {}", display_f64(stats.mean))
}

fn write_base_composition(
    formatter: &mut fmt::Formatter<'_>,
    bases: &BaseComposition,
) -> fmt::Result {
    writeln!(formatter, "  bases:")?;
    writeln!(formatter, "    total: {}", bases.total)?;
    writeln!(formatter, "    A: {}", bases.a)?;
    writeln!(formatter, "    C: {}", bases.c)?;
    writeln!(formatter, "    G: {}", bases.g)?;
    writeln!(formatter, "    T: {}", bases.t)?;
    writeln!(formatter, "    U: {}", bases.u)?;
    writeln!(formatter, "    N: {}", bases.n)?;
    writeln!(formatter, "    ambiguous: {}", bases.ambiguous)?;
    writeln!(formatter, "    gap: {}", bases.gap)?;
    writeln!(formatter, "    other: {}", bases.other)?;
    writeln!(formatter, "    GC: {}", bases.gc_count())?;
    writeln!(
        formatter,
        "    GC_fraction: {}",
        display_f64(bases.gc_fraction())
    )
}

fn write_quality_stats(formatter: &mut fmt::Formatter<'_>, quality: &QualityStats) -> fmt::Result {
    writeln!(formatter, "  qualities:")?;
    writeln!(formatter, "    read_count: {}", quality.read_count)?;
    writeln!(formatter, "    base_count: {}", quality.base_count)?;
    writeln!(
        formatter,
        "    min_phred: {}",
        quality
            .min_phred
            .map(|value| value.to_string())
            .unwrap_or_else(|| "-".to_string())
    )?;
    writeln!(
        formatter,
        "    max_phred: {}",
        quality
            .max_phred
            .map(|value| value.to_string())
            .unwrap_or_else(|| "-".to_string())
    )?;
    writeln!(
        formatter,
        "    mean_phred: {}",
        display_f64(quality.mean_phred)
    )?;
    writeln!(formatter, "    q20_bases: {}", quality.q20_bases)?;
    writeln!(formatter, "    q30_bases: {}", quality.q30_bases)?;
    writeln!(
        formatter,
        "    q20_fraction: {}",
        display_f64(quality.q20_fraction)
    )?;
    writeln!(
        formatter,
        "    q30_fraction: {}",
        display_f64(quality.q30_fraction)
    )?;
    write_float_stats(
        formatter,
        "per_read_mean_phred",
        &quality.per_read_mean_phred,
    )
}

fn write_sam_header_stats(
    formatter: &mut fmt::Formatter<'_>,
    header: &SamHeaderStats,
) -> fmt::Result {
    writeln!(formatter, "  header:")?;
    writeln!(formatter, "    header_records: {}", header.header_records)?;
    writeln!(formatter, "    HD: {}", header.hd_records)?;
    writeln!(formatter, "    SQ: {}", header.sq_records)?;
    writeln!(formatter, "    RG: {}", header.rg_records)?;
    writeln!(formatter, "    PG: {}", header.pg_records)?;
    writeln!(formatter, "    CO: {}", header.co_records)?;
    writeln!(
        formatter,
        "    version: {}",
        header.version.as_deref().unwrap_or("-")
    )?;
    writeln!(
        formatter,
        "    sort_order: {}",
        header.sort_order.as_deref().unwrap_or("-")
    )?;
    writeln!(
        formatter,
        "    group_order: {}",
        header.group_order.as_deref().unwrap_or("-")
    )?;
    writeln!(
        formatter,
        "    declared_reference_bases: {}",
        header.declared_reference_bases
    )?;
    write_references(formatter, &header.references)
}

fn write_bam_header_stats(
    formatter: &mut fmt::Formatter<'_>,
    header: &BamHeaderStats,
) -> fmt::Result {
    writeln!(formatter, "  header:")?;
    writeln!(
        formatter,
        "    header_text_bytes: {}",
        header.header_text_bytes
    )?;
    writeln!(
        formatter,
        "    binary_header_text_bytes: {}",
        header.binary_header_text_bytes
    )?;
    writeln!(formatter, "    header_lines: {}", header.header_lines)?;
    writeln!(formatter, "    reference_count: {}", header.reference_count)?;
    writeln!(
        formatter,
        "    declared_reference_bases: {}",
        header.declared_reference_bases
    )?;
    write_references(formatter, &header.references)
}

fn write_alignment_stats(
    formatter: &mut fmt::Formatter<'_>,
    stats: &AlignmentStats,
) -> fmt::Result {
    writeln!(formatter, "  alignments:")?;
    writeln!(formatter, "    records: {}", stats.records)?;
    writeln!(formatter, "    primary: {}", stats.primary)?;
    writeln!(formatter, "    secondary: {}", stats.secondary)?;
    writeln!(formatter, "    supplementary: {}", stats.supplementary)?;
    writeln!(formatter, "    mapped: {}", stats.mapped)?;
    writeln!(formatter, "    unmapped: {}", stats.unmapped)?;
    writeln!(formatter, "    paired: {}", stats.paired)?;
    writeln!(formatter, "    properly_paired: {}", stats.properly_paired)?;
    writeln!(formatter, "    mate_unmapped: {}", stats.mate_unmapped)?;
    writeln!(
        formatter,
        "    reverse_complemented: {}",
        stats.reverse_complemented
    )?;
    writeln!(
        formatter,
        "    mate_reverse_complemented: {}",
        stats.mate_reverse_complemented
    )?;
    writeln!(formatter, "    first_segment: {}", stats.first_segment)?;
    writeln!(formatter, "    last_segment: {}", stats.last_segment)?;
    writeln!(formatter, "    duplicate: {}", stats.duplicate)?;
    writeln!(formatter, "    qc_fail: {}", stats.qc_fail)?;
    writeln!(formatter, "    mapq_zero: {}", stats.mapq_zero)?;
    writeln!(
        formatter,
        "    mapq_unavailable: {}",
        stats.mapq_unavailable
    )?;
    write_number_stats(formatter, "mapq", &stats.mapq)?;
    write_number_stats(formatter, "template_length", &stats.template_length)?;
    write_length_stats(formatter, "query_lengths", &stats.query_lengths)?;
    write_length_stats(formatter, "reference_lengths", &stats.reference_lengths)?;
    write_map(formatter, "cigar_ops", &stats.cigar_ops)?;
    write_map(formatter, "optional_tag_counts", &stats.optional_tag_counts)?;
    write_map(
        formatter,
        "records_by_reference",
        &stats.records_by_reference,
    )
}

fn write_references(
    formatter: &mut fmt::Formatter<'_>,
    references: &[ReferenceSequenceStats],
) -> fmt::Result {
    writeln!(formatter, "    references:")?;
    if references.is_empty() {
        writeln!(formatter, "      -")
    } else {
        for reference in references {
            writeln!(formatter, "      {}: {}", reference.name, reference.length)?;
        }
        Ok(())
    }
}

fn write_pod5_channels(
    formatter: &mut fmt::Formatter<'_>,
    channels: &[Pod5ChannelStats],
) -> fmt::Result {
    writeln!(formatter, "  channels:")?;
    if channels.is_empty() {
        writeln!(formatter, "    -")
    } else {
        for channel in channels {
            writeln!(
                formatter,
                "    {}: reads={} samples={} duration_seconds={}",
                channel.channel,
                channel.read_count,
                channel.sample_count,
                display_f64(channel.duration_seconds)
            )?;
        }
        Ok(())
    }
}

fn write_pod5_runs(formatter: &mut fmt::Formatter<'_>, runs: &[Pod5RunStats]) -> fmt::Result {
    writeln!(formatter, "  runs:")?;
    if runs.is_empty() {
        writeln!(formatter, "    -")
    } else {
        for run in runs {
            writeln!(
                formatter,
                "    {}: reads={} samples={} duration_seconds={} sample_rate={} sample_id={} experiment={} flow_cell={} kit={} software={}",
                run.acquisition_id,
                run.read_count,
                run.sample_count,
                display_f64(run.duration_seconds),
                run.sample_rate,
                run.sample_id,
                run.experiment_name,
                run.flow_cell_id,
                run.sequencing_kit,
                run.software
            )?;
        }
        Ok(())
    }
}

fn write_map<K: fmt::Display, V: fmt::Display>(
    formatter: &mut fmt::Formatter<'_>,
    name: &str,
    map: &BTreeMap<K, V>,
) -> fmt::Result {
    writeln!(formatter, "  {name}:")?;
    if map.is_empty() {
        writeln!(formatter, "    -")
    } else {
        for (key, value) in map {
            writeln!(formatter, "    {key}: {value}")?;
        }
        Ok(())
    }
}

fn display_u64(value: Option<u64>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "-".to_string())
}

fn display_i64(value: Option<i64>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "-".to_string())
}

fn display_f64(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.4}"))
        .unwrap_or_else(|| "-".to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    const FASTA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../fasta/ace2_fragments.fasta");
    const FASTQ: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../fastq/UDP0057_sub100.fastq");
    const SAM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../sam/aligned.sam");
    const BAM: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../bam/aligned.bam");
    const POD5: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../pod5/A_100.pod5");

    #[test]
    fn fasta_fixture_stats_capture_sequence_rollups() {
        let stats = fasta_stats(FASTA).unwrap();

        assert_eq!(stats.records, 6);
        assert_eq!(stats.sequence_lengths.total, 13_500);
        assert_eq!(stats.sequence_lengths.n50, Some(2_250));
        assert_eq!(stats.records_with_description, 6);
        assert_eq!(stats.bases.gc_count(), 6_980);
        assert_eq!(stats.duplicate_id_records, 0);
    }

    #[test]
    fn fastq_fixture_stats_capture_quality_rollups() {
        let stats = fastq_stats(FASTQ).unwrap();

        assert_eq!(stats.reads, 100);
        assert_eq!(stats.read_lengths.total, 68_700);
        assert_eq!(stats.qualities.base_count, 68_700);
        assert_eq!(stats.qualities.q20_bases, 43_271);
        assert_eq!(stats.qualities.q30_bases, 36_286);
        assert_eq!(stats.reads_with_n, 0);
    }

    #[test]
    fn sam_and_bam_fixture_alignment_stats_match() {
        let sam = sam_stats(SAM).unwrap();
        let bam = bam_stats(BAM).unwrap();

        assert_eq!(sam.alignments.records, bam.alignments.records);
        assert_eq!(sam.alignments.mapped, bam.alignments.mapped);
        assert_eq!(sam.alignments.query_lengths, bam.alignments.query_lengths);
        assert_eq!(
            sam.alignments.reference_lengths,
            bam.alignments.reference_lengths
        );
        assert_eq!(sam.alignments.cigar_ops, bam.alignments.cigar_ops);
        assert_eq!(
            sam.alignments.optional_tag_counts,
            bam.alignments.optional_tag_counts
        );
    }

    #[test]
    fn pod5_fixture_stats_capture_metadata_without_signal_decompression() {
        let stats = pod5_stats(POD5).unwrap();

        assert_eq!(stats.read_count, 100);
        assert_eq!(stats.signal_count, 100);
        assert_eq!(stats.run_info_count, 1);
        assert_eq!(stats.total_samples, 1_126_116);
        assert_eq!(stats.runs.len(), 1);
        assert_eq!(stats.pore_types.get("not_set"), Some(&100));
        assert_eq!(stats.end_reasons.get("signal_positive"), Some(&100));
    }
}
