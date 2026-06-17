//! Shared error and result types for Brust crates.
//!
//! The format crates keep their existing `std::io::Result` APIs for now, but
//! parse and validation failures can carry this domain error as the inner error
//! of an `io::Error`. The public `brust` facade re-exports these types as
//! `brust::Error`, `brust::Diagnostic`, and `brust::Format`.

use std::fmt;
use std::io;

/// Format associated with a domain diagnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// FASTA text.
    Fasta,
    /// FASTQ text.
    Fastq,
    /// SAM text.
    Sam,
    /// BAM binary alignment data.
    Bam,
    /// POD5 signal container data.
    Pod5,
}

impl Format {
    /// Returns the conventional uppercase format name.
    pub fn name(self) -> &'static str {
        match self {
            Self::Fasta => "FASTA",
            Self::Fastq => "FASTQ",
            Self::Sam => "SAM",
            Self::Bam => "BAM",
            Self::Pod5 => "POD5",
        }
    }
}

/// Structured context for malformed format data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    /// Human-readable validation failure.
    pub message: String,
    /// One-based input line number, when the parser can identify it.
    pub line: Option<usize>,
    /// Field or column name associated with the failure, when known.
    pub field: Option<String>,
}

impl Diagnostic {
    /// Creates a diagnostic with a message and no source location.
    pub fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            line: None,
            field: None,
        }
    }

    /// Adds one-based line context.
    pub fn with_line(mut self, line: usize) -> Self {
        self.line = Some(line);
        self
    }

    /// Adds field context.
    pub fn with_field(mut self, field: impl Into<String>) -> Self {
        self.field = Some(field.into());
        self
    }
}

/// Brust-wide error type.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// A lower-level I/O error represented as text so the shared error remains
    /// cloneable and easy to carry in `io::Error` values.
    Io(String),
    /// Malformed FASTA data.
    InvalidFasta(Diagnostic),
    /// Malformed FASTQ data.
    InvalidFastq(Diagnostic),
    /// Malformed SAM data.
    InvalidSam(Diagnostic),
    /// Malformed BAM data.
    InvalidBam(Diagnostic),
    /// Malformed POD5 data.
    InvalidPod5(Diagnostic),
}

impl Error {
    /// Creates an invalid-data error for a specific format.
    pub fn invalid(format: Format, message: impl Into<String>) -> Self {
        let diagnostic = Diagnostic::new(message);
        match format {
            Format::Fasta => Self::InvalidFasta(diagnostic),
            Format::Fastq => Self::InvalidFastq(diagnostic),
            Format::Sam => Self::InvalidSam(diagnostic),
            Format::Bam => Self::InvalidBam(diagnostic),
            Format::Pod5 => Self::InvalidPod5(diagnostic),
        }
    }

    /// Adds one-based line context to a format diagnostic.
    pub fn with_line(mut self, line: usize) -> Self {
        if let Some(diagnostic) = self.diagnostic_mut() {
            diagnostic.line = Some(line);
        }
        self
    }

    /// Adds field context to a format diagnostic.
    pub fn with_field(mut self, field: impl Into<String>) -> Self {
        if let Some(diagnostic) = self.diagnostic_mut() {
            diagnostic.field = Some(field.into());
        }
        self
    }

    /// Returns the format associated with this error, when it is a domain
    /// parse/validation error.
    pub fn format(&self) -> Option<Format> {
        match self {
            Self::InvalidFasta(_) => Some(Format::Fasta),
            Self::InvalidFastq(_) => Some(Format::Fastq),
            Self::InvalidSam(_) => Some(Format::Sam),
            Self::InvalidBam(_) => Some(Format::Bam),
            Self::InvalidPod5(_) => Some(Format::Pod5),
            Self::Io(_) => None,
        }
    }

    /// Returns diagnostic details for malformed data errors.
    pub fn diagnostic(&self) -> Option<&Diagnostic> {
        match self {
            Self::InvalidFasta(diagnostic)
            | Self::InvalidFastq(diagnostic)
            | Self::InvalidSam(diagnostic)
            | Self::InvalidBam(diagnostic)
            | Self::InvalidPod5(diagnostic) => Some(diagnostic),
            Self::Io(_) => None,
        }
    }

    fn diagnostic_mut(&mut self) -> Option<&mut Diagnostic> {
        match self {
            Self::InvalidFasta(diagnostic)
            | Self::InvalidFastq(diagnostic)
            | Self::InvalidSam(diagnostic)
            | Self::InvalidBam(diagnostic)
            | Self::InvalidPod5(diagnostic) => Some(diagnostic),
            Self::Io(_) => None,
        }
    }
}

/// Brust result alias.
pub type Result<T> = std::result::Result<T, Error>;

impl From<io::Error> for Error {
    fn from(error: io::Error) -> Self {
        if let Some(error) = error
            .get_ref()
            .and_then(|inner| inner.downcast_ref::<Error>())
        {
            return error.clone();
        }

        Self::Io(error.to_string())
    }
}

impl From<Error> for io::Error {
    fn from(error: Error) -> Self {
        let kind = match error {
            Error::Io(_) => io::ErrorKind::Other,
            Error::InvalidFasta(_)
            | Error::InvalidFastq(_)
            | Error::InvalidSam(_)
            | Error::InvalidBam(_)
            | Error::InvalidPod5(_) => io::ErrorKind::InvalidData,
        };
        io::Error::new(kind, error)
    }
}

impl fmt::Display for Error {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(message) => formatter.write_str(message),
            Self::InvalidFasta(diagnostic) => {
                write_diagnostic(formatter, Format::Fasta, diagnostic)
            }
            Self::InvalidFastq(diagnostic) => {
                write_diagnostic(formatter, Format::Fastq, diagnostic)
            }
            Self::InvalidSam(diagnostic) => write_diagnostic(formatter, Format::Sam, diagnostic),
            Self::InvalidBam(diagnostic) => write_diagnostic(formatter, Format::Bam, diagnostic),
            Self::InvalidPod5(diagnostic) => write_diagnostic(formatter, Format::Pod5, diagnostic),
        }
    }
}

impl std::error::Error for Error {}

fn write_diagnostic(
    formatter: &mut fmt::Formatter<'_>,
    format: Format,
    diagnostic: &Diagnostic,
) -> fmt::Result {
    write!(formatter, "invalid {}", format.name())?;
    if let Some(line) = diagnostic.line {
        write!(formatter, " at line {line}")?;
    }
    if let Some(field) = &diagnostic.field {
        write!(formatter, " in {field}")?;
    }
    write!(formatter, ": {}", diagnostic.message)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_includes_format_line_and_field() {
        let error = Error::invalid(Format::Fastq, "quality length exceeds sequence length")
            .with_line(4)
            .with_field("QUAL");

        assert_eq!(
            error.to_string(),
            "invalid FASTQ at line 4 in QUAL: quality length exceeds sequence length"
        );
    }

    #[test]
    fn converts_to_invalid_data_io_error() {
        let error: io::Error = Error::invalid(Format::Fasta, "empty ID").into();

        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert_eq!(error.to_string(), "invalid FASTA: empty ID");
    }

    #[test]
    fn preserves_domain_error_when_converting_from_io_error() {
        let original = Error::invalid(Format::Sam, "bad CIGAR").with_line(12);
        let wrapped: io::Error = original.clone().into();
        let converted = Error::from(wrapped);

        assert_eq!(converted, original);
    }
}
