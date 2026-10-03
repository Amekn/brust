use clap::{Parser, Subcommand, ValueEnum};
use std::io::Write;
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(name = "brust", version, about, long_about = None)]
#[command(about = "Bioinformatics format processing toolkit", long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Stats {
        format: FormatArg,
        input: PathBuf,
    },
    Validate {
        format: FormatArg,
        input: PathBuf,
    },
    Convert {
        #[command(subcommand)]
        command: ConvertCommands,
    },
}

#[derive(Subcommand)]
enum ConvertCommands {
    FastqToFasta {
        input: PathBuf,
        output: PathBuf,
    },
    FastqToSam {
        input: PathBuf,
        output: PathBuf,
    },
    FastqToBam {
        input: PathBuf,
        output: PathBuf,
        /// BGZF compression threads
        #[arg(short = 't', long, default_value_t = 1, value_parser = clap::value_parser!(u32).range(1..))]
        threads: u32,
    },
    SamToBam {
        input: PathBuf,
        output: PathBuf,
        /// BGZF compression threads
        #[arg(short = 't', long, default_value_t = 1, value_parser = clap::value_parser!(u32).range(1..))]
        threads: u32,
    },
    BamToSam {
        input: PathBuf,
        output: PathBuf,
    },
    SamToFastq {
        input: PathBuf,
        output: PathBuf,
    },
    BamToFastq {
        input: PathBuf,
        output: PathBuf,
    },
}

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum FormatArg {
    Fasta,
    Fastq,
    Sam,
    Bam,
    Pod5,
}

impl From<FormatArg> for brust::Format {
    fn from(format: FormatArg) -> Self {
        match format {
            FormatArg::Fasta => Self::Fasta,
            FormatArg::Fastq => Self::Fastq,
            FormatArg::Sam => Self::Sam,
            FormatArg::Bam => Self::Bam,
            FormatArg::Pod5 => Self::Pod5,
        }
    }
}

impl ConvertCommands {
    fn into_parts(
        self,
    ) -> (
        brust::convert::Conversion,
        PathBuf,
        PathBuf,
        brust::ConvertOptions,
    ) {
        use brust::convert::Conversion;

        let options = brust::ConvertOptions::default();
        match self {
            Self::FastqToFasta { input, output } => {
                (Conversion::FastqToFasta, input, output, options)
            }
            Self::FastqToSam { input, output } => (Conversion::FastqToSam, input, output, options),
            Self::FastqToBam {
                input,
                output,
                threads,
            } => (
                Conversion::FastqToBam,
                input,
                output,
                options.threads(threads as usize),
            ),
            Self::SamToBam {
                input,
                output,
                threads,
            } => (
                Conversion::SamToBam,
                input,
                output,
                options.threads(threads as usize),
            ),
            Self::BamToSam { input, output } => (Conversion::BamToSam, input, output, options),
            Self::SamToFastq { input, output } => (Conversion::SamToFastq, input, output, options),
            Self::BamToFastq { input, output } => (Conversion::BamToFastq, input, output, options),
        }
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}

/// Writes one line to stdout.
///
/// A closed stdout (`brust stats ... | head -1`) means the reader stopped
/// listening, so it is not an error; `println!` would panic there instead.
fn print_line(text: impl std::fmt::Display) -> std::result::Result<(), String> {
    match writeln!(std::io::stdout(), "{text}") {
        Err(error) if error.kind() != std::io::ErrorKind::BrokenPipe => {
            Err(format!("Failed to write output: {error}"))
        }
        _ => Ok(()),
    }
}

fn run(cli: Cli) -> std::result::Result<(), String> {
    match cli.command {
        Commands::Stats { format, input } => {
            let stats = brust::stats::stats(format.into(), &input)
                .map_err(|error| format!("Stats failed for {}: {}", input.display(), error))?;
            print_line(stats.display())
        }
        Commands::Validate { format, input } => {
            brust::validate::validate(format.into(), &input)
                .map_err(|error| format!("Validation failed for {}: {}", input.display(), error))?;
            print_line(format_args!("The {} file is valid.", input.display()))
        }
        Commands::Convert { command } => {
            let (conversion, input, output, options) = command.into_parts();
            brust::convert::convert_with(conversion, &input, &output, &options).map_err(
                |error| {
                    format!(
                        "Conversion failed for {} ({} -> {}): {}",
                        conversion.name(),
                        input.display(),
                        output.display(),
                        error
                    )
                },
            )?;
            print_line(format_args!(
                "Conversion completed: {} -> {}",
                input.display(),
                output.display()
            ))
        }
    }
}
