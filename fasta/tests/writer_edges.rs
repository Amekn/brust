use brust_fasta::{Fasta, FastaRecord, FastaWriter};
use std::io;

fn write(sequence: &str, width: usize) -> io::Result<Vec<u8>> {
    let mut writer = FastaWriter::from_writer(Vec::new());
    writer.set_line_width(width);
    writer.write_record(&FastaRecord::new("seq1".into(), None, sequence.into()))?;
    Ok(writer.into_inner())
}

#[test]
fn written_sequences_read_back_unchanged_or_are_rejected() {
    // The reader takes a line starting with '>' as a header and trims trailing
    // whitespace, so the writer must refuse any line that would do either.
    for (sequence, width, accepted) in [
        (">ACGT", 0, false),
        (">ACGT", 60, false),
        ("AC>GT", 2, false),
        ("AC>GT", 0, true),
        ("AC>GT", 3, true),
        ("AC GT", 3, false),
        ("AC GT", 0, true),
        ("AC GT", 2, true),
        ("ACGT ", 0, false),
        ("ACGT ", 2, false),
        ("\tACGT", 0, true),
    ] {
        match write(sequence, width) {
            Ok(output) => {
                assert!(accepted, "{sequence:?} at width {width} was written");
                let fasta = Fasta::from_reader(&output[..]).unwrap();
                assert_eq!(fasta.records[0].sequence, sequence, "width {width}");
            }
            Err(error) => {
                assert!(!accepted, "{sequence:?} at width {width}: {error}");
                assert_eq!(error.kind(), io::ErrorKind::InvalidData);
            }
        }
    }
}

#[test]
fn wrapping_keeps_multibyte_characters_whole() {
    let output = write("AéTé", 2).unwrap();

    assert_eq!(
        String::from_utf8(output.clone()).unwrap(),
        ">seq1\nAé\nTé\n"
    );
    assert_eq!(
        Fasta::from_reader(&output[..]).unwrap().records[0].sequence,
        "AéTé"
    );
}
