//! Malformed BAM is reported as `InvalidData` carrying a Brust diagnostic for
//! the format at fault, as the crate documentation promises.

use brust_bam::bgzf::{EOF_BLOCK, compress_block};
use brust_bam::{BamReader, BamRecord, SamToBamConverter};
use brust_core::{Error, Format};
use sam::Sam;
use std::io;

/// The format of the Brust diagnostic inside `error`, if it has one.
fn diagnostic_format(error: &io::Error) -> Option<Format> {
    error
        .get_ref()
        .and_then(|inner| inner.downcast_ref::<Error>())
        .and_then(Error::format)
}

fn assert_bam_diagnostic(error: io::Error) {
    assert_eq!(error.kind(), io::ErrorKind::InvalidData, "{error}");
    assert_eq!(diagnostic_format(&error), Some(Format::Bam), "{error}");
}

#[test]
fn malformed_records_carry_a_bam_diagnostic() {
    // Too short for the fixed fields.
    assert_bam_diagnostic(BamRecord::new(31, vec![0; 31]).unwrap_err());
    // block_size disagrees with the data.
    assert_bam_diagnostic(BamRecord::new(40, vec![0; 36]).unwrap_err());
}

#[test]
fn corrupt_deflate_data_carries_a_bam_diagnostic() {
    let mut block = Vec::new();
    compress_block(b"BAM\x01", &mut block).unwrap();
    // Corrupt the deflate payload after the 18-byte BGZF header.
    let deflate_end = block.len() - 8;
    for byte in &mut block[18..deflate_end] {
        *byte = 0xff;
    }
    block.extend_from_slice(&EOF_BLOCK);

    let error = BamReader::from_reader(&block[..]).map(|_| ()).unwrap_err();
    assert_bam_diagnostic(error);
}

#[test]
fn bgzf_errors_keep_their_diagnostic_through_header_reads() {
    // The magic reads from a good block; l_text then hits a block with a bad
    // CRC, a BGZF error that reading the header must not flatten to text.
    let mut stream = Vec::new();
    compress_block(b"BAM\x01", &mut stream).unwrap();
    let mut second = Vec::new();
    compress_block(&[0; 8], &mut second).unwrap();
    let crc = second.len() - 8;
    second[crc] ^= 0xff;
    stream.extend_from_slice(&second);
    stream.extend_from_slice(&EOF_BLOCK);

    let error = BamReader::from_reader(&stream[..]).map(|_| ()).unwrap_err();
    assert_bam_diagnostic(error);
}

#[test]
fn sam_records_naming_undeclared_references_are_sam_errors() {
    let sam = Sam::from_reader(&b"@SQ\tSN:ref\tLN:10\n"[..]).unwrap();
    let converter = SamToBamConverter::new(&sam.header).unwrap();
    let record = Sam::from_reader(&b"r1\t0\tchrX\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"[..])
        .unwrap()
        .records
        .remove(0);

    let error = converter.convert_record(&record).unwrap_err();
    assert_eq!(diagnostic_format(&error), Some(Format::Sam), "{error}");
}
