//! BAM header text is kept byte for byte; reference names must be UTF-8.

use brust_bam::{BamHeader, BamReader, BamWriter, BgzfReader, BgzfWriter};
use std::io::{self, Read, Write};

fn header(text: &[u8]) -> BamHeader {
    BamHeader {
        magic: *b"BAM\x01",
        l_text: 0,
        text: text.to_vec(),
        n_ref: 0,
    }
}

fn write(header: &BamHeader) -> Vec<u8> {
    let mut writer = BamWriter::from_writer(Vec::new());
    writer.write_header(header, &[]).unwrap();
    writer.finish().unwrap()
}

#[test]
fn header_text_that_is_not_utf8_round_trips_unchanged() {
    // A Latin-1 comment; decoding it lossily used to rewrite the byte as
    // U+FFFD (three bytes), changing the text and its length.
    let text = b"@HD\tVN:1.6\n@CO\tcaf\xe9\n";
    let written = write(&header(text));

    let reader = BamReader::from_reader(&written[..]).unwrap();
    assert_eq!(reader.header.text, text);
    assert_eq!(write(&reader.header), written);
    assert!(reader.header.text_str().is_err());
}

#[test]
fn header_text_reads_as_utf8_when_it_is() {
    let written = write(&header(b"@CO\tcaf\xc3\xa9\n"));
    let reader = BamReader::from_reader(&written[..]).unwrap();

    assert_eq!(reader.header.text_str().unwrap(), "@CO\tcafé\n");
}

#[test]
fn reference_names_that_are_not_utf8_are_rejected() {
    // A reference name is a SAM RNAME; replacing bytes would silently rename it.
    let mut writer = BamWriter::from_writer(Vec::new());
    let mut with_ref = header(b"");
    with_ref.n_ref = 1;
    let reference = brust_bam::BamRef {
        l_name: 4,
        name: "ref".to_string(),
        l_seq: 10,
    };
    writer.write_header(&with_ref, &[reference]).unwrap();
    let compressed = writer.finish().unwrap();
    let mut raw = Vec::new();
    BgzfReader::new(&compressed[..])
        .read_to_end(&mut raw)
        .unwrap();
    let name = raw
        .windows(4)
        .position(|window| window == b"ref\0")
        .unwrap();
    raw[name] = 0xff;
    let mut bgzf = BgzfWriter::new(Vec::new());
    bgzf.write_all(&raw).unwrap();
    let patched = bgzf.finish().unwrap();

    let error = BamReader::from_reader(&patched[..])
        .map(|_| ())
        .unwrap_err();
    assert_eq!(error.kind(), io::ErrorKind::InvalidData);
}
