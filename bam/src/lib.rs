//! BAM reader and writer primitives.
//!
//! The crate exposes both streaming and materialized APIs:
//!
//! - [`BamReader`] reads a gzip-compatible BAM stream, parses the binary header
//!   and reference dictionary up front, then returns one [`BamRecord`] at a
//!   time.
//! - [`Bam`] owns the header, references, and every parsed record and can be
//!   cloned or written back out.
//!
//! Record fields follow the BAM binary layout. Alignment positions are
//! zero-based, unmapped positions are represented as `-1`, CIGAR and sequence
//! data remain packed in the record, and helper methods expose SAM-style
//! strings when needed. The writer serializes the stored BAM structures into
//! BGZF blocks and validates internal lengths before emitting records.
//!
//! The [`bgzf`] module holds the BGZF layer on its own: [`BgzfReader`],
//! [`BgzfWriter`], and [`bgzf::compress_block`]. [`BamWriter`] compresses on the
//! calling thread by default. [`BamWriter::from_path_with_threads`] and
//! [`BamWriter::from_writer_with_threads`] compress on worker threads instead,
//! with the same output bytes for any thread count.

use brust_core::{AtomicFile, Error, Format};
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Write};
use std::path::Path;

pub mod bgzf;

pub use bgzf::{BgzfReader, BgzfVirtualOffset, BgzfWriter};

/// A fully materialized BAM payload.
///
/// `Bam` owns its header, references, and all records, so cloning this type
/// performs a deep copy of the parsed BAM data.
#[derive(Debug, Clone, PartialEq)]
pub struct Bam {
    /// BAM header metadata and SAM header text.
    pub header: BamHeader,
    /// Reference sequence dictionary.
    pub refs: Vec<BamRef>,
    /// All alignment records loaded from the stream.
    pub records: Vec<BamRecord>,
}

/// Streaming BAM parser over a gzip-compatible BAM reader.
///
/// BAM files are normally BGZF, which is compatible with multi-member gzip
/// decoding. `BamReader` owns the active decoder and parser cursor, so it is
/// intentionally not cloneable. Use [`BamReader::read_all`] or
/// [`Bam::from_path`] when a cloneable, in-memory representation is needed.
pub struct BamReader<R: Read = File> {
    /// BAM header metadata and SAM header text.
    pub header: BamHeader,
    /// Reference sequence dictionary.
    pub refs: Vec<BamRef>,
    reader: BgzfReader<R>,
}

/// Streaming BAM writer over any writable byte stream.
///
/// `BamWriter` writes BAM binary data in BGZF blocks and appends the standard
/// BGZF EOF block from [`BamWriter::finish`]. Use [`BamWriter::write_all`] for
/// a materialized [`Bam`], or write a header and then stream records manually.
///
/// Compression runs on the calling thread unless the writer is built with
/// [`BamWriter::from_path_with_threads`] or
/// [`BamWriter::from_writer_with_threads`]. Either way the bytes written are
/// the same.
///
/// After an I/O error from the wrapped writer, or a panic partway through a
/// call, the BGZF stream is incomplete, so every later call on the writer
/// returns an error.
pub struct BamWriter<W: Write = File> {
    writer: BgzfWriter<W>,
    // Names and lengths of the references written in the header, once it is.
    references: Option<Vec<(String, u32)>>,
}

/// Converts SAM alignment records to BAM records using one SAM header.
///
/// `SamToBamConverter` owns the BAM header, reference dictionary, and reference
/// lookup derived from a [`sam::SamHeader`]. It is intended for streaming
/// conversions: construct it once from the parsed header, write
/// [`SamToBamConverter::header`] and [`SamToBamConverter::refs`] to a
/// [`BamWriter`], then call [`SamToBamConverter::convert_record`] for each SAM
/// alignment record as it is read.
#[derive(Debug, Clone)]
pub struct SamToBamConverter {
    header: BamHeader,
    refs: Vec<BamRef>,
    reference_ids: HashMap<String, i32>,
}

/// BAM record paired with the BGZF virtual offset where the record begins.
#[derive(Debug, Clone, PartialEq)]
pub struct PositionedBamRecord {
    /// Virtual offset of the record's `block_size` field.
    pub virtual_offset: BgzfVirtualOffset,
    /// Parsed BAM alignment record.
    pub record: BamRecord,
}

/// BAM header information from the binary file header.
///
/// The `text` field holds the SAM header text bytes exactly as stored, with
/// trailing NUL padding removed; [`BamHeader::text_str`] reads it as UTF-8.
/// The original `l_text` value is preserved and is used by the writer when it
/// is non-zero.
#[derive(Debug, Clone, PartialEq)]
pub struct BamHeader {
    /// BAM magic bytes, always `[0x42, 0x41, 0x4D, 0x01]`.
    pub magic: [u8; 4],
    /// Length of the SAM header text in bytes, including any NUL padding.
    pub l_text: u32,
    /// SAM header text bytes with trailing NUL padding removed. They are kept
    /// as stored, so text that is not UTF-8 still round-trips unchanged.
    pub text: Vec<u8>,
    /// Number of reference sequences in the BAM reference dictionary.
    pub n_ref: u32,
}

/// Reference sequence dictionary entry.
///
/// Contains a reference sequence name and its declared sequence length.
#[derive(Debug, Clone, PartialEq)]
pub struct BamRef {
    /// Length of the reference name plus one NUL terminator byte.
    pub l_name: u32,
    /// Reference sequence name without the trailing NUL terminator.
    pub name: String,
    /// Reference sequence length in bases.
    pub l_seq: u32,
}

/// A BAM alignment record.
#[derive(Debug, Clone, PartialEq)]
pub struct BamRecord {
    /// Fixed-width BAM record core fields.
    pub fixed: BamRecordFixed,
    /// Variable-width read name, CIGAR, sequence, and quality fields.
    pub variable: BamRecordVariable,
    /// Auxiliary tags following the core and variable-length record fields.
    pub auxiliary: Vec<BamRecordAuxiliary>,
}

/// Fixed-length core fields of a BAM alignment record.
#[derive(Debug, Clone, PartialEq)]
pub struct BamRecordFixed {
    /// Total alignment record length, excluding the `block_size` field itself.
    pub block_size: u32,
    /// Reference sequence ID, or `-1` for an unmapped read.
    pub ref_id: i32,
    /// Zero-based leftmost coordinate, equal to SAM `POS - 1`.
    pub pos: i32,
    /// Length of the read name field, including the trailing NUL byte.
    pub l_read_name: u8,
    /// Mapping quality (`MAPQ`).
    pub mapq: u8,
    /// BAI index bin.
    pub bin: u16,
    /// Number of packed CIGAR operations.
    pub n_cigar_op: u16,
    /// Bitwise SAM flags.
    pub flag: u16,
    /// Length of the decoded read sequence in bases.
    pub l_seq: u32,
    /// Reference ID of the next segment, or `-1` when absent.
    pub next_ref_id: i32,
    /// Zero-based leftmost coordinate of the next segment, or `-1` when absent.
    pub next_pos: i32,
    /// Template length (`TLEN`).
    pub tlen: i32,
}

impl BamRecordFixed {
    /// Parses the fixed 32-byte BAM record core.
    ///
    /// `block_size` is the record length read from the stream, excluding the
    /// four-byte `block_size` field itself. `data` begins at `refID` and must
    /// contain at least the fixed core bytes.
    pub fn new(block_size: u32, data: &[u8]) -> io::Result<Self> {
        if data.len() < 32 {
            return Err(invalid_data("data too short"));
        }

        let ref_id = match data.get(0..4) {
            Some(d) => i32::from_le_bytes(d.try_into().unwrap()),
            None => return Err(invalid_data("data too short")),
        };
        let pos = match data.get(4..8) {
            Some(d) => i32::from_le_bytes(d.try_into().unwrap()),
            None => return Err(invalid_data("data too short")),
        };
        let l_read_name = data[8];
        let mapq = data[9];
        let bin = match data.get(10..12) {
            Some(d) => u16::from_le_bytes(d.try_into().unwrap()),
            None => return Err(invalid_data("data too short")),
        };
        let n_cigar_op = match data.get(12..14) {
            Some(d) => u16::from_le_bytes(d.try_into().unwrap()),
            None => return Err(invalid_data("data too short")),
        };
        let flag = match data.get(14..16) {
            Some(d) => u16::from_le_bytes(d.try_into().unwrap()),
            None => return Err(invalid_data("data too short")),
        };
        let l_seq = match data.get(16..20) {
            Some(d) => u32::from_le_bytes(d.try_into().unwrap()),
            None => return Err(invalid_data("data too short")),
        };
        let next_ref_id = match data.get(20..24) {
            Some(d) => i32::from_le_bytes(d.try_into().unwrap()),
            None => return Err(invalid_data("data too short")),
        };
        let next_pos = match data.get(24..28) {
            Some(d) => i32::from_le_bytes(d.try_into().unwrap()),
            None => return Err(invalid_data("data too short")),
        };
        let tlen = match data.get(28..32) {
            Some(d) => i32::from_le_bytes(d.try_into().unwrap()),
            None => return Err(invalid_data("data too short")),
        };
        Ok(Self {
            block_size,
            ref_id,
            pos,
            l_read_name,
            mapq,
            bin,
            n_cigar_op,
            flag,
            l_seq,
            next_ref_id,
            next_pos,
            tlen,
        })
    }
}

/// Variable-length fields of a BAM alignment record.
#[derive(Debug, Clone, PartialEq)]
pub struct BamRecordVariable {
    /// Read name without the trailing NUL terminator.
    pub read_name: String,
    /// Packed BAM CIGAR operations (`op_len << 4 | op`).
    pub cigar: Vec<u32>,
    /// Packed 4-bit BAM sequence bytes.
    pub seq: Vec<u8>,
    /// Raw Phred-scaled base qualities, or `0xff` bytes when unavailable.
    pub qual: Vec<u8>,
}

impl BamRecordVariable {
    /// Parses variable-length fields from a BAM record payload.
    ///
    /// `fixed` contains offsets and lengths from the fixed record core. `data`
    /// begins at `refID` and does not include the four-byte `block_size` field.
    /// The returned offset points to the first auxiliary field.
    pub fn new(fixed: &BamRecordFixed, data: &[u8]) -> io::Result<(Self, usize)> {
        let mut offset = 32;

        let read_name = match data.get(offset..offset + fixed.l_read_name as usize) {
            Some(name) => String::from_utf8_lossy(name)
                .trim_end_matches('\0')
                .to_string(),
            None => {
                return Err(invalid_data("read_name cannot be extracted"));
            }
        };
        offset += fixed.l_read_name as usize;

        let cigar_len_bytes = (fixed.n_cigar_op as usize)
            .checked_mul(4)
            .ok_or_else(|| invalid_data("CIGAR byte length overflow"))?;

        let cigar_bytes = match data.get(offset..offset + cigar_len_bytes) {
            Some(slice) => slice,
            None => {
                return Err(invalid_data("CIGAR data out of bounds"));
            }
        };

        let mut cigar = Vec::with_capacity(fixed.n_cigar_op as usize);
        for chunk in cigar_bytes.chunks_exact(4) {
            let mut bytes = [0u8; 4];
            bytes.copy_from_slice(chunk);
            let packed = u32::from_le_bytes(bytes);
            if packed & 0x0f > 8 {
                return Err(invalid_data("CIGAR operation code out of bounds"));
            }
            cigar.push(packed);
        }
        offset += cigar_len_bytes;

        let seq_byte_length = fixed.l_seq.div_ceil(2) as usize;
        let seq = match data.get(offset..offset + seq_byte_length) {
            Some(slice) => slice.to_vec(),
            None => {
                return Err(invalid_data("Sequence data out of bounds"));
            }
        };
        offset += seq_byte_length;

        let qual = match data.get(offset..offset + fixed.l_seq as usize) {
            Some(slice) => slice.to_vec(),
            None => {
                return Err(invalid_data("Quality data out of bounds"));
            }
        };

        offset += fixed.l_seq as usize;
        Ok((
            Self {
                read_name,
                cigar,
                seq,
                qual,
            },
            offset,
        ))
    }
}

/// Auxiliary tag/value field from a BAM alignment record.
#[derive(Debug, Clone, PartialEq)]
pub struct BamRecordAuxiliary {
    /// Two-character SAM/BAM auxiliary tag, such as `NM`, `AS`, or `RG`.
    pub tag: String,
    /// Parsed auxiliary value.
    pub value: BamAuxValue,
}

impl BamRecordAuxiliary {
    /// Parses one auxiliary field and advances `offset` past it.
    ///
    /// The parser validates value type tags and byte lengths while preserving
    /// the two-byte auxiliary tag as stored.
    pub fn new(offset: &mut usize, data: &[u8]) -> io::Result<Self> {
        let tag = match data.get(*offset..*offset + 2) {
            Some(tag) => String::from_utf8_lossy(tag).to_string(),
            None => {
                return Err(invalid_data("auxiliary tag is incomplete"));
            }
        };
        *offset += 2;

        let val_type = match data.get(*offset) {
            Some(val_type) => *val_type,
            None => {
                return Err(invalid_data("auxiliary value type is missing"));
            }
        };
        *offset += 1;

        let value = match val_type {
            b'A' => {
                let v = match data.get(*offset) {
                    Some(v) => *v,
                    None => {
                        return Err(invalid_data("auxiliary A value is missing"));
                    }
                };
                *offset += 1;
                BamAuxValue::A(v)
            }

            b'c' => {
                let v = match data.get(*offset) {
                    Some(v) => *v as i8,
                    None => {
                        return Err(invalid_data("auxiliary c value is missing"));
                    }
                };
                *offset += 1;
                BamAuxValue::c(v)
            }

            b'C' => {
                let v = match data.get(*offset) {
                    Some(v) => *v,
                    None => {
                        return Err(invalid_data("auxiliary C value is missing"));
                    }
                };
                *offset += 1;
                BamAuxValue::C(v)
            }

            b's' => {
                let v = match data.get(*offset..*offset + 2) {
                    Some(v) => i16::from_le_bytes(v.try_into().unwrap()),
                    None => {
                        return Err(invalid_data("auxiliary s value is incomplete"));
                    }
                };
                *offset += 2;
                BamAuxValue::s(v)
            }

            b'S' => {
                let v = match data.get(*offset..*offset + 2) {
                    Some(v) => u16::from_le_bytes(v.try_into().unwrap()),
                    None => {
                        return Err(invalid_data("auxiliary S value is incomplete"));
                    }
                };
                *offset += 2;
                BamAuxValue::S(v)
            }

            b'i' => {
                let v = match data.get(*offset..*offset + 4) {
                    Some(v) => i32::from_le_bytes(v.try_into().unwrap()),
                    None => {
                        return Err(invalid_data("auxiliary i value is incomplete"));
                    }
                };
                *offset += 4;
                BamAuxValue::i(v)
            }

            b'I' => {
                let v = match data.get(*offset..*offset + 4) {
                    Some(v) => u32::from_le_bytes(v.try_into().unwrap()),
                    None => {
                        return Err(invalid_data("auxiliary I value is incomplete"));
                    }
                };
                *offset += 4;
                BamAuxValue::I(v)
            }

            b'f' => {
                let v = match data.get(*offset..*offset + 4) {
                    Some(v) => f32::from_le_bytes(v.try_into().unwrap()),
                    None => {
                        return Err(invalid_data("auxiliary f value is incomplete"));
                    }
                };
                *offset += 4;
                BamAuxValue::f(v)
            }

            b'Z' => {
                let start = *offset;

                while *offset < data.len() && data[*offset] != 0 {
                    *offset += 1;
                }

                if *offset >= data.len() {
                    return Err(invalid_data("auxiliary Z string is missing NUL terminator"));
                }

                let s = String::from_utf8_lossy(&data[start..*offset]).to_string();
                // Skip the NUL terminator stored in the BAM payload.
                *offset += 1;

                BamAuxValue::Z(s)
            }

            b'H' => {
                let start = *offset;

                while *offset < data.len() && data[*offset] != 0 {
                    *offset += 1;
                }

                if *offset >= data.len() {
                    return Err(invalid_data("auxiliary H string is missing NUL terminator"));
                }

                let s = String::from_utf8_lossy(&data[start..*offset]).to_string();
                // Skip the NUL terminator stored in the BAM payload.
                *offset += 1;

                BamAuxValue::H(s)
            }

            b'B' => {
                let subtype = match data.get(*offset) {
                    Some(subtype) => *subtype,
                    None => {
                        return Err(invalid_data("auxiliary B subtype is missing"));
                    }
                };
                *offset += 1;

                let count = match data.get(*offset..*offset + 4) {
                    Some(v) => u32::from_le_bytes(v.try_into().unwrap()) as usize,
                    None => {
                        return Err(invalid_data("auxiliary B count is incomplete"));
                    }
                };
                *offset += 4;
                // The count comes from the file: check the values are present
                // before reserving room for them.
                let value_size = match subtype {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    _ => 4,
                };
                if count
                    .checked_mul(value_size)
                    .is_none_or(|bytes| bytes > data.len().saturating_sub(*offset))
                {
                    return Err(invalid_data("auxiliary B array is longer than its record"));
                }

                match subtype {
                    b'c' => {
                        let mut values = Vec::with_capacity(count);

                        for _ in 0..count {
                            let v = match data.get(*offset) {
                                Some(v) => *v as i8,
                                None => {
                                    return Err(invalid_data(
                                        "auxiliary B,c array value is missing",
                                    ));
                                }
                            };
                            *offset += 1;
                            values.push(v);
                        }

                        BamAuxValue::B(BamAuxArray::c(values))
                    }

                    b'C' => {
                        let mut values = Vec::with_capacity(count);

                        for _ in 0..count {
                            let v = match data.get(*offset) {
                                Some(v) => *v,
                                None => {
                                    return Err(invalid_data(
                                        "auxiliary B,C array value is missing",
                                    ));
                                }
                            };
                            *offset += 1;
                            values.push(v);
                        }

                        BamAuxValue::B(BamAuxArray::C(values))
                    }

                    b's' => {
                        let mut values = Vec::with_capacity(count);

                        for _ in 0..count {
                            let v = match data.get(*offset..*offset + 2) {
                                Some(v) => i16::from_le_bytes(v.try_into().unwrap()),
                                None => {
                                    return Err(invalid_data(
                                        "auxiliary B,s array value is incomplete",
                                    ));
                                }
                            };
                            *offset += 2;
                            values.push(v);
                        }

                        BamAuxValue::B(BamAuxArray::s(values))
                    }

                    b'S' => {
                        let mut values = Vec::with_capacity(count);

                        for _ in 0..count {
                            let v = match data.get(*offset..*offset + 2) {
                                Some(v) => u16::from_le_bytes(v.try_into().unwrap()),
                                None => {
                                    return Err(invalid_data(
                                        "auxiliary B,S array value is incomplete",
                                    ));
                                }
                            };
                            *offset += 2;
                            values.push(v);
                        }

                        BamAuxValue::B(BamAuxArray::S(values))
                    }

                    b'i' => {
                        let mut values = Vec::with_capacity(count);

                        for _ in 0..count {
                            let v = match data.get(*offset..*offset + 4) {
                                Some(v) => i32::from_le_bytes(v.try_into().unwrap()),
                                None => {
                                    return Err(invalid_data(
                                        "auxiliary B,i array value is incomplete",
                                    ));
                                }
                            };
                            *offset += 4;
                            values.push(v);
                        }

                        BamAuxValue::B(BamAuxArray::i(values))
                    }

                    b'I' => {
                        let mut values = Vec::with_capacity(count);

                        for _ in 0..count {
                            let v = match data.get(*offset..*offset + 4) {
                                Some(v) => u32::from_le_bytes(v.try_into().unwrap()),
                                None => {
                                    return Err(invalid_data(
                                        "auxiliary B,I array value is incomplete",
                                    ));
                                }
                            };
                            *offset += 4;
                            values.push(v);
                        }

                        BamAuxValue::B(BamAuxArray::I(values))
                    }

                    b'f' => {
                        let mut values = Vec::with_capacity(count);

                        for _ in 0..count {
                            let v = match data.get(*offset..*offset + 4) {
                                Some(v) => f32::from_le_bytes(v.try_into().unwrap()),
                                None => {
                                    return Err(invalid_data(
                                        "auxiliary B,f array value is incomplete",
                                    ));
                                }
                            };
                            *offset += 4;
                            values.push(v);
                        }

                        BamAuxValue::B(BamAuxArray::f(values))
                    }

                    _ => {
                        return Err(invalid_data(format!(
                            "unknown auxiliary B array subtype: {}",
                            subtype as char
                        )));
                    }
                }
            }

            _ => {
                return Err(invalid_data(format!(
                    "unknown auxiliary value type: {}",
                    val_type as char
                )));
            }
        };

        Ok(Self { tag, value })
    }
}

/// BAM auxiliary tag value.
#[derive(Debug, Clone, PartialEq)]
pub enum BamAuxValue {
    /// Raw byte for an `A` auxiliary value.
    A(u8),
    /// Signed 8-bit integer (`c`).
    #[allow(non_camel_case_types)]
    c(i8),
    /// Unsigned 8-bit integer (`C`).
    C(u8),
    /// Signed 16-bit integer (`s`).
    #[allow(non_camel_case_types)]
    s(i16),
    /// Unsigned 16-bit integer (`S`).
    S(u16),
    /// Signed 32-bit integer (`i`).
    #[allow(non_camel_case_types)]
    i(i32),
    /// Unsigned 32-bit integer (`I`).
    I(u32),
    /// 32-bit floating point number (`f`).
    #[allow(non_camel_case_types)]
    f(f32),
    /// String value from a NUL-terminated `Z` payload, stored without the NUL.
    Z(String),
    /// Hex string from a NUL-terminated `H` payload, stored without the NUL.
    H(String),
    /// Homogeneous numeric array (`B`).
    B(BamAuxArray),
}

/// Homogeneous BAM auxiliary numeric array.
#[derive(Debug, Clone, PartialEq)]
pub enum BamAuxArray {
    /// Signed 8-bit integer array (`B,c`).
    #[allow(non_camel_case_types)]
    c(Vec<i8>),
    /// Unsigned 8-bit integer array (`B,C`).
    C(Vec<u8>),
    /// Signed 16-bit integer array (`B,s`).
    #[allow(non_camel_case_types)]
    s(Vec<i16>),
    /// Unsigned 16-bit integer array (`B,S`).
    S(Vec<u16>),
    /// Signed 32-bit integer array (`B,i`).
    #[allow(non_camel_case_types)]
    i(Vec<i32>),
    /// Unsigned 32-bit integer array (`B,I`).
    I(Vec<u32>),
    /// 32-bit floating point array (`B,f`).
    #[allow(non_camel_case_types)]
    f(Vec<f32>),
}

impl BamRecord {
    /// Parses a BAM alignment record from raw payload bytes.
    ///
    /// `block_size` is not included in `data`. `data` begins at `refID` and
    /// contains the full alignment record payload. The parser validates the
    /// fixed/variable field boundaries and parses auxiliary fields until the end
    /// of the payload.
    pub fn new(block_size: u32, data: Vec<u8>) -> io::Result<Self> {
        if block_size as usize != data.len() {
            return Err(invalid_data("block_size does not match record data length"));
        }

        if data.len() < 32 {
            return Err(invalid_data(
                "BAM record is smaller than the 32-byte fixed core (exclude block_size)",
            ));
        }

        let fixed = BamRecordFixed::new(block_size, &data)?;
        let (variable, mut var_offset) = BamRecordVariable::new(&fixed, &data)?;
        let mut auxiliary = Vec::new();
        while var_offset < data.len() {
            let aux = BamRecordAuxiliary::new(&mut var_offset, &data)?;
            auxiliary.push(aux);
        }

        let record = Self {
            fixed,
            variable,
            auxiliary,
        };
        record.check_query_length()?;
        Ok(record)
    }

    /// Appends this record in BAM binary form to `out`.
    ///
    /// The output is the `block_size` field (`u32`, little-endian) followed by
    /// the record payload, which is what [`BamWriter::write_record`] feeds into
    /// the BGZF stream, with the same record-layout checks; only the writer can
    /// also check reference IDs against its dictionary. `out` is left
    /// unchanged when it returns an error. A non-zero `fixed.block_size` must
    /// equal the encoded payload length.
    pub fn encode(&self, out: &mut Vec<u8>) -> io::Result<()> {
        let mut payload = Vec::new();
        encode_bam_record_payload(self, &mut payload)?;

        let block_size = payload.len() as u32;
        if self.fixed.block_size != 0 && self.fixed.block_size != block_size {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "BAM record block_size does not match encoded payload length",
            ));
        }

        out.extend_from_slice(&block_size.to_le_bytes());
        out.extend_from_slice(&payload);
        Ok(())
    }

    /// As htslib does: a mapped record with SEQ and a CIGAR must have them
    /// agree, using the real CIGAR (see [`BamRecord::cigar_ops`]).
    fn check_query_length(&self) -> io::Result<()> {
        if self.fixed.l_seq == 0 || self.fixed.flag & 0x4 != 0 {
            return Ok(());
        }
        let cigar = self.cigar_ops();
        if cigar.is_empty() {
            return Ok(());
        }
        let query_len = cigar
            .iter()
            .filter(|&&op| matches!(op & 0x0f, 0 | 1 | 4 | 7 | 8))
            .map(|&op| u64::from(op >> 4))
            .sum::<u64>();
        if query_len != u64::from(self.fixed.l_seq) {
            return Err(invalid_data("BAM CIGAR and query sequence lengths differ"));
        }
        Ok(())
    }

    /// Returns this record's read name.
    pub fn read_name(&self) -> &str {
        &self.variable.read_name
    }

    /// Returns the record's real CIGAR as packed BAM operations.
    ///
    /// An alignment with more than 65,535 operations is stored as a `kSmN`
    /// placeholder (`k` = `l_seq`) with the real CIGAR in a `CG:B:I` tag (SAM
    /// spec 4.2.2). For such a record this returns the tag's operations, using
    /// htslib's test: the record is placed, its first operation soft-clips the
    /// whole read, and `CG` is a `B,I` or `B,i` array at least as long as the
    /// stored CIGAR. Otherwise it returns `variable.cigar`.
    pub fn cigar_ops(&self) -> Cow<'_, [u32]> {
        match self.long_cigar_tag() {
            Some(index) => match &self.auxiliary[index].value {
                BamAuxValue::B(BamAuxArray::I(ops)) => Cow::Borrowed(ops),
                BamAuxValue::B(BamAuxArray::i(ops)) => {
                    Cow::Owned(ops.iter().map(|&op| op as u32).collect())
                }
                _ => unreachable!("long_cigar_tag checks the array type"),
            },
            None => Cow::Borrowed(&self.variable.cigar),
        }
    }

    /// Whether this record's real CIGAR is the one in its `CG` tag (see
    /// [`BamRecord::cigar_ops`]); that tag is then not part of the record's
    /// SAM form.
    pub fn stores_cigar_in_cg_tag(&self) -> bool {
        self.long_cigar_tag().is_some()
    }

    /// Index of the `CG` tag holding the real CIGAR, if any.
    fn long_cigar_tag(&self) -> Option<usize> {
        let first = *self.variable.cigar.first()?;
        if self.fixed.ref_id < 0
            || self.fixed.pos < 0
            || first & 0x0f != 4
            || first >> 4 != self.fixed.l_seq
        {
            return None;
        }
        let index = self.auxiliary.iter().position(|aux| aux.tag == "CG")?;
        let ops = match &self.auxiliary[index].value {
            BamAuxValue::B(BamAuxArray::I(ops)) => ops.len(),
            BamAuxValue::B(BamAuxArray::i(ops)) => ops.len(),
            _ => return None,
        };
        (ops >= self.variable.cigar.len()).then_some(index)
    }

    /// Returns the record's real CIGAR (see [`BamRecord::cigar_ops`]) as a
    /// SAM-style string.
    pub fn cigar_string(&self) -> String {
        const CIGAR_OPS: &[u8; 9] = b"MIDNSHP=X";
        let cigar = self.cigar_ops();
        if cigar.is_empty() {
            return "*".to_string();
        }

        let mut decoded = String::new();
        for op in cigar.iter() {
            let op_code = (op & 0x0f) as usize;
            let op_char = CIGAR_OPS.get(op_code).copied().unwrap_or(b'?') as char;
            decoded.push_str(&(op >> 4).to_string());
            decoded.push(op_char);
        }
        decoded
    }

    /// Decodes the packed BAM sequence bytes into a base string.
    pub fn sequence_string(&self) -> String {
        const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
        // l_seq and seq are public fields; never index past the packed data.
        let l_seq = (self.fixed.l_seq as usize).min(self.variable.seq.len() * 2);
        let mut decoded = String::with_capacity(l_seq);
        for i in 0..l_seq {
            let byte = self.variable.seq[i / 2];
            let code = if i % 2 == 0 { byte >> 4 } else { byte & 0x0f };
            decoded.push(BASES[code as usize] as char);
        }
        decoded
    }

    /// Decodes BAM quality bytes into a SAM-style quality string.
    ///
    /// If every stored quality byte is `0xff`, the unavailable-quality marker
    /// `*` is returned.
    pub fn quality_string(&self) -> String {
        if self.variable.qual.iter().all(|quality| *quality == 0xff) {
            return "*".to_string();
        }

        self.variable
            .qual
            .iter()
            .map(|quality| char::from(quality.saturating_add(33)))
            .collect()
    }

    /// Returns the sequence in its original sequencing orientation.
    ///
    /// A record on the reverse strand (flag `0x10`) stores its sequence reverse-complemented,
    /// so this reverse-complements [`sequence_string`](Self::sequence_string) back. Case and
    /// the gap symbols `-` and `.` are kept, `U` becomes `A`, and any other byte that is not an
    /// IUPAC nucleotide code (such as `=`) becomes `N`. See [`seq::complement`] for the exact
    /// rules. Any other record's sequence is returned unchanged, and a record with no sequence
    /// gives an empty string.
    pub fn original_sequence_string(&self) -> String {
        let sequence = self.sequence_string();
        if !self.is_reverse_complemented() {
            return sequence;
        }
        String::from_utf8(seq::reverse_complement(sequence.as_bytes()))
            .expect("reverse_complement output is ASCII")
    }

    /// Returns the quality string in its original sequencing orientation.
    ///
    /// A record on the reverse strand (flag `0x10`) stores its qualities reversed, so this
    /// reverses [`quality_string`](Self::quality_string) back. Any other record's qualities are
    /// returned unchanged, and the unavailable-quality marker `*` stays `*`.
    pub fn original_quality_string(&self) -> String {
        let quality = self.quality_string();
        if !self.is_reverse_complemented() || quality == "*" {
            return quality;
        }
        quality.chars().rev().collect()
    }

    /// Returns the first auxiliary tag value matching `tag`.
    pub fn aux(&self, tag: &str) -> Option<&BamAuxValue> {
        self.auxiliary
            .iter()
            .find(|auxiliary| auxiliary.tag == tag)
            .map(|auxiliary| &auxiliary.value)
    }

    /// Returns `true` when the multi-segment template flag (`0x1`) is set.
    pub fn has_multiple_segments(&self) -> bool {
        self.fixed.flag & sam::flags::MULTIPLE_SEGMENTS != 0
    }

    /// Returns `true` when the properly-aligned flag (`0x2`) is set.
    pub fn is_properly_aligned(&self) -> bool {
        self.fixed.flag & sam::flags::PROPERLY_ALIGNED != 0
    }

    /// Returns `true` when the segment-unmapped flag (`0x4`) is set.
    pub fn is_unmapped(&self) -> bool {
        self.fixed.flag & sam::flags::UNMAPPED != 0
    }

    /// Returns `true` when the mate/next-segment-unmapped flag (`0x8`) is set.
    pub fn is_next_unmapped(&self) -> bool {
        self.fixed.flag & sam::flags::NEXT_UNMAPPED != 0
    }

    /// Returns `true` when the reverse-complemented flag (`0x10`) is set.
    pub fn is_reverse_complemented(&self) -> bool {
        self.fixed.flag & sam::flags::REVERSE_COMPLEMENTED != 0
    }

    /// Returns `true` when the mate/next-segment reverse-complemented flag (`0x20`) is set.
    pub fn is_next_reverse_complemented(&self) -> bool {
        self.fixed.flag & sam::flags::NEXT_REVERSE_COMPLEMENTED != 0
    }

    /// Returns `true` when this is the first segment in a template (`0x40`).
    pub fn is_first_segment(&self) -> bool {
        self.fixed.flag & sam::flags::FIRST_SEGMENT != 0
    }

    /// Returns `true` when this is the last segment in a template (`0x80`).
    pub fn is_last_segment(&self) -> bool {
        self.fixed.flag & sam::flags::LAST_SEGMENT != 0
    }

    /// Returns `true` when this is a secondary alignment (`0x100`).
    pub fn is_secondary(&self) -> bool {
        self.fixed.flag & sam::flags::SECONDARY != 0
    }

    /// Returns `true` when this segment fails platform/vendor checks (`0x200`).
    pub fn is_filtered(&self) -> bool {
        self.fixed.flag & sam::flags::FILTERED != 0
    }

    /// Returns `true` when this segment is marked duplicate (`0x400`).
    pub fn is_duplicate(&self) -> bool {
        self.fixed.flag & sam::flags::DUPLICATE != 0
    }

    /// Returns `true` when this is a supplementary alignment (`0x800`).
    pub fn is_supplementary(&self) -> bool {
        self.fixed.flag & sam::flags::SUPPLEMENTARY != 0
    }

    /// Converts this BAM record to a SAM alignment record using the supplied
    /// reference dictionary.
    pub fn to_sam_record(&self, refs: &[BamRef]) -> io::Result<sam::SamRecord> {
        bam_record_to_sam_record(self, refs)
    }

    /// Formats this BAM record as a SAM alignment line using the supplied
    /// reference dictionary.
    pub fn to_sam_line(&self, refs: &[BamRef]) -> io::Result<String> {
        format_sam_line(&self.to_sam_record(refs)?)
    }
}

impl SamToBamConverter {
    /// Builds a converter and BAM header/reference dictionary from a SAM header.
    ///
    /// The SAM header is validated by formatting it into the BAM header text and
    /// extracting all `@SQ` records into the BAM reference dictionary. Records
    /// converted with this value must reference only these declared references,
    /// unless they are unmapped and use `*`.
    pub fn new(header: &sam::SamHeader) -> io::Result<Self> {
        let mut validation_sink = Vec::new();
        sam::Sam {
            header: header.clone(),
            records: Vec::new(),
        }
        .to_writer(&mut validation_sink)?;

        let refs = bam_refs_from_sam_header(header)?;
        let reference_ids = reference_ids_from_refs(&refs);
        let text = sam_header_text(header)?.into_bytes();
        let header = BamHeader {
            magic: [0x42, 0x41, 0x4D, 0x01],
            l_text: text.len() as u32,
            text,
            n_ref: refs.len() as u32,
        };

        Ok(Self {
            header,
            refs,
            reference_ids,
        })
    }

    /// Returns the BAM header derived from the SAM header.
    pub fn header(&self) -> &BamHeader {
        &self.header
    }

    /// Returns the BAM reference dictionary derived from `@SQ` header records.
    pub fn refs(&self) -> &[BamRef] {
        &self.refs
    }

    /// Converts one SAM alignment record to a BAM record.
    ///
    /// This validates the SAM record, resolves reference names against the
    /// converter's reference dictionary, packs CIGAR/sequence/quality fields,
    /// converts optional tags, and computes the BAM record `block_size`.
    pub fn convert_record(&self, record: &sam::SamRecord) -> io::Result<BamRecord> {
        bam_record_from_sam_record(record, &self.reference_ids)
    }
}

impl Bam {
    /// Opens a BAM file and materializes all records into memory.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        BamReader::from_path(path)?.read_all()
    }

    /// Materializes all records from a compressed BAM byte stream.
    pub fn from_reader<R: Read>(reader: R) -> io::Result<Self> {
        BamReader::from_reader(reader)?.read_all()
    }

    /// Writes this BAM payload to a filesystem path.
    pub fn to_path<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let writer = BamWriter::from_path(path)?;
        self.write_with(writer)
    }

    /// Writes this BAM payload to a filesystem path atomically.
    ///
    /// Nothing appears at `path` unless the whole payload is written and
    /// committed. See [`BamWriter::from_path_atomic`].
    pub fn to_path_atomic<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut writer = BamWriter::from_path_atomic(path)?;
        writer.write_all(self)?;
        writer.commit()
    }

    /// Writes this BAM payload to a writable byte stream.
    pub fn to_writer<W: Write>(&self, writer: W) -> io::Result<()> {
        let writer = BamWriter::from_writer(writer);
        self.write_with(writer)
    }

    /// Converts a supported materialized SAM payload into BAM.
    ///
    /// The conversion supports SAM records whose reference names are declared in
    /// `@SQ`, textual sequence and quality fields, standard CIGAR operations,
    /// and scalar/string/hex/numeric-array optional fields.
    pub fn from_sam(sam: &sam::Sam) -> io::Result<Self> {
        let converter = SamToBamConverter::new(&sam.header)?;
        let records = sam
            .records
            .iter()
            .map(|record| converter.convert_record(record))
            .collect::<io::Result<Vec<_>>>()?;

        Ok(Self {
            header: converter.header.clone(),
            refs: converter.refs.clone(),
            records,
        })
    }

    /// Converts all records in this BAM payload to SAM records.
    pub fn to_sam_records(&self) -> io::Result<Vec<sam::SamRecord>> {
        self.records
            .iter()
            .map(|record| record.to_sam_record(&self.refs))
            .collect()
    }

    fn write_with<W: Write>(&self, mut writer: BamWriter<W>) -> io::Result<()> {
        writer.write_all(self)?;
        writer.finish().map(|_| ())
    }
}

impl BamReader<File> {
    /// Opens a BAM file from a filesystem path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::open(path)?;
        Self::from_reader(file)
    }

    /// Opens a BAM file from a filesystem path.
    ///
    /// This is a convenience alias for [`BamReader::from_path`].
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_path(path)
    }
}

impl<R: Read> BamReader<R> {
    /// Creates a streaming parser from a compressed BAM byte stream.
    ///
    /// The BAM header and reference dictionary are parsed immediately and
    /// exposed through [`BamReader::header`] and [`BamReader::refs`].
    pub fn from_reader(reader: R) -> io::Result<Self> {
        let mut reader = BgzfReader::new(reader);
        let header = BamHeader::read_from(&mut reader)?;
        // n_ref comes from the file, so the list grows as entries are read.
        let mut refs = Vec::new();

        for _ in 0..header.n_ref {
            refs.push(BamRef::read_from(&mut reader)?);
        }

        Ok(Self {
            reader,
            header,
            refs,
        })
    }

    /// Sets whether the BAM stream must end with the standard BGZF EOF block.
    ///
    /// The default is `false`: a stream that ends without the marker is read
    /// as complete, as before. When `true`, reaching the end of the stream is
    /// an [`io::ErrorKind::InvalidData`] error unless the last BGZF block is
    /// byte-for-byte the 28-byte EOF block. A final empty block that differs
    /// in any byte fails too.
    ///
    /// The check happens at the end of the stream and is reported once, by
    /// [`BamReader::read_record`], [`BamReader::records`] or
    /// [`BamReader::read_all`]. After it the reader behaves as at end of
    /// stream. Set this before the reader reaches the end, because a lenient
    /// end is cached and is not re-checked. Empty blocks in the middle of the
    /// stream are still skipped.
    ///
    /// A block that fails to read, such as one cut short, is reported by its
    /// own error, not also as a missing EOF block.
    ///
    /// A joined stream cut exactly after an interior EOF marker ends with a
    /// valid marker, so this check cannot detect that truncation.
    pub fn set_require_eof_block(&mut self, require: bool) {
        self.reader.set_require_eof_block(require);
    }

    /// Returns the current BGZF virtual offset.
    pub fn virtual_offset(&self) -> BgzfVirtualOffset {
        self.reader.virtual_offset()
    }

    /// Reads the next BAM alignment record.
    ///
    /// Returns `Ok(None)` only when the stream is at EOF before a new record's
    /// `block_size` field begins.
    pub fn read_record(&mut self) -> io::Result<Option<BamRecord>> {
        Ok(self
            .read_record_with_virtual_offset()?
            .map(|positioned| positioned.record))
    }

    /// Reads the next BAM alignment record with its starting BGZF virtual offset.
    pub fn read_record_with_virtual_offset(&mut self) -> io::Result<Option<PositionedBamRecord>> {
        let virtual_offset = self.virtual_offset();
        let mut block_size = [0u8; 4];
        if !read_exact_or_eof(&mut self.reader, &mut block_size, "BAM record block_size")? {
            return Ok(None);
        }

        let block_size = u32::from_le_bytes(block_size);

        let record_data = read_exact_vec(&mut self.reader, block_size, "BAM record")?;
        let record = BamRecord::new(block_size, record_data)?;
        // As htslib does: reference IDs must be -1 or index the dictionary.
        let references = self.refs.len();
        let in_range = |id: i32| id == -1 || usize::try_from(id).is_ok_and(|id| id < references);
        if !in_range(record.fixed.ref_id) || !in_range(record.fixed.next_ref_id) {
            return Err(invalid_data(
                "BAM record reference ID is outside the reference dictionary",
            ));
        }
        Ok(Some(PositionedBamRecord {
            virtual_offset,
            record,
        }))
    }

    /// Reads the next BAM alignment record.
    ///
    /// This is a compatibility alias for [`BamReader::read_record`].
    pub fn read(&mut self) -> io::Result<Option<BamRecord>> {
        self.read_record()
    }

    /// Returns an iterator over records in the stream.
    pub fn records(&mut self) -> BamRecords<'_, R> {
        BamRecords { reader: self }
    }

    /// Consumes this reader and materializes the entire BAM stream.
    pub fn read_all(mut self) -> io::Result<Bam> {
        let mut records = Vec::new();
        while let Some(record) = self.read_record()? {
            records.push(record);
        }

        Ok(Bam {
            header: self.header,
            refs: self.refs,
            records,
        })
    }
}

impl BamWriter<File> {
    /// Creates or truncates a BAM file at a filesystem path.
    ///
    /// Use [`BamWriter::from_path_atomic`] to keep an existing file intact until
    /// the new one is complete.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self::from_writer(file))
    }

    /// Creates or truncates a BAM file at a filesystem path.
    ///
    /// This is a convenience alias for [`BamWriter::from_path`].
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Self::from_path(path)
    }

    /// Creates or truncates a BAM file at a filesystem path, compressing BGZF
    /// blocks on `threads` worker threads.
    ///
    /// A `threads` of 0 or 1 compresses on the calling thread, like
    /// [`BamWriter::from_path`]. Larger values start that many workers, capped
    /// at [`bgzf::MAX_THREADS`]. Given the same calls, the file is
    /// byte-for-byte the same as [`BamWriter::from_path`] writes, whatever the
    /// thread count.
    pub fn from_path_with_threads<P: AsRef<Path>>(path: P, threads: usize) -> io::Result<Self> {
        let file = File::create(path)?;
        Ok(Self::from_writer_with_threads(file, threads))
    }
}

impl BamWriter<AtomicFile> {
    /// Creates a BAM writer that publishes its output atomically.
    ///
    /// Nothing appears at `path` until [`commit`](Self::commit) renames the
    /// finished file into place; dropping the writer discards the output. See
    /// [`AtomicFile`] for the one error `commit` can return after the rename.
    pub fn from_path_atomic<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self::from_writer(AtomicFile::create(path)?))
    }

    /// Like [`BamWriter::from_path_atomic`], compressing BGZF blocks on
    /// `threads` worker threads.
    ///
    /// See [`BamWriter::from_path_with_threads`] for how `threads` is treated.
    /// Dropping the writer without committing stops the workers and discards
    /// the output.
    pub fn from_path_atomic_with_threads<P: AsRef<Path>>(
        path: P,
        threads: usize,
    ) -> io::Result<Self> {
        Ok(Self::from_writer_with_threads(
            AtomicFile::create(path)?,
            threads,
        ))
    }

    /// Finishes the stream, then flushes, syncs and renames the output into place.
    ///
    /// This writes the last block and the BGZF EOF block, and joins any worker
    /// threads. See [`AtomicFile::commit`] for the exact steps and errors.
    ///
    /// Returns an [`io::ErrorKind::InvalidInput`] error and publishes nothing
    /// if no header was written, since an EOF block alone is not a BAM.
    pub fn commit(self) -> io::Result<()> {
        if self.references.is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "BAM writer has no header to commit",
            ));
        }
        self.finish()?.commit()
    }
}

impl<W: Write> BamWriter<W> {
    /// Creates a BAM writer from a writable byte stream.
    pub fn from_writer(writer: W) -> Self {
        Self {
            writer: BgzfWriter::new(writer),
            references: None,
        }
    }

    /// Creates a BAM writer from a writable byte stream, compressing BGZF
    /// blocks on `threads` worker threads.
    ///
    /// A `threads` of 0 or 1 compresses on the calling thread, like
    /// [`BamWriter::from_writer`]. Larger values start that many workers,
    /// capped at [`bgzf::MAX_THREADS`]. Given the same calls, the bytes
    /// written are identical to [`BamWriter::from_writer`]'s, whatever the
    /// thread count. Workers only compress; all writes to `writer` happen on
    /// the calling thread, so `writer` need not be [`Send`].
    pub fn from_writer_with_threads(writer: W, threads: usize) -> Self {
        Self {
            writer: BgzfWriter::with_threads(writer, threads),
            references: None,
        }
    }

    /// Writes the BAM header and reference dictionary.
    ///
    /// The header can be written only once. The stored reference count must
    /// match the supplied reference slice.
    pub fn write_header(&mut self, header: &BamHeader, refs: &[BamRef]) -> io::Result<()> {
        if self.references.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "BAM header has already been written",
            ));
        }
        let mut data = Vec::new();
        write_bam_header(&mut data, header, refs)?;
        self.writer.write_all(&data)?;
        self.references = Some(
            refs.iter()
                .map(|reference| (reference.name.clone(), reference.l_seq))
                .collect(),
        );
        Ok(())
    }

    /// Writes one BAM alignment record.
    ///
    /// A header must have already been written, and the record's reference IDs
    /// must be -1 or index its reference dictionary, as readers require. Packed
    /// CIGAR, sequence, quality, and auxiliary fields are serialized from the
    /// record as stored.
    pub fn write_record(&mut self, record: &BamRecord) -> io::Result<()> {
        let Some(references) = &self.references else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "BAM header must be written before records",
            ));
        };
        let in_range =
            |id: i32| id == -1 || usize::try_from(id).is_ok_and(|id| id < references.len());
        if !in_range(record.fixed.ref_id) || !in_range(record.fixed.next_ref_id) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "BAM record reference ID is outside the written reference dictionary",
            ));
        }
        let mut data = Vec::new();
        record.encode(&mut data)?;
        self.writer.write_all(&data)
    }

    /// Writes one BAM alignment record.
    ///
    /// This is a compatibility alias for [`BamWriter::write_record`].
    pub fn write(&mut self, record: &BamRecord) -> io::Result<()> {
        self.write_record(record)
    }

    /// Writes a materialized BAM payload.
    ///
    /// The first call writes the payload's header. A later call appends the
    /// payload's records only if its reference dictionary (names and lengths)
    /// is the one already written, since record reference IDs index it;
    /// otherwise it returns an [`io::ErrorKind::InvalidInput`] error without
    /// writing anything.
    pub fn write_all(&mut self, bam: &Bam) -> io::Result<()> {
        match &self.references {
            None => self.write_header(&bam.header, &bam.refs)?,
            Some(references) => {
                let same = references.len() == bam.refs.len()
                    && references
                        .iter()
                        .zip(&bam.refs)
                        .all(|((name, length), reference)| {
                            *name == reference.name && *length == reference.l_seq
                        });
                if !same {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "BAM payload's reference dictionary differs from the one already written",
                    ));
                }
            }
        }
        for record in &bam.records {
            self.write_record(record)?;
        }
        Ok(())
    }

    /// Flushes pending BGZF blocks to the underlying writer.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }

    /// Finishes the BGZF stream, writes the EOF block, and returns the wrapped byte stream.
    pub fn finish(self) -> io::Result<W> {
        self.writer.finish()
    }
}

/// Iterator over records from a [`BamReader`].
pub struct BamRecords<'a, R: Read> {
    reader: &'a mut BamReader<R>,
}

impl<R: Read> Iterator for BamRecords<'_, R> {
    type Item = io::Result<BamRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_record() {
            Ok(Some(record)) => Some(Ok(record)),
            Ok(None) => None,
            Err(err) => Some(Err(err)),
        }
    }
}

impl BamHeader {
    /// Reads a BAM header from a decompressed BAM byte stream.
    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if magic != [0x42, 0x41, 0x4D, 0x01] {
            return Err(invalid_data("invalid BAM magic bytes"));
        }

        let l_text = read_u32_le(reader, "BAM header text length")?;
        let mut text = read_exact_vec(reader, l_text, "BAM header text")?;
        let end = text
            .iter()
            .rposition(|&byte| byte != 0)
            .map_or(0, |index| index + 1);
        text.truncate(end);

        let n_ref = read_u32_le(reader, "BAM reference count")?;

        Ok(Self {
            magic,
            l_text,
            text,
            n_ref,
        })
    }
}

impl BamHeader {
    /// Returns the header text as UTF-8.
    ///
    /// Returns an `InvalidData` error when the stored text is not UTF-8, which
    /// SAM header text must be.
    pub fn text_str(&self) -> io::Result<&str> {
        std::str::from_utf8(&self.text)
            .map_err(|error| invalid_data(format!("BAM header text is not UTF-8: {error}")))
    }
}

impl BamRef {
    /// Reads one reference dictionary entry from a decompressed BAM stream.
    ///
    /// The name must be UTF-8: it is a SAM reference name, and replacing bytes
    /// would silently rename the reference.
    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let l_name = read_u32_le(reader, "reference name length")?;
        let mut name = read_exact_vec(reader, l_name, "BAM reference name")?;
        let end = name
            .iter()
            .rposition(|&byte| byte != 0)
            .map_or(0, |index| index + 1);
        name.truncate(end);
        let name =
            String::from_utf8(name).map_err(|_| invalid_data("BAM reference name is not UTF-8"))?;

        let l_seq = read_u32_le(reader, "reference sequence length")?;

        Ok(Self {
            l_name,
            name,
            l_seq,
        })
    }
}

fn bam_record_to_sam_record(record: &BamRecord, refs: &[BamRef]) -> io::Result<sam::SamRecord> {
    let rname = reference_name(refs, record.fixed.ref_id)?;
    let rnext = mate_reference_name(refs, record.fixed.ref_id, record.fixed.next_ref_id)?;
    let pos = sam_position(record.fixed.pos)?;
    let pnext = sam_position(record.fixed.next_pos)?;
    let seq = if record.fixed.l_seq == 0 {
        "*".to_string()
    } else {
        record.sequence_string()
    };
    let long_cigar_tag = record.long_cigar_tag();
    if long_cigar_tag.is_some() && record.cigar_ops().iter().any(|op| op & 0x0f > 8) {
        return Err(invalid_data("BAM CG tag holds an invalid CIGAR operation"));
    }
    let optional = record
        .auxiliary
        .iter()
        .enumerate()
        .filter(|(index, _)| Some(*index) != long_cigar_tag)
        .map(|(_, auxiliary)| bam_aux_to_sam_optional(auxiliary))
        .collect::<io::Result<Vec<_>>>()?;

    Ok(sam::SamRecord::new(
        record.read_name().to_string(),
        record.fixed.flag,
        rname,
        pos,
        record.fixed.mapq,
        record.cigar_string(),
        rnext,
        pnext,
        record.fixed.tlen,
        seq,
        record.quality_string(),
        optional,
    ))
}

fn reference_name(refs: &[BamRef], ref_id: i32) -> io::Result<String> {
    if ref_id < 0 {
        return Ok("*".to_string());
    }

    refs.get(ref_id as usize)
        .map(|reference| reference.name.clone())
        .ok_or_else(|| invalid_data("BAM record reference ID is out of bounds"))
}

fn mate_reference_name(refs: &[BamRef], ref_id: i32, next_ref_id: i32) -> io::Result<String> {
    if next_ref_id < 0 {
        return Ok("*".to_string());
    }
    if next_ref_id == ref_id && ref_id >= 0 {
        return Ok("=".to_string());
    }

    reference_name(refs, next_ref_id)
}

fn sam_position(position: i32) -> io::Result<u32> {
    if position < 0 {
        return Ok(0);
    }

    u32::try_from(position)
        .ok()
        .and_then(|position| position.checked_add(1))
        .ok_or_else(|| invalid_data("BAM position cannot be represented as SAM POS"))
}

fn bam_aux_to_sam_optional(auxiliary: &BamRecordAuxiliary) -> io::Result<sam::SamOptionalField> {
    let value = match &auxiliary.value {
        BamAuxValue::A(value) => sam::SamOptionalValue::Character(char::from(*value)),
        BamAuxValue::c(value) => sam::SamOptionalValue::Integer(i64::from(*value)),
        BamAuxValue::C(value) => sam::SamOptionalValue::Integer(i64::from(*value)),
        BamAuxValue::s(value) => sam::SamOptionalValue::Integer(i64::from(*value)),
        BamAuxValue::S(value) => sam::SamOptionalValue::Integer(i64::from(*value)),
        BamAuxValue::i(value) => sam::SamOptionalValue::Integer(i64::from(*value)),
        BamAuxValue::I(value) => sam::SamOptionalValue::Integer(i64::from(*value)),
        BamAuxValue::f(value) => sam::SamOptionalValue::Float(*value),
        BamAuxValue::Z(value) => sam::SamOptionalValue::String(value.clone()),
        BamAuxValue::H(value) => sam::SamOptionalValue::Hex(hex_string_to_bytes(value)?),
        BamAuxValue::B(value) => sam::SamOptionalValue::Array(bam_array_to_sam_array(value)),
    };

    Ok(sam::SamOptionalField::new(auxiliary.tag.clone(), value))
}

fn bam_array_to_sam_array(array: &BamAuxArray) -> sam::SamOptionalArray {
    match array {
        BamAuxArray::c(values) => sam::SamOptionalArray::Int8(values.clone()),
        BamAuxArray::C(values) => sam::SamOptionalArray::UInt8(values.clone()),
        BamAuxArray::s(values) => sam::SamOptionalArray::Int16(values.clone()),
        BamAuxArray::S(values) => sam::SamOptionalArray::UInt16(values.clone()),
        BamAuxArray::i(values) => sam::SamOptionalArray::Int32(values.clone()),
        BamAuxArray::I(values) => sam::SamOptionalArray::UInt32(values.clone()),
        BamAuxArray::f(values) => sam::SamOptionalArray::Float(values.clone()),
    }
}

fn format_sam_line(record: &sam::SamRecord) -> io::Result<String> {
    record.to_sam_line()
}

fn bam_refs_from_sam_header(header: &sam::SamHeader) -> io::Result<Vec<BamRef>> {
    let mut refs = Vec::new();
    for record in header.records_of_type("SQ") {
        let name = record
            .value("SN")
            .ok_or_else(|| invalid_data("SAM @SQ record missing SN"))?;
        let length = record
            .value("LN")
            .ok_or_else(|| invalid_data("SAM @SQ record missing LN"))?
            .parse::<u32>()
            .map_err(|_| invalid_data("SAM @SQ LN is not a u32"))?;
        refs.push(BamRef {
            l_name: (name.len() + 1) as u32,
            name: name.to_string(),
            l_seq: length,
        });
    }

    Ok(refs)
}

fn reference_ids_from_refs(refs: &[BamRef]) -> HashMap<String, i32> {
    refs.iter()
        .enumerate()
        .map(|(index, reference)| (reference.name.clone(), index as i32))
        .collect()
}

fn sam_header_text(header: &sam::SamHeader) -> io::Result<String> {
    let mut output = String::new();
    for record in &header.records {
        output.push_str(&sam_header_record_line(record)?);
        output.push('\n');
    }
    Ok(output)
}

fn sam_header_record_line(record: &sam::SamHeaderRecord) -> io::Result<String> {
    if record.record_type == "CO" {
        return Ok(format!(
            "@CO\t{}",
            record.comment.as_deref().unwrap_or_default()
        ));
    }

    if record.record_type.len() != 2 {
        return Err(invalid_data(
            "SAM header record type must be two characters",
        ));
    }

    let mut line = format!("@{}", record.record_type);
    for field in &record.fields {
        line.push('\t');
        line.push_str(&field.tag);
        line.push(':');
        line.push_str(&field.value);
    }
    Ok(line)
}

fn bam_record_from_sam_record(
    record: &sam::SamRecord,
    reference_ids: &HashMap<String, i32>,
) -> io::Result<BamRecord> {
    record.validate()?;

    let ref_id = sam_reference_id(&record.rname, reference_ids)?;
    let next_ref_id = if record.rnext == "=" {
        ref_id
    } else {
        sam_reference_id(&record.rnext, reference_ids)?
    };
    let pos = bam_position(record.pos)?;
    let next_pos = bam_position(record.pnext)?;
    let mut cigar = pack_cigar(&record.cigar)?;
    let seq = pack_sequence(&record.seq)?;
    let l_seq = if record.seq == "*" {
        0
    } else {
        u32::try_from(record.seq.len()).map_err(|_| invalid_data("SAM SEQ length exceeds u32"))?
    };
    let qual = pack_quality(&record.qual, l_seq as usize)?;
    let read_name_len = record
        .qname
        .len()
        .checked_add(1)
        .ok_or_else(|| invalid_data("SAM QNAME length overflow"))?;
    let l_read_name =
        u8::try_from(read_name_len).map_err(|_| invalid_data("SAM QNAME is too long for BAM"))?;
    let bin = bam_bin(record, pos)?;
    let mut auxiliary = record
        .optional
        .iter()
        .map(sam_optional_to_bam_aux)
        .collect::<io::Result<Vec<_>>>()?;
    if cigar.len() > usize::from(u16::MAX) {
        // SAM spec 4.2.2: store a `kSmN` placeholder and the real CIGAR in CG.
        if auxiliary.iter().any(|aux| aux.tag == "CG") {
            return Err(invalid_data(
                "SAM record has more than 65535 CIGAR operations and its own CG tag",
            ));
        }
        let reference_len = cigar
            .iter()
            .filter(|&&op| matches!(op & 0x0f, 0 | 2 | 3 | 7 | 8))
            .try_fold(0u32, |total, op| total.checked_add(op >> 4))
            .filter(|&total| total <= MAX_CIGAR_OP_LEN)
            .ok_or_else(|| invalid_data("SAM CIGAR reference length exceeds BAM's limit"))?;
        if l_seq > MAX_CIGAR_OP_LEN {
            return Err(invalid_data(
                "SAM SEQ is too long for a BAM long-CIGAR record",
            ));
        }
        let real = std::mem::replace(&mut cigar, vec![l_seq << 4 | 4, reference_len << 4 | 3]);
        auxiliary.push(BamRecordAuxiliary {
            tag: "CG".to_string(),
            value: BamAuxValue::B(BamAuxArray::I(real)),
        });
    }
    let n_cigar_op = u16::try_from(cigar.len()).expect("long CIGARs use a 2-operation placeholder");

    let mut bam_record = BamRecord {
        fixed: BamRecordFixed {
            block_size: 0,
            ref_id,
            pos,
            l_read_name,
            mapq: record.mapq,
            bin,
            n_cigar_op,
            flag: record.flag,
            l_seq,
            next_ref_id,
            next_pos,
            tlen: record.tlen,
        },
        variable: BamRecordVariable {
            read_name: record.qname.clone(),
            cigar,
            seq,
            qual,
        },
        auxiliary,
    };

    let mut payload = Vec::new();
    encode_bam_record_payload(&bam_record, &mut payload)?;
    bam_record.fixed.block_size =
        u32::try_from(payload.len()).map_err(|_| invalid_data("BAM record is too large"))?;
    Ok(bam_record)
}

fn sam_reference_id(value: &str, reference_ids: &HashMap<String, i32>) -> io::Result<i32> {
    if value == "*" {
        return Ok(-1);
    }

    reference_ids
        .get(value)
        .copied()
        // The SAM input is at fault here, not a BAM.
        .ok_or_else(|| {
            Error::invalid(
                Format::Sam,
                format!("SAM reference {value} is missing from @SQ header"),
            )
            .into()
        })
}

fn bam_position(position: u32) -> io::Result<i32> {
    if position == 0 {
        return Ok(-1);
    }
    i32::try_from(position - 1).map_err(|_| invalid_data("SAM position exceeds BAM i32 range"))
}

fn pack_cigar(cigar: &str) -> io::Result<Vec<u32>> {
    if cigar == "*" {
        return Ok(Vec::new());
    }

    let mut packed = Vec::new();
    let mut len = 0u32;
    let mut saw_digit = false;
    for byte in cigar.bytes() {
        match byte {
            b'0'..=b'9' => {
                saw_digit = true;
                len = len
                    .checked_mul(10)
                    .and_then(|value| value.checked_add(u32::from(byte - b'0')))
                    .ok_or_else(|| invalid_data("SAM CIGAR operation length overflow"))?;
            }
            b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'=' | b'X' => {
                // Zero-length operations are allowed, as in SAMv1 and htslib.
                if !saw_digit {
                    return Err(invalid_data("invalid SAM CIGAR operation length"));
                }
                if len > MAX_CIGAR_OP_LEN {
                    return Err(invalid_data(
                        "SAM CIGAR operation length exceeds BAM's 28-bit limit",
                    ));
                }
                packed.push((len << 4) | cigar_op_code(byte));
                len = 0;
                saw_digit = false;
            }
            _ => return Err(invalid_data("invalid SAM CIGAR operation")),
        }
    }

    if saw_digit {
        return Err(invalid_data("SAM CIGAR ends with length but no operation"));
    }
    Ok(packed)
}

/// Largest CIGAR operation length a BAM record can store (28 bits).
const MAX_CIGAR_OP_LEN: u32 = (1 << 28) - 1;

fn cigar_op_code(op: u8) -> u32 {
    match op {
        b'M' => 0,
        b'I' => 1,
        b'D' => 2,
        b'N' => 3,
        b'S' => 4,
        b'H' => 5,
        b'P' => 6,
        b'=' => 7,
        b'X' => 8,
        _ => unreachable!("validated CIGAR op"),
    }
}

fn pack_sequence(seq: &str) -> io::Result<Vec<u8>> {
    if seq == "*" {
        return Ok(Vec::new());
    }

    let mut packed = Vec::with_capacity(seq.len().div_ceil(2));
    for chunk in seq.as_bytes().chunks(2) {
        let high = base_code(chunk[0]);
        let low = chunk.get(1).copied().map_or(0, base_code);
        packed.push((high << 4) | low);
    }
    Ok(packed)
}

/// Returns a base's 4-bit BAM code.
///
/// As SAM spec 4.2.3 says, characters without a code are stored as `N`;
/// like htslib, `U` is stored as `T`.
fn base_code(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'=' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' | b'U' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        _ => 15,
    }
}

fn pack_quality(qual: &str, len: usize) -> io::Result<Vec<u8>> {
    if qual == "*" {
        return Ok(vec![0xff; len]);
    }
    if qual.len() != len {
        return Err(invalid_data("SAM QUAL length does not match SEQ length"));
    }

    qual.bytes()
        .map(|byte| {
            if !(33..=126).contains(&byte) {
                return Err(invalid_data("SAM QUAL contains an invalid byte"));
            }
            Ok(byte - 33)
        })
        .collect()
}

/// Computes the BAM `bin` field as htslib does.
///
/// SAM spec 4.2.1 treats an unmapped record, or an alignment that consumes no
/// reference bases, as length one, so an unplaced read (`pos` -1) gets
/// `reg2bin(-1, 0)` = 4680.
fn bam_bin(record: &sam::SamRecord, pos: i32) -> io::Result<u16> {
    let reference_len = if record.is_unmapped() || record.cigar == "*" {
        0
    } else {
        record.reference_len_from_cigar()?
    };
    let beg = i64::from(pos);
    let end = i64::try_from(reference_len.max(1))
        .ok()
        .and_then(|len| beg.checked_add(len))
        .ok_or_else(|| invalid_data("SAM alignment end overflows"))?;
    Ok(reg2bin(beg, end))
}

/// htslib's `hts_reg2bin(beg, end, 14, 5)`, cut to BAM's 16-bit field as
/// htslib does: bins for positions past 2^29 don't fit the BAI scheme.
fn reg2bin(beg: i64, end: i64) -> u16 {
    let end = end - 1;
    let bin = if beg >> 14 == end >> 14 {
        4681 + (beg >> 14)
    } else if beg >> 17 == end >> 17 {
        585 + (beg >> 17)
    } else if beg >> 20 == end >> 20 {
        73 + (beg >> 20)
    } else if beg >> 23 == end >> 23 {
        9 + (beg >> 23)
    } else if beg >> 26 == end >> 26 {
        1 + (beg >> 26)
    } else {
        0
    };
    bin as u16
}

fn sam_optional_to_bam_aux(field: &sam::SamOptionalField) -> io::Result<BamRecordAuxiliary> {
    Ok(BamRecordAuxiliary {
        tag: field.tag.clone(),
        value: sam_optional_value_to_bam(&field.value)?,
    })
}

fn sam_optional_value_to_bam(value: &sam::SamOptionalValue) -> io::Result<BamAuxValue> {
    match value {
        sam::SamOptionalValue::Character(value) => {
            if !value.is_ascii() {
                return Err(invalid_data("SAM A optional value is not ASCII"));
            }
            Ok(BamAuxValue::A(*value as u8))
        }
        sam::SamOptionalValue::Integer(value) => {
            if let Ok(value) = i32::try_from(*value) {
                Ok(BamAuxValue::i(value))
            } else if let Ok(value) = u32::try_from(*value) {
                Ok(BamAuxValue::I(value))
            } else {
                Err(invalid_data(
                    "SAM integer optional value is out of BAM range",
                ))
            }
        }
        sam::SamOptionalValue::Float(value) => Ok(BamAuxValue::f(*value)),
        sam::SamOptionalValue::String(value) => Ok(BamAuxValue::Z(value.clone())),
        sam::SamOptionalValue::Hex(value) => Ok(BamAuxValue::H(bytes_to_hex(value))),
        sam::SamOptionalValue::Array(value) => Ok(BamAuxValue::B(sam_array_to_bam_array(value))),
    }
}

fn sam_array_to_bam_array(array: &sam::SamOptionalArray) -> BamAuxArray {
    match array {
        sam::SamOptionalArray::Int8(values) => BamAuxArray::c(values.clone()),
        sam::SamOptionalArray::UInt8(values) => BamAuxArray::C(values.clone()),
        sam::SamOptionalArray::Int16(values) => BamAuxArray::s(values.clone()),
        sam::SamOptionalArray::UInt16(values) => BamAuxArray::S(values.clone()),
        sam::SamOptionalArray::Int32(values) => BamAuxArray::i(values.clone()),
        sam::SamOptionalArray::UInt32(values) => BamAuxArray::I(values.clone()),
        sam::SamOptionalArray::Float(values) => BamAuxArray::f(values.clone()),
    }
}

fn bytes_to_hex(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789ABCDEF";
    let mut output = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        output.push(HEX[(byte >> 4) as usize] as char);
        output.push(HEX[(byte & 0x0f) as usize] as char);
    }
    output
}

fn hex_string_to_bytes(value: &str) -> io::Result<Vec<u8>> {
    if !value.len().is_multiple_of(2) {
        return Err(invalid_data("BAM H auxiliary value has odd length"));
    }

    // Decode byte pairs: slicing the string could split a character, and
    // u8::from_str_radix would accept a '+' sign.
    value
        .as_bytes()
        .chunks_exact(2)
        .map(|pair| Some(hex_digit(pair[0])? << 4 | hex_digit(pair[1])?))
        .collect::<Option<Vec<u8>>>()
        .ok_or_else(|| invalid_data("BAM H auxiliary value is not hex"))
}

fn hex_digit(byte: u8) -> Option<u8> {
    char::from(byte).to_digit(16).map(|digit| digit as u8)
}

/// Reads exactly `len` bytes, where `len` comes from the file.
///
/// The buffer grows as data arrives instead of being sized from `len` up
/// front, so a bogus length can't force a huge allocation.
fn read_exact_vec<R: Read>(reader: &mut R, len: u32, what: &str) -> io::Result<Vec<u8>> {
    const INITIAL_CAPACITY: usize = 64 * 1024;
    let mut data = Vec::with_capacity((len as usize).min(INITIAL_CAPACITY));
    reader.take(u64::from(len)).read_to_end(&mut data)?;
    if data.len() != len as usize {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!("{what} is truncated"),
        ));
    }
    Ok(data)
}

fn read_u32_le<R: Read>(reader: &mut R, field: &str) -> io::Result<u32> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes).map_err(|err| {
        // An error that already carries a Brust diagnostic (from BGZF, say)
        // is passed on as it is; only plain I/O errors get the field name.
        if err.get_ref().is_some_and(|inner| inner.is::<Error>()) {
            err
        } else {
            io::Error::new(err.kind(), format!("failed to read {field}: {err}"))
        }
    })?;
    Ok(u32::from_le_bytes(bytes))
}

/// Fills `buf`, returning `false` at a clean EOF before the first byte; `what`
/// names the field for a truncation error.
fn read_exact_or_eof<R: Read>(reader: &mut R, buf: &mut [u8], what: &str) -> io::Result<bool> {
    let mut bytes_read = 0;
    while bytes_read < buf.len() {
        match reader.read(&mut buf[bytes_read..]) {
            Ok(0) if bytes_read == 0 => return Ok(false),
            Ok(0) => {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!("unexpected EOF while reading {what}"),
                ));
            }
            Ok(n) => bytes_read += n,
            Err(err) if err.kind() == io::ErrorKind::Interrupted => {}
            Err(err) => return Err(err),
        }
    }

    Ok(true)
}

fn write_bam_header<W: Write>(
    writer: &mut W,
    header: &BamHeader,
    refs: &[BamRef],
) -> io::Result<()> {
    if header.magic != [0x42, 0x41, 0x4D, 0x01] {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid BAM magic bytes",
        ));
    }
    if header.n_ref as usize != refs.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "BAM header reference count does not match refs length",
        ));
    }

    let text = header.text.as_slice();
    if header.l_text as usize > 0 && (header.l_text as usize) < text.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "BAM header text length is smaller than text bytes",
        ));
    }
    let l_text = if header.l_text == 0 {
        text.len() as u32
    } else {
        header.l_text
    };

    writer.write_all(&header.magic)?;
    writer.write_all(&l_text.to_le_bytes())?;
    writer.write_all(text)?;
    for _ in text.len()..l_text as usize {
        writer.write_all(&[0])?;
    }
    writer.write_all(&header.n_ref.to_le_bytes())?;

    for reference in refs {
        write_bam_ref(writer, reference)?;
    }

    Ok(())
}

fn write_bam_ref<W: Write>(writer: &mut W, reference: &BamRef) -> io::Result<()> {
    let name = reference.name.as_bytes();
    if reference.l_name as usize > 0 && (reference.l_name as usize) < name.len() + 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "BAM reference name length is smaller than name plus NUL",
        ));
    }
    let l_name = if reference.l_name == 0 {
        (name.len() + 1) as u32
    } else {
        reference.l_name
    };

    writer.write_all(&l_name.to_le_bytes())?;
    writer.write_all(name)?;
    writer.write_all(&[0])?;
    for _ in name.len() + 1..l_name as usize {
        writer.write_all(&[0])?;
    }
    writer.write_all(&reference.l_seq.to_le_bytes())
}

fn encode_bam_record_payload(record: &BamRecord, payload: &mut Vec<u8>) -> io::Result<()> {
    let fixed = &record.fixed;
    let variable = &record.variable;
    let read_name = variable.read_name.as_bytes();

    if fixed.l_read_name as usize > 0 && (fixed.l_read_name as usize) < read_name.len() + 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "BAM read name length is smaller than name plus NUL",
        ));
    }
    if fixed.n_cigar_op as usize != variable.cigar.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "BAM n_cigar_op does not match CIGAR vector length",
        ));
    }
    if fixed.l_seq.div_ceil(2) as usize != variable.seq.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "BAM packed sequence length does not match l_seq",
        ));
    }
    if fixed.l_seq as usize != variable.qual.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "BAM quality length does not match l_seq",
        ));
    }

    let l_read_name = if fixed.l_read_name == 0 {
        u8::try_from(read_name.len() + 1).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "BAM read name is longer than 254 bytes",
            )
        })?
    } else {
        fixed.l_read_name
    };

    payload.extend_from_slice(&fixed.ref_id.to_le_bytes());
    payload.extend_from_slice(&fixed.pos.to_le_bytes());
    payload.push(l_read_name);
    payload.push(fixed.mapq);
    payload.extend_from_slice(&fixed.bin.to_le_bytes());
    payload.extend_from_slice(&fixed.n_cigar_op.to_le_bytes());
    payload.extend_from_slice(&fixed.flag.to_le_bytes());
    payload.extend_from_slice(&fixed.l_seq.to_le_bytes());
    payload.extend_from_slice(&fixed.next_ref_id.to_le_bytes());
    payload.extend_from_slice(&fixed.next_pos.to_le_bytes());
    payload.extend_from_slice(&fixed.tlen.to_le_bytes());
    payload.extend_from_slice(read_name);
    payload.push(0);
    for _ in read_name.len() + 1..l_read_name as usize {
        payload.push(0);
    }
    for cigar in &variable.cigar {
        payload.extend_from_slice(&cigar.to_le_bytes());
    }
    payload.extend_from_slice(&variable.seq);
    payload.extend_from_slice(&variable.qual);
    for auxiliary in &record.auxiliary {
        encode_bam_auxiliary(auxiliary, payload)?;
    }

    Ok(())
}

fn encode_bam_auxiliary(auxiliary: &BamRecordAuxiliary, payload: &mut Vec<u8>) -> io::Result<()> {
    let tag = auxiliary.tag.as_bytes();
    if tag.len() != 2 || !tag.iter().all(u8::is_ascii_alphanumeric) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "BAM auxiliary tag must contain two alphanumeric ASCII bytes",
        ));
    }
    payload.extend_from_slice(tag);

    match &auxiliary.value {
        BamAuxValue::A(value) => {
            payload.push(b'A');
            payload.push(*value);
        }
        BamAuxValue::c(value) => {
            payload.push(b'c');
            payload.push(*value as u8);
        }
        BamAuxValue::C(value) => {
            payload.push(b'C');
            payload.push(*value);
        }
        BamAuxValue::s(value) => {
            payload.push(b's');
            payload.extend_from_slice(&value.to_le_bytes());
        }
        BamAuxValue::S(value) => {
            payload.push(b'S');
            payload.extend_from_slice(&value.to_le_bytes());
        }
        BamAuxValue::i(value) => {
            payload.push(b'i');
            payload.extend_from_slice(&value.to_le_bytes());
        }
        BamAuxValue::I(value) => {
            payload.push(b'I');
            payload.extend_from_slice(&value.to_le_bytes());
        }
        BamAuxValue::f(value) => {
            payload.push(b'f');
            payload.extend_from_slice(&value.to_le_bytes());
        }
        BamAuxValue::Z(value) => encode_bam_aux_string(payload, b'Z', value)?,
        BamAuxValue::H(value) => encode_bam_aux_string(payload, b'H', value)?,
        BamAuxValue::B(value) => encode_bam_aux_array(payload, value)?,
    }

    Ok(())
}

fn encode_bam_aux_string(payload: &mut Vec<u8>, value_type: u8, value: &str) -> io::Result<()> {
    if value.as_bytes().contains(&0) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "BAM auxiliary string must not contain NUL bytes",
        ));
    }
    payload.push(value_type);
    payload.extend_from_slice(value.as_bytes());
    payload.push(0);
    Ok(())
}

fn encode_bam_aux_array(payload: &mut Vec<u8>, value: &BamAuxArray) -> io::Result<()> {
    payload.push(b'B');
    match value {
        BamAuxArray::c(values) => {
            encode_bam_array_header(payload, b'c', values.len())?;
            payload.extend(values.iter().map(|value| *value as u8));
        }
        BamAuxArray::C(values) => {
            encode_bam_array_header(payload, b'C', values.len())?;
            payload.extend_from_slice(values);
        }
        BamAuxArray::s(values) => {
            encode_bam_array_header(payload, b's', values.len())?;
            for value in values {
                payload.extend_from_slice(&value.to_le_bytes());
            }
        }
        BamAuxArray::S(values) => {
            encode_bam_array_header(payload, b'S', values.len())?;
            for value in values {
                payload.extend_from_slice(&value.to_le_bytes());
            }
        }
        BamAuxArray::i(values) => {
            encode_bam_array_header(payload, b'i', values.len())?;
            for value in values {
                payload.extend_from_slice(&value.to_le_bytes());
            }
        }
        BamAuxArray::I(values) => {
            encode_bam_array_header(payload, b'I', values.len())?;
            for value in values {
                payload.extend_from_slice(&value.to_le_bytes());
            }
        }
        BamAuxArray::f(values) => {
            encode_bam_array_header(payload, b'f', values.len())?;
            for value in values {
                payload.extend_from_slice(&value.to_le_bytes());
            }
        }
    }

    Ok(())
}

fn encode_bam_array_header(payload: &mut Vec<u8>, subtype: u8, len: usize) -> io::Result<()> {
    let len = u32::try_from(len).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "BAM auxiliary array length exceeds u32",
        )
    })?;
    payload.push(subtype);
    payload.extend_from_slice(&len.to_le_bytes());
    Ok(())
}

fn invalid_data(message: impl Into<String>) -> io::Error {
    Error::invalid(Format::Bam, message).into()
}

#[cfg(test)]
mod bam_tests {
    use super::*;
    use std::collections::HashMap;
    use std::fs::File;
    use std::io;
    use std::path::PathBuf;

    const ALIGNED_BAM: &str = "aligned.bam";
    const UNALIGNED_BAM: &str = "unaligned.bam";

    fn fixture_path(name: &str) -> String {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join(name)
            .to_string_lossy()
            .into_owned()
    }

    fn workspace_fixture_path(relative: &str) -> String {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("..")
            .join(relative)
            .to_string_lossy()
            .into_owned()
    }

    fn open_bam(name: &str) -> BamReader {
        let path = fixture_path(name);
        BamReader::new(&path).unwrap_or_else(|err| panic!("failed to open {name}: {err}"))
    }

    fn read_all_records(name: &str) -> Vec<BamRecord> {
        let mut bam = open_bam(name);
        let mut records = Vec::new();
        while let Some(record) = bam
            .read()
            .unwrap_or_else(|err| panic!("failed to read {name}: {err}"))
        {
            records.push(record);
        }
        assert!(
            bam.read()
                .expect("second EOF read should succeed")
                .is_none()
        );
        records
    }

    fn push_i16(data: &mut Vec<u8>, value: i16) {
        data.extend_from_slice(&value.to_le_bytes());
    }

    fn push_u16(data: &mut Vec<u8>, value: u16) {
        data.extend_from_slice(&value.to_le_bytes());
    }

    fn push_i32(data: &mut Vec<u8>, value: i32) {
        data.extend_from_slice(&value.to_le_bytes());
    }

    fn push_u32(data: &mut Vec<u8>, value: u32) {
        data.extend_from_slice(&value.to_le_bytes());
    }

    fn push_f32(data: &mut Vec<u8>, value: f32) {
        data.extend_from_slice(&value.to_le_bytes());
    }

    #[allow(clippy::too_many_arguments)]
    fn fixed_core(
        ref_id: i32,
        pos: i32,
        l_read_name: u8,
        mapq: u8,
        bin: u16,
        n_cigar_op: u16,
        flag: u16,
        l_seq: u32,
        next_ref_id: i32,
        next_pos: i32,
        tlen: i32,
    ) -> Vec<u8> {
        let mut data = Vec::new();
        push_i32(&mut data, ref_id);
        push_i32(&mut data, pos);
        data.push(l_read_name);
        data.push(mapq);
        push_u16(&mut data, bin);
        push_u16(&mut data, n_cigar_op);
        push_u16(&mut data, flag);
        push_u32(&mut data, l_seq);
        push_i32(&mut data, next_ref_id);
        push_i32(&mut data, next_pos);
        push_i32(&mut data, tlen);
        data
    }

    fn packed_cigar(len: u32, op: u32) -> u32 {
        (len << 4) | op
    }

    fn decode_cigar(cigar: &[u32]) -> String {
        const OPS: &[u8; 9] = b"MIDNSHP=X";
        let mut decoded = String::new();
        for op in cigar {
            decoded.push_str(&(op >> 4).to_string());
            decoded.push(OPS[(op & 0x0f) as usize] as char);
        }
        decoded
    }

    fn decode_seq(seq: &[u8], l_seq: u32) -> String {
        const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
        let mut decoded = String::with_capacity(l_seq as usize);
        for i in 0..l_seq as usize {
            let byte = seq[i / 2];
            let code = if i % 2 == 0 { byte >> 4 } else { byte & 0x0f };
            decoded.push(BASES[code as usize] as char);
        }
        decoded
    }

    fn parse_aux(data: Vec<u8>) -> BamRecordAuxiliary {
        let mut offset = 0;
        let aux = BamRecordAuxiliary::new(&mut offset, &data).expect("aux field should parse");
        assert_eq!(offset, data.len());
        aux
    }

    fn aux_value<'a>(record: &'a BamRecord, tag: &str) -> &'a BamAuxValue {
        record
            .aux(tag)
            .unwrap_or_else(|| panic!("missing aux tag {tag}"))
    }

    fn assert_aux_f32(value: &BamAuxValue, expected: f32) {
        match value {
            BamAuxValue::f(actual) => assert!((actual - expected).abs() < 0.0001),
            other => panic!("expected f32 aux value, got {other:?}"),
        }
    }

    fn assert_aux_int(value: &BamAuxValue, expected: i64) {
        let actual = match value {
            BamAuxValue::c(value) => *value as i64,
            BamAuxValue::C(value) => *value as i64,
            BamAuxValue::s(value) => *value as i64,
            BamAuxValue::S(value) => *value as i64,
            BamAuxValue::i(value) => *value as i64,
            BamAuxValue::I(value) => *value as i64,
            other => panic!("expected integer aux value, got {other:?}"),
        };
        assert_eq!(actual, expected);
    }

    fn assert_aux_f32_array(value: &BamAuxValue, expected: &[f32]) {
        match value {
            BamAuxValue::B(BamAuxArray::f(actual)) => {
                assert_eq!(actual.len(), expected.len());
                for (actual, expected) in actual.iter().zip(expected) {
                    assert!((actual - expected).abs() < 0.0001);
                }
            }
            other => panic!("expected f32 aux array, got {other:?}"),
        }
    }

    fn tag_counts(records: &[BamRecord]) -> HashMap<String, usize> {
        let mut counts = HashMap::new();
        for record in records {
            for aux in &record.auxiliary {
                *counts.entry(aux.tag.clone()).or_insert(0) += 1;
            }
        }
        counts
    }

    #[test]
    fn materialized_bam_can_be_deep_cloned() {
        let original = Bam::from_path(fixture_path(ALIGNED_BAM)).expect("Bam should materialize");
        assert_eq!(original.refs.len(), 1);
        assert_eq!(original.records.len(), 100);

        let mut cloned = original.clone();
        assert_eq!(cloned, original);

        cloned.header.text.extend_from_slice(b"@CO\tclone-only\n");
        cloned.refs[0].name.push_str("_clone");
        cloned.records[0].variable.read_name.push_str("_clone");

        assert_ne!(cloned, original);
        assert_eq!(original.refs[0].name, "fc_reference");
        assert_eq!(
            original.records[0].read_name(),
            "2f890ab0-63b9-45f6-ba8d-3c367ca26a63"
        );
        assert!(!original.header.text_str().unwrap().contains("clone-only"));
    }

    #[test]
    fn reader_can_materialize_all_records() {
        let reader = open_bam(UNALIGNED_BAM);
        let bam = reader.read_all().expect("reader should materialize");

        assert_eq!(bam.header.n_ref, 0);
        assert!(bam.refs.is_empty());
        assert_eq!(bam.records.len(), 100);
    }

    #[test]
    fn reader_records_iterator_streams_until_eof() {
        let mut reader = open_bam(ALIGNED_BAM);
        let records = reader
            .records()
            .collect::<io::Result<Vec<_>>>()
            .expect("records iterator should parse all records");

        assert_eq!(records.len(), 100);
        assert!(
            reader
                .read_record()
                .expect("post-iterator EOF should succeed")
                .is_none()
        );
    }

    #[test]
    fn bam_from_reader_materializes_compressed_stream() {
        let file = File::open(fixture_path(ALIGNED_BAM)).expect("fixture should open");
        let bam = Bam::from_reader(file).expect("Bam should materialize from reader");

        assert_eq!(bam.header.n_ref, 1);
        assert_eq!(bam.refs[0].name, "fc_reference");
        assert_eq!(bam.records.len(), 100);
    }

    #[test]
    fn writer_round_trips_materialized_aligned_bam() {
        let bam = Bam::from_path(fixture_path(ALIGNED_BAM)).expect("Bam should materialize");
        let mut output = Vec::new();
        bam.to_writer(&mut output).expect("Bam should write");
        let round_tripped = Bam::from_reader(&output[..]).expect("written Bam should parse");

        assert_eq!(&output[..4], b"\x1f\x8b\x08\x04");
        assert_eq!(&output[12..16], b"BC\x02\x00");
        assert!(output.ends_with(&bgzf::EOF_BLOCK));
        assert_eq!(round_tripped, bam);
    }

    #[test]
    fn record_encode_matches_writer_bytes() {
        use flate2::read::MultiGzDecoder;

        fn decompress(bytes: &[u8]) -> Vec<u8> {
            let mut decoded = Vec::new();
            MultiGzDecoder::new(bytes)
                .read_to_end(&mut decoded)
                .expect("BGZF stream should decompress");
            decoded
        }

        let bam = Bam::from_path(fixture_path(ALIGNED_BAM)).expect("Bam should materialize");

        let mut header_only = BamWriter::from_writer(Vec::new());
        header_only
            .write_header(&bam.header, &bam.refs)
            .expect("header should write");
        let header_len = decompress(&header_only.finish().expect("writer should finish")).len();

        let mut writer = BamWriter::from_writer(Vec::new());
        writer.write_all(&bam).expect("Bam should write");
        let stream = decompress(&writer.finish().expect("writer should finish"));

        let mut encoded = Vec::new();
        for record in &bam.records {
            record.encode(&mut encoded).expect("record should encode");
        }

        assert_eq!(bam.records.len(), 100);
        assert_eq!(&stream[header_len..], &encoded[..]);
    }

    #[test]
    fn record_encode_rejects_what_write_record_rejects() {
        let bam = Bam::from_path(fixture_path(ALIGNED_BAM)).expect("Bam should materialize");
        let valid = bam.records[0].clone();
        assert!(!valid.variable.qual.is_empty());
        assert!(!valid.auxiliary.is_empty());

        let mut short_name_len = valid.clone();
        short_name_len.fixed.l_read_name = 1;
        assert!(short_name_len.variable.read_name.len() > 1);

        let mut extra_cigar_op = valid.clone();
        extra_cigar_op.fixed.n_cigar_op += 1;

        let mut longer_l_seq = valid.clone();
        longer_l_seq.fixed.l_seq += 2;

        let mut short_qual = valid.clone();
        short_qual.variable.qual.pop();

        let mut bad_aux_tag = valid.clone();
        bad_aux_tag.auxiliary[0].tag = "X".to_string();

        let mut wrong_block_size = valid.clone();
        wrong_block_size.fixed.block_size = 1;

        // The unmodified record is accepted by both paths.
        let mut encoded = Vec::new();
        valid
            .encode(&mut encoded)
            .expect("valid record should encode");
        let mut writer = BamWriter::from_writer(Vec::new());
        writer
            .write_header(&bam.header, &bam.refs)
            .expect("header should write");
        writer
            .write_record(&valid)
            .expect("valid record should write");

        for (name, record) in [
            ("l_read_name", short_name_len),
            ("n_cigar_op", extra_cigar_op),
            ("l_seq", longer_l_seq),
            ("qual", short_qual),
            ("aux tag", bad_aux_tag),
            ("block_size", wrong_block_size),
        ] {
            let mut out = Vec::new();
            let encode_err = record
                .encode(&mut out)
                .expect_err(&format!("{name}: encode should fail"));
            assert!(out.is_empty(), "{name}: a failed encode must not write");

            // The header goes in first so write_record reaches record validation.
            let mut writer = BamWriter::from_writer(Vec::new());
            writer
                .write_header(&bam.header, &bam.refs)
                .expect("header should write");
            let write_err = writer
                .write_record(&record)
                .expect_err(&format!("{name}: write_record should fail"));

            assert_eq!(encode_err.to_string(), write_err.to_string(), "{name}");
        }
    }

    #[test]
    fn writer_round_trips_materialized_unaligned_bam() {
        let bam = Bam::from_path(fixture_path(UNALIGNED_BAM)).expect("Bam should materialize");
        let mut output = Vec::new();
        bam.to_writer(&mut output).expect("Bam should write");
        let round_tripped = Bam::from_reader(&output[..]).expect("written Bam should parse");

        assert_eq!(round_tripped, bam);
    }

    #[test]
    fn streaming_writer_round_trips_records() {
        let bam = Bam::from_path(fixture_path(ALIGNED_BAM)).expect("Bam should materialize");
        let mut writer = BamWriter::from_writer(Vec::new());
        writer
            .write_header(&bam.header, &bam.refs)
            .expect("header should write");
        for record in bam.records.iter().take(3) {
            writer.write_record(record).expect("record should write");
        }
        let output = writer.finish().expect("writer should finish");
        let round_tripped = Bam::from_reader(&output[..]).expect("written Bam should parse");

        assert_eq!(round_tripped.header, bam.header);
        assert_eq!(round_tripped.refs, bam.refs);
        assert_eq!(round_tripped.records, bam.records[..3]);
    }

    #[test]
    fn sam_to_bam_converter_streams_records_matching_materialized_conversion() {
        let sam_path = workspace_fixture_path("sam/aligned.sam");
        let mut reader = sam::SamReader::from_path(&sam_path).expect("SAM fixture should open");
        let converter =
            SamToBamConverter::new(&reader.header).expect("converter should derive BAM metadata");
        let mut writer = BamWriter::from_writer(Vec::new());

        writer
            .write_header(converter.header(), converter.refs())
            .expect("converted header should write");

        let mut streamed_record_count = 0;
        while let Some(record) = reader.read_record().expect("SAM record should parse") {
            let record = converter
                .convert_record(&record)
                .expect("SAM record should convert to BAM");
            writer
                .write_record(&record)
                .expect("BAM record should write");
            streamed_record_count += 1;
        }

        let output = writer.finish().expect("streaming BAM writer should finish");
        let streamed = Bam::from_reader(&output[..]).expect("streamed BAM should parse");
        let materialized_sam = sam::Sam::from_path(&sam_path).expect("SAM fixture should load");
        let materialized =
            Bam::from_sam(&materialized_sam).expect("materialized SAM should convert");

        assert_eq!(streamed_record_count, 100);
        assert_eq!(streamed, materialized);
    }

    const ORIENTATION_SAM: &[u8] = b"@HD\tVN:1.6\n\
@SQ\tSN:ref\tLN:100\n\
fwd\t0\tref\t1\t60\t9M\t*\t0\t0\tAACGRYKMN\t!#%')+-/1\n\
rev\t16\tref\t1\t60\t9M\t*\t0\t0\tAACGRYKMN\t!#%')+-/1\n\
revstar\t16\tref\t1\t60\t9M\t*\t0\t0\tAACGRYKMN\t*\n";

    fn orientation_records() -> (Vec<sam::SamRecord>, Vec<BamRecord>) {
        let sam = sam::Sam::from_reader(ORIENTATION_SAM).expect("orientation SAM should parse");
        let converter =
            SamToBamConverter::new(&sam.header).expect("converter should derive BAM metadata");
        let bam = sam
            .records
            .iter()
            .map(|record| converter.convert_record(record))
            .collect::<io::Result<Vec<_>>>()
            .expect("SAM records should convert to BAM");
        (sam.records, bam)
    }

    #[test]
    fn bam_original_orientation_matches_sam() {
        let (sam_records, bam_records) = orientation_records();

        assert_eq!(sam_records[0].flag, 0);
        assert_eq!(sam_records[1].flag, 16);
        for (sam_record, bam_record) in sam_records.iter().zip(&bam_records) {
            assert_eq!(
                bam_record.original_sequence_string(),
                sam_record.original_seq()
            );
            assert_eq!(
                bam_record.original_quality_string(),
                sam_record.original_qual()
            );
        }

        assert_eq!(bam_records[0].original_sequence_string(), "AACGRYKMN");
        assert_eq!(bam_records[0].original_quality_string(), "!#%')+-/1");
        assert_eq!(bam_records[1].original_sequence_string(), "NKMRYCGTT");
        assert_eq!(bam_records[1].original_quality_string(), "1/-+)'%#!");
    }

    #[test]
    fn bam_reverse_record_without_qualities_reports_star() {
        let (_, bam_records) = orientation_records();
        let record = &bam_records[2];

        assert!(record.is_reverse_complemented());
        assert_eq!(record.quality_string(), "*");
        assert_eq!(record.original_quality_string(), "*");
        assert_eq!(record.original_sequence_string(), "NKMRYCGTT");
    }

    #[test]
    fn fixed_core_parser_reads_all_fields() {
        let data = fixed_core(7, 42, 9, 60, 4681, 3, 16, 123, -1, -1, -250);
        let fixed = BamRecordFixed::new(321, &data).expect("fixed core should parse");

        assert_eq!(fixed.block_size, 321);
        assert_eq!(fixed.ref_id, 7);
        assert_eq!(fixed.pos, 42);
        assert_eq!(fixed.l_read_name, 9);
        assert_eq!(fixed.mapq, 60);
        assert_eq!(fixed.bin, 4681);
        assert_eq!(fixed.n_cigar_op, 3);
        assert_eq!(fixed.flag, 16);
        assert_eq!(fixed.l_seq, 123);
        assert_eq!(fixed.next_ref_id, -1);
        assert_eq!(fixed.next_pos, -1);
        assert_eq!(fixed.tlen, -250);
    }

    #[test]
    fn variable_fields_parse_offsets_and_packed_bytes() {
        let mut data = fixed_core(0, 10, 6, 20, 0, 2, 0, 5, -1, -1, 0);
        let fixed = BamRecordFixed::new(data.len() as u32, &data).expect("fixed core should parse");
        data.extend_from_slice(b"read1\0");
        push_u32(&mut data, packed_cigar(10, 0));
        push_u32(&mut data, packed_cigar(1, 1));
        data.extend_from_slice(&[0x12, 0x48, 0xf0]);
        data.extend_from_slice(&[30, 31, 32, 33, 34]);

        let (variable, offset) =
            BamRecordVariable::new(&fixed, &data).expect("variable fields should parse");

        assert_eq!(offset, 54);
        assert_eq!(variable.read_name, "read1");
        assert_eq!(decode_cigar(&variable.cigar), "10M1I");
        assert_eq!(decode_seq(&variable.seq, fixed.l_seq), "ACGTN");
        assert_eq!(variable.qual, vec![30, 31, 32, 33, 34]);
    }

    #[test]
    fn record_parser_combines_fixed_variable_and_auxiliary_fields() {
        let mut data = fixed_core(0, 10, 6, 42, 4681, 1, 0, 4, -1, -1, 0);
        data.extend_from_slice(b"read1\0");
        push_u32(&mut data, packed_cigar(4, 0));
        data.extend_from_slice(&[0x12, 0x48]);
        data.extend_from_slice(&[30, 31, 32, 33]);
        data.extend_from_slice(b"NM");
        data.push(b'i');
        push_i32(&mut data, 2);
        data.extend_from_slice(b"tp");
        data.push(b'A');
        data.push(b'P');
        data.extend_from_slice(b"RG");
        data.push(b'Z');
        data.extend_from_slice(b"group1\0");

        let record = BamRecord::new(data.len() as u32, data).expect("record should parse");

        assert_eq!(record.fixed.ref_id, 0);
        assert_eq!(record.variable.read_name, "read1");
        assert_eq!(decode_cigar(&record.variable.cigar), "4M");
        assert_eq!(decode_seq(&record.variable.seq, record.fixed.l_seq), "ACGT");
        assert_eq!(record.variable.qual, vec![30, 31, 32, 33]);
        assert_eq!(record.auxiliary.len(), 3);
        assert_eq!(*aux_value(&record, "NM"), BamAuxValue::i(2));
        assert_eq!(*aux_value(&record, "tp"), BamAuxValue::A(b'P'));
        assert_eq!(
            *aux_value(&record, "RG"),
            BamAuxValue::Z("group1".to_string())
        );
    }

    #[test]
    fn auxiliary_parser_accepts_all_scalar_types() {
        assert_eq!(parse_aux(b"taAZ".to_vec()).value, BamAuxValue::A(b'Z'));
        assert_eq!(
            parse_aux(vec![b't', b'c', b'c', 0xff]).value,
            BamAuxValue::c(-1)
        );
        assert_eq!(
            parse_aux(vec![b't', b'C', b'C', 255]).value,
            BamAuxValue::C(255)
        );

        let mut signed_short = b"tss".to_vec();
        push_i16(&mut signed_short, -1234);
        assert_eq!(parse_aux(signed_short).value, BamAuxValue::s(-1234));

        let mut unsigned_short = b"tSS".to_vec();
        push_u16(&mut unsigned_short, 50000);
        assert_eq!(parse_aux(unsigned_short).value, BamAuxValue::S(50000));

        let mut signed_int = b"tii".to_vec();
        push_i32(&mut signed_int, -123456);
        assert_eq!(parse_aux(signed_int).value, BamAuxValue::i(-123456));

        let mut unsigned_int = b"tII".to_vec();
        push_u32(&mut unsigned_int, 123456);
        assert_eq!(parse_aux(unsigned_int).value, BamAuxValue::I(123456));

        let mut float = b"tff".to_vec();
        push_f32(&mut float, 3.5);
        assert_aux_f32(&parse_aux(float).value, 3.5);

        assert_eq!(
            parse_aux(b"tzZhello\0".to_vec()).value,
            BamAuxValue::Z("hello".to_string())
        );
        assert_eq!(
            parse_aux(b"thH0A0B\0".to_vec()).value,
            BamAuxValue::H("0A0B".to_string())
        );
    }

    #[test]
    fn auxiliary_parser_accepts_all_array_types() {
        let mut b_c = b"bcBc".to_vec();
        push_u32(&mut b_c, 2);
        b_c.extend_from_slice(&[0xff, 2]);
        assert_eq!(
            parse_aux(b_c).value,
            BamAuxValue::B(BamAuxArray::c(vec![-1, 2]))
        );

        let mut b_c_unsigned = b"bCBC".to_vec();
        push_u32(&mut b_c_unsigned, 2);
        b_c_unsigned.extend_from_slice(&[1, 255]);
        assert_eq!(
            parse_aux(b_c_unsigned).value,
            BamAuxValue::B(BamAuxArray::C(vec![1, 255]))
        );

        let mut b_s = b"bsBs".to_vec();
        push_u32(&mut b_s, 2);
        push_i16(&mut b_s, -10);
        push_i16(&mut b_s, 20);
        assert_eq!(
            parse_aux(b_s).value,
            BamAuxValue::B(BamAuxArray::s(vec![-10, 20]))
        );

        let mut b_s_unsigned = b"bSBS".to_vec();
        push_u32(&mut b_s_unsigned, 2);
        push_u16(&mut b_s_unsigned, 10);
        push_u16(&mut b_s_unsigned, 65000);
        assert_eq!(
            parse_aux(b_s_unsigned).value,
            BamAuxValue::B(BamAuxArray::S(vec![10, 65000]))
        );

        let mut b_i = b"biBi".to_vec();
        push_u32(&mut b_i, 2);
        push_i32(&mut b_i, -1000);
        push_i32(&mut b_i, 2000);
        assert_eq!(
            parse_aux(b_i).value,
            BamAuxValue::B(BamAuxArray::i(vec![-1000, 2000]))
        );

        let mut b_i_unsigned = b"bIBI".to_vec();
        push_u32(&mut b_i_unsigned, 2);
        push_u32(&mut b_i_unsigned, 1000);
        push_u32(&mut b_i_unsigned, 2000);
        assert_eq!(
            parse_aux(b_i_unsigned).value,
            BamAuxValue::B(BamAuxArray::I(vec![1000, 2000]))
        );

        let mut b_f = b"bfBf".to_vec();
        push_u32(&mut b_f, 2);
        push_f32(&mut b_f, 1.25);
        push_f32(&mut b_f, 2.5);
        assert_aux_f32_array(&parse_aux(b_f).value, &[1.25, 2.5]);
    }

    #[test]
    fn malformed_records_and_fields_return_errors() {
        let data = fixed_core(0, 0, 1, 0, 0, 0, 0, 0, -1, -1, 0);
        assert!(BamRecordFixed::new(0, &[0; 8]).is_err());
        assert!(BamRecord::new(data.len() as u32 + 1, data.clone()).is_err());
        assert!(BamRecord::new(31, vec![0; 31]).is_err());

        let fixed = BamRecordFixed::new(0, &fixed_core(0, 0, 5, 0, 0, 0, 0, 0, -1, -1, 0))
            .expect("fixed core should parse");
        let mut truncated_name = fixed_core(0, 0, 5, 0, 0, 0, 0, 0, -1, -1, 0);
        truncated_name.extend_from_slice(b"abc");
        assert!(BamRecordVariable::new(&fixed, &truncated_name).is_err());

        let fixed = BamRecordFixed::new(0, &fixed_core(0, 0, 1, 0, 0, 1, 0, 0, -1, -1, 0))
            .expect("fixed core should parse");
        let mut truncated_cigar = fixed_core(0, 0, 1, 0, 0, 1, 0, 0, -1, -1, 0);
        truncated_cigar.extend_from_slice(b"\0");
        truncated_cigar.extend_from_slice(&[0, 0, 0]);
        assert!(BamRecordVariable::new(&fixed, &truncated_cigar).is_err());

        let fixed = BamRecordFixed::new(0, &fixed_core(0, 0, 1, 0, 0, 0, 0, 3, -1, -1, 0))
            .expect("fixed core should parse");
        let mut truncated_seq = fixed_core(0, 0, 1, 0, 0, 0, 0, 3, -1, -1, 0);
        truncated_seq.extend_from_slice(b"\0");
        truncated_seq.extend_from_slice(&[0x12]);
        assert!(BamRecordVariable::new(&fixed, &truncated_seq).is_err());

        let mut truncated_qual = fixed_core(0, 0, 1, 0, 0, 0, 0, 3, -1, -1, 0);
        truncated_qual.extend_from_slice(b"\0");
        truncated_qual.extend_from_slice(&[0x12, 0x30]);
        truncated_qual.extend_from_slice(&[30, 31]);
        assert!(BamRecordVariable::new(&fixed, &truncated_qual).is_err());
    }

    #[test]
    fn malformed_auxiliary_fields_return_errors() {
        let mut offset = 0;
        assert!(BamRecordAuxiliary::new(&mut offset, b"N").is_err());

        let mut offset = 0;
        assert!(BamRecordAuxiliary::new(&mut offset, b"NMq").is_err());

        let mut offset = 0;
        assert!(BamRecordAuxiliary::new(&mut offset, b"RGZunterminated").is_err());

        let mut unknown_array = b"aaBx".to_vec();
        push_u32(&mut unknown_array, 0);
        let mut offset = 0;
        assert!(BamRecordAuxiliary::new(&mut offset, &unknown_array).is_err());

        let mut truncated_array = b"aaBc".to_vec();
        push_u32(&mut truncated_array, 2);
        truncated_array.push(1);
        let mut offset = 0;
        assert!(BamRecordAuxiliary::new(&mut offset, &truncated_array).is_err());
    }

    #[test]
    fn aligned_header_and_reference_match_samtools_metadata() {
        let reader = open_bam(ALIGNED_BAM);
        let lines: Vec<&str> = reader.header.text_str().unwrap().lines().collect();

        assert_eq!(reader.header.magic, [0x42, 0x41, 0x4D, 0x01]);
        assert_eq!(reader.header.l_text, 568);
        assert_eq!(reader.header.n_ref, 1);
        assert_eq!(lines.len(), 6);
        assert_eq!(lines[0], "@HD\tVN:1.6\tSO:unsorted\tGO:query");
        assert_eq!(lines[1], "@SQ\tSN:fc_reference\tLN:687");
        assert!(lines[2].starts_with("@PG\tID:minimap2\tPN:minimap2\tVN:2.30-r1287"));
        assert_eq!(
            lines[5],
            "@PG\tID:samtools.2\tPN:samtools\tPP:samtools.1\tVN:1.19.2\tCL:samtools view -b -o aligned_test.bam -"
        );

        assert_eq!(reader.refs.len(), 1);
        assert_eq!(reader.refs[0].l_name, 13);
        assert_eq!(reader.refs[0].name, "fc_reference");
        assert_eq!(reader.refs[0].l_seq, 687);
    }

    #[test]
    fn unaligned_header_matches_samtools_metadata() {
        let reader = open_bam(UNALIGNED_BAM);
        let lines: Vec<&str> = reader.header.text_str().unwrap().lines().collect();

        assert_eq!(reader.header.magic, [0x42, 0x41, 0x4D, 0x01]);
        assert_eq!(reader.header.l_text, 960);
        assert_eq!(reader.header.n_ref, 0);
        assert!(reader.refs.is_empty());
        assert_eq!(lines.len(), 6);
        assert_eq!(lines[0], "@HD\tVN:1.6\tSO:unknown");
        assert!(lines[1].contains("ID:basecaller"));
        assert!(lines[5].starts_with("@RG\tID:169d620f-94f3-4892-a68a-bc35be9e241a"));
        assert!(lines[5].contains("BC:TCCATTCCCTCCGATAGATGAAAC"));
    }

    #[test]
    fn aligned_records_match_samtools_summary() {
        let records = read_all_records(ALIGNED_BAM);
        let lengths: Vec<usize> = records
            .iter()
            .map(|record| record.fixed.l_seq as usize)
            .collect();

        assert_eq!(records.len(), 100);
        assert_eq!(lengths.iter().sum::<usize>(), 68125);
        assert_eq!(lengths.iter().min(), Some(&615));
        assert_eq!(lengths.iter().max(), Some(&740));
        assert!(records.iter().all(|record| record.fixed.ref_id == 0));
        assert!(records.iter().all(|record| record.fixed.pos >= 0));
        assert!(records.iter().all(|record| record.fixed.mapq == 60));
        assert!(records.iter().all(|record| record.fixed.flag & 0x4 == 0));
        assert!(records.iter().all(|record| record.fixed.next_ref_id == -1));
        assert!(records.iter().all(|record| record.fixed.next_pos == -1));
        assert!(records.iter().all(|record| record.fixed.tlen == 0));
        assert!(
            records
                .iter()
                .all(|record| !record.variable.cigar.is_empty())
        );
        assert!(
            records
                .iter()
                .all(|record| record.variable.read_name.len() == 36)
        );
        assert!(records.iter().all(|record| record.auxiliary.len() == 10));

        let counts = tag_counts(&records);
        for tag in ["NM", "ms", "AS", "nn", "tp", "cm", "s1", "s2", "de", "rl"] {
            assert_eq!(
                counts.get(tag),
                Some(&100),
                "unexpected count for tag {tag}"
            );
        }
    }

    #[test]
    fn first_aligned_record_decodes_like_samtools_view() {
        let mut bam = open_bam(ALIGNED_BAM);
        let record = bam
            .read()
            .expect("read should succeed")
            .expect("record exists");

        assert_eq!(record.fixed.block_size, 1203);
        assert_eq!(record.fixed.ref_id, 0);
        assert_eq!(record.fixed.pos, 1);
        assert_eq!(record.fixed.l_read_name, 37);
        assert_eq!(record.fixed.mapq, 60);
        assert_eq!(record.fixed.n_cigar_op, 5);
        assert_eq!(record.fixed.flag, 16);
        assert_eq!(record.fixed.l_seq, 712);
        assert_eq!(record.read_name(), "2f890ab0-63b9-45f6-ba8d-3c367ca26a63");
        assert_eq!(record.cigar_string(), "7S333M1I334M37S");

        let seq = record.sequence_string();
        assert_eq!(seq.len(), 712);
        assert!(seq.starts_with("TCGCCAGCGCCATCCTGCTGTCATCCGCGTCTGTCCCTGCACCGTCCGGCGCTGGAGGAT"));
        assert!(seq.ends_with("GTAAGCTAAGATCAGGTGTC"));

        let qual = record.quality_string();
        assert_eq!(qual.len(), 712);
        assert!(qual.starts_with("####$$#)<?DOSMOSLQSO"));
        assert!(qual.ends_with("31/*)*))(''''$$%$$$#"));

        assert_aux_int(aux_value(&record, "NM"), 8);
        assert_aux_int(aux_value(&record, "ms"), 1286);
        assert_aux_int(aux_value(&record, "AS"), 1286);
        assert_aux_int(aux_value(&record, "nn"), 0);
        assert_eq!(*aux_value(&record, "tp"), BamAuxValue::A(b'P'));
        assert_aux_int(aux_value(&record, "cm"), 103);
        assert_aux_int(aux_value(&record, "s1"), 614);
        assert_aux_int(aux_value(&record, "s2"), 0);
        assert_aux_f32(aux_value(&record, "de"), 0.012);
        assert_aux_int(aux_value(&record, "rl"), 0);
    }

    #[test]
    fn unaligned_records_match_samtools_summary() {
        let records = read_all_records(UNALIGNED_BAM);
        let lengths: Vec<usize> = records
            .iter()
            .map(|record| record.fixed.l_seq as usize)
            .collect();

        assert_eq!(records.len(), 100);
        assert_eq!(lengths.iter().sum::<usize>(), 158745);
        assert_eq!(lengths.iter().min(), Some(&259));
        assert_eq!(lengths.iter().max(), Some(&5895));
        assert!(records.iter().all(|record| record.fixed.ref_id == -1));
        assert!(records.iter().all(|record| record.fixed.pos == -1));
        assert!(records.iter().all(|record| record.fixed.mapq == 0));
        assert!(records.iter().all(|record| record.fixed.flag & 0x4 != 0));
        assert!(records.iter().all(|record| record.fixed.n_cigar_op == 0));
        assert!(
            records
                .iter()
                .all(|record| record.variable.cigar.is_empty())
        );
        assert!(records.iter().all(|record| record.fixed.next_ref_id == -1));
        assert!(records.iter().all(|record| record.fixed.next_pos == -1));
        assert!(records.iter().all(|record| record.fixed.tlen == 0));
        assert!(
            records
                .iter()
                .all(|record| record.variable.read_name.len() == 36)
        );

        let counts = tag_counts(&records);
        for tag in [
            "BC", "bv", "bi", "qs", "du", "ns", "ts", "mx", "ch", "st", "rn", "fn", "sm", "sd",
            "sv", "dx", "RG", "po", "er", "me",
        ] {
            assert_eq!(
                counts.get(tag),
                Some(&100),
                "unexpected count for tag {tag}"
            );
        }
        assert_eq!(counts.get("pi"), Some(&5));
        assert_eq!(counts.get("sp"), Some(&5));
    }

    #[test]
    fn first_unaligned_record_decodes_like_samtools_view() {
        let mut bam = open_bam(UNALIGNED_BAM);
        let record = bam
            .read()
            .expect("read should succeed")
            .expect("record exists");

        assert_eq!(record.fixed.ref_id, -1);
        assert_eq!(record.fixed.pos, -1);
        assert_eq!(record.fixed.l_read_name, 37);
        assert_eq!(record.fixed.mapq, 0);
        assert_eq!(record.fixed.n_cigar_op, 0);
        assert_eq!(record.fixed.flag, 4);
        assert_eq!(record.fixed.l_seq, 798);
        assert_eq!(record.read_name(), "d3843747-8d64-429e-bc47-44763b006ad1");
        assert_eq!(record.cigar_string(), "*");

        let seq = record.sequence_string();
        assert_eq!(seq.len(), 798);
        assert!(seq.starts_with("AGCGTCAGATGTGTATAAGAGACAGTCTCGATCTTAGCTTACCATGCATGCATCAAAATC"));
        assert!(seq.ends_with("ATCTCCGAGCCCACGAGACA"));

        let qual = record.quality_string();
        assert_eq!(qual.len(), 798);
        assert!(qual.starts_with(";;<CAEIHHKMGIGGBA???"));
        assert!(qual.ends_with("99IIEG88EC??=FGEEFGH"));

        assert_eq!(
            *aux_value(&record, "BC"),
            BamAuxValue::Z("SQK-NBD114-24_barcode11".to_string())
        );
        assert_eq!(
            *aux_value(&record, "bv"),
            BamAuxValue::Z("var1".to_string())
        );
        assert_aux_f32_array(
            aux_value(&record, "bi"),
            &[1.0, 27.0, 45.0, 1.0, 914.0, 43.0, 1.0],
        );
        assert_aux_f32(aux_value(&record, "qs"), 21.506);
        assert_aux_f32(aux_value(&record, "du"), 2.1404);
        assert_aux_int(aux_value(&record, "ns"), 10702);
        assert_aux_int(aux_value(&record, "ts"), 1048);
        assert_aux_int(aux_value(&record, "mx"), 3);
        assert_aux_int(aux_value(&record, "ch"), 457);
        assert_eq!(
            *aux_value(&record, "st"),
            BamAuxValue::Z("2026-05-14T08:00:13.265000+00:00".to_string())
        );
        assert_aux_int(aux_value(&record, "rn"), 2355);
        assert_eq!(
            *aux_value(&record, "fn"),
            BamAuxValue::Z("14_05_2026.pod5".to_string())
        );
        assert_aux_f32(aux_value(&record, "sm"), 93.6924);
        assert_aux_f32(aux_value(&record, "sd"), 23.5067);
        assert_eq!(*aux_value(&record, "sv"), BamAuxValue::Z("pa".to_string()));
        assert_aux_int(aux_value(&record, "dx"), 0);
        assert_eq!(
            *aux_value(&record, "RG"),
            BamAuxValue::Z("169d620f-94f3-4892-a68a-bc35be9e241a_dna_r10.4.1_e8.2_400bps_sup@v5.2.0_SQK-NBD114-24_barcode11".to_string())
        );
        assert_eq!(
            *aux_value(&record, "po"),
            BamAuxValue::Z("not_set".to_string())
        );
        assert_eq!(
            *aux_value(&record, "er"),
            BamAuxValue::Z("signal_positive".to_string())
        );
        assert_aux_int(aux_value(&record, "me"), 1433);
    }
}
