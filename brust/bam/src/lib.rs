//! Minimal BAM parsing primitives.
//!
//! The crate exposes two ownership models:
//!
//! - [`BamReader`] is a streaming parser over a compressed BAM byte stream. It
//!   owns decoder state and is intentionally not [`Clone`].
//! - [`Bam`] is a fully materialized BAM payload containing the header,
//!   references, and records. It is [`Clone`] for deep-copy workflows.
//!
//! Positions in record fields follow the BAM binary format: alignment positions
//! are zero-based and unmapped positions are represented as `-1`.

use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufReader, Read};
use std::path::Path;

/// A fully materialized BAM file.
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

/// Streaming BAM parser over a BGZF/gzip-compressed reader.
///
/// `BamReader` owns an active decoder and file/reader cursor, so it is not
/// cloneable. Use [`BamReader::read_all`] or [`Bam::from_path`] when a deep
/// copyable, materialized representation is needed.
pub struct BamReader<R: Read = File> {
    /// BAM header metadata and SAM header text.
    pub header: BamHeader,
    /// Reference sequence dictionary.
    pub refs: Vec<BamRef>,
    reader: BufReader<MultiGzDecoder<R>>,
}

/// BAM header information.
///
/// Contains the BAM magic bytes, SAM header text length, SAM header text, and
/// reference count.
#[derive(Debug, Clone, PartialEq)]
pub struct BamHeader {
    /// BAM magic bytes, always `[0x42, 0x41, 0x4D, 0x01]`.
    pub magic: [u8; 4],
    /// Length of the SAM header text in bytes, including any NUL padding.
    pub l_text: u32,
    /// Plain SAM header text; not necessarily NUL-terminated in the file.
    pub text: String,
    /// Number of reference sequences in the BAM reference dictionary.
    pub n_ref: u32,
}

/// Reference sequence information.
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
    /// Auxiliary tags following the core record payload.
    pub auxiliary: Vec<BamRecordAuxiliary>,
}

/// Fixed-length fields of a BAM alignment record.
#[derive(Debug, Clone, PartialEq)]
pub struct BamRecordFixed {
    /// Total alignment record length, excluding the `block_size` field itself.
    pub block_size: u32,
    /// Reference sequence ID, or `-1` for a read without a mapping position.
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
    /// Zero-based leftmost coordinate of the next segment.
    pub next_pos: i32,
    /// Template length (`TLEN`).
    pub tlen: i32,
}

impl BamRecordFixed {
    /// Parses the fixed 32-byte BAM record core.
    ///
    /// `block_size` is the record length read from the stream, excluding the
    /// four-byte `block_size` field itself. `data` begins at `refID`.
    pub fn new(block_size: u32, data: &[u8]) -> io::Result<Self> {
        if data.len() < 32 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "data too short"));
        }

        let ref_id = match data.get(0..4) {
            Some(d) => i32::from_le_bytes(d.try_into().unwrap()),
            None => return Err(io::Error::new(io::ErrorKind::InvalidData, "data too short")),
        };
        let pos = match data.get(4..8) {
            Some(d) => i32::from_le_bytes(d.try_into().unwrap()),
            None => return Err(io::Error::new(io::ErrorKind::InvalidData, "data too short")),
        };
        let l_read_name = data[8];
        let mapq = data[9];
        let bin = match data.get(10..12) {
            Some(d) => u16::from_le_bytes(d.try_into().unwrap()),
            None => return Err(io::Error::new(io::ErrorKind::InvalidData, "data too short")),
        };
        let n_cigar_op = match data.get(12..14) {
            Some(d) => u16::from_le_bytes(d.try_into().unwrap()),
            None => return Err(io::Error::new(io::ErrorKind::InvalidData, "data too short")),
        };
        let flag = match data.get(14..16) {
            Some(d) => u16::from_le_bytes(d.try_into().unwrap()),
            None => return Err(io::Error::new(io::ErrorKind::InvalidData, "data too short")),
        };
        let l_seq = match data.get(16..20) {
            Some(d) => u32::from_le_bytes(d.try_into().unwrap()),
            None => return Err(io::Error::new(io::ErrorKind::InvalidData, "data too short")),
        };
        let next_ref_id = match data.get(20..24) {
            Some(d) => i32::from_le_bytes(d.try_into().unwrap()),
            None => return Err(io::Error::new(io::ErrorKind::InvalidData, "data too short")),
        };
        let next_pos = match data.get(24..28) {
            Some(d) => i32::from_le_bytes(d.try_into().unwrap()),
            None => return Err(io::Error::new(io::ErrorKind::InvalidData, "data too short")),
        };
        let tlen = match data.get(28..32) {
            Some(d) => i32::from_le_bytes(d.try_into().unwrap()),
            None => return Err(io::Error::new(io::ErrorKind::InvalidData, "data too short")),
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
    /// Raw Phred-scaled base qualities.
    pub qual: Vec<u8>,
}

impl BamRecordVariable {
    /// Parses variable-length fields from a BAM record payload.
    ///
    /// `fixed` contains offsets and lengths from the fixed record core. `data`
    /// begins at `refID` and does not include the four-byte `block_size` field.
    pub fn new(fixed: &BamRecordFixed, data: &[u8]) -> io::Result<(Self, usize)> {
        let mut offset = 32;

        let read_name = match data.get(offset..offset + fixed.l_read_name as usize) {
            Some(name) => String::from_utf8_lossy(name)
                .trim_end_matches('\0')
                .to_string(),
            None => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "read_name cannot be extracted",
                ));
            }
        };
        offset += fixed.l_read_name as usize;

        let cigar_len_bytes = (fixed.n_cigar_op as usize).checked_mul(4).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "CIGAR byte length overflow")
        })?;

        let cigar_bytes = match data.get(offset..offset + cigar_len_bytes) {
            Some(slice) => slice,
            None => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CIGAR data out of bounds",
                ));
            }
        };

        let mut cigar = Vec::with_capacity(fixed.n_cigar_op as usize);
        for chunk in cigar_bytes.chunks_exact(4) {
            let mut bytes = [0u8; 4];
            bytes.copy_from_slice(chunk);
            let packed = u32::from_le_bytes(bytes);
            if packed & 0x0f > 8 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CIGAR operation code out of bounds",
                ));
            }
            cigar.push(packed);
        }
        offset += cigar_len_bytes;

        let seq_byte_length = fixed.l_seq.div_ceil(2) as usize;
        let seq = match data.get(offset..offset + seq_byte_length) {
            Some(slice) => slice.to_vec(),
            None => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Sequence data out of bounds",
                ));
            }
        };
        offset += seq_byte_length;

        let qual = match data.get(offset..offset + fixed.l_seq as usize) {
            Some(slice) => slice.to_vec(),
            None => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Quality data out of bounds",
                ));
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
    pub fn new(offset: &mut usize, data: &[u8]) -> io::Result<Self> {
        let tag = match data.get(*offset..*offset + 2) {
            Some(tag) => String::from_utf8_lossy(tag).to_string(),
            None => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "auxiliary tag is incomplete",
                ));
            }
        };
        *offset += 2;

        let val_type = match data.get(*offset) {
            Some(val_type) => *val_type,
            None => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "auxiliary value type is missing",
                ));
            }
        };
        *offset += 1;

        let value = match val_type {
            b'A' => {
                let v = match data.get(*offset) {
                    Some(v) => *v,
                    None => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "auxiliary A value is missing",
                        ));
                    }
                };
                *offset += 1;
                BamAuxValue::A(v)
            }

            b'c' => {
                let v = match data.get(*offset) {
                    Some(v) => *v as i8,
                    None => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "auxiliary c value is missing",
                        ));
                    }
                };
                *offset += 1;
                BamAuxValue::c(v)
            }

            b'C' => {
                let v = match data.get(*offset) {
                    Some(v) => *v,
                    None => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "auxiliary C value is missing",
                        ));
                    }
                };
                *offset += 1;
                BamAuxValue::C(v)
            }

            b's' => {
                let v = match data.get(*offset..*offset + 2) {
                    Some(v) => i16::from_le_bytes(v.try_into().unwrap()),
                    None => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "auxiliary s value is incomplete",
                        ));
                    }
                };
                *offset += 2;
                BamAuxValue::s(v)
            }

            b'S' => {
                let v = match data.get(*offset..*offset + 2) {
                    Some(v) => u16::from_le_bytes(v.try_into().unwrap()),
                    None => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "auxiliary S value is incomplete",
                        ));
                    }
                };
                *offset += 2;
                BamAuxValue::S(v)
            }

            b'i' => {
                let v = match data.get(*offset..*offset + 4) {
                    Some(v) => i32::from_le_bytes(v.try_into().unwrap()),
                    None => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "auxiliary i value is incomplete",
                        ));
                    }
                };
                *offset += 4;
                BamAuxValue::i(v)
            }

            b'I' => {
                let v = match data.get(*offset..*offset + 4) {
                    Some(v) => u32::from_le_bytes(v.try_into().unwrap()),
                    None => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "auxiliary I value is incomplete",
                        ));
                    }
                };
                *offset += 4;
                BamAuxValue::I(v)
            }

            b'f' => {
                let v = match data.get(*offset..*offset + 4) {
                    Some(v) => f32::from_le_bytes(v.try_into().unwrap()),
                    None => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "auxiliary f value is incomplete",
                        ));
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
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "auxiliary Z string is missing NUL terminator",
                    ));
                }

                let s = String::from_utf8_lossy(&data[start..*offset]).to_string();
                *offset += 1; // consume NUL

                BamAuxValue::Z(s)
            }

            b'H' => {
                let start = *offset;

                while *offset < data.len() && data[*offset] != 0 {
                    *offset += 1;
                }

                if *offset >= data.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "auxiliary H string is missing NUL terminator",
                    ));
                }

                let s = String::from_utf8_lossy(&data[start..*offset]).to_string();
                *offset += 1; // consume NUL

                BamAuxValue::H(s)
            }

            b'B' => {
                let subtype = match data.get(*offset) {
                    Some(subtype) => *subtype,
                    None => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "auxiliary B subtype is missing",
                        ));
                    }
                };
                *offset += 1;

                let count = match data.get(*offset..*offset + 4) {
                    Some(v) => u32::from_le_bytes(v.try_into().unwrap()) as usize,
                    None => {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "auxiliary B count is incomplete",
                        ));
                    }
                };
                *offset += 4;

                match subtype {
                    b'c' => {
                        let mut values = Vec::with_capacity(count);

                        for _ in 0..count {
                            let v = match data.get(*offset) {
                                Some(v) => *v as i8,
                                None => {
                                    return Err(io::Error::new(
                                        io::ErrorKind::InvalidData,
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
                                    return Err(io::Error::new(
                                        io::ErrorKind::InvalidData,
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
                                    return Err(io::Error::new(
                                        io::ErrorKind::InvalidData,
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
                                    return Err(io::Error::new(
                                        io::ErrorKind::InvalidData,
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
                                    return Err(io::Error::new(
                                        io::ErrorKind::InvalidData,
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
                                    return Err(io::Error::new(
                                        io::ErrorKind::InvalidData,
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
                                    return Err(io::Error::new(
                                        io::ErrorKind::InvalidData,
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
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("unknown auxiliary B array subtype: {}", subtype as char),
                        ));
                    }
                }
            }

            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("unknown auxiliary value type: {}", val_type as char),
                ));
            }
        };

        Ok(Self { tag, value })
    }
}

/// BAM auxiliary tag value.
#[derive(Debug, Clone, PartialEq)]
pub enum BamAuxValue {
    /// Printable ASCII character (`A`).
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
    /// NUL-terminated printable string (`Z`).
    Z(String),
    /// NUL-terminated hexadecimal string (`H`).
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
    /// contains the full alignment record payload.
    pub fn new(block_size: u32, data: Vec<u8>) -> io::Result<Self> {
        if block_size as usize != data.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "block_size does not match record data length",
            ));
        }

        if data.len() < 32 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
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

        Ok(Self {
            fixed,
            variable,
            auxiliary,
        })
    }

    /// Returns this record's read name.
    pub fn read_name(&self) -> &str {
        &self.variable.read_name
    }

    /// Returns the record's CIGAR as a SAM-style string.
    pub fn cigar_string(&self) -> String {
        const CIGAR_OPS: &[u8; 9] = b"MIDNSHP=X";
        if self.variable.cigar.is_empty() {
            return "*".to_string();
        }

        let mut decoded = String::new();
        for op in &self.variable.cigar {
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
        let l_seq = self.fixed.l_seq as usize;
        let mut decoded = String::with_capacity(l_seq);
        for i in 0..l_seq {
            let byte = self.variable.seq[i / 2];
            let code = if i % 2 == 0 { byte >> 4 } else { byte & 0x0f };
            decoded.push(BASES[code as usize] as char);
        }
        decoded
    }

    /// Decodes BAM quality bytes into a SAM-style quality string.
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

    /// Returns the first auxiliary tag value matching `tag`.
    pub fn aux(&self, tag: &str) -> Option<&BamAuxValue> {
        self.auxiliary
            .iter()
            .find(|auxiliary| auxiliary.tag == tag)
            .map(|auxiliary| &auxiliary.value)
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
    pub fn from_reader(reader: R) -> io::Result<Self> {
        let mut reader = BufReader::new(MultiGzDecoder::new(reader));
        let header = BamHeader::read_from(&mut reader)?;
        let mut refs = Vec::with_capacity(header.n_ref as usize);

        for _ in 0..header.n_ref {
            refs.push(BamRef::read_from(&mut reader)?);
        }

        Ok(Self {
            reader,
            header,
            refs,
        })
    }

    /// Reads the next BAM alignment record.
    ///
    /// Returns `Ok(None)` only when the stream is at EOF before a new record's
    /// `block_size` field begins.
    pub fn read_record(&mut self) -> io::Result<Option<BamRecord>> {
        let mut block_size = [0u8; 4];
        if !read_exact_or_eof(&mut self.reader, &mut block_size)? {
            return Ok(None);
        }

        let block_size = u32::from_le_bytes(block_size);

        let mut record_data = vec![0u8; block_size as usize];
        self.reader.read_exact(&mut record_data)?;
        BamRecord::new(block_size, record_data).map(Some)
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
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid BAM magic bytes",
            ));
        }

        let l_text = read_u32_le(reader, "BAM header text length")?;
        let mut text = vec![0u8; l_text as usize];
        reader.read_exact(&mut text)?;
        let text = String::from_utf8_lossy(&text)
            .trim_end_matches('\0')
            .to_string();

        let n_ref = read_u32_le(reader, "BAM reference count")?;

        Ok(Self {
            magic,
            l_text,
            text,
            n_ref,
        })
    }
}

impl BamRef {
    /// Reads one reference dictionary entry from a decompressed BAM stream.
    pub fn read_from<R: Read>(reader: &mut R) -> io::Result<Self> {
        let l_name = read_u32_le(reader, "reference name length")?;
        let mut name = vec![0u8; l_name as usize];
        reader.read_exact(&mut name)?;
        let name = String::from_utf8_lossy(&name)
            .trim_end_matches('\0')
            .to_string();

        let l_seq = read_u32_le(reader, "reference sequence length")?;

        Ok(Self {
            l_name,
            name,
            l_seq,
        })
    }
}

fn read_u32_le<R: Read>(reader: &mut R, field: &str) -> io::Result<u32> {
    let mut bytes = [0u8; 4];
    reader
        .read_exact(&mut bytes)
        .map_err(|err| io::Error::new(err.kind(), format!("failed to read {field}: {err}")))?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_exact_or_eof<R: Read>(reader: &mut R, buf: &mut [u8]) -> io::Result<bool> {
    let mut bytes_read = 0;
    while bytes_read < buf.len() {
        match reader.read(&mut buf[bytes_read..]) {
            Ok(0) if bytes_read == 0 => return Ok(false),
            Ok(0) => {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "unexpected EOF while reading BAM record block_size",
                ));
            }
            Ok(n) => bytes_read += n,
            Err(err) if err.kind() == io::ErrorKind::Interrupted => {}
            Err(err) => return Err(err),
        }
    }

    Ok(true)
}

#[cfg(test)]
mod unit_tests {
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

        cloned.header.text.push_str("@CO\tclone-only\n");
        cloned.refs[0].name.push_str("_clone");
        cloned.records[0].variable.read_name.push_str("_clone");

        assert_ne!(cloned, original);
        assert_eq!(original.refs[0].name, "fc_reference");
        assert_eq!(
            original.records[0].read_name(),
            "2f890ab0-63b9-45f6-ba8d-3c367ca26a63"
        );
        assert!(!original.header.text.contains("clone-only"));
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
        let lines: Vec<&str> = reader.header.text.lines().collect();

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
        let lines: Vec<&str> = reader.header.text.lines().collect();

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
