/// Complement of one nucleotide byte, keeping case.
///
/// Every IUPAC code is handled: A and T swap, C and G swap, R and Y swap, K and M
/// swap, B and V swap, D and H swap, and S, W and N map to themselves. `U` maps to
/// `A`. Lowercase input gives lowercase output. The gap symbols `-` and `.` are
/// returned unchanged, and any other byte becomes `N`.
///
/// ```
/// assert_eq!(brust_seq::complement(b'A'), b'T');
/// assert_eq!(brust_seq::complement(b'r'), b'y');
/// assert_eq!(brust_seq::complement(b'-'), b'-');
/// assert_eq!(brust_seq::complement(b'*'), b'N');
/// ```
pub fn complement(base: u8) -> u8 {
    let upper = match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' | b'U' => b'A',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'S' => b'S',
        b'W' => b'W',
        b'N' => b'N',
        b'-' | b'.' => return base,
        _ => return b'N',
    };
    if base.is_ascii_lowercase() {
        upper.to_ascii_lowercase()
    } else {
        upper
    }
}

/// Reverse complement of a nucleotide sequence, keeping case.
///
/// The sequence is reversed and [`complement`] is applied to each byte, so
/// `U` becomes `A` and unknown bytes become `N`. Callers that want uppercase output
/// upper-case the result.
///
/// ```
/// assert_eq!(brust_seq::reverse_complement(b"AACG"), b"CGTT");
/// assert_eq!(brust_seq::reverse_complement(b"acgtRykm"), b"kmrYacgt");
/// ```
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&base| complement(base)).collect()
}

/// Whether a concrete base fits an IUPAC nucleotide code.
///
/// Both arguments are case-insensitive, and `U` counts as `T` in either. `base` must
/// be a concrete A, C, G or T: an ambiguous or unknown `base` never matches. `N` as
/// the `code` matches any concrete base, and a `code` that is not an IUPAC code
/// matches nothing.
///
/// ```
/// assert!(brust_seq::iupac_matches(b'R', b'g'));
/// assert!(brust_seq::iupac_matches(b'W', b'U'));
/// assert!(!brust_seq::iupac_matches(b'R', b'C'));
/// assert!(!brust_seq::iupac_matches(b'N', b'N'));
/// ```
pub fn iupac_matches(code: u8, base: u8) -> bool {
    let base = normalise(base);
    match normalise(code) {
        b'A' => base == b'A',
        b'C' => base == b'C',
        b'G' => base == b'G',
        b'T' => base == b'T',
        b'R' => matches!(base, b'A' | b'G'),
        b'Y' => matches!(base, b'C' | b'T'),
        b'S' => matches!(base, b'C' | b'G'),
        b'W' => matches!(base, b'A' | b'T'),
        b'K' => matches!(base, b'G' | b'T'),
        b'M' => matches!(base, b'A' | b'C'),
        b'B' => matches!(base, b'C' | b'G' | b'T'),
        b'D' => matches!(base, b'A' | b'G' | b'T'),
        b'H' => matches!(base, b'A' | b'C' | b'T'),
        b'V' => matches!(base, b'A' | b'C' | b'G'),
        b'N' => matches!(base, b'A' | b'C' | b'G' | b'T'),
        _ => false,
    }
}

/// Upper-case a byte and treat `U` as `T`.
fn normalise(byte: u8) -> u8 {
    match byte.to_ascii_uppercase() {
        b'U' => b'T',
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_pairs_every_iupac_code_in_both_cases() {
        let codes = b"ACGTRYKMBVDHSWN";
        let expected = b"TGCAYRMKVBHDSWN";
        for (code, want) in codes.iter().zip(expected) {
            assert_eq!(complement(*code), *want, "{}", *code as char);
            assert_eq!(
                complement(code.to_ascii_lowercase()),
                want.to_ascii_lowercase(),
                "{}",
                code.to_ascii_lowercase() as char
            );
        }
    }

    #[test]
    fn complement_maps_u_to_a_and_keeps_gaps() {
        assert_eq!(complement(b'U'), b'A');
        assert_eq!(complement(b'u'), b'a');
        assert_eq!(complement(b'-'), b'-');
        assert_eq!(complement(b'.'), b'.');
    }

    #[test]
    fn complement_turns_other_bytes_into_n() {
        for byte in [b'X', b'*', b'=', b'1', 0xC3] {
            assert_eq!(complement(byte), b'N', "{byte:#04x}");
        }
    }

    #[test]
    fn complement_turns_lowercase_unknown_bytes_into_uppercase_n() {
        // Unknown bytes become uppercase N, even when the input is lowercase.
        assert_eq!(complement(b'x'), b'N');
    }

    #[test]
    fn reverse_complement_keeps_case_and_reverses() {
        assert_eq!(reverse_complement(b"acgtRykm"), b"kmrYacgt");
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn reverse_complement_twice_is_identity() {
        let seq = b"ACGTRYKMBVDHSWNacgtrykmbvdhswn-.";
        assert_eq!(reverse_complement(&reverse_complement(seq)), seq);
    }

    #[test]
    fn iupac_memberships_match_nc_iub() {
        // NC-IUB Table 1: https://iubmb.qmul.ac.uk/misc/naseq.html#table1
        for (code, members) in [
            (b'A', "A"),
            (b'C', "C"),
            (b'G', "G"),
            (b'T', "T"),
            (b'R', "AG"),
            (b'Y', "CT"),
            (b'S', "CG"),
            (b'W', "AT"),
            (b'K', "GT"),
            (b'M', "AC"),
            (b'B', "CGT"),
            (b'D', "AGT"),
            (b'H', "ACT"),
            (b'V', "ACG"),
            (b'N', "ACGT"),
        ] {
            for base in b"ACGTN" {
                assert_eq!(
                    iupac_matches(code, *base),
                    members.as_bytes().contains(base)
                );
                assert_eq!(
                    iupac_matches(code.to_ascii_lowercase(), base.to_ascii_lowercase()),
                    members.as_bytes().contains(base)
                );
            }
        }
    }

    #[test]
    fn iupac_matches_treats_u_as_t() {
        assert!(iupac_matches(b'U', b'T'));
        assert!(iupac_matches(b'T', b'u'));
        assert!(iupac_matches(b'W', b'U'));
        assert!(!iupac_matches(b'N', b'N'));
        assert!(!iupac_matches(b'X', b'A'));
    }
}
