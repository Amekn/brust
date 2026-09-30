/// Translate one codon with the standard genetic code (NCBI table 1).
///
/// Input is case-insensitive, and `U` counts as `T`. Returns the one-letter amino
/// acid, or `*` for TAA, TAG and TGA. Returns `None` unless given exactly 3 bases
/// that are all A, C, G or T (or U).
///
/// ```
/// assert_eq!(brust_seq::translate_codon(b"atg"), Some(b'M'));
/// assert_eq!(brust_seq::translate_codon(b"UAA"), Some(b'*'));
/// assert_eq!(brust_seq::translate_codon(b"NNK"), None);
/// ```
pub fn translate_codon(codon: &[u8]) -> Option<u8> {
    let &[first, second, third] = codon else {
        return None;
    };
    let upper = [first, second, third].map(|base| match base.to_ascii_uppercase() {
        b'U' => b'T',
        other => other,
    });
    Some(match &upper {
        b"TTT" | b"TTC" => b'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
        b"ATT" | b"ATC" | b"ATA" => b'I',
        b"ATG" => b'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
        b"TAT" | b"TAC" => b'Y',
        b"TAA" | b"TAG" | b"TGA" => b'*',
        b"CAT" | b"CAC" => b'H',
        b"CAA" | b"CAG" => b'Q',
        b"AAT" | b"AAC" => b'N',
        b"AAA" | b"AAG" => b'K',
        b"GAT" | b"GAC" => b'D',
        b"GAA" | b"GAG" => b'E',
        b"TGT" | b"TGC" => b'C',
        b"TGG" => b'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
        _ => return None,
    })
}

/// Translate a sequence in reading frame 0 with the standard genetic code.
///
/// Each codon becomes [`translate_codon`], or `X` when that returns `None`. A
/// trailing 1 or 2 bases are ignored, and translation continues past `*`.
///
/// ```
/// assert_eq!(brust_seq::translate(b"ATGGCCTAAGGN"), b"MA*X");
/// assert_eq!(brust_seq::translate(b"ATGGC"), b"M");
/// ```
pub fn translate(seq: &[u8]) -> Vec<u8> {
    seq.chunks_exact(3)
        .map(|codon| translate_codon(codon).unwrap_or(b'X'))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_code_all_64_codons_and_20_amino_acid_identities() {
        // Synthetic test codons; independent NCBI Standard Code (transl_table=1),
        // verified 2026-09-05: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
        let expected = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
        let mut found = std::collections::BTreeSet::new();
        let mut index = 0;
        let mut stops = Vec::new();
        for a in b"TCAG" {
            for b in b"TCAG" {
                for c in b"TCAG" {
                    let codon = [*a, *b, *c];
                    assert_eq!(
                        translate_codon(&codon),
                        Some(expected[index]),
                        "{}",
                        String::from_utf8_lossy(&codon)
                    );
                    found.insert(translate_codon(&codon).unwrap());
                    if expected[index] == b'*' {
                        stops.push(codon);
                    }
                    index += 1;
                }
            }
        }
        assert_eq!(index, 64);
        assert_eq!(stops, [*b"TAA", *b"TAG", *b"TGA"]);
        found.remove(&b'*');
        assert_eq!(
            found.into_iter().collect::<Vec<_>>(),
            b"ACDEFGHIKLMNPQRSTVWY"
        );
        // All 20 identities independently audited against IUPAC-IUB Table 1:
        // https://iupac.qmul.ac.uk/AminoAcid/AA1n2.html
        for (name, symbol, codon) in [
            ("Alanine", b'A', b"GCT"),
            ("Cysteine", b'C', b"TGT"),
            ("Aspartic acid", b'D', b"GAT"),
            ("Glutamic acid", b'E', b"GAA"),
            ("Phenylalanine", b'F', b"TTT"),
            ("Glycine", b'G', b"GGT"),
            ("Histidine", b'H', b"CAT"),
            ("Isoleucine", b'I', b"ATT"),
            ("Lysine", b'K', b"AAA"),
            ("Leucine", b'L', b"TTA"),
            ("Methionine", b'M', b"ATG"),
            ("Asparagine", b'N', b"AAT"),
            ("Proline", b'P', b"CCT"),
            ("Glutamine", b'Q', b"CAA"),
            ("Arginine", b'R', b"CGT"),
            ("Serine", b'S', b"TCT"),
            ("Threonine", b'T', b"ACT"),
            ("Valine", b'V', b"GTT"),
            ("Tryptophan", b'W', b"TGG"),
            ("Tyrosine", b'Y', b"TAT"),
        ] {
            assert_eq!(translate_codon(codon), Some(symbol), "{name}");
        }
    }

    #[test]
    fn translate_codon_accepts_lowercase_and_u() {
        assert_eq!(translate_codon(b"atg"), Some(b'M'));
        assert_eq!(translate_codon(b"AUG"), Some(b'M'));
        assert_eq!(translate_codon(b"uaa"), Some(b'*'));
    }

    #[test]
    fn translate_codon_does_not_resolve_degenerate_codons() {
        // CTN is always leucine, but there is no degenerate-codon resolution.
        assert_eq!(translate_codon(b"CTN"), None);
    }

    #[test]
    fn translate_codon_rejects_ambiguous_and_wrong_length() {
        assert_eq!(translate_codon(b"NNK"), None);
        assert_eq!(translate_codon(b"AT"), None);
        assert_eq!(translate_codon(b"ATGA"), None);
        assert_eq!(translate_codon(b""), None);
    }

    #[test]
    fn translate_reads_frame_zero() {
        assert_eq!(translate(b"ATGGCCTAAGGN"), b"MA*X");
        assert_eq!(translate(b"ATGGC"), b"M");
        assert_eq!(translate(b""), b"");
    }
}
