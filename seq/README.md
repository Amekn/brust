# brust-seq

`brust-seq` provides sequence primitives for the Brust workspace: IUPAC
complement and reverse complement, IUPAC base matching, codon translation with
the standard genetic code, and Phred means taken over error probabilities. It
has no dependencies.

Most multi-format applications should depend on `brust` and use `brust::seq`.
Use `brust-seq` directly when you only need these helpers.

## Installation

```bash
cargo add brust-seq
```

```rust
use brust_seq::reverse_complement;
```

## API

Every function works on bytes, and every byte it returns is ASCII, so a result
converts to a `String` without loss.

- `complement(base)`: IUPAC complement of one byte, keeping case. `U` maps to
  `A`, the gap symbols `-` and `.` are unchanged, and any other byte becomes `N`.
- `reverse_complement(seq)`: reversed sequence with `complement` applied to each
  byte. Case is kept, so upper-case the result if you need uppercase.
- `iupac_matches(code, base)`: whether a concrete base (A, C, G or T) fits an
  IUPAC code. Case-insensitive, and `U` counts as `T`.
- `translate_codon(codon)`: one codon to its amino acid (NCBI table 1), `*` for a
  stop codon, or `None` unless given exactly 3 bases that are all A, C, G or T.
- `translate(seq)`: reading frame 0, with `X` for codons that cannot be
  translated. A trailing 1 or 2 bases are ignored, and translation continues
  past `*`.
- `PhredMean`: running Phred mean of finished Phred values (not Phred+33
  bytes), taken over error probabilities: `-10 * log10(mean(10^(-Q/10)))`. Use
  `add(phred)` and `mean()`, which is `None` until a value is added. The result
  is never below 0.0.
- `read_mean_phred(phred33)`: the same mean for one read's Phred+33 quality
  bytes, or `None` for an empty read. Bytes below 33 count as Phred 0.

## Example

```rust
use brust_seq::{iupac_matches, read_mean_phred, reverse_complement, translate};

fn main() {
    assert_eq!(reverse_complement(b"AACG"), b"CGTT");
    assert!(iupac_matches(b'R', b'g'));
    assert_eq!(translate(b"ATGGCCTAAGGN"), b"MA*X");
    // One Q0 base and one Q40 base: the arithmetic mean is 20.
    let qscore = read_mean_phred(b"!I").unwrap();
    assert!((qscore - 3.0098656839).abs() < 1e-9);
}
```
