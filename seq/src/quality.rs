use std::sync::LazyLock;

/// Error probability `10^(-p/10)` for each decoded Phred value, built once so
/// `read_mean_phred` does not call `powf` for every base.
static ERROR: LazyLock<[f64; 256]> = LazyLock::new(|| {
    let mut table = [0.0; 256];
    for (phred, error) in table.iter_mut().enumerate() {
        *error = 10f64.powf(-(phred as f64) / 10.0);
    }
    table
});

/// Phred value of a mean error probability, never below 0.0.
///
/// Written as a comparison rather than `.max(0.0)` so an exact-zero result is
/// `+0.0`, never `-0.0`.
fn phred_of_mean_error(error_sum: f64, count: u64) -> f64 {
    let phred = -10.0 * (error_sum / count as f64).log10();
    if phred <= 0.0 { 0.0 } else { phred }
}

/// Running Phred mean taken in error-probability space.
///
/// Averaging Phred values directly overstates quality, because a few bad bases
/// matter more than many good ones. `PhredMean` averages the error probabilities
/// `10^(-Q/10)` instead and converts the mean back to a Phred value:
/// `-10 * log10(mean(10^(-Q/10)))`. The result is clamped at 0.0, so it is never
/// `-0.0`.
///
/// ```
/// use brust_seq::PhredMean;
///
/// let mut mean = PhredMean::default();
/// assert_eq!(mean.mean(), None);
///
/// mean.add(0.0);
/// mean.add(40.0);
/// // The arithmetic mean would be 20.
/// assert!((mean.mean().unwrap() - 3.0098656839).abs() < 1e-9);
/// ```
#[derive(Debug, Clone, Default)]
pub struct PhredMean {
    count: u64,
    error_sum: f64,
}

impl PhredMean {
    /// Adds one Phred value, such as `30.0` for Q30.
    ///
    /// The value must be finite and non-negative. It is a Phred value, not a
    /// Phred+33 byte; use [`read_mean_phred`] for FASTQ quality strings.
    pub fn add(&mut self, phred: f64) {
        self.count += 1;
        self.error_sum += 10f64.powf(-phred / 10.0);
    }

    /// Phred value of the mean error probability, or `None` when nothing has
    /// been added. Never below 0.0.
    pub fn mean(&self) -> Option<f64> {
        (self.count > 0).then(|| phred_of_mean_error(self.error_sum, self.count))
    }
}

/// Phred value of the mean base error probability of one read, from its
/// Phred+33 quality bytes: `-10 * log10(mean(10^(-Q/10)))`.
///
/// Each byte is decoded with `saturating_sub(33)`, so a byte below 33 counts as
/// Phred 0. Returns `None` for an empty quality string. The result is never
/// below 0.0 and is never `-0.0`.
///
/// This gives the same value as [`PhredMean`] fed the decoded scores, but looks
/// error probabilities up in a table built once, so it stays fast on large
/// FASTQ files.
///
/// ```
/// use brust_seq::read_mean_phred;
///
/// assert_eq!(read_mean_phred(b""), None);
/// assert!((read_mean_phred(b"IIII").unwrap() - 40.0).abs() < 1e-9);
/// // One Q0 base and one Q40 base: the arithmetic mean would be 20.
/// assert!((read_mean_phred(b"!I").unwrap() - 3.0098656839).abs() < 1e-9);
/// ```
pub fn read_mean_phred(phred33: &[u8]) -> Option<f64> {
    if phred33.is_empty() {
        return None;
    }
    let table = &*ERROR;
    let error_sum: f64 = phred33
        .iter()
        .map(|&byte| table[usize::from(byte.saturating_sub(33))])
        .sum();
    Some(phred_of_mean_error(error_sum, phred33.len() as u64))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phred_mean_matches_worked_values() {
        // Worked values ported from NPTune's `PhredMean` test.
        for (scores, expected) in [
            (vec![], None),
            (vec![10.0], Some(10.0)),
            (vec![30.0; 5], Some(30.0)),
            (vec![93.0; 2], Some(93.0)),
            (vec![0.0; 2], Some(0.0)),
            (vec![0.0, 40.0], Some(3.0098656839)),
            (vec![10.0, 30.0], Some(12.9670862188)),
        ] {
            let mut mean = PhredMean::default();
            for &q in &scores {
                mean.add(q);
            }
            match expected {
                None => assert_eq!(mean.mean(), None),
                Some(q) => assert!((mean.mean().unwrap() - q).abs() < 1e-9, "{scores:?}"),
            }
        }
    }

    #[test]
    fn read_mean_phred_equals_phred_mean_for_every_score() {
        for p in 0..=93u8 {
            let mut mean = PhredMean::default();
            mean.add(f64::from(p));
            mean.add(f64::from(93 - p));
            let table = read_mean_phred(&[33 + p, 33 + (93 - p)]).unwrap();
            assert!((table - mean.mean().unwrap()).abs() < 1e-12, "p = {p}");
        }
    }

    #[test]
    fn read_mean_phred_decodes_like_the_stats_code() {
        assert_eq!(read_mean_phred(b""), None);
        // Bytes below 33 decode as Phred 0.
        assert_eq!(read_mean_phred(b"\x20!"), Some(0.0));
        // The clamp must give +0.0, never -0.0.
        let zero = read_mean_phred(b"!!").unwrap();
        assert_eq!(zero, 0.0);
        assert!(zero.is_sign_positive());
    }

    #[test]
    fn phred_mean_never_returns_negative_zero() {
        let mut mean = PhredMean::default();
        mean.add(0.0);
        mean.add(0.0);
        let zero = mean.mean().unwrap();
        assert_eq!(zero, 0.0);
        assert!(zero.is_sign_positive());
    }

    #[test]
    fn a_long_read_keeps_precision() {
        // 1,000,000 bases of Q20 ('5').
        let read = vec![b'5'; 1_000_000];
        assert!((read_mean_phred(&read).unwrap() - 20.0).abs() < 1e-9);
    }
}
