//! Sampling utilities inspired by `falcon.py`.
//!
//! These helpers are **not** constant-time and are meant for experimentation
//! and interoperability testing, not production use.

use core::f64::consts::PI;

/// Minimal RNG abstraction; when the `rand` feature is enabled, a `rand_core`
/// adapter is provided, otherwise a deterministic xorshift is available for
/// tests and examples.
pub trait RandomSource {
    fn next_u64(&mut self) -> u64;
}

/// A tiny xorshift PRNG; convenient for deterministic tests.
#[derive(Debug, Clone)]
pub struct XorShift64 {
    state: u64,
}

impl XorShift64 {
    pub fn new(seed: u64) -> Self {
        let seed = if seed == 0 { 1 } else { seed };
        Self { state: seed }
    }
}

impl RandomSource for XorShift64 {
    fn next_u64(&mut self) -> u64 {
        // xorshift64* with constants from Vigna.
        let mut x = self.state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.state = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
}

/// Adapter over `rand_core` RNGs when the `rand` feature is enabled.
#[cfg(feature = "rand")]
#[derive(Debug)]
pub struct RandCoreSource<R: rand_core::RngCore>(pub R);

#[cfg(feature = "rand")]
impl<R: rand_core::RngCore> RandomSource for RandCoreSource<R> {
    fn next_u64(&mut self) -> u64 {
        self.0.next_u64()
    }
}

#[inline]
fn uniform01<R: RandomSource>(rng: &mut R) -> f64 {
    // Take the top 53 bits to fill an f64 mantissa.
    let bits = rng.next_u64() >> 11;
    bits as f64 / (1u64 << 53) as f64
}

/// Box-Muller transform producing a standard normal variate.
pub fn gaussian_box_muller<R: RandomSource>(rng: &mut R) -> f64 {
    let u1 = uniform01(rng).max(f64::MIN_POSITIVE);
    let u2 = uniform01(rng);
    let r = (-2.0 * u1.ln()).sqrt();
    let theta = 2.0 * PI * u2;
    r * theta.cos()
}

/// Sample an integer from a discrete Gaussian centered at `center` with
/// standard deviation `sigma` (rounded Box-Muller).
pub fn discrete_gaussian<R: RandomSource>(
    rng: &mut R,
    sigma: f64,
    center: f64,
) -> i64 {
    let z = gaussian_box_muller(rng) * sigma + center;
    z.round() as i64
}

/// Bernoulli trial with success probability `p` (0 <= p <= 1).
pub fn bernoulli<R: RandomSource>(rng: &mut R, p: f64) -> bool {
    assert!((0.0..=1.0).contains(&p));
    uniform01(rng) < p
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_stream() {
        let mut rng = XorShift64::new(42);
        assert_eq!(rng.next_u64(), 6_255_019_084_209_693_600);
        assert_eq!(rng.next_u64(), 14_430_073_426_741_505_498);
    }

    #[test]
    fn gaussian_is_centered() {
        let mut rng = XorShift64::new(1234);
        let mut acc = 0.0;
        let samples = 2000;
        for _ in 0..samples {
            acc += discrete_gaussian(&mut rng, 2.5, 0.0) as f64;
        }
        let mean = acc / samples as f64;
        assert!(mean.abs() < 0.3, "mean drifted too far: {mean}");
    }

    #[test]
    fn gaussian_with_offset() {
        let mut rng = XorShift64::new(999);
        let mut acc = 0.0;
        let samples = 2000;
        for _ in 0..samples {
            acc += discrete_gaussian(&mut rng, 1.5, 5.0) as f64;
        }
        let mean = acc / samples as f64;
        assert!((mean - 5.0).abs() < 0.3, "mean drifted too far: {mean}");
    }

    #[test]
    fn bernoulli_frequency_matches_probability() {
        let mut rng = XorShift64::new(2024);
        let p = 0.3;
        let trials = 5000;
        let mut successes = 0;
        for _ in 0..trials {
            if bernoulli(&mut rng, p) {
                successes += 1;
            }
        }
        let freq = successes as f64 / trials as f64;
        assert!((freq - p).abs() < 0.03, "frequency {freq} deviates too far");
    }
}
