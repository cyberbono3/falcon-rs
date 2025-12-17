//! Falcon key generation scaffolding.
//!
//! This mirrors the reference flow at a high level: sample short polynomials
//! `f` and `g` from a discrete Gaussian, derive a public representative `h`,
//! and expose deterministic vs. RNG-backed paths. The exact Falcon arithmetic
//! (e.g., `h = g / f mod q` with proper norm checks) is not yet implemented;
//! this module focuses on a clean API and reproducible sampling.

use crate::math::field::modq_from_i64;
use crate::math::params::{MAX_LOGN, MIN_LOGN, ParameterSet, Parameters};
use crate::math::sampling::{RandomSource, XorShift64, discrete_gaussian};
use crate::math::{ModQ, Poly};
use crate::poly;

/// Placeholder public key container.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PublicKey {
    // TODO add field element representation
    pub h: Poly<ModQ>,
    pub params: Parameters,
}

impl PublicKey {
    pub fn new(h: Poly<ModQ>, params: Parameters) -> Self {
        Self { h, params }
    }

    /// Polynomial degree (`n`).
    pub fn degree(&self) -> usize {
        self.params.degree
    }

    /// `log2(n)`.
    pub fn logn(&self) -> u8 {
        self.params.log_degree
    }

    /// Access the public polynomial coefficients.
    pub fn coeffs(&self) -> &[crate::math::ModQ] {
        self.h.coeffs()
    }
}

/// Placeholder secret key container.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SecretKey {
    // TODO make all fields private
    pub f: Vec<i64>,
    pub g: Vec<i64>,
    pub params: Parameters,
}

impl SecretKey {
    pub(crate) fn new(f: Vec<i64>, g: Vec<i64>, params: Parameters) -> Self {
        Self { f, g, params }
    }

    /// Polynomial degree (`n`).
    pub fn degree(&self) -> usize {
        self.params.degree
    }

    /// `log2(n)`.
    pub fn logn(&self) -> u8 {
        self.params.log_degree
    }

    /// Accessors for the short secret polynomials.
    pub fn f(&self) -> &[i64] {
        &self.f
    }

    pub fn g(&self) -> &[i64] {
        &self.g
    }
}

/// Keypair tying public and secret parts together.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Keypair {
    pub public: PublicKey,
    pub secret: SecretKey,
}

/// Errors that can occur during key generation.
#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum KeygenError {
    #[error("unsupported logn {0}, expected between {MIN_LOGN} and {MAX_LOGN}")]
    UnsupportedLogN(u8),
    #[error("seed must be non-empty for deterministic key generation")]
    EmptySeed,
}

/// Deterministic key generation from a seed (useful for tests/examples).
pub fn keypair_from_seed(
    logn: u8,
    seed: &[u8],
) -> Result<Keypair, KeygenError> {
    let seed64 = derive_seed(seed)?;
    let mut rng = XorShift64::new(seed64);
    keypair_with_rng(logn, &mut rng)
}

/// Randomized key generation using any RNG that implements [`RandomSource`].
pub fn keypair_with_rng<R: RandomSource>(
    logn: u8,
    rng: &mut R,
) -> Result<Keypair, KeygenError> {
    let params = params_for_logn(logn)?;
    Ok(generate_keypair(params, rng))
}

fn ensure_supported_logn(logn: u8) -> Result<(), KeygenError> {
    if !(MIN_LOGN..=MAX_LOGN).contains(&logn) {
        return Err(KeygenError::UnsupportedLogN(logn));
    }
    Ok(())
}

fn derive_seed(seed: &[u8]) -> Result<u64, KeygenError> {
    if seed.is_empty() {
        return Err(KeygenError::EmptySeed);
    }
    let mut acc = 0u64;
    for (i, &b) in seed.iter().enumerate() {
        let shift = (i % 8) * 8;
        acc ^= (b as u64) << shift;
        acc = acc.rotate_left(5).wrapping_add(0x9E3779B97F4A7C15);
    }
    Ok(acc)
}

fn params_for_logn(logn: u8) -> Result<Parameters, KeygenError> {
    ensure_supported_logn(logn)?;
    Ok(ParameterSet::new_falcon(1 << logn as usize)
        .expect("logn validated")
        .params())
}

fn generate_keypair<R: RandomSource>(
    params: Parameters,
    rng: &mut R,
) -> Keypair {
    // Sample short polynomials f and g; in a full implementation we'd run
    // norm checks and compute h = g / f mod q.
    let f = sample_gaussian_vector(params.degree, rng);
    let g = sample_gaussian_vector(params.degree, rng);
    let h = build_placeholder_h(&g);

    Keypair {
        public: PublicKey::new(h, params),
        secret: SecretKey::new(f, g, params),
    }
}

fn sample_gaussian_vector<R: RandomSource>(
    degree: usize,
    rng: &mut R,
) -> Vec<i64> {
    let mut out = Vec::with_capacity(degree);
    for _ in 0..degree {
        out.push(discrete_gaussian(rng, 2.5, 0.0));
    }
    out
}

fn build_placeholder_h(g: &[i64]) -> Poly<ModQ> {
    let coeffs: Vec<_> = g.iter().map(|&x| modq_from_i64(x)).collect();
    poly!(coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::field::ModQExt;

    #[test]
    fn deterministic_seed_reproducible() {
        let seed = b"deterministic-seed";
        let k1 = keypair_from_seed(9, seed).unwrap();
        let k2 = keypair_from_seed(9, seed).unwrap();
        assert_eq!(k1, k2);
    }

    #[test]
    fn different_seeds_yield_different_keys() {
        let k1 = keypair_from_seed(9, b"a").unwrap();
        let k2 = keypair_from_seed(9, b"b").unwrap();
        assert_ne!(k1.public.h.coeffs(), k2.public.h.coeffs());
        assert_ne!(k1.secret.f, k2.secret.f);
    }

    #[test]
    fn rejects_unsupported_logn() {
        assert!(matches!(
            keypair_from_seed(8, b"x"),
            Err(KeygenError::UnsupportedLogN(8))
        ));
    }

    #[test]
    fn random_path_uses_rng() {
        let mut rng = XorShift64::new(123);
        let kp = keypair_with_rng(9, &mut rng).unwrap();
        assert_eq!(kp.public.h.coeffs().len(), 1 << 9);
        assert_eq!(kp.secret.f.len(), 1 << 9);
        assert!(kp.public.h.coeffs().iter().any(|c| c.value() != 0));
    }

    #[test]
    fn public_secret_params_align() {
        let kp = keypair_from_seed(10, b"params").unwrap();
        assert_eq!(kp.public.params, kp.secret.params);
        assert_eq!(kp.public.degree(), kp.secret.degree());
        assert_eq!(kp.public.logn(), kp.secret.logn());
    }
}
