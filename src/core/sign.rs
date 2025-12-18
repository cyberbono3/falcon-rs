//! Falcon signing scaffolding.
//!
//! This module implements a **non-production** signature generation pipeline
//! that mirrors the structure of `tprest/falcon.py`:
//! - hash-to-point (message → polynomial mod `q`)
//! - (floating-point) LDL-style decomposition (ffLDL)
//! - a Gaussian sampler that consumes the decomposition
//!
//! Important: the current crate's key generation is still a placeholder and
//! does not derive the full NTRU basis required for real Falcon. The signing
//! logic here is therefore best viewed as an API/flow skeleton for further
//! work, not a correct Falcon implementation.

use crate::core::keygen::SecretKey;
use crate::math::Parameters;
use crate::math::params::FALCON_Q;
use crate::math::sampling::{RandomSource, discrete_gaussian};

pub const NONCE_LEN: usize = 40;

/// Errors that can occur during signing.
#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum SignError {
    #[error("secret key parameters unsupported (logn = {0})")]
    UnsupportedLogN(u8),
}

/// Signing configuration knobs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SignConfig {
    /// Base Gaussian standard deviation for sampling.
    pub sigma: f64,
    /// Maximum number of rejection-sampling attempts.
    pub max_attempts: usize,
}

impl Default for SignConfig {
    fn default() -> Self {
        Self {
            sigma: 1.55,
            max_attempts: 64,
        }
    }
}

/// A detached signature (nonce + `s2` vector).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Signature {
    nonce: [u8; NONCE_LEN],
    s2: Vec<i16>,
    params: Parameters,
}

impl Signature {
    pub fn nonce(&self) -> &[u8; NONCE_LEN] {
        &self.nonce
    }

    pub fn s2(&self) -> &[i16] {
        &self.s2
    }

    pub fn params(&self) -> Parameters {
        self.params
    }
}

/// Sign `message` with the given secret key and RNG.
pub fn sign_with_rng<R: RandomSource>(
    secret: &SecretKey,
    message: &[u8],
    rng: &mut R,
) -> Result<Signature, SignError> {
    sign_with_rng_and_config(secret, message, rng, SignConfig::default())
}

pub fn sign_with_rng_and_config<R: RandomSource>(
    secret: &SecretKey,
    message: &[u8],
    rng: &mut R,
    config: SignConfig,
) -> Result<Signature, SignError> {
    let params = secret.params();
    if !params.is_supported() {
        return Err(SignError::UnsupportedLogN(params.log_degree));
    }

    let nonce = random_nonce(rng);
    let c = hash_to_point(params, &nonce, message);

    // Target vector (t0, t1) for sampling in the (2n)-dimensional lattice.
    // In Falcon: t0 is the hash point, t1 is zero.
    let t0: Vec<f64> = c.iter().map(|&x| x as f64).collect();
    let t1 = vec![0.0f64; params.degree];

    let gram = toy_gram_matrix(secret);
    let ldl = ffldl(&gram);

    // Sample and return only the `s2` part (Falcon signatures are encoded
    // as (nonce, s2); verification recomputes `s1` from public data).
    let (_s1, s2) = sample_short_vector(&ldl, &t0, &t1, rng, config);
    Ok(Signature { nonce, s2, params })
}

fn random_nonce<R: RandomSource>(rng: &mut R) -> [u8; NONCE_LEN] {
    let mut out = [0u8; NONCE_LEN];
    let mut i = 0usize;
    while i < NONCE_LEN {
        let w = rng.next_u64().to_le_bytes();
        let take = (NONCE_LEN - i).min(w.len());
        out[i..i + take].copy_from_slice(&w[..take]);
        i += take;
    }
    out
}

/// Hash-to-point scaffolding.
///
/// The reference uses SHAKE256(nonce || message) with rejection sampling into
/// `[0, q)`. To keep dependencies minimal (and avoid pulling in a SHAKE
/// implementation), we use a deterministic xorshift-like mixing routine.
pub fn hash_to_point(
    params: Parameters,
    nonce: &[u8; NONCE_LEN],
    message: &[u8],
) -> Vec<u16> {
    let seed = mix_seed64(0xF4_6A_C9_12_3B_77_E1_19, nonce, message);
    let mut x = seed;
    let mut out = Vec::with_capacity(params.degree);
    while out.len() < params.degree {
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        x = x.wrapping_mul(0x2545F4914F6CDD1D);
        let v = (x % (FALCON_Q as u64)) as u16;
        out.push(v);
    }
    out
}

fn mix_seed64(domain: u64, nonce: &[u8; NONCE_LEN], message: &[u8]) -> u64 {
    // A tiny non-cryptographic mixer to provide deterministic behavior.
    let mut acc = domain;
    for &b in nonce.iter().chain(message) {
        acc ^= b as u64;
        acc = acc.rotate_left(5).wrapping_add(0x9E3779B97F4A7C15);
    }
    acc
}

#[derive(Clone, Debug)]
struct GramMatrix {
    g00: Vec<f64>,
    g01: Vec<f64>,
    g11: Vec<f64>,
}

/// Build a toy Gram matrix from the placeholder secret key.
///
/// Real Falcon uses the NTRU basis and FFT-domain convolution; this builds a
/// per-coordinate SPD matrix to exercise the ffLDL + sampler pipeline.
fn toy_gram_matrix(secret: &SecretKey) -> GramMatrix {
    let f = secret.f();
    let g = secret.g();
    let n = secret.degree();
    debug_assert_eq!(f.len(), n);
    debug_assert_eq!(g.len(), n);

    let mut g00 = Vec::with_capacity(n);
    let mut g01 = Vec::with_capacity(n);
    let mut g11 = Vec::with_capacity(n);

    for (&fi, &gi) in f.iter().zip(g.iter()) {
        let fi = fi as f64;
        let gi = gi as f64;
        let a = fi * fi + gi * gi + 1e-9;
        let b = fi * gi;
        let c = fi * fi + 1e-9;
        g00.push(a);
        g01.push(b);
        g11.push(c);
    }

    GramMatrix { g00, g01, g11 }
}

/// A simplified "ffLDL" decomposition performed coefficient-wise.
#[derive(Clone, Debug)]
struct FfLdl {
    d00: Vec<f64>,
    d11: Vec<f64>,
    l10: Vec<f64>,
}

fn ffldl(gram: &GramMatrix) -> FfLdl {
    let n = gram.g00.len();
    debug_assert_eq!(gram.g01.len(), n);
    debug_assert_eq!(gram.g11.len(), n);

    let mut d00 = Vec::with_capacity(n);
    let mut d11 = Vec::with_capacity(n);
    let mut l10 = Vec::with_capacity(n);

    for ((&a0, &b), &c) in gram.g00.iter().zip(&gram.g01).zip(&gram.g11) {
        let a = a0.max(1e-12);
        let li = b / a;
        let di = (c - li * b).max(1e-12);
        d00.push(a);
        l10.push(li);
        d11.push(di);
    }

    FfLdl { d00, d11, l10 }
}

fn sample_short_vector<R: RandomSource>(
    ldl: &FfLdl,
    t0: &[f64],
    t1: &[f64],
    rng: &mut R,
    config: SignConfig,
) -> (Vec<i16>, Vec<i16>) {
    debug_assert_eq!(ldl.d00.len(), t0.len());
    debug_assert_eq!(ldl.d11.len(), t1.len());
    debug_assert_eq!(ldl.l10.len(), t0.len());

    let n = t0.len();
    let mut s1 = vec![0i16; n];
    let mut s2 = vec![0i16; n];

    for _ in 0..config.max_attempts.max(1) {
        for i in 0..n {
            let sigma1 = config.sigma * ldl.d11[i].sqrt();
            let z1 = discrete_gaussian(rng, sigma1, t1[i]);

            let sigma0 = config.sigma * ldl.d00[i].sqrt();
            let center0 = t0[i] - ldl.l10[i] * (z1 as f64);
            let z0 = discrete_gaussian(rng, sigma0, center0);

            s1[i] = z0.clamp(i16::MIN as i64, i16::MAX as i64) as i16;
            s2[i] = z1.clamp(i16::MIN as i64, i16::MAX as i64) as i16;
        }

        // A very lightweight rejection predicate to bound the signature.
        // Real Falcon uses a scheme-specific norm bound.
        let l2: i64 = s1
            .iter()
            .zip(&s2)
            .map(|(&a, &b)| (a as i64) * (a as i64) + (b as i64) * (b as i64))
            .sum();
        let bound = (n as i64) * 10_000;
        if l2 <= bound {
            return (s1, s2);
        }
    }

    (s1, s2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::keygen::Keypair;
    use crate::math::sampling::XorShift64;

    #[test]
    fn sign_is_deterministic_given_rng_stream() {
        let kp = Keypair::from_seed(9, b"sign-key").unwrap();
        let msg = b"hello";

        let mut rng1 = XorShift64::new(123);
        let mut rng2 = XorShift64::new(123);
        let sig1 = sign_with_rng(kp.secret(), msg, &mut rng1).unwrap();
        let sig2 = sign_with_rng(kp.secret(), msg, &mut rng2).unwrap();
        assert_eq!(sig1, sig2);
        assert_eq!(sig1.s2().len(), 1 << 9);
        assert_eq!(sig1.nonce().len(), NONCE_LEN);
    }

    #[test]
    fn hash_to_point_in_range() {
        let kp = Keypair::from_seed(9, b"hash-key").unwrap();
        let nonce = [0u8; NONCE_LEN];
        let params = kp.secret().params();
        let c = hash_to_point(params, &nonce, b"msg");
        assert_eq!(c.len(), kp.secret().degree());
        assert!(c.iter().all(|&x| x < FALCON_Q as u16));
    }
}
