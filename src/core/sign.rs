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
const DEFAULT_HASH_DOMAIN: u64 = 0xF4_6A_C9_12_3B_77_E1_19;
const MIN_VARIANCE: f64 = 1e-12;
const MIN_GRAM_BIAS: f64 = 1e-9;

/// Errors that can occur during signing.
#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum SignError {
    #[error("secret key parameters unsupported (logn = {0})")]
    UnsupportedLogN(u8),
    #[error("rejection sampling failed after {0} attempts")]
    RejectionSamplingFailed(usize),
}

/// Signing configuration knobs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SignConfig {
    /// Base Gaussian standard deviation for sampling.
    pub sigma: f64,
    /// Maximum number of rejection-sampling attempts.
    pub max_attempts: usize,
    /// Toy rejection bound factor (bound = `degree * factor`).
    pub l2_bound_factor: i64,
}

impl Default for SignConfig {
    fn default() -> Self {
        Self {
            sigma: 1.55,
            max_attempts: 64,
            l2_bound_factor: 10_000,
        }
    }
}

impl SignConfig {
    fn attempts(self) -> usize {
        self.max_attempts.max(1)
    }

    fn l2_bound(self, degree: usize) -> i64 {
        (degree as i64) * self.l2_bound_factor
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

/// Hash-to-point strategy.
pub trait HashToPoint {
    fn hash_to_point(
        &self,
        params: Parameters,
        nonce: &[u8; NONCE_LEN],
        message: &[u8],
    ) -> Vec<u16>;
}

/// Default, dependency-free hash-to-point scaffolding (non-cryptographic).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct XorHashToPoint {
    domain: u64,
}

impl Default for XorHashToPoint {
    fn default() -> Self {
        Self {
            domain: DEFAULT_HASH_DOMAIN,
        }
    }
}

impl XorHashToPoint {
    pub const fn new(domain: u64) -> Self {
        Self { domain }
    }
}

impl HashToPoint for XorHashToPoint {
    fn hash_to_point(
        &self,
        params: Parameters,
        nonce: &[u8; NONCE_LEN],
        message: &[u8],
    ) -> Vec<u16> {
        let seed = mix_seed64(self.domain, nonce, message);
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
}

/// Stateful signer, holding strategies for hashing and sampling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Signer<H = XorHashToPoint> {
    config: SignConfig,
    hasher: H,
}

impl Default for Signer {
    fn default() -> Self {
        Self {
            config: SignConfig::default(),
            hasher: XorHashToPoint::default(),
        }
    }
}

impl<H> Signer<H> {
    #[must_use]
    pub fn with_config(mut self, config: SignConfig) -> Self {
        self.config = config;
        self
    }

    #[must_use]
    pub fn with_hasher<NH: HashToPoint>(self, hasher: NH) -> Signer<NH> {
        Signer {
            config: self.config,
            hasher,
        }
    }
}

impl<H: HashToPoint> Signer<H> {
    /// Sign `message` with the given secret key and RNG.
    pub fn sign_with_rng<R: RandomSource>(
        &self,
        secret: &SecretKey,
        message: &[u8],
        rng: &mut R,
    ) -> Result<Signature, SignError> {
        let params = secret.params();
        if !params.is_supported() {
            return Err(SignError::UnsupportedLogN(params.log_degree));
        }

        let nonce = random_nonce(rng);
        let c = self.hasher.hash_to_point(params, &nonce, message);

        // Target vector (t0, t1) for sampling in the (2n)-dimensional lattice.
        // In Falcon: t = (c, 0), but the *output* still depends on `c` through
        // the secret basis. Since `SecretKey` is still a placeholder in this
        // crate, we inject the hash point into both components to keep this
        // scaffolding signature dependent on the message.
        //
        // NOTE: With real Falcon keys, the sampler uses the secret basis to
        // obtain a *short* lattice solution for a target modulo `q`. Since
        // `SecretKey` is still a placeholder in this crate, we scale the target
        // into a small centered interval to keep this scaffolding sampler
        // well-behaved.
        let t0: Vec<f64> = c.iter().map(|&x| point_to_target(x)).collect();
        let t1: Vec<f64> = t0.iter().map(|&x| -x).collect();

        let gram = toy_gram_matrix(secret);
        let ldl = ffldl(&gram);

        // Sample and return only the `s2` part (Falcon signatures are encoded
        // as (nonce, s2); verification recomputes `s1` from public data).
        let (_s1, s2) = sample_short_vector(&ldl, &t0, &t1, rng, self.config)?;
        Ok(Signature { nonce, s2, params })
    }
}

pub fn sign_with_rng<R: RandomSource>(
    secret: &SecretKey,
    message: &[u8],
    rng: &mut R,
) -> Result<Signature, SignError> {
    Signer::default().sign_with_rng(secret, message, rng)
}

pub fn sign_with_rng_and_config<R: RandomSource>(
    secret: &SecretKey,
    message: &[u8],
    rng: &mut R,
    config: SignConfig,
) -> Result<Signature, SignError> {
    Signer::default()
        .with_config(config)
        .sign_with_rng(secret, message, rng)
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
    XorHashToPoint::default().hash_to_point(params, nonce, message)
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

fn point_to_target(x: u16) -> f64 {
    let q = FALCON_Q as i32;
    let half_q = q / 2;
    let mut v = x as i32;
    if v > half_q {
        v -= q;
    }
    v as f64 / (FALCON_Q as f64)
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
        let a = fi * fi + gi * gi + MIN_GRAM_BIAS;
        let b = fi * gi;
        let c = fi * fi + MIN_GRAM_BIAS;
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
        let a = a0.max(MIN_VARIANCE);
        let li = b / a;
        let di = (c - li * b).max(MIN_VARIANCE);
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
) -> Result<(Vec<i16>, Vec<i16>), SignError> {
    debug_assert_eq!(ldl.d00.len(), t0.len());
    debug_assert_eq!(ldl.d11.len(), t1.len());
    debug_assert_eq!(ldl.l10.len(), t0.len());

    let n = t0.len();
    let mut s1 = vec![0i16; n];
    let mut s2 = vec![0i16; n];

    for _ in 0..config.attempts() {
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
        let bound = config.l2_bound(n);
        if l2 <= bound {
            return Ok((s1, s2));
        }
    }

    Err(SignError::RejectionSamplingFailed(config.attempts()))
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
    fn sign_differs_for_different_messages_with_same_rng_seed() {
        let kp = Keypair::from_seed(9, b"sign-key").unwrap();
        let msg1 = b"hello";
        let msg2 = b"world";

        let mut rng1 = XorShift64::new(777);
        let mut rng2 = XorShift64::new(777);
        let sig1 = sign_with_rng(kp.secret(), msg1, &mut rng1).unwrap();
        let sig2 = sign_with_rng(kp.secret(), msg2, &mut rng2).unwrap();
        assert_ne!(sig1, sig2);
        // Nonce is generated solely from RNG, so with identical RNG streams it
        // should match; the signature difference should come from `s2`.
        assert_eq!(sig1.nonce(), sig2.nonce());
        assert_ne!(sig1.s2(), sig2.s2());
    }

    #[test]
    fn signer_domain_affects_signature() {
        let kp = Keypair::from_seed(9, b"sign-key").unwrap();
        let msg = b"hello";

        let signer_a = Signer::default().with_hasher(XorHashToPoint::new(1));
        let signer_b = Signer::default().with_hasher(XorHashToPoint::new(2));

        let mut rng1 = XorShift64::new(555);
        let mut rng2 = XorShift64::new(555);
        let sig_a =
            signer_a.sign_with_rng(kp.secret(), msg, &mut rng1).unwrap();
        let sig_b =
            signer_b.sign_with_rng(kp.secret(), msg, &mut rng2).unwrap();
        assert_ne!(sig_a.s2(), sig_b.s2());
        assert_eq!(sig_a.nonce(), sig_b.nonce());
    }

    #[test]
    fn signing_twice_advances_nonce() {
        let kp = Keypair::from_seed(9, b"sign-key").unwrap();
        let msg = b"hello";

        let mut rng = XorShift64::new(999);
        let sig1 = sign_with_rng(kp.secret(), msg, &mut rng).unwrap();
        let sig2 = sign_with_rng(kp.secret(), msg, &mut rng).unwrap();
        assert_ne!(sig1.nonce(), sig2.nonce());
    }

    #[test]
    fn rejects_unsupported_secret_key_params() {
        let params = Parameters {
            name: "Falcon-256 (unsupported)",
            degree: 256,
            log_degree: 8,
            modulus: FALCON_Q,
            max_sig_len: 0,
            pk_len: 0,
            sk_len: 0,
        };
        let sk = SecretKey::new(
            vec![0; params.degree],
            vec![0; params.degree],
            params,
        );
        let mut rng = XorShift64::new(1);
        assert!(matches!(
            sign_with_rng(&sk, b"msg", &mut rng),
            Err(SignError::UnsupportedLogN(8))
        ));
    }

    #[test]
    fn rejection_failure_is_reported() {
        let kp = Keypair::from_seed(9, b"sign-key").unwrap();
        let msg = b"hello";

        let config = SignConfig {
            max_attempts: 3,
            l2_bound_factor: -1,
            ..SignConfig::default()
        };
        let mut rng = XorShift64::new(123);
        assert!(matches!(
            sign_with_rng_and_config(kp.secret(), msg, &mut rng, config),
            Err(SignError::RejectionSamplingFailed(3))
        ));
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

    #[test]
    fn hash_to_point_is_deterministic_for_same_inputs() {
        let kp = Keypair::from_seed(9, b"hash-key").unwrap();
        let params = kp.secret().params();
        let nonce = [7u8; NONCE_LEN];
        let c1 = hash_to_point(params, &nonce, b"msg");
        let c2 = hash_to_point(params, &nonce, b"msg");
        assert_eq!(c1, c2);
    }
}
