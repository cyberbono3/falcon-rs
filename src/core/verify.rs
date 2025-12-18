//! Falcon verification scaffolding.
//!
//! This follows the high-level Falcon verification logic:
//! - recompute the hash-to-point `c` from `(nonce, message)`
//! - compute `s1 = c - s2*h (mod q)` in the Falcon ring
//! - accept if the squared norm of `(s1, s2)` is below a bound
//!
//! Note: this crate's `core::sign` and `core::keygen` are still scaffolding;
//! bounds and hashing are therefore configurable, and the default hash-to-point
//! currently matches `core::sign` (not SHAKE256).

use crate::core::keygen::PublicKey;
use crate::core::sign::Signature;
use crate::math::field::{ModQExt, modq_from_i64};
use crate::math::{ModQ, Parameters, Poly};
use ark_ff::Zero;

/// Verification failures.
#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum VerifyError {
    #[error("public key parameters unsupported (logn = {0})")]
    UnsupportedLogN(u8),
    #[error("signature parameters do not match public key")]
    ParameterMismatch,
    #[error("signature s2 length mismatch (expected {expected}, got {actual})")]
    InvalidS2Length { expected: usize, actual: usize },
    #[error("signature norm too large")]
    NormTooLarge,
}

/// Verification configuration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct VerifyConfig {
    /// Upper bound on the squared norm: `sum_i (s1[i]^2 + s2[i]^2)`.
    pub l2_bound: i128,
}

impl VerifyConfig {
    pub fn for_params(params: Parameters) -> Self {
        // Scaffolding-only bound; real Falcon uses parameter-specific constants.
        Self {
            l2_bound: (params.degree as i128) * 10_000,
        }
    }
}

/// Verify `signature` against `message` and the `public` key.
pub fn verify(
    public: &PublicKey,
    message: &[u8],
    signature: &Signature,
) -> Result<(), VerifyError> {
    verify_with_config(
        public,
        message,
        signature,
        VerifyConfig::for_params(public.params()),
    )
}

pub fn verify_with_config(
    public: &PublicKey,
    message: &[u8],
    signature: &Signature,
    config: VerifyConfig,
) -> Result<(), VerifyError> {
    let params = public.params();
    if !params.is_supported() {
        return Err(VerifyError::UnsupportedLogN(params.log_degree));
    }
    if signature.params() != params {
        return Err(VerifyError::ParameterMismatch);
    }

    let n = params.degree;
    let s2 = signature.s2();
    if s2.len() != n {
        return Err(VerifyError::InvalidS2Length {
            expected: n,
            actual: s2.len(),
        });
    }

    let c =
        crate::core::sign::hash_to_point(params, signature.nonce(), message);
    let c_modq: Vec<ModQ> =
        c.into_iter().map(|x| modq_from_i64(x as i64)).collect();
    let s2_modq: Vec<ModQ> =
        s2.iter().map(|&x| modq_from_i64(x as i64)).collect();

    let prod = negacyclic_product_modq(&s2_modq, public.h(), n);
    let s1: Vec<ModQ> = c_modq
        .into_iter()
        .zip(prod)
        .map(|(ci, pi)| ci - pi)
        .collect();

    let norm = squared_norm(&s1, s2);
    if norm > config.l2_bound {
        return Err(VerifyError::NormTooLarge);
    }
    Ok(())
}

fn squared_norm(s1: &[ModQ], s2: &[i16]) -> i128 {
    s1.iter()
        .zip(s2)
        .map(|(a, &b)| {
            let a = a.centered() as i128;
            let b = b as i128;
            a * a + b * b
        })
        .sum()
}

fn negacyclic_product_modq(
    lhs: &[ModQ],
    rhs: &Poly<ModQ>,
    n: usize,
) -> Vec<ModQ> {
    let mut a = Poly::from(lhs.to_vec());
    a.resize(n);

    let mut b = rhs.clone();
    b.resize(n);

    let prod = a.ntt_product(&b);
    let mut out = vec![ModQ::zero(); n];
    for (k, &ck) in prod.coeffs().iter().enumerate() {
        if k < n {
            out[k] += ck;
        } else {
            out[k - n] -= ck;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::keygen::Keypair;
    use crate::core::sign::sign_with_rng;
    use crate::math::sampling::XorShift64;

    #[test]
    fn rejects_param_mismatch() {
        let kp = Keypair::from_seed(9, b"verify-key").unwrap();
        let mut rng = XorShift64::new(1);
        let sig = sign_with_rng(kp.secret(), b"msg", &mut rng).unwrap();
        let forged = Signature::new(
            *sig.nonce(),
            sig.s2().to_vec(),
            crate::math::params::ParameterSet::Falcon1024.params(),
        );
        assert!(matches!(
            verify(kp.public(), b"msg", &forged),
            Err(VerifyError::ParameterMismatch)
        ));
    }

    #[test]
    fn rejects_s2_length_mismatch() {
        let kp = Keypair::from_seed(9, b"verify-key").unwrap();
        let params = kp.public().params();
        let sig = Signature::new(
            [0u8; crate::core::sign::NONCE_LEN],
            vec![0i16; params.degree - 1],
            params,
        );
        assert!(matches!(
            verify(kp.public(), b"msg", &sig),
            Err(VerifyError::InvalidS2Length { .. })
        ));
    }

    #[test]
    fn accepts_with_permissive_bound() {
        let kp = Keypair::from_seed(9, b"verify-key").unwrap();
        let mut rng = XorShift64::new(123);
        let sig = sign_with_rng(kp.secret(), b"msg", &mut rng).unwrap();

        let config = VerifyConfig {
            l2_bound: i128::MAX,
        };
        verify_with_config(kp.public(), b"msg", &sig, config).unwrap();
    }

    #[test]
    fn rejects_with_tight_bound() {
        let kp = Keypair::from_seed(9, b"verify-key").unwrap();
        let mut rng = XorShift64::new(123);
        let sig = sign_with_rng(kp.secret(), b"msg", &mut rng).unwrap();

        let config = VerifyConfig { l2_bound: 0 };
        assert!(matches!(
            verify_with_config(kp.public(), b"msg", &sig, config),
            Err(VerifyError::NormTooLarge)
        ));
    }
}
