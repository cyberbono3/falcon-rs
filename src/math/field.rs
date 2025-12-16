//! Finite-field helpers for `Z/qZ` arithmetic.
//!
//! Falcon works over the prime modulus `q = 12289`. We adopt Arkworks' field
//! traits and Montgomery backend to ensure correctness while exposing a small,
//! ergonomic API for the rest of the crate.

#![allow(non_local_definitions)]

use ark_ff::{Field, Fp64, MontBackend, MontConfig, PrimeField};

use super::params::FALCON_Q;

#[derive(MontConfig)]
#[modulus = "12289"]
#[generator = "11"]
pub struct ModQConfig;

/// Field element modulo `12289` using Arkworks' Montgomery backend.
pub type ModQ = Fp64<MontBackend<ModQConfig, 1>>;

/// Convenience helpers layered on top of Arkworks' `PrimeField`.
pub trait ModQExt {
    /// Canonical representative in the range `[0, q)`.
    fn value(&self) -> u32;
    /// Centered representative in `[-q/2, q/2]`.
    fn centered(&self) -> i32;
    /// Exponentiation by a small integer.
    fn pow_small(self, exp: u32) -> Self;
    /// Multiplicative inverse, panicking on zero.
    fn inv_unchecked(self) -> Self;
}

/// Create a field element from an integer, reducing into the canonical range.
pub fn modq_from_i64(value: i64) -> ModQ {
    let modulus = FALCON_Q as i64;
    let mut v = value % modulus;
    if v < 0 {
        v += modulus;
    }
    ModQ::from(v as u64)
}

impl ModQExt for ModQ {
    #[inline]
    fn value(&self) -> u32 {
        // The first limb holds the reduced value for this small modulus.
        (self.into_bigint().0[0] % FALCON_Q as u64) as u32
    }

    #[inline]
    fn centered(&self) -> i32 {
        let q = FALCON_Q as i32;
        let v = self.value() as i32;
        if v > q / 2 { v - q } else { v }
    }

    #[inline]
    fn pow_small(self, exp: u32) -> Self {
        // Arkworks `pow` takes little-endian limbs.
        Field::pow(&self, [exp as u64])
    }

    #[inline]
    fn inv_unchecked(self) -> Self {
        self.inverse()
            .expect("attempted to invert zero in ModQ::inv_unchecked")
    }
}

/// Convenient macro for constructing [`ModQ`] values.
#[macro_export]
macro_rules! fe {
    ($value:expr) => {
        $crate::math::field::modq_from_i64($value as i64)
    };
}

/// Create a vector of [`ModQ`]s.
#[macro_export]
macro_rules! fe_vec {
    ($val:expr; $count:expr) => {
        vec![$crate::math::field::modq_from_i64($val as i64); $count]
    };
    ($($val:expr),* $(,)?) => {
        vec![$($crate::math::field::modq_from_i64($val as i64)),*]
    };
}

/// Create an array of [`ModQ`]s.
#[macro_export]
macro_rules! fe_array {
    ($val:expr; $count:expr) => {
        [$crate::math::field::modq_from_i64($val as i64); $count]
    };
    ($($val:expr),* $(,)?) => {
        [$($crate::math::field::modq_from_i64($val as i64)),*]
    };
}

#[cfg(test)]
mod tests {
    use super::{ModQ, ModQExt, modq_from_i64};
    use crate::math::params::FALCON_Q;
    use ark_ff::{One, Zero};

    #[test]
    fn add_and_sub_wrap() {
        let a = modq_from_i64(FALCON_Q as i64 - 1);
        let b = ModQ::one();
        assert_eq!((a + b).value(), 0);
        assert_eq!((ModQ::zero() - b).value(), FALCON_Q - 1);
    }

    #[test]
    fn mul_and_inv() {
        let a = modq_from_i64(7);
        let inv_a = a.inv_unchecked();
        assert_eq!((a * inv_a).value(), 1);
        let b = modq_from_i64(1_000);
        assert_eq!(((a * b) / a).value(), b.value());
    }

    #[test]
    fn pow_matches_iterated_mul() {
        let a = modq_from_i64(5);
        assert_eq!(a.pow_small(3).value(), (a * a * a).value());
        assert_eq!(a.pow_small(0), ModQ::one());
    }

    #[test]
    fn centered_representation() {
        let half = (FALCON_Q / 2) as i32;
        assert_eq!(modq_from_i64(half as i64).centered(), half);
        assert_eq!(
            modq_from_i64((half + 1) as i64).centered(),
            -((FALCON_Q as i32) - (half + 1))
        );
    }
}
