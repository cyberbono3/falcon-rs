//! Optional Arkworks compatibility layer.
//!
//! This module is only compiled when the `ark-compat` feature is enabled. It
//! exposes an Arkworks-native field type (`ArkModQ`) and simple conversions to
//! and from the crate's [`ModQ`] type. This allows reusing Arkworks gadgets
//! (e.g., MSM or constraint systems) while keeping a lightweight `ModQ`
//! internally.

#![cfg(feature = "ark-compat")]

use ark_ff::{Fp64, MontBackend, MontConfig};

use super::field::{ModQ, ModQExt, modq_from_i64};
use super::params::FALCON_Q;

#[derive(MontConfig)]
#[modulus = "12289"]
#[generator = "11"]
pub struct ModQConfig;
/// Arkworks-backed representation of `Z/12289Z`.
pub type ArkModQ = Fp64<MontBackend<ModQConfig, 1>>;

impl From<ModQ> for ArkModQ {
    fn from(value: ModQ) -> Self {
        ArkModQ::from(value.value() as u64)
    }
}

impl From<ArkModQ> for ModQ {
    fn from(value: ArkModQ) -> Self {
        let limbs = value.into_bigint().0;
        let v = limbs[0] % FALCON_Q as u64;
        modq_from_i64(v as i64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_modq_conversion() {
        let original = modq_from_i64(42);
        let ark: ArkModQ = original.into();
        let back: ModQ = ark.into();
        assert_eq!(back, original);
    }
}
