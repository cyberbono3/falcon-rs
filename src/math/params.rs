//! Falcon parameter set definitions and common constants.
//!
//! The values mirror the reference implementation (falcon.py / PQClean).
//! They are kept in one place to make it easy to share across the math,
//! encoding, keygen, and signing modules.

/// Prime modulus used by Falcon polynomials.
pub const FALCON_Q: u32 = 12_289;
/// `log2(FALCON_Q)`, used for bit sizing in encoding routines.
pub const LOGQ: u8 = 14;
/// Minimum degree supported by the scheme (`n = 2^logn`).
pub const MIN_LOGN: u8 = 9; // 2^9 = 512
/// Maximum degree supported by the scheme (`n = 2^logn`).
pub const MAX_LOGN: u8 = 10; // 2^10 = 1024

/// Supported Falcon parameter sets.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum ParameterSet {
    Falcon512,
    Falcon1024,
}

/// High-level description of a Falcon parameter set.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Parameters {
    /// Human-friendly name (e.g. "Falcon-512").
    pub name: &'static str,
    /// Degree of the lattice dimension (n).
    pub degree: usize,
    /// `log2(degree)`, convenient for table lookups.
    pub log_degree: u8,
    /// Public modulus `q`.
    pub modulus: u32,
    /// Maximum encoded signature length in bytes.
    pub max_sig_len: usize,
    /// Encoded public key length in bytes.
    pub pk_len: usize,
    /// Encoded secret key length in bytes.
    pub sk_len: usize,
}

impl Parameters {
    /// Returns whether `self` uses a supported `logn`.
    pub const fn is_supported(&self) -> bool {
        self.log_degree >= MIN_LOGN && self.log_degree <= MAX_LOGN
    }
}

impl ParameterSet {
    /// Construct a parameter set from a Falcon degree (512 or 1024).
    pub const fn new_falcon(degree: usize) -> Option<Self> {
        match degree {
            512 => Some(ParameterSet::Falcon512),
            1024 => Some(ParameterSet::Falcon1024),
            _ => None,
        }
    }

    /// Return all supported parameter sets.
    pub const fn all() -> &'static [ParameterSet] {
        &[ParameterSet::Falcon512, ParameterSet::Falcon1024]
    }

    /// Return the strongly typed parameters for this set.
    pub const fn params(self) -> Parameters {
        match self {
            ParameterSet::Falcon512 => Parameters {
                name: "Falcon-512",
                degree: 512,
                log_degree: 9,
                modulus: FALCON_Q,
                max_sig_len: 690,
                pk_len: 897,
                sk_len: 1_281,
            },
            ParameterSet::Falcon1024 => Parameters {
                name: "Falcon-1024",
                degree: 1024,
                log_degree: 10,
                modulus: FALCON_Q,
                max_sig_len: 1_330,
                pk_len: 1_793,
                sk_len: 2_305,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn params_are_supported() {
        for set in [ParameterSet::Falcon512, ParameterSet::Falcon1024] {
            assert!(set.params().is_supported());
        }
    }

    #[test]
    fn params_switch() {
        let p512 = ParameterSet::Falcon512.params();
        let p1024 = ParameterSet::Falcon1024.params();

        assert_eq!(p512.degree, 512);
        assert_eq!(p1024.max_sig_len, 1_330);
    }

    #[test]
    fn new_falcon_from_degree() {
        assert_eq!(
            ParameterSet::new_falcon(512),
            Some(ParameterSet::Falcon512)
        );
        assert_eq!(
            ParameterSet::new_falcon(1024),
            Some(ParameterSet::Falcon1024)
        );
        assert_eq!(ParameterSet::new_falcon(256), None);
    }

    #[test]
    fn all_returns_every_supported_set_once() {
        let all = ParameterSet::all();
        assert_eq!(all.len(), 2);
        assert!(all.contains(&ParameterSet::Falcon512));
        assert!(all.contains(&ParameterSet::Falcon1024));
    }

    #[test]
    fn table_entries_match_enum_order() {
        // Ensure params match the expected enum discriminants.
        assert_eq!(
            ParameterSet::Falcon512.params().degree,
            ParameterSet::Falcon512 as usize * 512 + 512
        );
        assert_eq!(ParameterSet::Falcon1024.params().degree, 1024);
    }
}
