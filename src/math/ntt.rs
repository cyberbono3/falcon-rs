//! Number Theoretic Transform (NTT) over `Z/qZ`.
//!
//! The constants mirror the Falcon reference implementation (`q = 12289`,
//! `n <= 1024`). Twiddle factors are derived from a primitive 1024-th root of
//! unity.

use ark_ff::One;

use super::field::{ModQ, ModQExt, modq_from_i64};
use super::util;

/// Maximum supported `n = 2^logn`.
pub const MAX_LOGN: usize = 10;
/// Maximum transform size.
pub const MAX_N: usize = 1 << MAX_LOGN;
/// Primitive 1024-th root of unity mod `q` (value before reduction).
pub const PRIMITIVE_NTT_ROOT_VALUE: u32 = 10302;

/// Return the primitive NTT root as a field element.
pub fn primitive_ntt_root() -> ModQ {
    modq_from_i64(PRIMITIVE_NTT_ROOT_VALUE as i64)
}

#[inline]
fn assert_supported(n: usize) {
    assert!(n.is_power_of_two(), "NTT size must be a power of two");
    assert!(n <= MAX_N, "NTT size exceeds Falcon bound (n <= {MAX_N})");
}

#[inline]
fn primitive_root_for(n: usize) -> ModQ {
    // PRIMITIVE_NTT_ROOT has order MAX_N, so raising it to (MAX_N / n)
    // yields a primitive n-th root.
    primitive_ntt_root().pow_small((MAX_N / n) as u32)
}

/// In-place forward NTT (cooley-tukey, breadth-first).
pub fn forward(values: &mut [ModQ]) {
    let n = values.len();
    assert_supported(n);

    let omega = primitive_root_for(n);
    util::bit_reverse(values);

    let mut len = 2;
    while len <= n {
        let wlen = omega.pow_small((n / len) as u32);
        for chunk in values.chunks_mut(len) {
            let mut w = ModQ::one();
            for i in 0..len / 2 {
                let u = chunk[i];
                let v = chunk[i + len / 2] * w;
                chunk[i] = u + v;
                chunk[i + len / 2] = u - v;
                w *= wlen;
            }
        }
        len <<= 1;
    }
}

/// In-place inverse NTT.
pub fn inverse(values: &mut [ModQ]) {
    let n = values.len();
    assert_supported(n);

    let omega_inv = primitive_root_for(n).inv_unchecked();
    util::bit_reverse(values);

    let mut len = 2;
    while len <= n {
        let wlen = omega_inv.pow_small((n / len) as u32);
        for chunk in values.chunks_mut(len) {
            let mut w = ModQ::one();
            for i in 0..len / 2 {
                let u = chunk[i];
                let v = chunk[i + len / 2] * w;
                chunk[i] = u + v;
                chunk[i + len / 2] = u - v;
                w *= wlen;
            }
        }
        len <<= 1;
    }

    let inv_n = modq_from_i64(n as i64).inv_unchecked();
    for c in values.iter_mut() {
        *c *= inv_n;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fe;

    #[test]
    fn roundtrip_identity() {
        let mut data = [
            fe!(1),
            fe!(2),
            fe!(3),
            fe!(4),
            fe!(5),
            fe!(6),
            fe!(7),
            fe!(8),
        ];
        forward(&mut data);
        inverse(&mut data);
        assert_eq!(
            data,
            [
                fe!(1),
                fe!(2),
                fe!(3),
                fe!(4),
                fe!(5),
                fe!(6),
                fe!(7),
                fe!(8)
            ]
        );
    }

    #[test]
    #[should_panic]
    fn rejects_non_power_of_two() {
        let mut data = [ModQ::one(); 10];
        forward(&mut data);
    }
}
