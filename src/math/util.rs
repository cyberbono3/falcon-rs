//! Shared helpers for transform-style algorithms.
//!
//! Keeping the bit-reversal permutation in one place avoids duplication between
//! the NTT and complex FFT implementations.

/// Apply an in-place bit-reversal permutation to `values`.
///
/// The length must be a power of two; callers should enforce their own size
/// checks to avoid surprising panics.
pub fn bit_reverse<T>(values: &mut [T]) {
    let n = values.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            values.swap(i, j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::bit_reverse;

    #[test]
    fn bit_reverse_permutation() {
        let mut data = [0, 1, 2, 3, 4, 5, 6, 7];
        bit_reverse(&mut data);
        assert_eq!(data, [0, 4, 2, 6, 1, 5, 3, 7]);
    }
}
