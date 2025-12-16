//! Small matrix helpers for Falcon's LDL decomposition and Gram computations.
//!
//! The reference implementation mostly needs 2×2 matrices; we keep the type
//! generic so it works for both modular and complex arithmetic.

use core::ops::{Add, Div, Mul, Sub};

/// Generic 2×2 matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mat2<T> {
    pub data: [[T; 2]; 2],
}

impl<T> Mat2<T> {
    pub const fn new(a00: T, a01: T, a10: T, a11: T) -> Self {
        Self {
            data: [[a00, a01], [a10, a11]],
        }
    }
}

impl<T> Mat2<T>
where
    T: Copy
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>,
{
    /// LDL decomposition for a symmetric matrix.
    ///
    /// Returns `(d00, d11, l10)` such that:
    /// ```text
    /// [d00, 0  ]       [1, 0  ] [d00, 0  ] [1, l10]
    /// [0,   d11] = L * [l10, 1] [0,   d11] [0,  1 ]
    /// ```
    pub fn ldl(&self) -> (T, T, T) {
        let a00 = self.data[0][0];
        let a01 = self.data[0][1];
        let a11 = self.data[1][1];
        let d0 = a00;
        let l10 = a01 / d0;
        let d1 = a11 - l10 * a01;
        (d0, d1, l10)
    }

    /// Multiply by a 2-vector.
    pub fn mul_vector(&self, v: [T; 2]) -> [T; 2] {
        [
            self.data[0][0] * v[0] + self.data[0][1] * v[1],
            self.data[1][0] * v[0] + self.data[1][1] * v[1],
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fe;
    use crate::math::fft::Complex;
    use num_traits::One;

    #[test]
    fn ldl_modq() {
        let m = Mat2::new(fe!(5), fe!(2), fe!(2), fe!(9));
        let (d0, d1, l10) = m.ldl();
        assert_eq!(d0, fe!(5));
        assert_eq!(l10, fe!(2) / fe!(5));
        assert_eq!(d1, fe!(9) - l10 * fe!(2));
        // Recompose the matrix to ensure the decomposition is sane.
        let recon = Mat2::new(d0, d0 * l10, d0 * l10, d1 + d0 * l10 * l10);
        assert_eq!(m, recon);
    }

    #[test]
    fn matrix_vector_complex() {
        let m = Mat2::new(
            Complex::new(1.0, 1.0),
            Complex::new(0.0, 1.0),
            Complex::new(2.0, 0.0),
            Complex::new(1.0, -1.0),
        );
        let out = m.mul_vector([Complex::one(), Complex::new(0.0, 1.0)]);
        // (1+i)*1 + i*i = (1+i) + (-1) = i
        assert_eq!(out[0], Complex::new(0.0, 1.0));
        // Second component sanity: (2)*1 + (1 - i)*i = 2 + i - i^2 = 3 + i
        assert_eq!(out[1], Complex::new(3.0, 1.0));
    }
}
