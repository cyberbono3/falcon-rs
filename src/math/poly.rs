//! Polynomial helpers mirroring the structure of `falcon.py`.
//!
//! The implementation keeps things lightweight but exposes the same operations
//! the Python reference relies on: coefficient-wise arithmetic, evaluation,
//! and NTT-backed convolution on `Z/qZ` polynomials.

use core::fmt;
use core::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

use ark_ff::{PrimeField, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, Polynomial as ArkPolynomial};

use super::fft::Complex;
use super::field::{ModQ, ModQExt};
use super::ntt;

/// A univariate polynomial with coefficients in a finite field.
#[derive(Clone, PartialEq, Eq)]
pub struct Poly<T: PrimeField> {
    coeffs: Vec<T>,
}

impl<T: PrimeField> fmt::Debug for Poly<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_tuple("Poly").field(&self.coeffs).finish()
    }
}

impl<T: PrimeField> Poly<T> {
    /// Construct a polynomial, trimming high-order zero coefficients.
    pub fn new(coeffs: impl Into<Vec<T>>) -> Self {
        let mut coeffs = coeffs.into();
        Self::normalize(&mut coeffs);
        Self { coeffs }
    }

    /// Degree of the polynomial (`-1` for the zero polynomial).
    pub fn degree(&self) -> isize {
        if self.coeffs.is_empty() {
            -1
        } else {
            (self.coeffs.len() - 1) as isize
        }
    }

    /// Coefficient slice.
    pub fn coeffs(&self) -> &[T] {
        &self.coeffs
    }

    /// Mutable coefficient slice.
    pub fn coeffs_mut(&mut self) -> &mut [T] {
        &mut self.coeffs
    }

    /// Resize to the requested number of coefficients.
    pub fn resize(&mut self, n: usize) {
        self.coeffs.resize(n, T::zero());
    }

    /// Evaluate at a point using Arkworks' dense polynomial logic.
    pub fn evaluate(&self, x: T) -> T {
        self.to_dense().evaluate(&x)
    }

    /// In-place addition with another polynomial.
    pub fn add_assign_poly(&mut self, rhs: &Self) {
        let lhs = self.to_dense();
        let rhs_dense = rhs.to_dense();
        let res = &lhs + &rhs_dense;
        self.coeffs = res.coeffs;
        Self::normalize(&mut self.coeffs);
    }

    /// In-place subtraction with another polynomial.
    pub fn sub_assign_poly(&mut self, rhs: &Self) {
        let lhs = self.to_dense();
        let rhs_dense = rhs.to_dense();
        let res = &lhs - &rhs_dense;
        self.coeffs = res.coeffs;
        Self::normalize(&mut self.coeffs);
    }

    /// Multiply by a scalar.
    pub fn scale(&mut self, scalar: T) {
        let mut dense = self.to_dense();
        dense.coeffs.iter_mut().for_each(|c| *c *= scalar);
        self.coeffs = dense.coeffs;
        Self::normalize(&mut self.coeffs);
    }

    fn normalize(coeffs: &mut Vec<T>) {
        while matches!(coeffs.last(), Some(c) if *c == T::zero()) {
            coeffs.pop();
        }
    }

    fn to_dense(&self) -> DensePolynomial<T> {
        DensePolynomial::from_coefficients_vec(self.coeffs.clone())
    }
}

impl<T: PrimeField> Add for Poly<T> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        self.add_assign_poly(&rhs);
        self
    }
}

impl<T: PrimeField> Sub for Poly<T> {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self::Output {
        self.sub_assign_poly(&rhs);
        self
    }
}

impl<T: PrimeField> AddAssign for Poly<T> {
    fn add_assign(&mut self, rhs: Self) {
        self.add_assign_poly(&rhs);
    }
}

impl<T: PrimeField> SubAssign for Poly<T> {
    fn sub_assign(&mut self, rhs: Self) {
        self.sub_assign_poly(&rhs);
    }
}

impl<T: PrimeField> Mul<T> for Poly<T> {
    type Output = Self;
    fn mul(mut self, rhs: T) -> Self::Output {
        self.scale(rhs);
        self
    }
}

impl<T: PrimeField> MulAssign<T> for Poly<T> {
    fn mul_assign(&mut self, rhs: T) {
        self.scale(rhs);
    }
}

impl<T: PrimeField> From<Vec<T>> for Poly<T> {
    fn from(value: Vec<T>) -> Self {
        Self::new(value)
    }
}

impl<T: PrimeField, const N: usize> From<[T; N]> for Poly<T> {
    fn from(value: [T; N]) -> Self {
        Self::new(Vec::from(value))
    }
}

/// `Z/qZ`-specific helpers that lean on the NTT.
impl Poly<ModQ> {
    /// Convert coefficients to complex numbers (useful for floating-point FFT
    /// experiments mirroring `falcon.py`).
    pub fn to_complex(&self) -> Vec<Complex> {
        self.coeffs
            .iter()
            .map(|c| Complex::new(c.value() as f64, 0.0))
            .collect()
    }

    /// Multiply two polynomials using the NTT (degree must fit in `n <= 1024`).
    pub fn ntt_product(&self, rhs: &Self) -> Self {
        let n = (self.coeffs.len().max(rhs.coeffs.len()) * 2)
            .next_power_of_two()
            .max(1);
        assert!(n <= ntt::MAX_N, "degree too large for Falcon NTT (n = {n})");

        let mut a = self.coeffs.clone();
        let mut b = rhs.coeffs.clone();
        a.resize(n, ModQ::zero());
        b.resize(n, ModQ::zero());

        ntt::forward(&mut a);
        ntt::forward(&mut b);
        for (ai, bi) in a.iter_mut().zip(&b) {
            *ai *= *bi;
        }
        ntt::inverse(&mut a);

        Poly::new(a)
    }

    /// Naive convolution (mostly for tests).
    pub fn schoolbook_product(&self, rhs: &Self) -> Self {
        let lhs = self.to_dense();
        let rhs_dense = rhs.to_dense();
        let res = &lhs * &rhs_dense;
        Poly::new(res.coeffs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fe;

    #[test]
    fn add_and_subtract() {
        let mut a: Poly<ModQ> = Poly::from([fe!(1), fe!(2), fe!(3)]);
        let b: Poly<ModQ> = Poly::from([fe!(5), fe!(6)]);
        a.add_assign_poly(&b);
        assert_eq!(a.coeffs(), &[fe!(6), fe!(8), fe!(3)]);
        a.sub_assign_poly(&b);
        assert_eq!(a.coeffs(), &[fe!(1), fe!(2), fe!(3)]);
    }

    #[test]
    fn evaluate_horner() {
        let p: Poly<ModQ> = Poly::from([fe!(3), fe!(0), fe!(2)]); // 2x^2 + 3
        assert_eq!(p.evaluate(fe!(2)).value(), fe!(11).value());
    }

    #[test]
    fn ntt_matches_schoolbook() {
        let p: Poly<ModQ> = Poly::from([fe!(1), fe!(2), fe!(3)]);
        let q: Poly<ModQ> = Poly::from([fe!(4), fe!(5)]);
        let prod_ntt = p.ntt_product(&q);
        let prod_naive = p.schoolbook_product(&q);
        assert_eq!(prod_ntt, prod_naive);
    }
}
