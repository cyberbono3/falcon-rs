//! Minimal complex FFT used by Falcon's floating-point routines.
//!
//! The implementation is deliberately small but mirrors the `falcon.py`
//! behavior: iterative Cooley-Tukey with bit-reversal permutation.

use core::f64::consts::PI;
use core::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign,
};

use num_traits::{One, Zero};

use super::util;

/// Lightweight complex number type (to avoid pulling external dependencies).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    pub fn abs2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    pub fn exp_i(theta: f64) -> Self {
        Self::new(theta.cos(), theta.sin())
    }
}

impl Zero for Complex {
    fn zero() -> Self {
        Complex::new(0.0, 0.0)
    }

    fn is_zero(&self) -> bool {
        self.re == 0.0 && self.im == 0.0
    }
}

impl One for Complex {
    fn one() -> Self {
        Complex::new(1.0, 0.0)
    }

    fn is_one(&self) -> bool {
        self.re == 1.0 && self.im == 0.0
    }
}

impl Add for Complex {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl Sub for Complex {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl Neg for Complex {
    type Output = Self;
    fn neg(self) -> Self::Output {
        Self::new(-self.re, -self.im)
    }
}

impl Mul for Complex {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl Div for Complex {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let denom = rhs.abs2();
        debug_assert!(denom != 0.0);
        self * rhs.conj() * (1.0 / denom)
    }
}

impl AddAssign for Complex {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for Complex {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl MulAssign for Complex {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl DivAssign<f64> for Complex {
    fn div_assign(&mut self, rhs: f64) {
        self.re /= rhs;
        self.im /= rhs;
    }
}

impl Mul<f64> for Complex {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.re * rhs, self.im * rhs)
    }
}

impl MulAssign<f64> for Complex {
    fn mul_assign(&mut self, rhs: f64) {
        self.re *= rhs;
        self.im *= rhs;
    }
}

fn bit_reverse(values: &mut [Complex]) {
    util::bit_reverse(values);
}

fn transform(values: &mut [Complex], invert: bool) {
    let n = values.len();
    assert!(n.is_power_of_two(), "FFT size must be a power of two");
    bit_reverse(values);

    let mut len = 2;
    while len <= n {
        let angle = 2.0 * PI / (len as f64) * if invert { -1.0 } else { 1.0 };
        let wlen = Complex::exp_i(angle);
        for chunk in values.chunks_mut(len) {
            let mut w = Complex::one();
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

    if invert {
        let inv_n = 1.0 / (n as f64);
        for v in values.iter_mut() {
            *v *= inv_n;
        }
    }
}

/// Forward complex FFT.
pub fn fft(values: &mut [Complex]) {
    transform(values, false);
}

/// Inverse complex FFT.
pub fn ifft(values: &mut [Complex]) {
    transform(values, true);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_fft() {
        let mut data = [
            Complex::new(1.0, 0.0),
            Complex::new(2.0, -1.0),
            Complex::new(3.0, 1.0),
            Complex::new(0.0, 0.0),
        ];
        let original = data;
        fft(&mut data);
        ifft(&mut data);
        for (a, b) in data.iter().zip(original.iter()) {
            assert!((a.re - b.re).abs() < 1e-9);
            assert!((a.im - b.im).abs() < 1e-9);
        }
    }
}
