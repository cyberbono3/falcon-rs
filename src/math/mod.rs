#[cfg(feature = "ark-compat")]
pub mod ark;
pub mod fft;
pub mod field;
pub mod mat;
pub mod ntt;
pub mod params;
pub mod poly;
pub mod sampling;
pub(crate) mod util;

#[cfg(feature = "ark-compat")]
pub use ark::ArkModQ;
pub use ark_ff::PrimeField;
pub use fft::Complex;
pub use field::ModQ;
pub use mat::Mat2;
pub use ntt::{MAX_LOGN, MAX_N, PRIMITIVE_NTT_ROOT_VALUE, primitive_ntt_root};
pub use params::{FALCON_Q, LOGQ, ParameterSet, Parameters};
pub use poly::Poly;
