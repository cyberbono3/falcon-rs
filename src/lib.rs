//! Falcon implementation scaffolding.
//!
//! The crate mirrors the layout of `tprest/falcon.py` while leaning on Rust
//! traits for modular arithmetic and transform-heavy helpers. Higher-level
//! signing and verification routines can be built on top of these math
//! building blocks.

pub mod core;
pub mod math;
