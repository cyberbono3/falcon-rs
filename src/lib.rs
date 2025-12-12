//! Falcon implementation scaffolding.
//!
//! For now, we expose parameter definitions; cryptographic primitives will
//! follow in subsequent modules.

pub mod params;

pub use params::{ParameterSet, Parameters, FALCON_Q, LOGQ};
