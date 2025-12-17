# falcon-rs

Rust rewrite of the math layer from [tprest/falcon.py](https://github.com/tprest/falcon.py). The goal is to keep the reference structure—`ntt`, `fft`, `poly`, `mat`—while leaning on Rust traits for field operations and clear separation of concerns. Everything lives directly under `src/` in this crate.

## What’s included

- `field`: trait-based modular arithmetic over `Z/12289Z` (`ModQ`) inspired by the `dilithium_threshold` `field_element.rs` layout.
- `poly`: generic polynomial container with `Z/qZ` helpers, NTT-backed convolution, and Horner evaluation.
- `ntt`: radix-2 NTT/INTT for sizes up to 1024 using the Falcon modulus.
- `fft`: dependency-free complex FFT mirroring the reference floating-point routines.
- `mat`: tiny 2×2 matrix type with LDL decomposition for Gram matrix manipulations.
- `sampling`: deterministic test-friendly Box–Muller and Bernoulli utilities (not constant-time).

Higher-level Falcon keygen/sign/verify routines can build on these blocks while keeping close to the Python reference’s naming and flow.

## Usage

```rust
use falcon_rs::{fe, ModQ, Poly};

let a: Poly<ModQ> = poly!(fe!(1), fe!(2));
let b: Poly<ModQ> = poly!(fe!(3), fe!(4));
let c = a.ntt_product(&b); // convolution via NTT
assert_eq!(c.coeffs(), &[fe!(3), fe!(10), fe!(8)]);

// Using the helper macro for convenience
let p: Poly<ModQ> = poly!(1, 2, 3);
assert_eq!(p.coeffs(), &[fe!(1), fe!(2), fe!(3)]);
```

## Notes

- The sampling helpers are meant for testing and prototyping, not for production randomness or constant-time guarantees.
- Module boundaries intentionally mirror `falcon.py`, making it easier to compare against the reference codebase.
- Enable the `rand` feature to use a `rand_core`-backed RNG adapter (`RandCoreSource`); otherwise a deterministic xorshift RNG is provided for tests/examples.
