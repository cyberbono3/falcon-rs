use super::keygen::{Keypair, PublicKey};
use super::sign::{Signature, sign_with_rng};
use super::verify::{VerifyConfig, VerifyError, verify, verify_with_config};
use crate::math::field::modq_from_i64;
use crate::math::params::FALCON_Q;
use crate::math::params::ParameterSet;
use crate::math::sampling::XorShift64;
use crate::math::{ModQ, Poly};

fn end_to_end(degree: usize) {
    let logn = ParameterSet::new_falcon(degree)
        .expect("test only covers Falcon-512/1024")
        .params()
        .log_degree;
    let kp = Keypair::from_seed(logn, b"e2e-key").unwrap();
    let params = kp.public().params();
    let n = params.degree;
    assert_eq!(n, degree);

    let mut rng = XorShift64::new(42);
    let message = b"end-to-end message";
    let sig = sign_with_rng(kp.secret(), message, &mut rng).unwrap();
    assert_eq!(sig.s2().len(), n);

    verify_with_config(
        kp.public(),
        message,
        &sig,
        VerifyConfig {
            l2_bound: i128::MAX,
        },
    )
    .unwrap();

    assert!(matches!(
        verify_with_config(
            kp.public(),
            message,
            &sig,
            VerifyConfig { l2_bound: 0 },
        ),
        Err(VerifyError::NormTooLarge)
    ));

    let mismatched_params = match params.log_degree {
        9 => ParameterSet::Falcon1024.params(),
        10 => ParameterSet::Falcon512.params(),
        _ => unreachable!("logn is validated by keygen"),
    };
    let forged =
        Signature::new(*sig.nonce(), sig.s2().to_vec(), mismatched_params);
    assert!(matches!(
        verify_with_config(
            kp.public(),
            message,
            &forged,
            VerifyConfig {
                l2_bound: i128::MAX
            },
        ),
        Err(VerifyError::ParameterMismatch)
    ));

    let wrong_len =
        Signature::new(*sig.nonce(), vec![0i16; n.saturating_sub(1)], params);
    assert!(matches!(
        verify_with_config(
            kp.public(),
            message,
            &wrong_len,
            VerifyConfig {
                l2_bound: i128::MAX
            },
        ),
        Err(VerifyError::InvalidS2Length { .. })
    ));
}

#[test]
fn falcon512_end_to_end() {
    end_to_end(512);
}

#[test]
fn falcon1024_end_to_end() {
    end_to_end(1024);
}

#[test]
fn verify_rejects_unsupported_public_key_params() {
    let params = crate::math::params::Parameters {
        name: "Falcon-256 (unsupported)",
        degree: 256,
        log_degree: 8,
        modulus: FALCON_Q,
        max_sig_len: 0,
        pk_len: 0,
        sk_len: 0,
    };
    let mut h: Poly<ModQ> = Poly::from(vec![modq_from_i64(1)]);
    h.resize(params.degree);
    let pk = PublicKey::new(h, params);

    let sig = Signature::new(
        [0u8; super::sign::NONCE_LEN],
        vec![0i16; params.degree],
        params,
    );
    assert!(matches!(
        verify(&pk, b"msg", &sig),
        Err(VerifyError::UnsupportedLogN(8))
    ));
}

#[test]
fn verify_rejects_oversized_signature_under_default_bound_falcon512() {
    let kp = Keypair::from_seed(9, b"e2e-key").unwrap();
    let params = kp.public().params();
    let mut rng = XorShift64::new(42);
    let message = b"end-to-end message";
    let sig = sign_with_rng(kp.secret(), message, &mut rng).unwrap();

    let oversized =
        Signature::new(*sig.nonce(), vec![i16::MAX; params.degree], params);
    assert!(matches!(
        verify(kp.public(), message, &oversized),
        Err(VerifyError::NormTooLarge)
    ));
}
