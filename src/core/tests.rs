use super::keygen::{Keypair, PublicKey};
use super::sign::{
    HashToPoint, Signature, XorHashToPoint, hash_to_point, sign_with_rng,
};
use super::verify::{VerifyConfig, VerifyError, verify, verify_with_config};
use crate::math::field::modq_from_i64;
use crate::math::params::FALCON_Q;
use crate::math::params::ParameterSet;
use crate::math::sampling::{RandomSource, XorShift64};
use crate::math::{ModQ, Poly};

fn random_bytes(rng: &mut XorShift64, len: usize) -> Vec<u8> {
    let mut out = vec![0u8; len];
    let mut i = 0usize;
    while i < len {
        let w = rng.next_u64().to_le_bytes();
        let take = (len - i).min(w.len());
        out[i..i + take].copy_from_slice(&w[..take]);
        i += take;
    }
    out
}

fn random_nonce(rng: &mut XorShift64) -> [u8; super::sign::NONCE_LEN] {
    let mut out = [0u8; super::sign::NONCE_LEN];
    let mut i = 0usize;
    while i < out.len() {
        let w = rng.next_u64().to_le_bytes();
        let take = (out.len() - i).min(w.len());
        out[i..i + take].copy_from_slice(&w[..take]);
        i += take;
    }
    out
}

fn random_s2(rng: &mut XorShift64, len: usize) -> Vec<i16> {
    (0..len)
        .map(|_| (rng.next_u64() as i64 % 512 - 256) as i16)
        .collect()
}

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

#[test]
fn prop_hash_to_point_invariants_falcon512_and_falcon1024() {
    let mut rng = XorShift64::new(0xC0FFEE);
    for &degree in &[512usize, 1024usize] {
        let params = ParameterSet::new_falcon(degree).unwrap().params();
        let cases = if degree == 512 { 20 } else { 10 };
        for _ in 0..cases {
            let nonce = random_nonce(&mut rng);
            let msg_len = (rng.next_u64() % 200) as usize;
            let msg = random_bytes(&mut rng, msg_len);

            let c1 = hash_to_point(params, &nonce, &msg);
            let c2 = hash_to_point(params, &nonce, &msg);
            assert_eq!(c1, c2);
            assert_eq!(c1.len(), degree);
            assert!(c1.iter().all(|&x| x < FALCON_Q as u16));

            let h1 = XorHashToPoint::new(1).hash_to_point(params, &nonce, &msg);
            let h2 = XorHashToPoint::new(2).hash_to_point(params, &nonce, &msg);
            assert_ne!(h1, h2);
        }
    }
}

#[test]
fn prop_sign_verify_roundtrip_permissive_bound_falcon512_and_falcon1024() {
    for &degree in &[512usize, 1024usize] {
        let logn = ParameterSet::new_falcon(degree)
            .unwrap()
            .params()
            .log_degree;
        let kp = Keypair::from_seed(logn, b"prop-key").unwrap();
        let cases = if degree == 512 { 8 } else { 3 };

        for i in 0..cases {
            let mut rng = XorShift64::new(10_000 + i as u64);
            let msg = format!("property-message-{degree}-{i}").into_bytes();
            let sig = sign_with_rng(kp.secret(), &msg, &mut rng).unwrap();
            verify_with_config(
                kp.public(),
                &msg,
                &sig,
                VerifyConfig {
                    l2_bound: i128::MAX,
                },
            )
            .unwrap();
        }
    }
}

#[test]
fn prop_verify_rejects_invalid_s2_length() {
    let mut rng = XorShift64::new(0xBADC0DE);
    for &degree in &[512usize, 1024usize] {
        let logn = ParameterSet::new_falcon(degree)
            .unwrap()
            .params()
            .log_degree;
        let kp = Keypair::from_seed(logn, b"prop-key").unwrap();
        let params = kp.public().params();

        let cases = if degree == 512 { 8 } else { 3 };
        for _ in 0..cases {
            let nonce = random_nonce(&mut rng);
            let msg = random_bytes(&mut rng, 32);

            let mut bad_len = (rng.next_u64() as usize) % (degree + 10);
            if bad_len == degree {
                bad_len = degree - 1;
            }

            let sig =
                Signature::new(nonce, random_s2(&mut rng, bad_len), params);
            assert!(matches!(
                verify(kp.public(), &msg, &sig),
                Err(VerifyError::InvalidS2Length { .. })
            ));
        }
    }
}

#[test]
fn prop_verify_rejects_param_mismatch() {
    let mut rng = XorShift64::new(0xFACEFEED);
    for &(degree, other) in &[(512usize, 1024usize), (1024usize, 512usize)] {
        let logn = ParameterSet::new_falcon(degree)
            .unwrap()
            .params()
            .log_degree;
        let kp = Keypair::from_seed(logn, b"prop-key").unwrap();
        let good_params = kp.public().params();
        let bad_params = ParameterSet::new_falcon(other).unwrap().params();

        let nonce = random_nonce(&mut rng);
        let msg = random_bytes(&mut rng, 48);
        let sig = Signature::new(
            nonce,
            random_s2(&mut rng, good_params.degree),
            bad_params,
        );
        assert!(matches!(
            verify(kp.public(), &msg, &sig),
            Err(VerifyError::ParameterMismatch)
        ));
    }
}
