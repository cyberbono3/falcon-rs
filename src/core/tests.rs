use super::keygen::Keypair;
use super::sign::{Signature, sign_with_rng};
use super::verify::{VerifyConfig, VerifyError, verify_with_config};
use crate::math::params::ParameterSet;
use crate::math::sampling::XorShift64;

fn end_to_end(logn: u8) {
    let kp = Keypair::from_seed(logn, b"e2e-key").unwrap();
    let params = kp.public().params();
    let n = params.degree;

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
    end_to_end(9);
}

#[test]
fn falcon1024_end_to_end() {
    end_to_end(10);
}
