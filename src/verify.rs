use crate::field::Field;
use crate::proof::Proof;
use crate::r1cs::QAP;
use crate::setup::CRS;

pub fn verify(crs: &CRS, qap: &QAP, public_inputs: &[Field], proof: &Proof) -> bool {
    let (lhs, rhs) = verify_verbose(crs, qap, public_inputs, proof);
    lhs == rhs
}

pub fn verify_verbose(
    crs: &CRS,
    qap: &QAP,
    public_inputs: &[Field],
    proof: &Proof,
) -> (Field, Field) {
    assert_eq!(public_inputs.len(), qap.num_pub);

    let mut public_sum = Field::ZERO;
    for i in 0..qap.num_pub {
        public_sum = public_sum + public_inputs[i] * crs.public_terms[i];
    }

    let lhs = proof.a * proof.b;
    let rhs = crs.alpha * crs.beta + public_sum * crs.gamma + proof.c * crs.delta;

    (lhs, rhs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::proof::prove;
    use crate::r1cs::R1CS;
    use crate::setup::setup;
    use crate::setup::test_helpers::toxic_waste;

    fn build_proof_for_x(x: u64) -> (CRS, QAP, Proof, Vec<Field>) {
        let r1cs = R1CS::circuit();
        let w = R1CS::build_witness(x);
        let domain = vec![Field::new(5), Field::new(7)];
        let qap = r1cs.to_qap(&domain);

        let (alpha, beta, gamma, delta, tau) = toxic_waste();
        let crs = setup(&qap, alpha, beta, gamma, delta, tau);
        let proof = prove(&crs, &qap, &w, Field::new(3), Field::new(6));

        let public_inputs = vec![w[0], w[1]];
        (crs, qap, proof, public_inputs)
    }

    #[test]
    fn verification_matches_8_equals_8() {
        let (crs, qap, proof, public_inputs) = build_proof_for_x(4);
        assert!(verify(&crs, &qap, &public_inputs, &proof));
    }

    #[test]
    fn show_lhs_equals_rhs_directly() {
        let (crs, qap, proof, public_inputs) = build_proof_for_x(4);
        let (lhs, rhs) = verify_verbose(&crs, &qap, &public_inputs, &proof);

        println!("LHS = {}, RHS = {}", lhs, rhs);
        assert_eq!(lhs, rhs);
        assert_eq!(lhs, Field::new(8));
    }

    #[test]
    fn rejects_wrong_public_input() {
        let (crs, qap, proof, _) = build_proof_for_x(4);
        let wrong_public = vec![Field::new(1), Field::new(4)];
        assert!(!verify(&crs, &qap, &wrong_public, &proof));
    }

    #[test]
    fn works_end_to_end_for_other_witnesses() {
        for x in [0u64, 1, 2, 6, 11] {
            let (crs, qap, proof, public_inputs) = build_proof_for_x(x);
            assert!(
                verify(&crs, &qap, &public_inputs, &proof),
                "failed for x={}",
                x
            );
        }
    }
}
