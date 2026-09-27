use crate::field::Field;
use crate::r1cs::QAP;
use crate::setup::{CRS, linear_combination};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Proof {
    pub a: Field,
    pub b: Field,
    pub c: Field,
}

pub fn prove(crs: &CRS, qap: &QAP, witness: &[Field], r: Field, s: Field) -> Proof {
    let num_pub = qap.num_pub;
    let num_vars = qap.num_vars;

    let mut u_sum = Field::ZERO;
    let mut v_sum = Field::ZERO;
    let mut c_sum = Field::ZERO;

    for i in 0..num_vars {
        let ui = linear_combination(&qap.u[i], &crs.tau_powers);
        let vi = linear_combination(&qap.v[i], &crs.tau_powers);
        u_sum = u_sum + ui * witness[i];
        v_sum = v_sum + vi * witness[i];

        if i >= num_pub {
            c_sum = c_sum + witness[i] * crs.private_terms[i - num_pub];
        }
    }

    let h_poly = qap.compute_h(witness);
    let h_term = linear_combination(&h_poly, &crs.h_powers);

    let a = crs.alpha + u_sum + r * crs.delta;
    let b = crs.beta + v_sum + s * crs.delta;
    let c = c_sum + h_term + a * s + r * b - r * s * crs.delta;

    Proof { a, b, c }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r1cs::R1CS;
    use crate::setup::setup;
    use crate::setup::test_helpers::toxic_waste;

    #[test]
    fn proof_works() {
        let r1cs = R1CS::circuit();
        let w = R1CS::build_witness(4);
        let domain = vec![Field::new(5), Field::new(7)];
        let qap = r1cs.to_qap(&domain);

        let (alpha, beta, gamma, delta, tau) = toxic_waste();
        let crs = setup(&qap, alpha, beta, gamma, delta, tau);

        let r = Field::new(3);
        let s = Field::new(6);
        let proof = prove(&crs, &qap, &w, r, s);

        assert_eq!(proof.a, Field::new(11));
        assert_eq!(proof.b, Field::new(9));
        assert_eq!(proof.c, Field::new(11));
    }

    #[test]
    fn shows_pi_works() {
        // ([11]_1, [11]_1, [9]_2)
        let r1cs = R1CS::circuit();
        let w = R1CS::build_witness(4);
        let domain = vec![Field::new(5), Field::new(7)];
        let qap = r1cs.to_qap(&domain);

        let (alpha, beta, gamma, delta, tau) = toxic_waste();
        let crs = setup(&qap, alpha, beta, gamma, delta, tau);
        let proof = prove(&crs, &qap, &w, Field::new(3), Field::new(6));

        let pi = (proof.a, proof.c, proof.b);
        println!("pi = ({}, {}, {})", pi.0, pi.1, pi.2);

        assert_eq!(pi, (Field::new(11), Field::new(11), Field::new(9)));
    }
}
