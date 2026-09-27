use crate::field::Field;
use crate::poly::Poly;
use crate::r1cs::QAP;

#[derive(Debug, Clone)]
pub struct CRS {
    pub alpha: Field,
    pub beta: Field,
    pub gamma: Field,
    pub delta: Field,

    pub tau_powers: Vec<Field>,
    pub public_terms: Vec<Field>,
    pub private_terms: Vec<Field>,
    pub h_powers: Vec<Field>,
}

/// Evaluates a polynomial "in the exponent" by dotting its coefficients
/// against a matching list of precomputed powers/terms -- shared by
/// setup() (building the CRS) and prove() (using it).
pub fn linear_combination(poly: &Poly, terms: &[Field]) -> Field {
    let mut result = Field::ZERO;
    for (i, c) in poly.coeffs().iter().enumerate() {
        if i < terms.len() {
            result = result + *c * terms[i];
        }
    }
    result
}

pub fn setup(qap: &QAP, alpha: Field, beta: Field, gamma: Field, delta: Field, tau: Field) -> CRS {
    let n = qap.t.degree();

    let mut tau_powers = Vec::with_capacity(n);
    let mut pow = Field::ONE;
    for _ in 0..n {
        tau_powers.push(pow);
        pow = pow * tau;
    }

    let t_at_tau = qap.t.eval(tau);

    let mut public_terms = Vec::with_capacity(qap.num_pub);
    for i in 0..qap.num_pub {
        let ui = qap.u[i].eval(tau);
        let vi = qap.v[i].eval(tau);
        let wi = qap.w[i].eval(tau);
        let combined = (beta * ui) + (alpha * vi) + wi;
        public_terms.push(combined / gamma);
    }

    let mut private_terms = Vec::with_capacity(qap.num_vars - qap.num_pub);
    for i in qap.num_pub..qap.num_vars {
        let ui = qap.u[i].eval(tau);
        let vi = qap.v[i].eval(tau);
        let wi = qap.w[i].eval(tau);
        let combined = (beta * ui) + (alpha * vi) + wi;
        private_terms.push(combined / delta);
    }

    let h_len = if n >= 1 { n - 1 } else { 0 };
    let mut h_powers = Vec::with_capacity(h_len);
    let mut pow_tau = Field::ONE;
    for _ in 0..h_len {
        h_powers.push(pow_tau * t_at_tau / delta);
        pow_tau = pow_tau * tau;
    }

    CRS {
        alpha,
        beta,
        gamma,
        delta,
        tau_powers,
        public_terms,
        private_terms,
        h_powers,
    }
}

#[cfg(test)]
pub(crate) mod test_helpers {
    use crate::field::Field;

    pub fn toxic_waste() -> (Field, Field, Field, Field, Field) {
        // alpha=2, beta=4, gamma=8, delta=9, tau=10(x)
        (
            Field::new(2),
            Field::new(4),
            Field::new(8),
            Field::new(9),
            Field::new(10),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::r1cs::R1CS;
    use test_helpers::toxic_waste;

    #[test]
    fn tau_powers_and_h_powers_have_expected_lengths() {
        let r1cs = R1CS::circuit();
        let domain = vec![Field::new(5), Field::new(7)];
        let qap = r1cs.to_qap(&domain);
        let (alpha, beta, gamma, delta, tau) = toxic_waste();
        let crs = setup(&qap, alpha, beta, gamma, delta, tau);

        assert_eq!(crs.tau_powers.len(), 2);
        assert_eq!(crs.h_powers.len(), 1);
        assert_eq!(crs.public_terms.len(), 2);
        assert_eq!(crs.private_terms.len(), 2);
    }

    #[test]
    fn h_powers_first_term() {
        let r1cs = R1CS::circuit();
        let domain = vec![Field::new(5), Field::new(7)];
        let qap = r1cs.to_qap(&domain);
        let (alpha, beta, gamma, delta, tau) = toxic_waste();
        let crs = setup(&qap, alpha, beta, gamma, delta, tau);

        assert_eq!(crs.h_powers[0], Field::new(6));
    }
}
