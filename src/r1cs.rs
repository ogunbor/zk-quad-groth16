use crate::field::Field;
use crate::poly::Poly;

#[derive(Debug, Clone)]
pub struct R1CS {
    pub num_vars: usize,
    pub num_pub: usize,
    pub num_constraints: usize,
    pub a: Vec<Vec<Field>>,
    pub b: Vec<Vec<Field>>,
    pub c: Vec<Vec<Field>>,
}

#[derive(Debug, Clone)]
pub struct QAP {
    pub num_vars: usize,
    pub num_pub: usize,
    pub u: Vec<Poly>,
    pub v: Vec<Poly>,
    pub w: Vec<Poly>,
    pub t: Poly,
}

impl R1CS {
    pub fn new(num_vars: usize, num_pub: usize) -> Self {
        R1CS {
            num_vars,
            num_pub,
            num_constraints: 0,
            a: Vec::new(),
            b: Vec::new(),
            c: Vec::new(),
        }
    }

    pub fn add_constraint(&mut self, a: Vec<Field>, b: Vec<Field>, c: Vec<Field>) {
        assert_eq!(a.len(), self.num_vars);
        assert_eq!(b.len(), self.num_vars);
        assert_eq!(c.len(), self.num_vars);
        self.a.push(a);
        self.b.push(b);
        self.c.push(c);
        self.num_constraints += 1;
    }

    /// The circuit: y = x^2 + x + 9.
    /// Witness order: [one, y, x, s1], s1 = x*x.
    pub fn circuit() -> Self {
        let mut r1cs = R1CS::new(4, 2); // 4 vars, first 2 (one, y) public
        r1cs.add_constraint(
            vec![Field::new(0), Field::new(0), Field::new(1), Field::new(0)],
            vec![Field::new(0), Field::new(0), Field::new(1), Field::new(0)],
            vec![Field::new(0), Field::new(0), Field::new(0), Field::new(1)],
        );
        r1cs.add_constraint(
            vec![Field::new(9), Field::new(0), Field::new(1), Field::new(1)],
            vec![Field::new(1), Field::new(0), Field::new(0), Field::new(0)],
            vec![Field::new(0), Field::new(1), Field::new(0), Field::new(0)],
        );
        r1cs
    }

    /// Checks (A_i . w) * (B_i . w) == (C_i . w) for every constraint row.
    pub fn is_satisfied(&self, w: &[Field]) -> bool {
        for i in 0..self.num_constraints {
            let av = dot(&self.a[i], w);
            let bv = dot(&self.b[i], w);
            let cv = dot(&self.c[i], w);
            if av * bv != cv {
                return false;
            }
        }
        true
    }

    /// Builds [one, y, x, s1] for a given x, computing y = x^2 + x + 9.
    pub fn build_witness(x: u64) -> Vec<Field> {
        let one = Field::ONE;
        let x = Field::new(x);
        let s1 = x * x;
        let y = s1 + x + Field::new(9);
        vec![one, y, x, s1]
    }

    /// Converts this R1CS into a QAP over the given evaluation domain
    /// (one point per constraint row). randoms used: {5, 7}.
    pub fn to_qap(&self, domain: &[Field]) -> QAP {
        assert_eq!(domain.len(), self.num_constraints);
        let n = self.num_constraints;
        let m = self.num_vars;

        let mut u = Vec::with_capacity(m);
        let mut v = Vec::with_capacity(m);
        let mut w = Vec::with_capacity(m);

        for var_idx in 0..m {
            let a_pts: Vec<(Field, Field)> =
                (0..n).map(|i| (domain[i], self.a[i][var_idx])).collect();
            let b_pts: Vec<(Field, Field)> =
                (0..n).map(|i| (domain[i], self.b[i][var_idx])).collect();
            let c_pts: Vec<(Field, Field)> =
                (0..n).map(|i| (domain[i], self.c[i][var_idx])).collect();

            u.push(Poly::lagrange_interpolate(&a_pts));
            v.push(Poly::lagrange_interpolate(&b_pts));
            w.push(Poly::lagrange_interpolate(&c_pts));
        }

        let mut t = Poly::one();
        for &d in domain {
            t = t * Poly(vec![-d, Field::ONE]);
        }

        QAP {
            num_vars: m,
            num_pub: self.num_pub,
            u,
            v,
            w,
            t,
        }
    }
}

impl QAP {
    /// h(x) = (A(x)*B(x) - C(x)) / t(x), where A/B/C are the witness-weighted
    /// combinations of u/v/w. Panics if the remainder isn't zero, i.e. if the
    /// witness doesn't actually satisfy the circuit.
    pub fn compute_h(&self, witness: &[Field]) -> Poly {
        let mut u_sum = Poly::zero();
        let mut v_sum = Poly::zero();
        let mut w_sum = Poly::zero();

        for i in 0..self.num_vars {
            u_sum = u_sum + self.u[i].clone() * Poly::constant(witness[i]);
            v_sum = v_sum + self.v[i].clone() * Poly::constant(witness[i]);
            w_sum = w_sum + self.w[i].clone() * Poly::constant(witness[i]);
        }

        let lhs = u_sum * v_sum;
        let numerator = lhs - w_sum;
        let (h, remainder) = numerator.poly_div(&self.t);

        assert!(
            remainder.is_zero(),
            "QAP divisibility check failed: remainder = {}",
            remainder
        );
        h
    }
}

fn dot(row: &[Field], w: &[Field]) -> Field {
    let mut acc = Field::ZERO;
    for (c, wi) in row.iter().zip(w.iter()) {
        acc = acc + *c * *wi;
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn witness_for_x_4_matches_paper_and_satisfies_circuit() {
        let r1cs = R1CS::circuit();
        let w = R1CS::build_witness(4);

        assert_eq!(
            w,
            vec![Field::new(1), Field::new(3), Field::new(4), Field::new(3)]
        );
        assert!(r1cs.is_satisfied(&w));
    }

    #[test]
    fn qap_matches_on_abc_and_h() {
        let r1cs = R1CS::circuit();
        let w = R1CS::build_witness(4);
        let domain = vec![Field::new(5), Field::new(7)];
        let qap = r1cs.to_qap(&domain);

        let h = qap.compute_h(&w);
        // h(x) = 4 (constant).
        assert_eq!(h.degree(), 0);
        assert_eq!(h.0[0], Field::new(4));
    }

    #[test]
    fn sweep_of_x_values_all_satisfy_and_have_valid_qap() {
        let r1cs = R1CS::circuit();
        let domain = vec![Field::new(5), Field::new(7)];
        let qap = r1cs.to_qap(&domain);
        for x in 0..13 {
            let w = R1CS::build_witness(x);
            assert!(r1cs.is_satisfied(&w), "x={} should satisfy", x);
            let _h = qap.compute_h(&w); // panics on bad witness
        }
    }
}
