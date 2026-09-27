use crate::field::Field;
use std::fmt;
use std::ops::{Add, AddAssign, Mul, Sub};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Poly(pub Vec<Field>);

impl Poly {
    pub fn is_zero(&self) -> bool {
        self.degree() == 0 && self.0[0].is_zero()
    }

    pub fn zero() -> Self {
        Poly(vec![Field::ZERO])
    }

    pub fn one() -> Self {
        Poly(vec![Field::ONE])
    }

    pub fn constant(c: Field) -> Self {
        Poly(vec![c])
    }

    pub fn coeffs(&self) -> &[Field] {
        &self.0
    }

    pub fn degree(&self) -> usize {
        if self.0.is_empty() {
            return 0;
        }
        let mut d = self.0.len() - 1;
        while d > 0 && self.0[d].is_zero() {
            d -= 1;
        }
        d
    }

    pub fn trim(mut self) -> Self {
        let d = self.degree();
        self.0.truncate(d + 1);
        self
    }

    /// Horner's method: ((c_n*x + c_{n-1})*x + ... )*x + c0.
    pub fn eval(&self, x: Field) -> Field {
        let mut result = Field::ZERO;
        for coeff in self.0.iter().rev() {
            result = result * x + *coeff;
        }
        result
    }

    /// Standard polynomial long division: returns (quotient, remainder).
    pub fn poly_div(&self, divisor: &Poly) -> (Poly, Poly) {
        let a_deg = self.degree();
        let b_deg = divisor.degree();

        if a_deg < b_deg || self.is_zero() {
            return (Poly::zero(), self.clone());
        }

        let mut quotient = vec![Field::ZERO; a_deg - b_deg + 1];
        let mut remainder = self.0.clone();
        remainder.resize(a_deg + 1, Field::ZERO);

        let leading_inv = divisor.0[b_deg].inv();

        for i in (0..=a_deg - b_deg).rev() {
            let lead_idx = i + b_deg;
            let lead = remainder[lead_idx];
            if lead.is_zero() {
                continue;
            }
            let factor = lead * leading_inv;
            quotient[i] = factor;
            for j in 0..=b_deg {
                remainder[i + j] = remainder[i + j] - factor * divisor.0[j];
            }
        }

        (Poly(quotient).trim(), Poly(remainder).trim())
    }

    /// Lagrange interpolation through the given (x, y) points.
    pub fn lagrange_interpolate(points: &[(Field, Field)]) -> Self {
        let n = points.len();
        let mut result = Poly::zero();

        for i in 0..n {
            let (xi, yi) = points[i];
            let mut numerator = Poly::one();
            let mut denominator = Field::ONE;

            for j in 0..n {
                if i == j {
                    continue;
                }
                let (xj, _) = points[j];
                numerator = numerator * Poly(vec![-xj, Field::ONE]);
                denominator = denominator * (xi - xj);
            }

            let term = numerator * Poly::constant(yi * denominator.inv());
            result = result + term;
        }

        result.trim()
    }
}

impl Add for Poly {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let len = std::cmp::max(self.0.len(), rhs.0.len());
        let mut coeffs = vec![Field::ZERO; len];
        for (i, c) in self.0.iter().enumerate() {
            coeffs[i] = coeffs[i] + *c;
        }
        for (i, c) in rhs.0.iter().enumerate() {
            coeffs[i] = coeffs[i] + *c;
        }
        Poly(coeffs).trim()
    }
}

impl AddAssign for Poly {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}

impl Sub for Poly {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let len = std::cmp::max(self.0.len(), rhs.0.len());
        let mut coeffs = vec![Field::ZERO; len];
        for (i, c) in self.0.iter().enumerate() {
            coeffs[i] = coeffs[i] + *c;
        }
        for (i, c) in rhs.0.iter().enumerate() {
            coeffs[i] = coeffs[i] - *c;
        }
        Poly(coeffs).trim()
    }
}

impl Mul for Poly {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let len = self.0.len() + rhs.0.len();
        let mut coeffs = vec![Field::ZERO; len];
        for (i, a) in self.0.iter().enumerate() {
            for (j, b) in rhs.0.iter().enumerate() {
                coeffs[i + j] = coeffs[i + j] + *a * *b;
            }
        }
        Poly(coeffs).trim()
    }
}

impl Mul<Field> for Poly {
    type Output = Self;
    fn mul(self, rhs: Field) -> Self {
        Poly(self.0.iter().map(|c| *c * rhs).collect()).trim()
    }
}

impl fmt::Display for Poly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let d = self.degree();
        if d == 0 {
            return write!(f, "{}", self.0[0]);
        }
        let mut first = true;
        for i in (0..=d).rev() {
            let c = self.0[i];
            if c.is_zero() {
                continue;
            }
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            if i == 0 {
                write!(f, "{}", c)?;
            } else {
                if c != Field::ONE {
                    write!(f, "{}", c)?;
                }
                write!(f, "x")?;
                if i > 1 {
                    write!(f, "^{}", i)?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eval_uses_horners_method_correctly() {
        let p = Poly(vec![Field::new(1), Field::new(2)]);
        assert_eq!(p.eval(Field::new(3)), Field::new(7));
    }

    #[test]
    fn lagrange_reproduces_paper_basis_polys() {
        // L2(x) through (5,0),(7,1) (7x + 4 mod 13) should give L2(0)=4, L2(5)=0, L2(7)=1.
        let l2 = Poly::lagrange_interpolate(&[
            (Field::new(5), Field::ZERO),
            (Field::new(7), Field::ONE),
        ]);
        assert_eq!(l2.eval(Field::new(0)), Field::new(4));
        assert_eq!(l2.0[1], Field::new(7));
    }

    #[test]
    fn division_recovers_exact_quotient() {
        // (x-5)(x-7) divided by (x-5) should give exactly (x-7).
        let a = Poly(vec![-Field::new(5), Field::ONE]);
        let b = Poly(vec![-Field::new(7), Field::ONE]);
        let product = a.clone() * b.clone();
        let (q, r) = product.poly_div(&a);
        assert_eq!(q, b);
        assert!(r.is_zero());
    }
}
