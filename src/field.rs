use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Field(u64);

const MOD: u64 = 13;

impl Field {
    pub const ONE: Field = Field(1);
    pub const ZERO: Field = Field(0);

    pub fn new(value: u64) -> Self {
        Field(value % MOD)
    }

    pub fn is_zero(&self) -> bool {
        self.0 == 0
    }

    /// Fermat's little theorem: a^(p-2) = a^-1 (mod p).
    pub fn inv(&self) -> Self {
        assert!(!self.is_zero(), "cannot invert zero");
        self.pow(MOD - 2)
    }

    pub fn pow(&self, mut exp: u64) -> Self {
        let mut base = *self;
        let mut result = Field::ONE;
        while exp > 0 {
            if exp & 1 == 1 {
                result = result * base;
            }
            base = base * base;
            exp >>= 1;
        }
        result
    }
}

impl Add for Field {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let s = self.0 + rhs.0;
        Field(if s >= MOD { s - MOD } else { s })
    }
}

impl AddAssign for Field {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for Field {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        if self.0 >= rhs.0 {
            Field(self.0 - rhs.0)
        } else {
            Field(self.0 + MOD - rhs.0)
        }
    }
}

impl SubAssign for Field {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl Mul for Field {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Field((self.0 * rhs.0) % MOD)
    }
}

impl MulAssign for Field {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Neg for Field {
    type Output = Self;
    fn neg(self) -> Self {
        if self.is_zero() {
            Field::ZERO
        } else {
            Field(MOD - self.0)
        }
    }
}

impl Div for Field {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        self * rhs.inv()
    }
}

impl fmt::Display for Field {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_wraps_at_modulus() {
        assert_eq!(Field::new(9) + Field::new(9), Field::new(5)); // 18 mod 13
    }

    #[test]
    fn neg_and_sub_are_consistent() {
        let a = Field::new(3);
        assert_eq!(a + (-a), Field::ZERO);
        assert_eq!(Field::new(3) - Field::new(5), Field::new(11)); // -2 mod 13
    }
}
