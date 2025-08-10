use std::cmp::Ordering;

use num::Integer;

/// A ratio of two integers corresponds to a point on the projective line
/// This ratio is used to parametrize a point on the unit circle c(u:t) = ((u²-t²)/(u²+t²), 2ut/(u²+t²))
#[derive(Debug, Clone)]
pub struct Ratio(pub i64, pub i64);

impl PartialEq for Ratio {
    fn eq(&self, other: &Self) -> bool {
        self.0 * other.1 - self.1 * other.0 == 0
    }
}

impl Eq for Ratio {}

impl Ord for Ratio {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.eq(other) {
            return Ordering::Equal;
        } else if self.0 == 0 && other.0 * other.1 >= 0 {
            return Ordering::Greater;
        } else if self.0 == 0 && other.0 * other.1 < 0 {
            return Ordering::Less;
        } else if other.0 == 0 && self.0 * self.1 >= 0 {
            return Ordering::Less;
        } else if other.0 == 0 && self.0 * self.1 < 0 {
            return Ordering::Greater;
        } else {
            let q1 = self.1 as f64 / self.0 as f64;
            let q2 = other.1 as f64 / other.0 as f64;
            // Safety: division by 0 is not possible since it is already checked in earlier branches
            return q1.partial_cmp(&q2).expect("Can neither be NaN nor ∞");
        }
    }
}

impl PartialOrd for Ratio {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ratio {
    pub fn conjugate(&self) -> Self {
        Self(self.0, -self.1)
    }

    pub fn zero() -> Self {
        Self(1, 0)
    }

    pub fn is_zero(&self) -> bool {
        self.eq(&Self::zero())
    }

    pub fn circular_mul(&self, other: &Self) -> Self {
        Self(
            self.0 * other.0 - self.1 * other.1,
            self.0 * other.1 + self.1 * other.0,
        )
    }

    pub fn circular_pow(&self, mut n: i32) -> Self {
        let mut p = if n >= 0 {
            self.clone()
        } else {
            self.clone().conjugate()
        };

        if n == 0 {
            return Ratio(1, 0);
        }

        if n > 0 {
            while n > 1 {
                p = p.circular_mul(self);
                n -= 1;
            }
        } else {
            let c = self.conjugate();
            while n < -1 {
                p = p.circular_mul(&c);
                n += 1;
            }
        }

        p
    }

    /// The remainder when self is divided by other
    /// This is equivalent to the rotation needed to rotate from other to self
    pub fn circular_rem(&self, other: &Self) -> Self {
        self.circular_mul(&other.conjugate())
    }

    /// Given two ratios (u:t) and (x:y), finds the smallest natural number n such that the circular remainder (v:w)
    /// of the ratio (r:s) = (u:t) ⊗ (u:t) ⊗ ... ⊗ (u:t) (n-times) with respect to (x:y)
    /// satisfies (1:0) <= (v:w) < (u:t)
    /// Interpretation: If (u:t) parametrizes a point on the circle c(u:t) = ((u²-t²)/(u²+t²), (2ut)/(u²+t²)), then n is the smallest possible
    /// number of multiplications such that c(r:s) does not cross or is equal to the point c(x:y)
    /// Convention is that if `self.is_zero()`, then self is considered to be a full turn and not zero
    pub fn circular_div(&self, other: &Self) -> u32 {
        if other > self && !self.is_zero() {
            return 0;
        }

        let mut c = other.clone();
        let mut rem = self.circular_rem(&c);
        let mut n: u32 = 1;
        let o = Ratio(1, 0);
        while !(o <= rem && rem < *other) {
            c = c.circular_mul(&other);
            c.reduce();
            rem = self.circular_rem(&c);
            n += 1;
        }

        n
    }

    pub fn reduce(&mut self) {
        if self.0 == 0 {
            self.1 = 1;
            return;
        }

        if self.1 == 0 {
            self.0 = 1;
            return;
        }

        if self.0 == self.1 {
            self.0 = 1;
            self.1 = 1;
            return;
        }

        let g = self.0.gcd(&self.1);

        self.0 /= g;
        self.1 /= g;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ord_1() {
        let z = Ratio::zero();
        let r1 = Ratio(0, 1);
        let r2 = Ratio(1, 1);
        let r3 = Ratio(3, 10);
        let r4 = Ratio(6, 20);
        let r5 = Ratio(1, -1);
        assert_eq!(r1.cmp(&r2), Ordering::Greater);
        assert_eq!(r2.cmp(&r3), Ordering::Less);
        assert_eq!(r3.cmp(&r4), Ordering::Equal);
        assert_eq!(r1.cmp(&r5), Ordering::Less);
        assert_eq!(r1.cmp(&r2), Ordering::Greater);
        assert_eq!(r1.cmp(&z), Ordering::Greater);
        assert_eq!(r2.cmp(&z), Ordering::Greater);
    }

    #[test]
    fn test_circular_pow_1() {
        let r1 = Ratio(1, 1);
        let p1 = r1.circular_pow(-2);
        assert_eq!(p1, Ratio(0, 1));
    }

    #[test]
    fn test_circular_rem_1() {
        let o = Ratio(1, 0);
        let r1 = Ratio(1, 1);
        let rem1 = o.circular_rem(&r1);
        assert_eq!(rem1, Ratio(1, -1));

        let r2 = Ratio(0, 1);
        let rem2 = o.circular_rem(&r2);
        assert_eq!(rem2, Ratio(0, 1));

        let rem3 = r2.circular_rem(&r1);
        assert_eq!(rem3, Ratio(1, 1));
    }

    #[test]
    fn test_circular_div_1() {
        let o = Ratio(1, 0);
        let r1 = Ratio(3, 4);
        let n1 = o.circular_div(&r1);
        assert_eq!(n1, 3);

        let r2 = Ratio(1, 2);
        let n2 = o.circular_div(&r2);
        assert_eq!(n2, 2);

        let r3 = Ratio(0, 1);
        let n3 = o.circular_div(&r3);
        assert_eq!(n3, 2);

        let r4 = Ratio(1, 1);
        let n4 = o.circular_div(&r4);
        assert_eq!(n4, 4);

        let r5 = Ratio(100, 63);
        let r6 = Ratio(100, 9);
        let n5 = r5.circular_div(&r6);
        assert_eq!(n5, 6);

        let n6 = r6.circular_div(&r5);
        assert_eq!(n6, 0);

        let r7 = Ratio(1, -1);
        let n7 = r7.circular_div(&r3);
        assert_eq!(n7, 1);
    }

    #[test]
    fn test_reduce_1() {
        let mut r1 = Ratio(125, 25);
        r1.reduce();
        assert_eq!(r1, Ratio(5, 1));
    }
}
