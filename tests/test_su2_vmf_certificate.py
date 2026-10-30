"""Independent rational checks for the SU(2) vMF interval certificate."""

import unittest
from fractions import Fraction
from math import factorial

from experiments.su2_vmf_certificate import (
    Interval,
    block_certificate,
    coefficient,
    derivative_enclosure,
    derivative_tables,
    evaluate,
)


def contains(interval, value):
    value = Fraction(value)
    return Fraction(interval.lo, 1 << 120) <= value <= Fraction(interval.hi, 1 << 120)


class IntervalArithmeticTests(unittest.TestCase):
    def test_signed_operations_contain_exact_rational_samples(self):
        samples = [
            (Fraction(-7, 3), Fraction(-2, 5)),
            (Fraction(-5, 4), Fraction(3, 7)),
            (Fraction(1, 9), Fraction(11, 6)),
        ]
        for a, b in samples:
            ia, ib = Interval.rational(a), Interval.rational(b)
            self.assertTrue(contains(ia + ib, a + b))
            self.assertTrue(contains(ia - ib, a - b))
            self.assertTrue(contains(ia * ib, a * b))
            if b:
                self.assertTrue(contains(ia / ib, a / b))

    def test_negative_denominator_reciprocal_is_ordered_and_outward(self):
        denominator = Interval.hull(Fraction(-7, 3), Fraction(-5, 4))
        quotient = Interval.rational(2) / denominator
        for d in (Fraction(-7, 3), Fraction(-9, 5), Fraction(-5, 4)):
            self.assertTrue(contains(quotient, 2 / d))
        with self.assertRaises(ValueError):
            Interval.rational(1) / Interval.hull(-1, 1)


class BesselSeriesTests(unittest.TestCase):
    def test_derivative_coefficients_from_independent_factorial_formula(self):
        # Differentiate a_n s^n directly, where
        # a_n = 1 / (4^n n! (n+1)!).
        for j in range(4):
            for m in range(8):
                n = m + j
                direct = Fraction(factorial(n), factorial(m)) * Fraction(
                    1, 4**n * factorial(n) * factorial(n + 1)
                )
                self.assertEqual(coefficient(j, m), direct)

    def test_geometric_remainder_contains_long_exact_series(self):
        radius = Fraction(21, 10)
        tables = derivative_tables(radius, 24)
        points = [Fraction(0), Fraction(1, 8), radius * radius]
        for j, table in enumerate(tables):
            for s in points:
                enclosure = evaluate(table, Interval.rational(s))
                # Exact partial sum plus an independent geometric upper bound
                # on all terms after degree 100 encloses the positive series.
                partial = sum((coefficient(j, m) * s**m for m in range(101)), Fraction(0))
                first = coefficient(j, 101) * s**101
                ratio = s / (4 * 102 * (101 + j + 2))
                upper = partial + first / (1 - ratio)
                self.assertLessEqual(Fraction(enclosure.lo, 1 << 120), partial)
                self.assertGreaterEqual(Fraction(enclosure.hi, 1 << 120), upper)


class CertificateDomainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.enclosure = derivative_enclosure(cells=4096)

    def test_invalid_series_domain_and_coupling_domain_fail_closed(self):
        with self.assertRaises(ValueError):
            derivative_tables(Fraction(10), 1)
        with self.assertRaises(ValueError):
            derivative_enclosure(radius=0)
        with self.assertRaises(ValueError):
            block_certificate(Fraction(351, 1000), self.enclosure)

    def test_source_normalization_and_monotone_baseline_threshold(self):
        # The pinned source uses beta_source = beta_W/4. These are the exact
        # substitutions in w, C, K, and alpha, independently of the JSON.
        row = block_certificate(Fraction(69, 200), self.enclosure, improved=False)
        exact = row["exact"]
        beta_w = Fraction(69, 200)
        beta_source = beta_w / 4
        tau = Fraction(1, 4)
        lam = Fraction(self.enclosure["lambda_lower"])
        third = Fraction(self.enclosure["third_upper"])
        w = 8 * beta_source**2 * (tau - lam)
        self.assertEqual(Fraction(exact["w_upper"]), w)
        # The baseline branch uses the source's sqrt(tau) <= 1/2 and M factor 1.
        # Check s and C using the returned exponential upper bound.
        from experiments.su2_vmf_certificate import exp_upper
        ew = exp_upper(w)
        s = tau + (ew - 1) * ew / 2
        c = 20 * beta_source * (4 * beta_source * tau * s + 16 * beta_source**2 * third)
        k = (1 + beta_source) * s / (1 - c)
        self.assertEqual(Fraction(exact["s_upper"]), s)
        self.assertEqual(Fraction(exact["C_upper"]), c)
        self.assertEqual(Fraction(exact["K_upper"]), k)
        self.assertEqual(Fraction(exact["alpha_upper"]), 40 * beta_source * k)

        below = block_certificate(Fraction(69, 200), self.enclosure, improved=False)
        above = block_certificate(Fraction(87, 250), self.enclosure, improved=False)
        self.assertTrue(below["certified"])
        self.assertFalse(above["certified"])
        self.assertLess(
            Fraction(below["exact"]["alpha_upper"]),
            Fraction(above["exact"]["alpha_upper"]),
        )

    def test_poincare_sharpening_matches_its_exact_rational_inputs(self):
        from experiments.su2_vmf_certificate import exp_upper

        beta_w = Fraction(87, 250)
        beta_source = beta_w / 4
        tau = Fraction(1, 4)
        lam = Fraction(self.enclosure["lambda_lower"])
        third = Fraction(self.enclosure["third_upper"])
        w = 8 * beta_source**2 * (tau - lam)
        ew = exp_upper(w)
        # 289/1000 > sqrt(tau/3), and Holley--Stroock plus gap 3
        # gives the residual covariance factor e^w/3.
        s = tau + (ew - 1) * ew * Fraction(289, 1000)
        c = 20 * beta_source * (
            4 * beta_source * tau * s + 16 * beta_source**2 * third * ew / 3
        )
        k = (1 + beta_source) * s / (1 - c)
        row = block_certificate(beta_w, self.enclosure, improved=True)
        exact = row["exact"]
        self.assertEqual(Fraction(exact["w_upper"]), w)
        self.assertEqual(Fraction(exact["s_upper"]), s)
        self.assertEqual(Fraction(exact["C_upper"]), c)
        self.assertEqual(Fraction(exact["K_upper"]), k)
        self.assertEqual(Fraction(exact["alpha_upper"]), 40 * beta_source * k)
        self.assertTrue(row["certified"])

        # It is a pointwise certificate with a monotone formula: the next
        # tested coupling lies beyond the certified alpha<1 range.
        next_row = block_certificate(Fraction(7, 20), self.enclosure, improved=True)
        self.assertFalse(next_row["certified"])
        self.assertLess(
            Fraction(exact["alpha_upper"]),
            Fraction(next_row["exact"]["alpha_upper"]),
        )


if __name__ == "__main__":
    unittest.main()
