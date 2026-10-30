"""Focused checks for the rational one-color SU(2) decimation certificate."""

from fractions import Fraction as F
import unittest

from experiments.su2_decimation_certificate import (
    calculate,
    decimation_certificate,
    mean_radius_upper,
)

def independent_mean_interval(radius, terms=64):
    """Independent positive-series enclosure for A(r)=2r Z'(r^2)/Z(r^2).

    Uses the coefficient ratio of Z(s)=sum s^m/(4^m m!(m+1)!), rather
    than the certificate's interval-polynomial/Horner implementation.
    """
    radius = F(radius)
    if radius == 0:
        return F(0), F(0)
    s = radius * radius
    z_terms = [F(1)]
    for m in range(terms):
        z_terms.append(z_terms[-1] * s / (4 * (m + 1) * (m + 2)))
    z_lo = sum(z_terms)
    z_next = z_terms[-1] * s / (4 * (terms + 1) * (terms + 2))
    z_ratio = s / (4 * (terms + 2) * (terms + 3))
    z_hi = z_lo + z_next / (1 - z_ratio)

    zp_terms = [F(m) * z_terms[m] / s for m in range(1, terms + 1)]
    zp_lo = sum(zp_terms)
    zp_next = F(terms + 1) * z_next / s
    zp_ratio = s / (4 * terms * (terms + 2))
    zp_hi = zp_lo + zp_next / (1 - zp_ratio)
    return 2 * radius * zp_lo / z_hi, 2 * radius * zp_hi / z_lo


def qmul(q, r):
    a, b, c, d = q
    e, f, g, h = r
    return (a*e-b*f-c*g-d*h,
            a*f+b*e+c*h-d*g,
            a*g-b*h+c*e+d*f,
            a*h+b*g-c*f+d*e)


def qconj(q):
    return (q[0], -q[1], -q[2], -q[3])


def qdot(q, r):
    return sum(a*b for a, b in zip(q, r))


def qnorm2(q):
    return qdot(q, q)


def qscale(a, q):
    return tuple(a*x for x in q)


class SU2DecimationCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = calculate(cells=4096)
        cls.enclosure = cls.report["enclosure"]
        cls.rows = {F(row["beta_W"]): row for row in cls.report["certificates"]}

    def test_mean_enclosure_contains_independent_positive_series(self):
        for beta in (F(0), F(1, 3), F(77, 200), F(11, 25)):
            radius = 5 * beta
            independent_lo, independent_hi = independent_mean_interval(radius)
            upper = mean_radius_upper(radius, self.enclosure)
            self.assertLessEqual(independent_lo, independent_hi)
            self.assertLessEqual(independent_hi, upper)

    def test_mean_domain_includes_zero_and_rejects_outside(self):
        self.assertEqual(mean_radius_upper(F(0), self.enclosure), 0)
        with self.assertRaises(ValueError):
            mean_radius_upper(F(-1, 10), self.enclosure)
        with self.assertRaises(ValueError):
            mean_radius_upper(F(12, 5) + F(1, 100), self.enclosure)

    def test_certificate_domain_and_endpoint_status_use_exact_fractions(self):
        at_domain_edge = decimation_certificate(F(2, 5), self.enclosure)
        self.assertFalse(at_domain_edge["certified"])
        with self.assertRaises(ValueError):
            decimation_certificate(F(401, 1000), self.enclosure)

        row_385 = self.rows[F(77, 200)]
        row_390 = self.rows[F(39, 100)]
        margin_385 = F(row_385["exact"]["margin_lower"])
        margin_390 = F(row_390["exact"]["margin_lower"])
        self.assertGreater(margin_385, 0)
        self.assertLess(margin_390, 0)
        self.assertTrue(row_385["certified"])
        self.assertFalse(row_390["certified"])

    def test_occurrence_formula_matches_saved_exact_components(self):
        row = self.rows[F(77, 200)]
        exact = row["exact"]
        same = F(exact["same_plaquette_occurrence_upper"])
        cross = F(exact["different_plaquette_occurrence_upper"])
        alpha = F(exact["column_sum_upper"])
        self.assertEqual(row["same_plaquette_occurrences"], 6)
        self.assertEqual(row["different_plaquette_occurrences"], 60)
        self.assertEqual(alpha, 6 * same + 60 * cross)
        self.assertEqual(F(exact["margin_lower"]), 1 - alpha)

    def test_exact_nonnegative_margin_and_formula_floor_at_beta_044(self):
        margin = F(self.rows[F(77, 200)]["exact"]["margin_lower"])
        self.assertGreater(margin.numerator, 0)
        lower_mean, _ = independent_mean_interval(F(11, 5))
        beta = F(11, 25)
        self.assertGreater(lower_mean, F(11, 25))
        # This is a lower floor for the proposed upper-bound formula when
        # delta=M=0, not a lower bound on the actual Dobrushin coefficient.
        floor = (6 * beta * F(11, 25) + 60 * beta**2 * F(1, 4)) / 4
        self.assertEqual(floor, F(2541, 2500))
        self.assertGreater(floor, 1)

    def test_quaternion_conformality_and_same_face_cancellation(self):
        one = (F(1), F(0), F(0), F(0))
        y0 = (F(3, 5), F(4, 5), F(0), F(0))
        y1 = (F(0), F(3, 5), F(4, 5), F(0))
        x = (F(0), F(0), F(3, 5), F(4, 5))
        self.assertEqual(qnorm2(y0), 1)
        self.assertEqual(qnorm2(y1), 1)
        self.assertEqual(qnorm2(x), 1)

        for step in (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)):
            y = tuple((1-step)*a + step*b for a, b in zip(y0, y1))
            # Quaternion multiplication is conformal even when y is the
            # nonunit point on the endpoint chord.
            self.assertEqual(qnorm2(qmul(y, x)), qnorm2(y) * qnorm2(x))
            self.assertEqual(qnorm2(qmul(x, y)), qnorm2(x) * qnorm2(y))
            for v in (y0, y1, qscale(F(1, 2), y0)):
                self.assertEqual(qdot(qmul(y, x), qmul(v, x)),
                                 qnorm2(x) * qdot(y, v))
                self.assertEqual(qdot(qmul(x, y), qmul(x, v)),
                                 qnorm2(x) * qdot(y, v))
                cy, cv = qconj(y), qconj(v)
                self.assertEqual(qdot(qmul(cy, x), qmul(cv, x)),
                                 qnorm2(x) * qdot(y, v))
        self.assertEqual(qmul(one, x), x)


if __name__ == "__main__":
    unittest.main()
