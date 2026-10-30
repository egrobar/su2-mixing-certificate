"""Outward integer intervals for SU(2) vMF derivatives and a local block bound.

This is a finite rational certificate, not a Lean proof or a physical gap.
The analytic response theorem and gauge geometry are separate obligations.
No numerical libraries or upstream certificate programs are used.
"""

from fractions import Fraction
from math import factorial
from pathlib import Path
import argparse
import json

SCALE = 1 << 120


def ceil_div(a, b):
    return -((-a) // b)


class Interval:
    """Endpoints are integer multiples of 2**-120; every operation rounds out."""

    def __init__(self, lo, hi=None):
        self.lo = int(lo)
        self.hi = int(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError("reversed interval")

    @classmethod
    def rational(cls, value):
        q = Fraction(value)
        return cls(q.numerator * SCALE // q.denominator,
                   ceil_div(q.numerator * SCALE, q.denominator))

    @classmethod
    def hull(cls, lo, hi):
        return cls(cls.rational(lo).lo, cls.rational(hi).hi)

    def __add__(self, other):
        if not isinstance(other, Interval):
            other = Interval.rational(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other if isinstance(other, Interval) else -Fraction(other))

    def __mul__(self, other):
        if not isinstance(other, Interval):
            other = Interval.rational(other)
        ends = [self.lo * other.lo, self.lo * other.hi,
                self.hi * other.lo, self.hi * other.hi]
        return Interval(min(ends) // SCALE, ceil_div(max(ends), SCALE))

    __rmul__ = __mul__

    def __truediv__(self, other):
        if not isinstance(other, Interval):
            other = Interval.rational(other)
        if other.lo <= 0 <= other.hi:
            raise ValueError("division across zero")
        inv = Interval((SCALE * SCALE) // other.hi,
                       ceil_div(SCALE * SCALE, other.lo))
        return self * inv

    def power(self, n):
        if n < 0:
            raise ValueError("negative power")
        out = Interval.rational(1)
        for _ in range(n):
            out = out * self
        return out

    def abs_upper(self):
        return max(abs(self.lo), abs(self.hi))

    def bounds(self):
        return [str(Fraction(self.lo, SCALE)), str(Fraction(self.hi, SCALE))]


def coefficient(j, m):
    # d^j/ds^j Z(s), Z(s)=2 I_1(sqrt(s))/sqrt(s).
    return Fraction(1, 4 ** (m + j) * factorial(m) * factorial(m + j + 1))


def derivative_tables(radius, degree):
    smax = Fraction(radius) ** 2
    tables = []
    for j in range(4):
        coeffs = [Interval.rational(coefficient(j, m)) for m in range(degree + 1)]
        first = coefficient(j, degree + 1) * smax ** (degree + 1)
        ratio = smax / (4 * (degree + 2) * (degree + j + 3))
        if not 0 <= ratio < 1:
            raise ValueError("geometric remainder unavailable")
        tail = Interval.hull(0, first / (1 - ratio))
        tables.append((coeffs, tail))
    return tables


def evaluate(table, s):
    coeffs, tail = table
    out = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        out = out * s + c
    return out + tail


def derivative_enclosure(radius=Fraction(21, 10), cells=4096, degree=24):
    radius = Fraction(radius)
    if radius <= 0 or cells <= 0 or degree < 1:
        raise ValueError("invalid certificate domain")
    tables = derivative_tables(radius, degree)
    lam_lower = SCALE
    tau_upper = 0
    third_upper = 0
    critical = {}
    for k in range(cells):
        r = Interval.hull(radius * k / cells, radius * (k + 1) / cells)
        s = r * r
        z0, z1, z2, z3 = [evaluate(t, s) for t in tables]
        p1 = z1 / z0
        p2 = z2 / z0 - p1 * p1
        p3 = z3 / z0 - 3 * z1 * z2 / (z0 * z0) + 2 * p1.power(3)
        transverse = 2 * p1
        radial = transverse + 4 * s * p2
        low = min(transverse.lo, radial.lo)
        high = max(transverse.hi, radial.hi)
        if low < lam_lower:
            lam_lower = low
            critical["hessian_lower_cell"] = k
        if high > tau_upper:
            tau_upper = high
            critical["hessian_upper_cell"] = k
        a = 12 * r * p2
        gamma = 8 * r.power(3) * p3
        # For At+Gamma*t^3, endpoints and any interior critical points
        # give an UPPER bound max(|A+Gamma|,2|A|/3). It is not equality.
        third = max((a + gamma).abs_upper(), ceil_div(2 * a.abs_upper(), 3))
        if third > third_upper:
            third_upper = third
            critical["third_upper_cell"] = k
    return {
        "radius": str(radius), "cells": cells, "polynomial_degree": degree,
        "rounding_bits": 120,
        "lambda_lower": str(Fraction(lam_lower, SCALE)),
        "tau_interval_upper": str(Fraction(tau_upper, SCALE)),
        "third_upper": str(Fraction(third_upper, SCALE)),
        "critical_cells": critical,
    }


def exp_upper(x, degree=32):
    x = Fraction(x)
    if x < 0:
        raise ValueError("nonnegative argument required")
    total = term = Fraction(1)
    for n in range(1, degree + 1):
        term *= x / n
        total += term
    first = term * x / (degree + 1)
    ratio = x / (degree + 2)
    if ratio >= 1:
        raise ValueError("exponential tail ratio")
    return total + first / (1 - ratio)


def block_certificate(beta, enclosure, improved=True):
    """Safe rational majorants; uses separately proved tau<=1/4 globally."""
    beta = Fraction(beta)
    if beta < 0 or 6 * beta > Fraction(enclosure["radius"]):
        raise ValueError("coupling outside certified field domain")
    lam = Fraction(enclosure["lambda_lower"])
    third = Fraction(enclosure["third_upper"])
    tau = Fraction(1, 4)
    if not 0 < lam <= tau:
        raise ValueError("invalid Hessian lower bound")
    w = beta * beta * (tau - lam) / 2
    ew = exp_upper(w)
    # 289/1000 > sqrt(1/12), verified by squaring rationally.
    sd_product = Fraction(289, 1000) if improved else Fraction(1, 2)
    s = tau + (ew - 1) * ew * sd_product
    residual_factor = ew / 3 if improved else Fraction(1)
    c = 5 * beta * (beta * tau * s + beta * beta * third * residual_factor)
    if c >= 1:
        return {"beta_W": str(beta), "certified": False, "reason": "C >= 1"}
    k = (1 + beta / 4) * s / (1 - c)
    alpha = 10 * beta * k
    vals = {"w_upper": w, "s_upper": s, "C_upper": c,
            "K_upper": k, "alpha_upper": alpha, "margin_lower": 1 - alpha}
    return {"beta_W": str(beta), "improved_tilt": improved,
            "certified": alpha < 1,
            "exact": {name: str(value) for name, value in vals.items()},
            "display_only": {name: float(value) for name, value in vals.items()}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=int, default=4096)
    parser.add_argument("--radius", default="21/10")
    parser.add_argument("--output", default=".artifacts/su2-vmf-certificate.json")
    args = parser.parse_args()
    enclosure = derivative_enclosure(Fraction(args.radius), args.cells)
    report = {
        "scope": "Uniform local three-link response, conditional on audited analytic lemmas; not a physical gap or Lean proof",
        "enclosure": enclosure,
        "enclosure_display_only": {x: float(Fraction(enclosure[x])) for x in
                                    ("lambda_lower", "tau_interval_upper", "third_upper")},
        "certificates": [block_certificate(b, enclosure, improved) for b in
                         ("1/3", "69/200", "87/250", "7/20") for improved in (False, True)],
    }
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps({"enclosure": report["enclosure_display_only"], "rows": [
        {"beta": row["beta_W"], "improved": row["improved_tilt"],
         "alpha_upper": row["display_only"]["alpha_upper"], "certified": row["certified"]}
        for row in report["certificates"]]}, indent=2))


if __name__ == "__main__":
    main()
