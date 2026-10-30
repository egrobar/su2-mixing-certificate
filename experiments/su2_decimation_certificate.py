"""Rational certificate for one exact color elimination in the forest gauge.

The geometric coloring, conformal quaternion identity, and projected Gibbs
specification are analytic/combinatorial inputs documented separately.
This is not a certificate for a cube with frozen original color-2 links.
"""

from fractions import Fraction
from pathlib import Path
import argparse
import json

from experiments.su2_vmf_certificate import (
    Interval, SCALE, derivative_enclosure, derivative_tables, evaluate, exp_upper,
)


def mean_radius_upper(radius, enclosure):
    """Bound A(r)=E_{nu_(r e1)}x1 = 2r Z'(r^2)/Z(r^2), including r=0."""
    radius = Fraction(radius)
    if not 0 <= radius <= Fraction(enclosure["radius"]):
        raise ValueError("mean parameter outside enclosure domain")
    tables = derivative_tables(
        Fraction(enclosure["radius"]), enclosure["polynomial_degree"]
    )
    squared = Interval.rational(radius * radius)
    z0, z1 = [evaluate(t, squared) for t in tables[:2]]
    return Fraction((2 * radius * z1 / z0).hi, SCALE)


def decimation_certificate(beta, enclosure):
    """Endpoint majorant, uniform over all retained-link boundary values.

    Monotonicity of the true A(r) (A'=Var x1 >=0) extends a successful
    endpoint to all smaller beta. Pointwise numerical monotonicity of an
    interval implementation is not assumed.
    """
    beta = Fraction(beta)
    if beta < 0 or 6 * beta > Fraction(enclosure["radius"]):
        raise ValueError("coupling outside certified field domain")
    lam = Fraction(enclosure["lambda_lower"])
    third = Fraction(enclosure["third_upper"])
    tau = Fraction(1, 4)  # separately proved global vMF Hessian upper bound
    if not 0 < lam <= tau or third < 0:
        raise ValueError("invalid derivative enclosure")
    delta = tau - lam
    mean = mean_radius_upper(5 * beta, enclosure)
    w = 3 * beta * beta * delta
    ew = exp_upper(w)
    sensitivity = tau + (ew - 1) * ew * Fraction(289, 1000)
    same = sensitivity * beta * mean + ew * 2 * beta * beta * delta / 3
    cross = sensitivity * beta * beta * tau + ew * beta**3 * third / 3
    alpha = 6 * same + 60 * cross
    exact = {
        "mean_upper": mean, "hessian_spread_upper": delta, "w_upper": w,
        "tilted_sensitivity_upper": sensitivity,
        "same_plaquette_occurrence_upper": same,
        "different_plaquette_occurrence_upper": cross,
        "column_sum_upper": alpha, "margin_lower": 1 - alpha,
    }
    return {
        "beta_W": str(beta), "certified": alpha < 1,
        "same_plaquette_occurrences": 6,
        "different_plaquette_occurrences": 60,
        "exact": {key: str(value) for key, value in exact.items()},
        "display_only": {key: float(value) for key, value in exact.items()},
    }


def calculate(cells=4096):
    enclosure = derivative_enclosure(radius=Fraction(12, 5), cells=cells)
    return {
        "scope": "Projected bulk forest specification after one exact color elimination; fixed lattice, no continuum gap or Lean proof",
        "excluded_boundary_claim": "Frozen original color-2 exterior may leave bare plaquette terms; this is not their conditional-law certificate",
        "enclosure": enclosure,
        "certificates": [decimation_certificate(beta, enclosure) for beta in
                         ("0", "1/3", "87/250", "19/50", "77/200", "39/100", "2/5")],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=int, default=4096)
    parser.add_argument("--output", default=".artifacts/su2-decimation-certificate.json")
    args = parser.parse_args()
    report = calculate(args.cells)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps([{
        "beta_W": row["beta_W"], "certified": row["certified"],
        "column_sum_upper_display": row["display_only"]["column_sum_upper"],
    } for row in report["certificates"]], indent=2))


if __name__ == "__main__":
    main()
