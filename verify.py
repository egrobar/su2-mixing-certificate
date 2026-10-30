"""Recompute the finite arithmetic/geometry certificate using only Python stdlib.

The response theorem and vMF analytic inputs remain mathematical obligations;
this command neither builds Lean nor certifies a continuum mass gap.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json
import platform

from experiments.su2_decimation_certificate import calculate
from experiments.su2_color_geometry import build_report, COLOR_ROWS


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def unwrapped_midpoint_geometry():
    """Check all 16 eliminated-link parity types without periodic wrapping.

    A link midpoint is represented by twice its integer coordinates. The
    metric is lattice l1 distance, so doubled distance <=4 means distance <=2.
    Local incidence has period two; these types cover every infinite star.
    """
    def shift(x, d, sign=1):
        return tuple(v + (sign if i == d else 0) for i, v in enumerate(x))

    def edges(x, u, v):
        return ((x, u), (shift(x, u), v), (shift(x, v), u), (x, v))

    def twice_midpoint(edge):
        x, d = edge
        return tuple(2*v + (i == d) for i, v in enumerate(x))

    maximum = 0
    count = 0
    for x in product(range(2), repeat=4):
        for mu in range(4):
            # Fixed links have color -1 in the literal table, so use the table
            # directly to select color 2 and avoid calling color_of on them.
            if COLOR_ROWS[x][mu] != 2:
                continue
            eliminated = (x, mu)
            members = set()
            for nu in range(4):
                if nu == mu:
                    continue
                for base in (x, shift(x, nu, -1)):
                    face = edges(base, min(mu, nu), max(mu, nu))
                    require(eliminated in face, "incorrect unwrapped incidence")
                    labels = [COLOR_ROWS[tuple(v % 2 for v in e[0])][e[1]] for e in face]
                    require(sorted(labels) == [-1, 0, 1, 2], "unwrapped face not rainbow")
                    members.update(e for e, c in zip(face, labels) if c in (0, 1))
            require(len(members) == 12, "unwrapped eliminated star arity")
            midpoints = [twice_midpoint(e) for e in members]
            maximum = max(maximum, max(sum(abs(a-b) for a, b in zip(p, q))
                                      for p in midpoints for q in midpoints))
            count += 1
    require(count == 16 and maximum <= 4, "midpoint range bound failed")
    return {"eliminated_parity_types": count,
            "max_doubled_midpoint_l1_distance": maximum,
            "midpoint_l1_range_upper": "2"}


def verify():
    report = calculate(cells=4096)
    geometry = build_report(period=6)
    enc = report["enclosure"]
    endpoint = next(row for row in report["certificates"] if F(row["beta_W"]) == F(77, 200))
    alpha = F(endpoint["exact"]["column_sum_upper"])
    target = F(990737, 1000000)
    require(endpoint["certified"] and alpha < target < 1, "exact endpoint inequality failed")
    require(F(endpoint["exact"]["margin_lower"]) == 1-alpha, "inconsistent margin")
    require(F(enc["radius"]) == F(12, 5) and 6*F(77,200) < F(enc["radius"]), "field domain guard")
    require(F(enc["lambda_lower"]) > F(1393407619,10**10), "Hessian lower bound")
    require(F(enc["third_upper"]) < F(617748750,10**10), "third derivative bound")
    require(F(289,1000)**2 > F(1,12), "square-root rational majorant")
    require(F(endpoint["exact"]["hessian_spread_upper"]) >= 0, "negative Hessian spread")
    require(F(endpoint["exact"]["column_sum_upper"]) ==
            6*F(endpoint["exact"]["same_plaquette_occurrence_upper"])
            +60*F(endpoint["exact"]["different_plaquette_occurrence_upper"]), "occurrence formula")
    incidence = geometry["factor_incidence"]
    for key, expected in (("retained_factors_per_link_values",[6]),
            ("same_face_row_occurrences_per_retained_link",[6]),
            ("same_face_column_occurrences_per_retained_link",[6]),
            ("cross_factor_row_occurrences_per_retained_link",[60]),
            ("cross_factor_column_occurrences_per_retained_link",[60])):
        require(incidence[key] == expected, "geometry count: "+key)
    require(geometry["color_table"]["face_types_checked"] == 96, "parity face coverage")
    require(incidence["factor_arity_histogram"] == {12: 1296}, "factor arity")
    require(incidence["distinct_retained_neighbors_histogram"] == {"52":2592}, "neighbor count")
    midpoint = unwrapped_midpoint_geometry()
    # A compact outward rational endpoint is saved instead of the 5393-character
    # exact internal fraction. Its digest permits comparison of the full value.
    decimal_scale = 10**15
    rounded = F(-((-alpha.numerator*decimal_scale)//alpha.denominator), decimal_scale)
    require(alpha <= rounded < target, "outward summary rounding")
    sources = ["experiments/su2_vmf_certificate.py", "experiments/su2_decimation_certificate.py",
               "experiments/su2_forest_geometry.py", "experiments/su2_color_geometry.py"]
    return {
        "certificate": {"beta_W_interval": ["0", "77/200"],
            "column_and_row_sum_upper_rational": str(rounded),
            "strict_display_bound_rational": str(target),
            "positive_margin_lower_rational": str(1-rounded),
            "exact_internal_alpha_sha256": hashlib.sha256(str(alpha).encode("ascii")).hexdigest(),
            "exact_fraction_decision_passed": True,
            "extension_to_smaller_beta": "Analytic A'(r)=Var(x_1)>=0 and positivity of all majorant factors, using fixed domain constants; no beta grid."},
        "parameters": enc,
        "geometry": {"parity_face_types":96,"period":6,"factor_arity":12,
            "factors_per_retained_link":6,"same_occurrences_row_and_column":6,
            "cross_occurrences_row_and_column":60,"distinct_neighbors":52,
            "unwrapped_midpoint_check":midpoint},
        "scope": "Finite arithmetic and local geometry certificate for projected bulk/free-box SU(2) forest specification; analytic response inputs are separate. No Lean build or continuum gap claim.",
        "runtime": {"implementation":platform.python_implementation(),
                    "python_version":platform.python_version(),"third_party_dependencies":[]},
        "source_sha256":{p:hashlib.sha256(Path(__file__).parent.joinpath(p).read_bytes()).hexdigest() for p in sources},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",default="verification/certificate-summary.json")
    args = parser.parse_args()
    summary = verify()
    destination = Path(args.output)
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("PASS: exact alpha < 990737/1000000 < 1 at beta_W=77/200")
    print("PASS: 96 parity faces; row/column occurrences 6+60; unwrapped midpoint range <=2")
    print("Summary: "+destination.as_posix())


if __name__ == "__main__":
    main()
