#!/usr/bin/env python3
"""Exact period-two coloring and one-color factor-incidence audit.

This imports the pinned forest geometry rather than duplicating its active
link or plaquette construction. It proves the displayed parity coloring is
rainbow on every plaquette type and counts the induced factors on a finite
period-six torus. It is a geometry check, not a Gibbs-response proof.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

from experiments import su2_forest_geometry as geometry


# Entry is the color of a positively oriented link at the listed residue;
# -1 marks the fixed link. Directions are columns 0,1,2,3.
COLOR_ROWS = {
    (0, 0, 0, 0): (-1, 0, 1, 2),
    (0, 0, 0, 1): (1, -1, 0, 2),
    (0, 0, 1, 0): (0, 2, 1, -1),
    (0, 0, 1, 1): (2, 1, 0, -1),
    (0, 1, 0, 0): (2, 0, -1, 1),
    (0, 1, 0, 1): (0, -1, 2, 1),
    (0, 1, 1, 0): (1, 2, -1, 0),
    (0, 1, 1, 1): (-1, 1, 2, 0),
    (1, 0, 0, 0): (-1, 1, 2, 0),
    (1, 0, 0, 1): (1, 2, -1, 0),
    (1, 0, 1, 0): (0, -1, 2, 1),
    (1, 0, 1, 1): (2, 0, -1, 1),
    (1, 1, 0, 0): (2, 1, 0, -1),
    (1, 1, 0, 1): (0, 2, 1, -1),
    (1, 1, 1, 0): (1, -1, 0, 2),
    (1, 1, 1, 1): (-1, 0, 1, 2),
}


def color_of(edge, rows=COLOR_ROWS):
    x, direction = edge
    color = rows[tuple(value % 2 for value in x)][direction]
    if color < 0:
        raise ValueError(f"fixed link has no active color: {edge}")
    return color


def check_color_table(rows=COLOR_ROWS):
    """Verify the 16-row assignment and all 96 period-two face types."""
    if set(rows) != set(geometry.SIGMA_ROWS):
        raise ValueError("color table must contain exactly the 16 geometry residues")

    class_counts = Counter()
    face_count = 0
    for parity, sigma_row in geometry.SIGMA_ROWS.items():
        colors = rows[parity]
        if len(colors) != 4 or colors[sigma_row[0]] != -1:
            raise ValueError(f"fixed-link color mismatch at residue {parity}")
        active_colors = [colors[d] for d in range(4) if d != sigma_row[0]]
        if sorted(active_colors) != [0, 1, 2]:
            raise ValueError(f"active links are not singly colored at residue {parity}")
        class_counts.update(active_colors)

        for mu in range(4):
            for nu in range(mu + 1, 4):
                edges = geometry.plaquette_edges(parity, mu, nu, p=2)
                active = [edge for edge in edges
                          if edge[1] != geometry.sig(edge[0])]
                face_colors = [color_of(edge, rows) for edge in active]
                if len(active) != 3 or sorted(face_colors) != [0, 1, 2]:
                    raise ValueError(
                        f"non-rainbow face at residue {parity}, directions {(mu, nu)}"
                    )
                face_count += 1

    if class_counts != Counter({0: 16, 1: 16, 2: 16}):
        raise ValueError("color classes are not balanced in the period-two cell")
    return {
        "parity_rows": 16,
        "face_types_checked": face_count,
        "all_faces_rainbow": True,
        "links_per_color_per_period_two_cell": dict(sorted(class_counts.items())),
    }


def factor_incidence_report(period=6, rows=COLOR_ROWS):
    """Enumerate eliminated stars and ordered response incidences."""
    check_color_table(rows)
    geometry_report = geometry.check_period(period)
    sites = tuple(geometry.product(range(period), repeat=4))
    active_links = {
        (x, direction)
        for x in sites for direction in range(4)
        if direction != geometry.sig(x)
    }
    link_colors = {edge: color_of(edge, rows) for edge in active_links}
    retained_links = {edge for edge in active_links if link_colors[edge] in (0, 1)}
    eliminated_links = {edge for edge in active_links if link_colors[edge] == 2}

    faces_by_eliminated = defaultdict(list)
    for x, mu, nu in geometry.all_plaquettes(period):
        face = (x, mu, nu)
        active = [edge for edge in geometry.plaquette_edges(x, mu, nu, period)
                  if edge in active_links]
        if len(active) != 3 or {link_colors[edge] for edge in active} != {0, 1, 2}:
            raise ValueError(f"period-{period} face failed the coloring: {face}")
        eliminated = next(edge for edge in active if link_colors[edge] == 2)
        retained = tuple(edge for edge in active if link_colors[edge] != 2)
        faces_by_eliminated[eliminated].append((face, retained))

    if set(faces_by_eliminated) != eliminated_links:
        raise ValueError("an eliminated link has no incident plaquette star")

    factors = {}
    factor_face_members = {}
    factor_arity_histogram = Counter()
    same_row = Counter()
    same_column = Counter()
    cross_row = Counter()
    cross_column = Counter()
    retained_factor_count = Counter()

    for eliminated, incidences in faces_by_eliminated.items():
        if len(incidences) != 6:
            raise ValueError(f"eliminated link has {len(incidences)} incident faces")
        by_face = []
        members = set()
        for face, retained in incidences:
            if len(retained) != 2 or {link_colors[e] for e in retained} != {0, 1}:
                raise ValueError(f"bad retained pair in face {face}")
            by_face.append(retained)
            members.update(retained)
            left, right = retained
            same_row[left] += 1
            same_row[right] += 1
            same_column[left] += 1
            same_column[right] += 1
        if len(members) != 12:
            raise ValueError(f"factor has arity {len(members)}, expected 12")
        factors[eliminated] = members
        factor_face_members[eliminated] = by_face
        factor_arity_histogram[len(members)] += 1

        for left_face, pair_left in enumerate(by_face):
            for right_face, pair_right in enumerate(by_face):
                if left_face == right_face:
                    continue
                for target in pair_left:
                    for source in pair_right:
                        cross_row[target] += 1
                        cross_column[source] += 1

        for edge in members:
            retained_factor_count[edge] += 1

    neighbors = {edge: set() for edge in retained_links}
    for members in factors.values():
        for edge in members:
            neighbors[edge].update(members - {edge})

    def values(counter, domain):
        return {counter.get(edge, 0) for edge in domain}

    row_same_values = values(same_row, retained_links)
    column_same_values = values(same_column, retained_links)
    row_cross_values = values(cross_row, retained_links)
    column_cross_values = values(cross_column, retained_links)
    factor_count_values = values(retained_factor_count, retained_links)
    if (row_same_values != {6} or column_same_values != {6}
            or row_cross_values != {60} or column_cross_values != {60}
            or factor_count_values != {6}):
        raise ValueError("row/column occurrence or factor-count control failed")

    return {
        "period": period,
        "positive_links": geometry_report["positive_links"],
        "active_links": len(active_links),
        "retained_links": len(retained_links),
        "eliminated_links": len(eliminated_links),
        "plaquettes": geometry_report["plaquettes"],
        "all_plaquettes_rainbow": True,
        "eliminated_factor_count": len(factors),
        "factor_arity_histogram": dict(sorted(factor_arity_histogram.items())),
        "retained_factors_per_link_values": sorted(factor_count_values),
        "same_face_row_occurrences_per_retained_link": sorted(row_same_values),
        "same_face_column_occurrences_per_retained_link": sorted(column_same_values),
        "cross_factor_row_occurrences_per_retained_link": sorted(row_cross_values),
        "cross_factor_column_occurrences_per_retained_link": sorted(column_cross_values),
        "distinct_retained_neighbors_histogram": {
            str(k): v for k, v in sorted(Counter(
                len(neighbors[edge]) for edge in retained_links
            ).items())
        },
        "max_distinct_retained_neighbors": max(map(len, neighbors.values())),
        "ordered_same_face_occurrences": sum(same_row.values()),
        "ordered_cross_factor_occurrences": sum(cross_row.values()),
        "scope": f"exact period-{period} incidence counts; projected bulk factors only",
    }


def build_report(period=6):
    return {
        "source_geometry": "experiments.su2_forest_geometry",
        "color_table": check_color_table(),
        "factor_incidence": factor_incidence_report(period),
        "finite_boundary_scope": (
            "Integrating a finite subset T0 leaves bare beta_W plaquette terms "
            "for plaquettes whose color-2 link is outside T0."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--period", type=int, default=6)
    parser.add_argument("--output", default=".artifacts/su2-color-geometry.json")
    args = parser.parse_args()
    report = build_report(args.period)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
