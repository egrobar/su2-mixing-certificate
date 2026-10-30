"""Focused checks for the pinned-source forest-gauge geometry audit."""

import unittest
from collections import defaultdict
from itertools import product

from experiments.su2_forest_geometry import (
    SIGMA_ROWS,
    all_plaquettes,
    check_finite_free_boxes,
    check_period,
    is_block_plaquette,
    plaquette_edges,
    sig,
)


class ForestGeometryTests(unittest.TestCase):
    def test_period_six_independent_owner_and_incidence_counts(self):
        p = 6
        sites = list(product(range(p), repeat=4))
        active = {(x, mu) for x in sites for mu in range(4) if mu != sig(x)}
        owners = defaultdict(list)
        block_faces = []
        for x, mu, nu in all_plaquettes(p):
            edges = plaquette_edges(x, mu, nu, p)
            # Recount the fixed-edge incidence directly from its table label.
            fixed_count = sum(edge[1] == sig(edge[0]) for edge in edges)
            self.assertEqual(fixed_count, 1)
            if is_block_plaquette(x, mu, nu):
                face = (x, mu, nu)
                block_faces.append(face)
                for edge in edges:
                    if edge in active:
                        owners[edge].append(face)

        self.assertEqual(len(active), 3 * p**4)
        self.assertEqual(len(block_faces), p**4)
        self.assertEqual(set(owners), active)
        self.assertTrue(all(len(v) == 1 for v in owners.values()))

        # For every non-block face the owner labels are all different. This
        # checks the partition/interaction relation, not just the global totals.
        for x, mu, nu in all_plaquettes(p):
            if is_block_plaquette(x, mu, nu):
                continue
            face_owners = [owners[e][0] for e in plaquette_edges(x, mu, nu, p)
                           if e in active]
            self.assertEqual(len(face_owners), 3)
            self.assertEqual(len(set(face_owners)), 3)

        report = check_period(p)
        self.assertEqual(report["external_partners_per_active_link"], [10])
        self.assertEqual(report["total_active_partners_per_link"], [12])

    def test_period_two_table_parity_identity(self):
        self.assertEqual(len(SIGMA_ROWS), 16)
        for parity, values in SIGMA_ROWS.items():
            fixed, *pi = values
            self.assertIn(fixed, range(4))
            other = set(range(4)) - {fixed}
            self.assertEqual({pi[d] for d in other}, other)
            for mu, nu in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
                # The three summands encode: fixed direction lies in the
                # face, or one of the two oriented arrows selects that face.
                self.assertEqual(
                    int(fixed in (mu, nu)) + int(pi[mu] == nu) + int(pi[nu] == mu),
                    1,
                    (parity, mu, nu),
                )

    def test_free_box_path_audit_is_independent_of_torus_period(self):
        report = check_finite_free_boxes()
        self.assertEqual(report["cases"], 4096)
        self.assertTrue(report["all_components_are_path_segments"])


if __name__ == "__main__":
    unittest.main()
