"""Focused exact checks for period-two SU(2) color elimination geometry."""

import unittest

from experiments.su2_color_geometry import (
    COLOR_ROWS,
    check_color_table,
    factor_incidence_report,
)


class SU2ColorGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.periodic = check_color_table()
        cls.period_six = factor_incidence_report(6)

    def test_period_two_table_is_balanced_and_rainbow(self):
        self.assertEqual(self.periodic["parity_rows"], 16)
        self.assertEqual(self.periodic["face_types_checked"], 96)
        self.assertTrue(self.periodic["all_faces_rainbow"])
        self.assertEqual(
            self.periodic["links_per_color_per_period_two_cell"],
            {0: 16, 1: 16, 2: 16},
        )

    def test_color_table_corruption_is_detected(self):
        corrupted = dict(COLOR_ROWS)
        corrupted[(0, 0, 0, 0)] = (-1, 1, 1, 2)
        with self.assertRaisesRegex(ValueError, "singly colored"):
            check_color_table(corrupted)

    def test_period_six_factor_arity_occurrences_and_neighbors(self):
        report = self.period_six
        self.assertEqual(report["period"], 6)
        self.assertEqual(report["active_links"], 3 * 6**4)
        self.assertEqual(report["retained_links"], 2 * 6**4)
        self.assertEqual(report["eliminated_links"], 6**4)
        self.assertEqual(report["factor_arity_histogram"], {12: 6**4})
        self.assertEqual(report["retained_factors_per_link_values"], [6])
        self.assertEqual(report["same_face_row_occurrences_per_retained_link"], [6])
        self.assertEqual(report["same_face_column_occurrences_per_retained_link"], [6])
        self.assertEqual(report["cross_factor_row_occurrences_per_retained_link"], [60])
        self.assertEqual(report["cross_factor_column_occurrences_per_retained_link"], [60])
        self.assertEqual(
            report["distinct_retained_neighbors_histogram"],
            {"52": 2 * 6**4},
        )
        self.assertEqual(report["max_distinct_retained_neighbors"], 52)
        self.assertEqual(report["ordered_same_face_occurrences"], 6 * 2 * 6**4)
        self.assertEqual(report["ordered_cross_factor_occurrences"], 60 * 2 * 6**4)
        self.assertTrue(report["all_plaquettes_rainbow"])


if __name__ == "__main__":
    unittest.main()
