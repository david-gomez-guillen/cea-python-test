"""Sanity checks of the model.

Each test states a property any correct cohort model must have. Run them from the
project folder with::

    python -m unittest discover tests

(or with ``pytest`` if it is installed).
"""

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import example  # noqa: E402
import model  # noqa: E402


class TestBuildingBlocks(unittest.TestCase):

    def test_rate_and_probability_are_inverse(self):
        for p in (0.0, 0.01, 0.1, 0.5, 0.9):
            self.assertAlmostEqual(model.rate_to_probability(model.probability_to_rate(p)), p)

    def test_hazard_ratio_of_one_changes_nothing(self):
        self.assertAlmostEqual(model.apply_hazard_ratio(0.1, 1.0), 0.1)

    def test_hazard_ratio_keeps_probability_below_one(self):
        # Multiplying the probability by the HR would give 1.8 here.
        self.assertLess(model.apply_hazard_ratio(0.6, 3.0), 1.0)

    def test_hazard_ratio_on_rate(self):
        # Surviving the cycle without the event: (1 - p) ** HR.
        self.assertAlmostEqual(model.apply_hazard_ratio(0.1, 0.5), 1 - 0.9 ** 0.5)

    def test_discount_factor_after_one_year(self):
        self.assertAlmostEqual(model.discount_factor(12, 0.03), 1 / 1.03)
        self.assertEqual(model.discount_factor(0, 0.03), 1.0)

    def test_rows_of_transition_matrix_add_up_to_one(self):
        matrix = model.transition_matrix(0.1, 0.01, 0.07)
        np.testing.assert_allclose(matrix.sum(axis=1), 1.0)

    def test_invalid_probabilities_are_rejected(self):
        with self.assertRaises(ValueError):
            model.transition_matrix(0.995, 0.01, 0.07)

    def test_value_by_year(self):
        by_year = [0.1 * (year + 1) for year in range(model.HORIZON_YEARS)]
        self.assertEqual(model.value_in_cycle(by_year, 0), 0.1)
        self.assertEqual(model.value_in_cycle(by_year, 12), 0.2)
        self.assertEqual(model.value_in_cycle(0.3, 100), 0.3)


class TestCohort(unittest.TestCase):

    def setUp(self):
        self.results = model.simulate(model.STRATEGIES, **model.base_case())

    def test_cohort_always_adds_up_to_one(self):
        for trace in self.results["cohort_info"].values():
            np.testing.assert_allclose(trace.sum(axis=1), 1.0)

    def test_nobody_comes_back_from_the_dead(self):
        for trace in self.results["cohort_info"].values():
            self.assertTrue(np.all(np.diff(trace[:, 2]) >= 0))

    def test_survival_curves_never_go_up(self):
        for curve in ("pfs", "os"):
            for values in self.results[curve].values():
                self.assertTrue(np.all(np.diff(values) <= 0))

    def test_pfs_is_below_os(self):
        for strategy in model.STRATEGIES:
            self.assertTrue(np.all(np.array(self.results["pfs"][strategy])
                                   <= np.array(self.results["os"][strategy])))

    def test_immunotherapy_delays_progression(self):
        pfs = self.results["pfs"]
        self.assertGreater(pfs["immunotherapy"][0], pfs["chemotherapy"][0])

    def test_costs_add_up(self):
        outcomes = self.results["outcomes"]
        parts = outcomes[["cost_drug", "cost_care", "cost_end_of_life", "cost_test"]].sum(axis=1)
        np.testing.assert_allclose(parts, outcomes["cost_total"])
        np.testing.assert_allclose(outcomes["cost_total"], self.results["summary"]["C"])


class TestStrategies(unittest.TestCase):

    def summary(self, **changes):
        params = model.base_case() | changes
        return model.simulate(model.STRATEGIES, **params)["summary"].set_index("strategy")

    def test_everyone_positive_means_testing_changes_only_the_cost(self):
        summary = self.summary(p_biomarker_positive=1.0)
        test = model.base_case()["cost_biomarker_test"]
        self.assertAlmostEqual(summary.loc["biomarker_immunotherapy", "E"],
                               summary.loc["immunotherapy", "E"])
        self.assertAlmostEqual(summary.loc["biomarker_immunotherapy", "C"],
                               summary.loc["immunotherapy", "C"] + test)

    def test_everyone_negative_means_testing_leads_to_chemotherapy(self):
        summary = self.summary(p_biomarker_positive=0.0)
        self.assertAlmostEqual(summary.loc["biomarker_immunotherapy", "E"],
                               summary.loc["chemotherapy", "E"])

    def test_no_discounting_gives_quality_adjusted_life_years(self):
        params = model.base_case() | {"discount": 0.0, "utility_pfs": 1.0,
                                      "utility_progressed": 1.0}
        outcomes = model.simulate(["chemotherapy"], **params)["outcomes"]
        self.assertAlmostEqual(outcomes.loc[0, "qalys"], outcomes.loc[0, "life_years"])

    def test_one_value_per_year_equals_one_value(self):
        p = model.base_case()["p_progression"]
        by_year = model.simulate(["chemotherapy"], p_progression=[p] * model.HORIZON_YEARS)
        constant = model.simulate(["chemotherapy"], p_progression=p)
        pd.testing.assert_frame_equal(by_year["summary"], constant["summary"])

    def test_only_probabilities_vary_by_year(self):
        with self.assertRaises(ValueError):
            model.simulate(["chemotherapy"], cost_chemo=[2500] * model.HORIZON_YEARS)

    def test_unknown_strategy_is_rejected(self):
        with self.assertRaises(ValueError):
            model.simulate(["surgery"])


class TestIncrementalAnalysis(unittest.TestCase):

    def test_dominance(self):
        summary = pd.DataFrame({"strategy": ["A", "B", "C", "D"],
                                "C": [0, 100, 150, 400],
                                "E": [1.0, 0.9, 2.0, 2.5]})
        table = example.incremental_analysis(summary).set_index("strategy")
        # B costs more than A and gives less: strongly dominated.
        self.assertEqual(table.loc["B", "status"], "dominated")
        self.assertAlmostEqual(table.loc["C", "ICER"], 150)
        self.assertAlmostEqual(table.loc["D", "ICER"], 500)

    def test_extended_dominance(self):
        summary = pd.DataFrame({"strategy": ["A", "B", "C"],
                                "C": [0, 900, 1000],
                                "E": [1.0, 1.5, 2.0]})
        table = example.incremental_analysis(summary).set_index("strategy")
        # B: 1800 per QALY, then C: 200 per QALY. Going straight to C is better.
        self.assertEqual(table.loc["B", "status"], "extendedly dominated")
        self.assertAlmostEqual(table.loc["C", "ICER"], 1000)


if __name__ == "__main__":
    unittest.main()
