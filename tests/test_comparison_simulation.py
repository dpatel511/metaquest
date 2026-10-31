"""Small validation fixtures for the development Monte Carlo harness."""
import copy
import json
from pathlib import Path
import runpy
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.stats import binomtest


ROOT = Path(__file__).resolve().parents[1]
benchmark = SimpleNamespace(**runpy.run_path(str(ROOT / "tools/benchmark_comparison.py")))
PROTOCOL = json.loads((ROOT / "benchmarks/comparison/protocol.json").read_text())


@pytest.mark.parametrize("successes", [0, 25, 250, 500])
def test_wilson_interval_agrees_with_independent_scipy_calculation(successes):
    expected = binomtest(successes, 500).proportion_ci(confidence_level=.95, method="wilson")
    np.testing.assert_allclose(benchmark.wilson(successes, 500), [expected.low, expected.high], atol=1e-15)


def test_zero_observed_fdp_does_not_claim_zero_uncertainty():
    interval = benchmark.mean_bounds([0] * 500)
    assert interval["mean"] == 0
    assert .05 < interval["upper_one_sided_95"] < .06


@pytest.mark.parametrize("scenario", PROTOCOL["scenarios"], ids=lambda scenario: scenario["name"])
def test_frozen_generators_are_reproducible_valid_compositions(scenario):
    first, samples = benchmark.generate(PROTOCOL, scenario, 0)
    second, _ = benchmark.generate(PROTOCOL, scenario, 0)
    other, _ = benchmark.generate(PROTOCOL, scenario, 1)
    np.testing.assert_array_equal(first.to_numpy(), second.to_numpy())
    assert not np.array_equal(first.to_numpy(), other.to_numpy())
    assert first.index.equals(samples.index)
    assert first.shape == (scenario["n_reference"] + scenario["n_comparison"], PROTOCOL["features"])
    assert np.isfinite(first.to_numpy()).all() and (first.to_numpy() >= 0).all()
    np.testing.assert_allclose(first.sum(axis=1), 1)
    assert bool((first.to_numpy() == 0).any()) == (scenario["name"] == "sparse_null")


def records():
    rows = []
    for scenario in PROTOCOL["scenarios"]:
        for method in PROTOCOL["parameters"]["fdr_methods"]:
            for replicate in range(PROTOCOL["replicates_per_scenario"]):
                rows.append({"scenario": scenario["name"], "method": method, "setting": "primary",
                             "replicate": replicate, "discoveries": 0, "fdp": 0,
                             "target_recall": 1, "both_target_directions_correct": True})
    return rows


def test_acceptance_keeps_failures_and_ignores_descriptive_sensitivity_rows():
    rows = records()
    assert all(item["passed"] for item in benchmark.summarize(rows, PROTOCOL))
    failed = copy.deepcopy(rows)
    for row in failed:
        if row["scenario"] == "balanced_clr_effect":
            row["target_recall"] = 0
    summary = benchmark.summarize(failed, PROTOCOL)
    assert all(not item["passed"] for item in summary if item["scenario"] == "balanced_clr_effect")
    rows.extend({**row, "setting": "sensitivity", "discoveries": 24} for row in failed[:50])
    assert all(item["passed"] for item in benchmark.summarize(rows, PROTOCOL))


def test_missing_simulation_studies_are_rejected():
    with pytest.raises(ValueError, match="Missing or duplicated"):
        benchmark.summarize(records()[1:], PROTOCOL)



def test_unsupported_tied_population_mass_is_not_reported_as_zero_fdr():
    tied = SimpleNamespace(**runpy.run_path(str(ROOT / "tools/benchmark_comparison_ties.py")))
    protocol = {"n_comparison_reference": [[3, 3]], "population_probabilities": [.5],
                "settings": {"methods": ["by", "bh"], "pseudocount": 1e-6,
                             "min_prevalence": .1, "alpha": .05}}
    result = tied.calculate(protocol)
    assert result["count_pair_states"] == 16
    assert result["blocked_state_method_evaluations"] == 28
    for row in result["population_evaluations"]:
        assert row["unsupported_probability_mass"] == pytest.approx(31 / 32)
        assert row["exact_global_null_fdr"] is None
        assert row["passes_nominal_alpha"] is None



def test_tied_diagnostic_uses_the_protocol_significance_level():
    tied = SimpleNamespace(**runpy.run_path(str(ROOT / "tools/benchmark_comparison_ties.py")))
    protocol = {"n_comparison_reference": [[10, 10]], "population_probabilities": [.5],
                "settings": {"methods": ["bh"], "pseudocount": 1e-6,
                             "min_prevalence": .1, "alpha": .05}}
    ordinary = tied.calculate(protocol)["population_evaluations"][0]
    protocol["settings"]["alpha"] = 1e-8
    stringent = tied.calculate(protocol)["population_evaluations"][0]
    assert ordinary["exact_global_null_fdr"] > 0
    assert stringent["exact_global_null_fdr"] == 0
    assert stringent["passes_nominal_alpha"]
