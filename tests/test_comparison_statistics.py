"""Independent arithmetic oracles; these are synthetic studies, not biological validation."""
from itertools import combinations
import math
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from metaquest.core.comparison import ComparisonError, group_tests


def observed(log_coordinates, pseudocount=1e-6):
    values = np.asarray(log_coordinates, dtype=float)
    weights = np.exp(values - values.max(axis=1, keepdims=True))
    weights /= weights.sum(axis=1, keepdims=True)
    proportions = (1 + values.shape[1] * pseudocount) * weights - pseudocount
    return pd.DataFrame(proportions, columns=[f"f{i}" for i in range(values.shape[1])])


def compare(values, size, **options):
    samples = pd.DataFrame({"group": ["case"] * size + ["control"] * (len(values) - size)}, index=values.index)
    settings = dict(group_column="group", reference_group="control", pseudocount=1e-6,
                    fdr_method="by", alpha=.05, min_prevalence=.1)
    settings.update(options)
    return group_tests(values, samples, **settings)


def u_statistic(x, y):
    return sum((first > second) + .5 * (first == second) for first in x for second in y)


def permutation_p(x, y):
    pooled = list(x) + list(y)
    center = len(x) * len(y) / 2
    observed_u = abs(u_statistic(x, y) - center)
    statistics = []
    for chosen in combinations(range(len(pooled)), len(x)):
        left = [pooled[index] for index in chosen]
        right = [pooled[index] for index in range(len(pooled)) if index not in chosen]
        statistics.append(abs(u_statistic(left, right) - center))
    return sum(value >= observed_u for value in statistics) / len(statistics)


def test_documented_small_untied_u_example_and_effect_direction():
    # SciPy documentation: males U=17, n=5 vs 4, exact two-sided p=1/9.
    first, second = [19, 22, 16, 29, 24], [20, 11, 17, 12]
    values = observed([[value / 100, -value / 100] for value in first + second])
    result, info = compare(values, 5)
    result = result.set_index("feature_id")
    assert u_statistic(first, second) == 17
    assert permutation_p(first, second) == pytest.approx(1 / 9)
    assert result.loc["f0", "u_statistic"] == 17
    assert result.loc["f0", "p_value"] == pytest.approx(1 / 9)
    assert result.loc["f0", "rank_biserial_effect"] == pytest.approx(.7)
    assert result.loc["f0", "median_clr_difference"] == pytest.approx(.075)
    assert result.loc["f1", "rank_biserial_effect"] == pytest.approx(-.7)
    assert info["small_tied_features"] == 0


@pytest.mark.parametrize("first,second", [([0, 2, 4], [1, 3, 5]), ([0, 1, 2], [3, 4, 5])])
def test_exact_small_sample_p_values_match_exhaustive_label_permutations(first, second):
    result, _ = compare(observed([[value / 10, -value / 10] for value in first + second]), len(first))
    np.testing.assert_allclose(result["p_value"], permutation_p(first, second))


def test_small_tied_profiles_require_a_supported_permutation_method():
    first, second = [0, 0, 1], [1, 1, 1]
    variance = 3 * 3 / 12 * (7 - ((2**3 - 2) + (4**3 - 4)) / (6 * 5))
    z = (abs(u_statistic(first, second) - 4.5) - .5) / math.sqrt(variance)
    approximate = math.erfc(z / math.sqrt(2))
    assert not np.isclose(approximate, permutation_p(first, second))
    with pytest.raises(ComparisonError, match="permutation.*10"):
        compare(observed([[value, -value] for value in first + second]), 3)


def test_small_unequal_tied_null_configuration_is_rejected():
    # Supplemental exact null calculation found BH FDR 0.08207 for n=3 vs 30.
    coordinates = [[value, -value] for value in [1, 0, 0] + [0] * 30]
    with pytest.raises(ComparisonError, match="permutation.*10"):
        compare(observed(coordinates), 3, fdr_method="bh")


def test_sufficiently_large_tied_samples_match_independent_normal_calculation():
    first, second = [0] * 8 + [1] * 2, [0] * 2 + [1] * 8
    result, info = compare(observed([[value, -value] for value in first + second]), 10)
    variance = 10 * 10 / 12 * (21 - 2 * (10**3 - 10) / (20 * 19))
    z = (abs(u_statistic(first, second) - 50) - .5) / math.sqrt(variance)
    np.testing.assert_allclose(result["p_value"], math.erfc(z / math.sqrt(2)))
    assert info["small_tied_features"] == 0


@pytest.mark.parametrize("method,multiplier", [("bh", 1), ("by", sum(1 / i for i in range(1, 6)))])
def test_fdr_correction_includes_all_tested_features(monkeypatch, method, multiplier):
    import scipy.stats

    probabilities = iter([.001, .01, .02, .06, .5])
    monkeypatch.setattr(scipy.stats, "mannwhitneyu", lambda *a, **k: SimpleNamespace(statistic=4.5, pvalue=next(probabilities)))
    coordinates = np.random.default_rng(613).uniform(-1, 1, size=(6, 5))
    result, info = compare(observed(coordinates), 3, fdr_method=method)
    result = result.set_index("feature_id").loc[[f"f{i}" for i in range(5)]]
    expected = np.minimum(np.array([.005, .025, 1 / 30, .075, .5]) * multiplier, 1)
    np.testing.assert_allclose(result["q_value"], expected)
    assert info["features_tested"] == 5


def test_prevalence_filter_changes_clr_reference_but_not_pseudocount_scale():
    values = pd.DataFrame({"a": [.1, .2, .3, .6, .65, .7], "b": [.2, .8, .7, .4, .35, .3],
                           "rare": [.7, 0, 0, 0, 0, 0]})
    result, info = compare(values, 3, reference_group="case", pseudocount=.01, min_prevalence=.5)
    result = result.set_index("feature_id")
    expected_clr = .5 * np.log((values["a"] + .01) / (values["b"] + .01))
    expected = np.median(expected_clr.iloc[3:]) - np.median(expected_clr.iloc[:3])
    assert result.loc["a", "median_clr_difference"] == pytest.approx(expected)
    assert result.loc["b", "median_clr_difference"] == pytest.approx(-expected)
    assert info["features_tested"] == 2 and info["features_excluded"] == 1


def test_group_reversal_and_row_column_order_preserve_statistical_relationships():
    coordinates = np.random.default_rng(231).uniform(-1, 1, size=(10, 4))
    values = observed(coordinates)
    first, _ = compare(values, 5)
    reverse, _ = compare(values, 5, reference_group="case")
    samples = pd.DataFrame({"group": ["case"] * 5 + ["control"] * 5})
    reordered, _ = group_tests(values.iloc[::-1, ::-1], samples, group_column="group", reference_group="control",
                               pseudocount=1e-6, fdr_method="by", alpha=.05, min_prevalence=.1)
    first, reverse, reordered = [frame.set_index("feature_id").sort_index() for frame in (first, reverse, reordered)]
    for name in ("p_value", "q_value", "median_clr_difference", "rank_biserial_effect"):
        np.testing.assert_allclose(first[name], reordered[name], atol=1e-14)
    np.testing.assert_allclose(first["p_value"], reverse["p_value"])
    np.testing.assert_allclose(first["q_value"], reverse["q_value"])
    np.testing.assert_allclose(first["median_clr_difference"], -reverse["median_clr_difference"])
    np.testing.assert_allclose(first["rank_biserial_effect"], -reverse["rank_biserial_effect"])
    np.testing.assert_allclose(first["u_statistic"] + reverse["u_statistic"], 25)
