import json

import numpy as np
import pandas as pd
import pytest

from metaquest.core.comparison import ComparisonError, group_tests, load_samples, run_comparison


def study(root, *, count=6, null=False):
    rows = []
    for number in range(count):
        name = f"{number:03d}"
        directory = root / name
        directory.mkdir()
        group = "control" if number < count // 2 else "case"
        high = 50 if null else 20 + number if group == "control" else 80 + number
        (directory / "bracken_report.tsv").write_text(
            "name\ttaxonomy_id\ttaxonomy_lvl\tnew_est_reads\n"
            f"First\t1\tS\t{high}\nSecond\t2\tS\t{100 - high}\n")
        summary = {"taxonomy": {"count_unit": "fragments", "total_input": 120,
                                  "classification_rate": 1, "top_taxa": [{"name": "First"}]}}
        metadata = {"status": "completed", "input_files": [{"sha256": f"input-{number}"}],
                    "effective_config": {"classification": {"taxonomic_level": "S"}},
                    "tool_provenance": {"kraken2": "2.1.6", "bracken": "3.1",
                        "kraken_database": {"manifest": {"release": "test-v1", "archive_md5": "fixture"}}}}
        (directory / "analysis_summary.json").write_text(json.dumps(summary))
        (directory / "analysis_metadata.json").write_text(json.dumps(metadata))
        rows.append({"sample_id": name, "result_dir": name, "group": group})
    sheet = root / "samples.tsv"
    pd.DataFrame(rows).to_csv(sheet, sep="\t", index=False)
    return sheet


def test_comparison_joins_full_profiles_and_preserves_sample_ids(tmp_path):
    sheet = study(tmp_path)
    samples = load_samples(sheet, "group")
    assert samples.index[0] == "000"
    output = tmp_path / "comparison"
    path = run_comparison(sheet, output, group_column="group")
    result = json.loads(path.read_text())
    assert result["feature_count"] == 2  # not the single taxon in the displayed top-N list
    matrix = pd.read_csv(output / "feature_values.tsv", sep="\t", dtype={"sample_id": str})
    assert matrix["sample_id"].tolist() == samples.index.tolist()
    assert (output / "report.html").is_file()
    distances = pd.read_csv(output / "bray_curtis.tsv", sep="\t", index_col=0).to_numpy()
    assert np.allclose(distances, distances.T) and np.allclose(np.diag(distances), 0)


def test_formal_test_requires_explicit_design_parameters(tmp_path):
    sheet = study(tmp_path)
    with pytest.raises(ComparisonError, match="require.*--group-column"):
        run_comparison(sheet, tmp_path / "out", test="clr-mannwhitney")


def test_group_test_effect_direction_and_fdr():
    samples = pd.DataFrame({"group": ["control"] * 5 + ["case"] * 5}, index=[f"s{i}" for i in range(10)])
    relative = pd.DataFrame({"one": [.15, .16, .17, .18, .19, .75, .76, .77, .78, .79]}, index=samples.index)
    relative["two"] = 1 - relative["one"]
    result, info = group_tests(relative, samples, group_column="group", reference_group="control",
                                pseudocount=1e-6, fdr_method="by", alpha=.05, min_prevalence=.1)
    result = result.set_index("feature_id")
    assert result.loc["one", "median_clr_difference"] > 0
    assert result.loc["two", "median_clr_difference"] < 0
    assert result["p_value"].between(0, 1, inclusive="neither").all()
    assert (result["q_value"] >= result["p_value"]).all()
    assert result["passes_fdr"].all()
    assert info["fdr_method"] == "by"


def test_identical_profiles_do_not_produce_significance(tmp_path):
    sheet = study(tmp_path, null=True)
    run_comparison(sheet, tmp_path / "out", group_column="group", test="clr-mannwhitney",
                   reference_group="control", pseudocount=1e-6)
    tests = pd.read_csv(tmp_path / "out/group_tests.tsv", sep="\t")
    assert (tests["p_value"] == 1).all() and not tests["passes_fdr"].any()
    assert (tests["u_statistic"] == 4.5).all()
    assert (tests["rank_biserial_effect"] == 0).all()


def test_undefined_nonconstant_p_values_are_rejected(monkeypatch):
    from types import SimpleNamespace
    import scipy.stats

    monkeypatch.setattr(scipy.stats, "mannwhitneyu", lambda *a, **k: SimpleNamespace(statistic=0., pvalue=float("nan")))
    samples = pd.DataFrame({"group": ["a"] * 3 + ["b"] * 3})
    values = pd.DataFrame({"one": [.1, .2, .3, .6, .7, .8], "two": [.9, .8, .7, .4, .3, .2]})
    with pytest.raises(ComparisonError, match="Undefined Mann-Whitney"):
        group_tests(values, samples, group_column="group", reference_group="a", pseudocount=1e-6,
                    fdr_method="by", alpha=.05, min_prevalence=.1)


def test_duplicate_samples_and_provenance_mismatch_are_rejected(tmp_path):
    sheet = study(tmp_path)
    frame = pd.read_csv(sheet, sep="\t", dtype=str)
    frame.loc[1, "result_dir"] = frame.loc[0, "result_dir"]
    frame.to_csv(sheet, sep="\t", index=False)
    with pytest.raises(ComparisonError, match="same result directory"):
        load_samples(sheet)
    frame.loc[1, "result_dir"] = "001"
    frame.to_csv(sheet, sep="\t", index=False)
    path = tmp_path / "001/analysis_metadata.json"
    metadata = json.loads(path.read_text())
    metadata["tool_provenance"]["kraken_database"]["manifest"]["release"] = "different"
    path.write_text(json.dumps(metadata))
    with pytest.raises(ComparisonError, match="Incompatible"):
        run_comparison(sheet, tmp_path / "out")


def test_technical_replicates_are_not_independent(tmp_path):
    sheet = study(tmp_path)
    source = tmp_path / "000/analysis_metadata.json"
    repeated = tmp_path / "001/analysis_metadata.json"
    repeated.write_text(source.read_text())
    with pytest.raises(ComparisonError, match="independent biological replicates"):
        run_comparison(sheet, tmp_path / "out", group_column="group", test="clr-mannwhitney",
                       reference_group="control", pseudocount=1e-6)


def test_repeated_subjects_and_small_groups_are_rejected():
    samples = pd.DataFrame({"group": ["a", "a", "b", "b"], "subject_id": ["x", "y", "x", "y"]})
    values = pd.DataFrame({"one": [.2, .3, .6, .7], "two": [.8, .7, .4, .3]})
    options = dict(group_column="group", reference_group="a", pseudocount=1e-6,
                   fdr_method="by", alpha=.05, min_prevalence=.1)
    with pytest.raises(ComparisonError, match="Repeated subject"):
        group_tests(values, samples, **options)
    with pytest.raises(ComparisonError, match="at least three"):
        group_tests(values, samples.drop(columns="subject_id"), **options)


def test_comparison_allows_thread_differences_but_not_preprocessing_changes(tmp_path):
    sheet = study(tmp_path)
    for index in range(6):
        path = tmp_path / f"{index:03d}/analysis_metadata.json"
        metadata = json.loads(path.read_text())
        metadata["effective_config"]["classification"]["threads"] = index + 1
        path.write_text(json.dumps(metadata))
    run_comparison(sheet, tmp_path / "allowed")
    path = tmp_path / "001/analysis_metadata.json"
    metadata = json.loads(path.read_text())
    metadata["effective_config"]["preprocessing"] = {"qualified_quality_phred": 30}
    path.write_text(json.dumps(metadata))
    with pytest.raises(ComparisonError, match="Incompatible"):
        run_comparison(sheet, tmp_path / "rejected")


def test_comparison_rejects_overwriting_results_or_sample_parent(tmp_path):
    sheet = study(tmp_path)
    with pytest.raises(ComparisonError, match="must not replace"):
        run_comparison(sheet, tmp_path)
    output = tmp_path / "comparison"
    output.mkdir()
    (output / "feature_values.tsv").write_text("original")
    with pytest.raises(ComparisonError, match="already contains"):
        run_comparison(sheet, output)
    assert (output / "feature_values.tsv").read_text() == "original"


def test_functional_comparison_joins_namespace_terms_and_checks_attribution(tmp_path):
    sheet = study(tmp_path, count=2)
    for number in range(2):
        directory = tmp_path / f"{number:03d}"
        (directory / "functional.tsv").write_text(
            f"namespace\tterm\ttpm\nKO\tK00001\t{10 + number}\nKO\tK00002\t20\nCOG\tA\t100\n")
        summary = json.loads((directory / "analysis_summary.json").read_text())
        summary["abundance"] = {"functional_abundance": "functional.tsv", "tool_version": "39.26",
                                "min_identity": .95, "ambiguous": "random", "count_unit": "fragments",
                                "functional": {"term_attribution": "full"}}
        (directory / "analysis_summary.json").write_text(json.dumps(summary))
        metadata = json.loads((directory / "analysis_metadata.json").read_text())
        metadata["effective_config"]["annotation"] = {"evalue": 1e-6, "threads": number + 1}
        metadata["tool_provenance"].update({"eggnog_mapper": "2.1.15", "pyrodigal": "3.7.1", "megahit": "1.2.9"})
        metadata["tool_provenance"]["functional_database"] = {
            "manifest": {"release": "5.0.2", "artifacts": [{"sha256": "test-db"}]}}
        (directory / "analysis_metadata.json").write_text(json.dumps(metadata))
    output = tmp_path / "KO"
    run_comparison(sheet, output, feature="KO")
    values = pd.read_csv(output / "feature_values.tsv", sep="\t")
    assert values.columns.tolist() == ["sample_id", "K00001", "K00002"]
    assert values["K00001"].tolist() == [10, 11]
    path = tmp_path / "001/analysis_summary.json"
    summary = json.loads(path.read_text())
    summary["abundance"]["functional"]["term_attribution"] = "split"
    path.write_text(json.dumps(summary))
    with pytest.raises(ComparisonError, match="Incompatible"):
        run_comparison(sheet, tmp_path / "mixed", feature="KO")
