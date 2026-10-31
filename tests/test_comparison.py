import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

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
                                  "classification_rate": 110 / 120, "kraken_classified": 110,
                                  "species_assigned": 100, "top_taxa": [{"name": "First"}]}}
        metadata = {"status": "completed", "input_files": [{"sha256": hashlib.sha256(f"input-{number}".encode()).hexdigest()}],
                    "effective_config": {"classification": {"taxonomic_level": "S"},
                        "preprocessing": {"enabled": True, "qualified_quality_phred": 20},
                        "assembly": {"assembler": "megahit", "threads": 4}},
                    "tool_provenance": {"kraken2": "2.1.6", "bracken": "3.1",
                        "fastp": "1.3.6", "bracken_model_read_length": 150,
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
    metadata["effective_config"]["preprocessing"]["qualified_quality_phred"] = 30
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


# Checkpoint 4: malformed studies must fail before exporting a comparison.
def update_json(path, update):
    value = json.loads(path.read_text())
    update(value)
    path.write_text(json.dumps(value))


def formal(sheet, output):
    return run_comparison(sheet, output, group_column="group", test="clr-mannwhitney",
                          reference_group="control", pseudocount=1e-6)


def test_subject_whitespace_cannot_hide_repeated_subjects(tmp_path):
    sheet = study(tmp_path)
    samples = pd.read_csv(sheet, sep="\t", dtype=str)
    samples["subject_id"] = [" person01", "person01 ", "p3", "p4", "p5", "p6"]
    samples.to_csv(sheet, sep="\t", index=False)
    with pytest.raises(ComparisonError, match="Repeated subject"):
        formal(sheet, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_shared_mate_cannot_be_treated_as_an_independent_sample(tmp_path):
    sheet = study(tmp_path)
    shared = json.loads((tmp_path / "000/analysis_metadata.json").read_text())["input_files"][0]
    update_json(tmp_path / "001/analysis_metadata.json", lambda d: d["input_files"].append(shared))
    with pytest.raises(ComparisonError, match="independent biological replicates"):
        formal(sheet, tmp_path / "out")


@pytest.mark.parametrize("value", [None, "", "not-a-sha256", 123])
def test_formal_input_checksums_must_be_complete_sha256(tmp_path, value):
    sheet = study(tmp_path)
    update_json(tmp_path / "001/analysis_metadata.json", lambda d: d["input_files"][0].update(sha256=value))
    with pytest.raises(ComparisonError, match="checksums"):
        formal(sheet, tmp_path / "out")


@pytest.mark.parametrize("change", [
    lambda d: d["tool_provenance"].update(bracken_model_read_length=100),
    lambda d: d["tool_provenance"].update(fastp="different"),
])
def test_effective_model_and_preprocessor_versions_must_match(tmp_path, change):
    sheet = study(tmp_path)
    update_json(tmp_path / "001/analysis_metadata.json", change)
    with pytest.raises(ComparisonError, match="Incompatible"):
        run_comparison(sheet, tmp_path / "out")


@pytest.mark.parametrize("key", ["bracken_model_read_length", "fastp"])
def test_missing_required_comparison_provenance_is_rejected(tmp_path, key):
    sheet = study(tmp_path)
    for number in range(6):
        update_json(tmp_path / f"{number:03d}/analysis_metadata.json", lambda d: d["tool_provenance"].pop(key))
    with pytest.raises(ComparisonError, match="provenance.*incomplete"):
        run_comparison(sheet, tmp_path / "out")


@pytest.mark.parametrize("value", ["NaN", "inf", "-1", "1.5", "1e308"])
def test_invalid_taxonomy_quantities_rejected_before_output(tmp_path, value):
    sheet = study(tmp_path)
    (tmp_path / "000/bracken_report.tsv").write_text(
        f"name\ttaxonomy_id\ttaxonomy_lvl\tnew_est_reads\nFirst\t1\tS\t{value}\nSecond\t2\tS\t{value}\n")
    with pytest.raises(ComparisonError, match="finite|counts|total"):
        run_comparison(sheet, tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("taxid", ["", " ", "-1", "0", "1.2", "not-a-taxid"])
def test_invalid_taxonomy_identifiers_are_rejected(tmp_path, taxid):
    sheet = study(tmp_path)
    path = tmp_path / "000/bracken_report.tsv"
    path.write_text(path.read_text().replace("First\t1\t", f"First\t{taxid}\t"))
    with pytest.raises(ComparisonError, match="feature IDs|taxonomy IDs"):
        run_comparison(sheet, tmp_path / "out")


def test_taxonomy_profile_total_must_reconcile_with_summary(tmp_path):
    sheet = study(tmp_path)
    update_json(tmp_path / "000/analysis_summary.json", lambda d: d["taxonomy"].update(species_assigned=99))
    with pytest.raises(ComparisonError, match="counts.*summary"):
        run_comparison(sheet, tmp_path / "out")


def test_taxonomy_rank_must_match_recorded_settings(tmp_path):
    sheet = study(tmp_path)
    path = tmp_path / "000/bracken_report.tsv"
    path.write_text(path.read_text().replace("\tS\t", "\tG\t"))
    with pytest.raises(ComparisonError, match="rank"):
        run_comparison(sheet, tmp_path / "out")


@pytest.mark.parametrize("name,content", [
    ("analysis_metadata.json", "[]"),
    ("bracken_report.tsv", "name\ttaxonomy_id\tnew_est_reads\tnew_est_reads\nFirst\t1\t100\t9\n"),
])
def test_malformed_profile_contracts_raise_actionable_errors(tmp_path, name, content):
    sheet = study(tmp_path)
    (tmp_path / "000" / name).write_text(content)
    with pytest.raises(ComparisonError, match="JSON object|Duplicate columns"):
        run_comparison(sheet, tmp_path / "out")


def test_missing_profile_is_an_actionable_comparison_error(tmp_path):
    sheet = study(tmp_path)
    (tmp_path / "000/bracken_report.tsv").unlink()
    with pytest.raises(ComparisonError, match="Cannot read"):
        run_comparison(sheet, tmp_path / "out")


def test_duplicate_sample_sheet_headers_are_rejected(tmp_path):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("sample_id\tresult_dir\tgroup\tgroup\na\ta\tx\ty\nb\tb\tx\ty\n")
    with pytest.raises(ComparisonError, match="Duplicate columns"):
        load_samples(sheet, "group")


def test_conflicting_feature_labels_are_not_silently_overwritten(tmp_path):
    sheet = study(tmp_path)
    path = tmp_path / "001/bracken_report.tsv"
    path.write_text(path.read_text().replace("First", "Different species"))
    with pytest.raises(ComparisonError, match="Conflicting.*label"):
        run_comparison(sheet, tmp_path / "out")


def test_comparison_rejects_nonempty_unrelated_output(tmp_path):
    sheet = study(tmp_path)
    output = tmp_path / "out"
    output.mkdir()
    (output / "report.html").write_text("retain original")
    with pytest.raises(ComparisonError, match="already contains"):
        run_comparison(sheet, output)
    assert (output / "report.html").read_text() == "retain original"


def test_text_report_records_the_actual_transformation_and_contrast(tmp_path):
    sheet = study(tmp_path)
    formal(sheet, tmp_path / "out")
    report = (tmp_path / "out/comparison_report.txt").read_text()
    for value in ("Pseudocount: 1e-06", "Prevalence threshold: 0.1", "case - control", "natural log", "no covariate adjustment"):
        assert value in report


@pytest.mark.parametrize("kind", ["extra_sample", "missing_sample", "duplicate_sample", "duplicate_feature", "negative", "nan", "not_proportions"])
def test_group_test_input_alignment_and_quantities_are_validated(kind):
    samples = pd.DataFrame({"group": ["a"] * 3 + ["b"] * 3}, index=[f"s{i}" for i in range(6)])
    values = pd.DataFrame({"one": [.1, .2, .3, .6, .7, .8], "two": [.9, .8, .7, .4, .3, .2]}, index=samples.index)
    if kind == "extra_sample":
        values.loc["extra"] = [.5, .5]
    elif kind == "missing_sample":
        values = values.iloc[:-1]
    elif kind == "duplicate_sample":
        values = pd.concat([values, values.iloc[:1]])
    elif kind == "duplicate_feature":
        values.columns = ["same", "same"]
    elif kind == "negative":
        values.iloc[0] = [-.1, 1.1]
    elif kind == "nan":
        values.iloc[0, 0] = np.nan
    else:
        values = values * 2
    with pytest.raises(ComparisonError, match="sample|feature|proportion|finite"):
        group_tests(values, samples, group_column="group", reference_group="a", pseudocount=1e-6,
                    fdr_method="by", alpha=.05, min_prevalence=.1)


def test_partial_profile_overlap_is_zero_filled_by_id_without_losing_samples(tmp_path):
    sheet = study(tmp_path, count=2)
    (tmp_path / "001/bracken_report.tsv").write_text(
        "name\ttaxonomy_id\ttaxonomy_lvl\tnew_est_reads\nThird\t3\tS\t40\nFirst\t1\tS\t60\n")
    run_comparison(sheet, tmp_path / "out")
    matrix = pd.read_csv(tmp_path / "out/feature_values.tsv", sep="\t", dtype={"sample_id": str})
    assert matrix.columns.tolist() == ["sample_id", "1", "2", "3"]
    assert matrix["sample_id"].tolist() == ["000", "001"]
    np.testing.assert_allclose(matrix.iloc[:, 1:].to_numpy(), [[20, 80, 0], [60, 0, 40]])
    relative = pd.read_csv(tmp_path / "out/relative_abundance.tsv", sep="\t")
    np.testing.assert_allclose(relative.iloc[:, 1:].to_numpy(), [[.2, .8, 0], [.6, 0, .4]])


def test_row_scaling_does_not_change_relative_abundances_or_formal_results(tmp_path):
    sheet = study(tmp_path)
    formal(sheet, tmp_path / "original")
    for number in range(6):
        path = tmp_path / f"{number:03d}/bracken_report.tsv"
        frame = pd.read_csv(path, sep="\t")
        frame["new_est_reads"] *= number + 1
        frame.to_csv(path, sep="\t", index=False)
        update_json(path.parent / "analysis_summary.json", lambda d: d["taxonomy"].update(
            species_assigned=100 * (number + 1), kraken_classified=110 * (number + 1), total_input=120 * (number + 1)))
    formal(sheet, tmp_path / "scaled")
    for name in ("relative_abundance.tsv", "group_tests.tsv"):
        pd.testing.assert_frame_equal(pd.read_csv(tmp_path / "original" / name, sep="\t"),
                                      pd.read_csv(tmp_path / "scaled" / name, sep="\t"))


@pytest.mark.parametrize("value", ["nan", "inf", "-1"])
def test_invalid_functional_tpm_is_rejected(tmp_path, value):
    sheet = study(tmp_path, count=2)
    # Reuse complete functional provenance from the supported fixture path.
    for number in range(2):
        directory = tmp_path / f"{number:03d}"
        (directory / "functional.tsv").write_text(f"namespace\tterm\ttpm\nKO\tK00001\t{value}\nKO\tK00002\t20\n")
        update_json(directory / "analysis_summary.json", lambda d: d.update(abundance={
            "functional_abundance": "functional.tsv", "tool_version": "39.26", "min_identity": .95,
            "ambiguous": "random", "count_unit": "fragments", "functional": {"term_attribution": "full"}}))
        update_json(directory / "analysis_metadata.json", lambda d: (
            d["effective_config"].update(annotation={"evalue": 1e-6}),
            d["tool_provenance"].update(eggnog_mapper="2.1.15", pyrodigal="3.7.1", megahit="1.2.9",
                functional_database={"manifest": {"release": "5.0.2", "artifacts": [{"sha256": "test-db"}]}})))
    with pytest.raises(ComparisonError, match="finite nonnegative"):
        run_comparison(sheet, tmp_path / "out", feature="KO")


@pytest.mark.parametrize("force", [False, True])
def test_compare_cli_accepts_its_own_log_and_preserves_forced_backups(tmp_path, force):
    sheet = study(tmp_path)
    output = tmp_path / "cli-output"
    if force:
        output.mkdir()
        (output / "report.html").write_text("retain old report")
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(Path(__file__).resolve().parents[1] / "src")
    command = [sys.executable, "-c", "from metaquest.cli import main; main()", "--no-color", "compare",
               "--sample-sheet", str(sheet), "--output", str(output), "--group-column", "group",
               "--test", "clr-mannwhitney", "--reference-group", "control", "--pseudocount", "1e-6"]
    if force:
        command.append("--force")
    result = subprocess.run(command, cwd=tmp_path, env=environment, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert result.stdout == "" and "\x1b" not in result.stderr
    assert (output / "report.html").is_file() and (output / "metaquest.log").is_file()
    if force:
        backups = list(tmp_path.glob("cli-output.metaquest-backup-*"))
        assert len(backups) == 1
        assert (backups[0] / "report.html").read_text() == "retain old report"


def test_comparison_api_preserves_an_existing_log_only_directory(tmp_path):
    sheet = study(tmp_path)
    output = tmp_path / "out"
    output.mkdir()
    (output / "metaquest.log").write_text("keep diagnostic log")
    run_comparison(sheet, output)
    assert (output / "metaquest.log").read_text() == "keep diagnostic log"


@pytest.mark.parametrize("quiet", [False, True])
def test_compare_cli_rejects_unsupported_small_ties_without_statistics(tmp_path, quiet):
    sheet = study(tmp_path)
    for number in range(6):
        path = tmp_path / f"{number:03d}/bracken_report.tsv"
        frame = pd.read_csv(path, sep="\t")
        frame["new_est_reads"] = [20, 80] if number < 3 else [80, 20]
        frame.to_csv(path, sep="\t", index=False)
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(Path(__file__).resolve().parents[1] / "src")
    output = tmp_path / "invalid-study"
    command = [sys.executable, "-c", "from metaquest.cli import main; main()", "--no-color"]
    if quiet:
        command.append("--quiet")
    command += ["compare", "--sample-sheet", str(sheet), "--output", str(output), "--group-column", "group",
                "--test", "clr-mannwhitney", "--reference-group", "control", "--pseudocount", "1e-6"]
    result = subprocess.run(command, cwd=tmp_path, env=environment, capture_output=True, text=True, timeout=60)
    assert result.returncode == 1 and result.stdout == ""
    assert "permutation-based" in result.stderr and "at least 10" in result.stderr
    assert not (output / "group_tests.tsv").exists()
    assert not (output / "comparison_summary.json").exists()
