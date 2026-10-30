"""Small calculation/contract fixtures; no biological accuracy claims."""
import copy
import csv
import json
from pathlib import Path
import runpy
import subprocess
import sys
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools/benchmark_taxonomy.py"
benchmark = SimpleNamespace(**runpy.run_path(str(SCRIPT)))
PROTOCOL = json.loads((ROOT / "benchmarks/zymo-even/protocol.json").read_text())


@pytest.fixture
def run(tmp_path):
    protocol = copy.deepcopy(PROTOCOL)
    path = tmp_path / "protocol.json"
    path.write_text(json.dumps(protocol))
    reports = [{"md5_checksum": item["md5"], "total_sequences": 1000000,
                "sampling_limited": True} for item in protocol["inputs"]]
    metadata = {"status": "completed", "workflow": {"read_mode": "paired"},
                "validation": {"status": "passed", "bypassed": False, "reports": reports},
                "tool_provenance": {"kraken_database": {"manifest": {
                    "name": "Kraken2 Standard-8 with Bracken", "release": "2026-06-26"}}}}
    (tmp_path / "analysis_metadata.json").write_text(json.dumps(metadata))
    (tmp_path / "analysis_summary.json").write_text(json.dumps({"taxonomy": {
        "total_input": 1000, "kraken_classified": 800, "unclassified_reads": 200, "species_assigned": 800}}))
    (tmp_path / "preprocessing").mkdir()
    (tmp_path / "preprocessing/fastp.json").write_text(json.dumps({"summary": {
        "before_filtering": {"total_reads": protocol["raw_read_pairs"] * 2},
        "after_filtering": {"total_reads": 2000}}}))
    (tmp_path / "kraken_report.txt").write_text("20\t200\t200\tU\t0\tunclassified\n80\t800\t0\tR\t1\troot\n")
    with (tmp_path / "bracken_report.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["taxonomy_id", "taxonomy_lvl", "name", "new_est_reads"])
        for member in protocol["expected_members"]:
            if member["bacterial"]:
                writer.writerow([member["taxids"][0], "S", member["name"], 100])
    return tmp_path, path


def alter(path, change):
    value = json.loads(path.read_text())
    change(value)
    path.write_text(json.dumps(value))


def test_full_input_is_not_confused_with_sampled_validation(run):
    root, protocol = run
    result = benchmark.evaluate(root, protocol)
    assert result["classification_input_fragments"] == 1000
    assert result["unclassified_fraction"] == .2
    assert result["metrics"][0]["bacterial_recall"] == 1
    assert result["metrics"][0]["community_recall"] == .8
    assert result["bacterial_baseline_gate_passed"]
    assert not result["whole_checkpoint_passed"]
    assert len(result["artifacts_sha256"]) == 5
    assert result["nominal_dna_l1_distance_not_read_accuracy"] == pytest.approx(.4)


def test_threshold_uses_all_fragments_not_only_classified():
    truth = [{"name": "a", "taxids": [1], "bacterial": True}]
    result = benchmark.score({1: {"count": 1, "name": "a"}}, 1000, truth, .01)
    assert result["bacterial_recall"] == 0
    assert result["community_precision"] is None
    assert result["expected"][0]["fraction_of_input"] == .001


def test_threshold_boundary_hybrid_group_and_low_mass_false_positives():
    truth = [{"name": "a", "taxids": [1], "bacterial": True},
             {"name": "hybrid", "taxids": [2, 3], "bacterial": False}]
    counts = {key: {"count": count, "name": str(key)}
              for key, count in ((1, 10), (2, 6), (3, 4), (4, 10), (5, 1))}
    result = benchmark.score(counts, 1000, truth, .01)
    assert result["community_true_positives"] == 2
    assert result["unexpected_reported_species"] == 1
    assert result["community_precision"] == pytest.approx(2 / 3)
    assert result["unexpected_fraction_of_input_all_abundances"] == .011
    assert len(result["unexpected"]) == 2


@pytest.mark.parametrize("value", [float("nan"), float("inf"), 0, -1, 1.1])
def test_invalid_threshold_rejected(value):
    with pytest.raises(ValueError, match="threshold"):
        benchmark.score({}, 1000, [{"name": "a", "taxids": [1], "bacterial": True}], value)


def test_overlapping_truth_rejected():
    truth = [{"name": "a", "taxids": [1], "bacterial": True},
             {"name": "b", "taxids": [1], "bacterial": False}]
    with pytest.raises(ValueError, match="Overlapping"):
        benchmark.score({}, 1000, truth, .001)


def test_zero_calls_have_undefined_precision_and_zero_recall():
    truth = [{"name": "a", "taxids": [1], "bacterial": True}]
    result = benchmark.score({}, 1000, truth, .001)
    assert result["community_precision"] is None
    assert result["community_recall"] == 0


def test_counts_cannot_exceed_input():
    with pytest.raises(ValueError, match="exceed"):
        benchmark.score({1: {"count": 1001}}, 1000,
                        [{"name": "a", "taxids": [1], "bacterial": True}], .001)


@pytest.mark.parametrize("rows, message", [
    (["1\tS\ta\t1", "1\tS\ta\t2"], "Duplicate"),
    (["1\tG\ta\t1"], "species-level"),
    (["1\tS\ta\t-1"], "Invalid"),
    (["1\tS\ta\tnan"], "invalid literal"),
    (["1\tS\ta\t1.5"], "invalid literal"),
])
def test_bracken_invalid_rows_rejected(tmp_path, rows, message):
    path = tmp_path / "bracken.tsv"
    path.write_text("taxonomy_id\ttaxonomy_lvl\tname\tnew_est_reads\n" + "\n".join(rows))
    with pytest.raises(ValueError, match=message):
        benchmark.read_counts(path)


@pytest.mark.parametrize("file, change, message", [
    ("analysis_metadata.json", lambda d: d.update(status="failed"), "not completed"),
    ("analysis_metadata.json", lambda d: d["workflow"].update(read_mode="single"), "paired"),
    ("analysis_metadata.json", lambda d: d["validation"].update(bypassed=True), "bypassed"),
    ("analysis_metadata.json", lambda d: d["validation"]["reports"][0].update(md5_checksum="wrong"), "MD5"),
    ("analysis_metadata.json", lambda d: d["tool_provenance"]["kraken_database"]["manifest"].update(release="old"), "release"),
    ("analysis_summary.json", lambda d: d["taxonomy"].update(species_assigned=1), "Bracken"),
    ("analysis_summary.json", lambda d: d["taxonomy"].update(kraken_classified=1), "Classified"),
    ("preprocessing/fastp.json", lambda d: d["summary"]["after_filtering"].update(total_reads=1000), "denominator"),
    ("preprocessing/fastp.json", lambda d: d["summary"]["before_filtering"].update(total_reads=2), "raw"),
])
def test_mismatched_provenance_and_accounting_rejected(run, file, change, message):
    root, protocol = run
    alter(root / file, change)
    with pytest.raises(ValueError, match=message):
        benchmark.evaluate(root, protocol)


def test_failure_gate_retains_results_and_nonzero_exit(run):
    root, protocol = run
    with (root / "bracken_report.tsv").open("a") as handle:
        handle.write("999\tS\tunexpected\t100\n")
    alter(root / "analysis_summary.json", lambda d: d["taxonomy"].update(species_assigned=900))
    (root / "kraken_report.txt").write_text("10\t100\t100\tU\t0\tu\n90\t900\t0\tR\t1\tr\n")
    alter(root / "analysis_summary.json", lambda d: d["taxonomy"].update(unclassified_reads=100))
    # Two unexpected detections cause precision 8/10, below the frozen 0.9 gate.
    with (root / "bracken_report.tsv").open("a") as handle:
        handle.write("998\tS\tunexpected2\t1\n")
    alter(root / "analysis_summary.json", lambda d: d["taxonomy"].update(species_assigned=901))
    (root / "kraken_report.txt").write_text("9.9\t99\t99\tU\t0\tu\n90.1\t901\t0\tR\t1\tr\n")
    alter(root / "analysis_summary.json", lambda d: d["taxonomy"].update(unclassified_reads=99, kraken_classified=901))
    result = subprocess.run([sys.executable, str(SCRIPT), str(root), "--protocol", str(protocol)],
                            capture_output=True, text=True)
    assert result.returncode == 1
    assert not json.loads(result.stdout)["bacterial_baseline_gate_passed"]
    assert result.stderr == ""


def test_invalid_run_cli_is_stderr(run):
    root, protocol = run
    (root / "kraken_report.txt").write_text("malformed\n")
    result = subprocess.run([sys.executable, str(SCRIPT), str(root), "--protocol", str(protocol)],
                            capture_output=True, text=True)
    assert result.returncode == 2
    assert result.stdout == ""
    assert "Malformed Kraken" in result.stderr


def test_unsampled_validation_must_have_full_count(run):
    root, protocol = run
    alter(root / "analysis_metadata.json", lambda d: d["validation"]["reports"][0].update(sampling_limited=False))
    with pytest.raises(ValueError, match="Unsampled"):
        benchmark.evaluate(root, protocol)


def test_bracken_cannot_exceed_classified(run):
    root, protocol = run
    with (root / "bracken_report.tsv").open("a") as handle:
        handle.write("999\tS\tunexpected\t1\n")
    with pytest.raises(ValueError, match="classified"):
        benchmark.evaluate(root, protocol)
