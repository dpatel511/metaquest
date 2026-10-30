import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import sys
from types import SimpleNamespace

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "tools/benchmark_reference.py"
benchmark = SimpleNamespace(**runpy.run_path(str(SCRIPT)))


@pytest.fixture
def quast(tmp_path):
    (tmp_path / "genome_stats").mkdir()
    (tmp_path / "report.tsv").write_text("Assembly\tcontigs\nGenome fraction (%)\t77.551\n")
    # 40/80 bs bases and 80/100 bases in each other member: 760/980 overall.
    lines = ["reference chromosomes:"]
    for prefix in benchmark.MEMBERS:
        non_n, covered = (80, 40) if prefix == "bs" else (100, 80)
        lines.append(f"\t{prefix}_contig1 (total length: 100 bp, total length without N's: {non_n} bp, maximal covered length: {covered} bp)")
    (tmp_path / "genome_stats/genome_info.txt").write_text("\n".join(lines) + "\n")
    reference = "".join(f">{prefix}_contig1\n" + ("A" * 80 + "N" * 20 if prefix == "bs" else "A" * 100) + "\n"
                        for prefix in benchmark.MEMBERS)
    (tmp_path / "reference.fa").write_text(reference)
    (tmp_path / "quast_corrected_input").mkdir()
    (tmp_path / "quast_corrected_input/reference.fa").write_text(reference)
    return tmp_path


def test_coverage_is_base_weighted_and_excludes_reference_ns(quast):
    result = benchmark.aggregate(quast, quast / "reference.fa")
    assert result["per_member"]["bs"]["coverage_percent"] == 50
    assert result["per_member"]["ef"]["coverage_percent"] == 80
    assert result["total_covered_non_N_bases"] == 760
    assert result["total_reference_non_N_bases"] == 980
    assert not result["whole_checkpoint_passed"]


@pytest.mark.parametrize("old, new, message", [
    ("bs_contig1", "unknown_contig1", "Unknown"),
    ("maximal covered length: 40", "maximal covered length: 81", "Invalid"),
    ("total length without N's: 80", "total length without N's: 101", "Invalid"),
    ("ec_contig1", "ef_contig1", "Duplicate"),
])
def test_invalid_reference_rows_rejected(quast, old, new, message):
    path = quast / "genome_stats/genome_info.txt"
    path.write_text(path.read_text().replace(old, new))
    with pytest.raises(ValueError, match=message):
        benchmark.aggregate(quast, quast / "reference.fa")


def test_multi_assembly_maxima_cannot_be_treated_as_one_assembly(quast):
    (quast / "report.tsv").write_text("Assembly\ta\tb\nGenome fraction (%)\t77.551\t80\n")
    with pytest.raises(ValueError, match="Exactly one"):
        benchmark.aggregate(quast, quast / "reference.fa")


def test_aggregate_disagreement_rejected(quast):
    (quast / "report.tsv").write_text("Assembly\tcontigs\nGenome fraction (%)\t80.000\n")
    with pytest.raises(ValueError, match="disagrees"):
        benchmark.aggregate(quast, quast / "reference.fa")


def test_reference_cli_checks_checksum_and_keeps_diagnostics_on_stderr(quast):
    reference = quast / "reference.fa"
    protocol = quast / "protocol.json"
    command = [sys.executable, str(SCRIPT), str(quast), "--reference", str(reference), "--protocol", str(protocol)]
    protocol.write_text(json.dumps({"reference": {"md5": "incorrect"}}))
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 2 and result.stdout == ""
    assert "Reference MD5" in result.stderr
    protocol.write_text(json.dumps({"reference": {"md5": hashlib.md5(reference.read_bytes()).hexdigest()}}))
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0 and result.stderr == ""
    assert json.loads(result.stdout)["total_reference_non_N_bases"] == 980


def test_reference_cli_rejects_malformed_protocol_without_traceback(quast):
    reference = quast / "reference.fa"
    protocol = quast / "protocol.json"
    protocol.write_text("[]")
    result = subprocess.run([sys.executable, str(SCRIPT), str(quast), "--reference", str(reference),
                             "--protocol", str(protocol)], capture_output=True, text=True)
    assert result.returncode == 2 and result.stdout == ""
    assert "Reference assessment rejected" in result.stderr and "Traceback" not in result.stderr


def test_quast_output_must_use_the_same_reference_content(quast):
    corrected = quast / "quast_corrected_input/reference.fa"
    corrected.write_text(corrected.read_text().replace("AAAA", "CCCC", 1))
    with pytest.raises(ValueError, match="reference content"):
        benchmark.aggregate(quast, quast / "reference.fa")
