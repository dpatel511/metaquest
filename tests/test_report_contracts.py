"""Regression contracts for ambiguous labels and missing stage evidence."""

import json

import pandas as pd

from metaquest.pipeline.context import PipelineContext, ClassificationResult
from metaquest.reporting.stable_reporter import generate_stable_reports
from metaquest.reporting.steps import write_step_reports
from metaquest.settings import load_config


def test_taxonomy_plot_keeps_same_named_taxa_separate(tmp_path):
    table = tmp_path / "bracken_report.tsv"
    table.write_text(
        "name\ttaxonomy_id\tnew_est_reads\nUnknown species\t10\t70\nUnknown species\t20\t30\n"
    )
    ctx = PipelineContext(load_config(db_dir=tmp_path / "db"), [], tmp_path)
    ctx.classification = ClassificationResult(
        tmp_path / "kraken",
        table,
        total_reads=100,
        kraken_classified_reads=100,
        species_assigned_reads=100,
    )
    generate_stable_reports(ctx)
    summary = json.loads((tmp_path / "analysis_summary.json").read_text())
    assert [r["taxonomy_id"] for r in summary["taxonomy"]["top_taxa"]] == ["10", "20"]
    data = pd.read_csv(tmp_path / "plots/plot_data/taxonomy_top.tsv", sep="\t")
    assert data["taxon"].is_unique
    assert data["percent_all_input"].tolist() == [70, 30]


def test_comparison_composition_has_unambiguous_ids(tmp_path):
    from metaquest.reporting.comparison_report import write_comparison_report

    relative = pd.DataFrame({"10": [0.7, 0.2], "20": [0.3, 0.8]}, index=["001", "002"])
    distances = pd.DataFrame([[0, 0.5], [0.5, 0]], index=relative.index, columns=relative.index)
    summary = {
        "feature": "taxonomy",
        "sample_count": 2,
        "feature_count": 2,
        "input_unit": "counts",
        "statistics": {"method": "descriptive only"},
        "relative_denominator": "reported total",
        "limitations": [],
    }
    write_comparison_report(
        tmp_path,
        relative,
        {"10": "Same", "20": "Same"},
        pd.DataFrame({"sample": ["001", "002"]}),
        distances,
        None,
        summary,
    )
    by_id = pd.read_csv(
        tmp_path / "plots/composition_by_id.tsv", sep="\t", dtype={"sample_id": str}
    ).set_index("sample_id")
    pd.testing.assert_frame_equal(by_id[["10", "20"]], relative, check_names=False)
    labeled = pd.read_csv(tmp_path / "plots/composition.tsv", sep="\t")
    assert {"Same [10]", "Same [20]"}.issubset(labeled.columns)
    heatmap = (tmp_path / "plots/feature_heatmap.svg").read_text()
    assert "Same [10]" in heatmap and "Same [20]" in heatmap


def test_stage_payload_preserves_null_zero_and_escaped_text(tmp_path):
    summary = {
        "validation": {"status": "skipped", "bypassed": True},
        "annotation": {"functional_status": "skipped", "annotated_genes": None},
        "screening": {
            "status": "completed",
            "protein_matches": 0,
            "taxon_candidates": 0,
            "taxon_table": "taxa.tsv",
            "protein_table": "proteins.tsv",
            "label": "<script>bad</script>",
        },
    }
    for name in ("taxa.tsv", "proteins.tsv"):
        (tmp_path / name).write_text("gene_id\tlabel\n")
    reports = write_step_reports(tmp_path, summary)
    assert reports[0]["status"] == "skipped" and reports[5]["status"] == "skipped"
    assert reports[3]["status"] == "not_run" and reports[7]["status"] == "completed"
    payload = json.loads((tmp_path / "reports/07_screening.json").read_text())
    assert payload["metrics"]["protein_matches"] == 0
    page = (tmp_path / "reports/07_screening.html").read_text()
    assert "<script>bad</script>" not in page and "&lt;script&gt;" in page
