"""Independent count fixtures for the taxonomy evidence/report contract."""
import csv
from dataclasses import replace
import json

import pytest

from metaquest.pipeline.context import ClassificationResult, PipelineContext
from metaquest.reporting.stable_reporter import generate_stable_reports
from metaquest.settings import load_config


@pytest.fixture
def context(tmp_path, monkeypatch):
    monkeypatch.setattr("metaquest.reporting.plots.generate_plots", lambda *args: [])
    config = load_config(db_dir=tmp_path / "db")
    config = replace(config, reporting=replace(config.reporting, top_taxa=1))
    ctx = PipelineContext(config, [], tmp_path)
    kraken = tmp_path / "kraken_report.txt"
    kraken.write_text("40\t80\t80\tU\t0\tunclassified\n"
                      "60\t120\t0\tR\t1\troot\n"
                      "30\t60\t20\tS\t10\tSame name\n"
                      "20\t40\t40\tS1\t11\t  Child\n"
                      "10\t20\t20\tS\t20\tSame name\n"
                      "20\t40\t40\tG\t30\tOther genus\n")
    bracken = tmp_path / "bracken_report.tsv"
    bracken.write_text("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n"
                       "Same name\t20\tS\t20\t0\t20\t0.22222\n"
                       "Same name\t10\tS\t60\t10\t70\t0.77778\n")
    ctx.classification = ClassificationResult(kraken, bracken, species_count=2,
        species_assigned_reads=90, total_reads=200, kraken_classified_reads=120,
        unclassified_reads=80, observed_read_length=151, bracken_read_length=150)
    (tmp_path / "analysis_metadata.json").write_text(json.dumps({
        "effective_config": {"classification": {"taxonomic_level": "S", "min_hit_groups": 5,
                                                  "bracken_threshold": 12}},
        "tool_provenance": {"kraken2": "Kraken fixture", "bracken": "Bracken fixture",
            "kraken_database": {"configured_release": "fixture release", "manifest": {
                "release": "fixture release", "archive_md5": "a" * 32}}}}))
    return ctx


def evidence(ctx):
    with (ctx.output_dir / "taxonomy_evidence.tsv").open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.mark.parametrize("read_mode,unit", [("paired", "fragments"), ("interleaved", "fragments"), ("single", "reads")])
def test_taxonomy_evidence_joins_ids_not_labels_and_preserves_denominators(context, read_mode, unit):
    context.read_mode = read_mode
    generate_stable_reports(context)
    rows = {row["taxonomy_id"]: row for row in evidence(context)}
    assert set(rows) == {"10", "20"}  # Full table, although top_taxa is limited to one.
    assert rows["10"]["kraken_direct_count"] == "20"
    assert rows["10"]["kraken_clade_count"] == "60"
    assert rows["10"]["bracken_assigned_count"] == "60"  # Not the direct-assignment column.
    assert rows["10"]["bracken_added_count"] == "10"
    assert rows["10"]["bracken_estimated_count"] == "70"
    assert float(rows["10"]["fraction_all_input"]) == 0.35
    assert float(rows["10"]["fraction_bracken_profile"]) == pytest.approx(70 / 90)
    assert rows["20"]["bracken_added_count"] == "0"
    assert {row["count_unit"] for row in rows.values()} == {unit}
    summary = json.loads((context.output_dir / "analysis_summary.json").read_text())["taxonomy"]
    assert summary["species_assigned"] == summary["rank_assigned"] == 90
    assert summary["classified_above_species"] == summary["classified_not_in_bracken"] == 30
    assert summary["top_taxa"][0]["estimated_count"] == 70
    assert summary["evidence_schema_version"] == 1
    methods = json.loads((context.output_dir / "taxonomy_methods.json").read_text())
    assert methods["schema_version"] == 1
    assert methods["effective_settings"]["min_hit_groups"] == 5
    assert methods["effective_settings"]["bracken_threshold"] == 12
    assert methods["tool_versions"] == {"kraken2": "Kraken fixture", "bracken": "Bracken fixture"}
    assert methods["database"]["manifest"]["archive_md5"] == "a" * 32
    assert len(methods["sources"]["bracken_report"]["sha256"]) == 64
    for name in ("report.html", "reports/02_taxonomy.html"):
        page = (context.output_dir / name).read_text()
        assert "Kraken direct" in page and "Bracken estimate" in page
        assert "taxonomy_evidence.tsv" in page and "taxonomy_methods.json" in page
        assert "Classified above species" not in page


def test_missing_kraken_rows_are_unavailable_not_zero(context):
    path = context.classification.kraken_report
    path.write_text("\n".join(line for line in path.read_text().splitlines() if "\t20\tSame name" not in line) + "\n")
    generate_stable_reports(context)
    rows = {row["taxonomy_id"]: row for row in evidence(context)}
    assert rows["20"]["kraken_direct_count"] == rows["20"]["kraken_clade_count"] == ""
    assert rows["20"]["bracken_added_count"] == "0"
    methods = json.loads((context.output_dir / "taxonomy_methods.json").read_text())
    assert methods["taxa_without_kraken_row"] == 1


def test_legacy_missing_sources_and_html_names(context):
    context.classification.kraken_report.unlink()
    context.classification.bracken_file.write_text("name\tnew_est_reads\n<script>taxon</script>\t90\n")
    (context.output_dir / "analysis_metadata.json").unlink()
    generate_stable_reports(context)
    row, = evidence(context)
    assert row["taxonomy_id"] == row["kraken_direct_count"] == row["bracken_added_count"] == ""
    methods = json.loads((context.output_dir / "taxonomy_methods.json").read_text())
    assert methods["tool_versions"] == {"kraken2": None, "bracken": None}
    assert methods["database"] is None
    assert "<script>taxon</script>" not in (context.output_dir / "report.html").read_text()
    assert "&lt;script&gt;taxon&lt;/script&gt;" in (context.output_dir / "reports/02_taxonomy.html").read_text()


def test_literal_na_names_and_recorded_provenance_are_preserved(context):
    path = context.classification.bracken_file
    path.write_text(path.read_text().replace("Same name", "NA"))
    metadata = context.output_dir / "analysis_metadata.json"
    before = metadata.read_bytes()
    generate_stable_reports(context)
    assert {row["name"] for row in evidence(context)} == {"NA"}
    assert metadata.read_bytes() == before


@pytest.mark.parametrize("total", [0, 200])
def test_empty_profile_has_headers_and_explicit_zero_denominator(context, total):
    context.classification.bracken_file.write_text(context.classification.bracken_file.read_text().splitlines()[0] + "\n")
    context.classification.kraken_report.write_text(f"100\t{total}\t{total}\tU\t0\tunclassified\n0\t0\t0\tR\t1\troot\n")
    context.classification.total_reads = context.classification.unclassified_reads = total
    context.classification.kraken_classified_reads = context.classification.species_assigned_reads = 0
    generate_stable_reports(context)
    assert evidence(context) == []
    assert "bracken_estimated_count" in (context.output_dir / "taxonomy_evidence.tsv").read_text()
    assert "NaN" not in (context.output_dir / "analysis_summary.json").read_text()
    if not total:
        taxonomy = json.loads((context.output_dir / "analysis_summary.json").read_text())["taxonomy"]
        assert taxonomy["classified_fraction"] is taxonomy["unclassified_fraction"] is None


@pytest.mark.parametrize("rank,label", [("G", "genus"), ("F", "family")])
def test_non_species_rank_is_labelled_correctly(context, rank, label):
    context.config = replace(context.config, classification=replace(context.config.classification, taxonomic_level=rank))
    path = context.classification.bracken_file
    path.write_text(path.read_text().replace("\tS\t", f"\t{rank}\t"))
    metadata = context.output_dir / "analysis_metadata.json"
    payload = json.loads(metadata.read_text())
    payload["effective_config"]["classification"]["taxonomic_level"] = rank
    metadata.write_text(json.dumps(payload))
    generate_stable_reports(context)
    summary = json.loads((context.output_dir / "analysis_summary.json").read_text())["taxonomy"]
    assert summary["taxonomic_level"] == rank
    assert {row["rank"] for row in evidence(context)} == {rank}
    assert f"Bracken-assigned {label}" in (context.output_dir / "report.html").read_text()


@pytest.mark.parametrize("value", ["-1", "1.5", "NaN", "inf"])
def test_invalid_bracken_estimates_fail_instead_of_truncating(context, value):
    path = context.classification.bracken_file
    path.write_text(path.read_text().replace("\t70\t0.77778", f"\t{value}\t0.77778"))
    with pytest.raises(ValueError):
        generate_stable_reports(context)


def test_duplicate_taxonomy_ids_are_rejected(context):
    path = context.classification.bracken_file
    path.write_text(path.read_text().replace("Same name\t20\tS", "Same name\t10\tS"))
    with pytest.raises(ValueError, match="[Dd]uplicate"):
        generate_stable_reports(context)


def test_duplicate_kraken_ids_are_rejected(context):
    path = context.classification.kraken_report
    path.write_text(path.read_text() + "30\t60\t20\tS\t10\tSame name\n")
    with pytest.raises(ValueError, match="[Dd]uplicate"):
        generate_stable_reports(context)


@pytest.mark.parametrize("content", ["", "not a report\n", "60\t-120\t0\tR\t1\troot\n",
                                    "60\t120\t121\tR\t1\troot\n"])
def test_invalid_kraken_sources_are_rejected(context, content):
    context.classification.kraken_report.write_text(content)
    with pytest.raises(ValueError):
        generate_stable_reports(context)


def test_bracken_support_counts_must_reconcile(context):
    path = context.classification.bracken_file
    path.write_text(path.read_text().replace("\t60\t10\t70\t", "\t60\t9\t70\t"))
    with pytest.raises(ValueError, match="[Cc]ount|reconcile"):
        generate_stable_reports(context)


@pytest.mark.parametrize("field,value", [("species_assigned_reads", 91), ("total_reads", 201), ("kraken_classified_reads", 89)])
def test_invalid_accounting_fails_instead_of_clamping(context, field, value):
    setattr(context.classification, field, value)
    with pytest.raises(ValueError, match="[Cc]ount|denominator|reconcile"):
        generate_stable_reports(context)


@pytest.mark.parametrize("change", ["count", "fraction", "source_hash", "id"])
def test_independent_acceptance_checker_detects_tampered_evidence(context, change):
    from tools.verify_run import verify_taxonomy_evidence

    generate_stable_reports(context)
    taxonomy = json.loads((context.output_dir / "analysis_summary.json").read_text())["taxonomy"]
    assert verify_taxonomy_evidence(context.output_dir, taxonomy) == 2
    path = context.output_dir / "taxonomy_evidence.tsv"
    if change == "count":
        path.write_text(path.read_text().replace("\t70\t", "\t71\t"))
    elif change == "fraction":
        path.write_text(path.read_text().replace("\t0.35\t", "\t0.36\t"))
    elif change == "id":
        path.write_text(path.read_text().replace("10\tSame name", "999\tSame name"))
    else:
        path = context.output_dir / "taxonomy_methods.json"
        methods = json.loads(path.read_text())
        methods["sources"]["kraken_report"]["sha256"] = "0" * 64
        path.write_text(json.dumps(methods))
    with pytest.raises(ValueError, match="[Tt]axonomy"):
        verify_taxonomy_evidence(context.output_dir, taxonomy)
