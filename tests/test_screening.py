import csv
import json
import os
import shutil
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest

from metaquest.core import screening
from metaquest.core.screening import ScreeningError, read_protein_evidence, run_screening, screen_taxa
from metaquest.pipeline.runner import build_default_pipeline
from metaquest.settings import load_config


def fixture_files(root):
    proteins = root / "genes.faa"
    proteins.write_text(">g1\nMAAA\n>g2\nMCCC\n")
    gff = root / "genes.gff3"
    gff.write_text("c1\tfixture\tCDS\t1\t12\t.\t+\t0\tID=g1;partial=00\n"
                   "c2\tfixture\tCDS\t20\t31\t.\t-\t0\tID=g2;partial=10\n")
    bracken = root / "bracken.tsv"
    bracken.write_text("name\ttaxonomy_id\ttaxonomy_lvl\tnew_est_reads\nEscherichia coli\t562\tS\t8\nOther\t999\tS\t2\n")
    database = root / "db"
    database.mkdir()
    for name in screening.AMR_FILES:
        (database / name).write_text("fixture-v1\n")
    return proteins, gff, bracken, database


def raw_report(path, gene="g1"):
    path.write_text(
        "Protein id\tElement symbol\tElement name\tType\tSubtype\tScope\tMethod\t% Identity to reference\t% Coverage of reference\n"
        f"{gene}\tsymbol\tCurated family\tAMR\tAMR\tcore\tHMM\tNA\tNA\n"
    )


def test_taxon_screen_uses_exact_ids_not_name_or_genus(tmp_path):
    _, _, bracken, _ = fixture_files(tmp_path)
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text("taxid\tname\tsource\tnote\n562\tEscherichia coli\thttps://example.test\tCandidate only\n")
    hits, size = screen_taxa(bracken, catalog, 10)
    assert size == 1
    assert hits[0]["estimated_count"] == 8
    assert hits[0]["fraction_all_input"] == .8
    assert len(hits) == 1
    bracken.write_text("name\ttaxonomy_id\ttaxonomy_lvl\tnew_est_reads\nEscherichia coli\t999\tS\t8\n")
    assert screen_taxa(bracken, catalog, 10)[0] == []


def test_catalog_duplicates_and_genus_level_are_rejected(tmp_path):
    _, _, bracken, _ = fixture_files(tmp_path)
    catalog = tmp_path / "catalog.tsv"
    catalog.write_text("taxid\tname\tsource\tnote\n562\tE coli\tref\tNote\n562\tDuplicate\tref\tNote\n")
    with pytest.raises(ScreeningError, match="duplicate"):
        screen_taxa(bracken, catalog, 10)
    bracken.write_text(bracken.read_text().replace("\tS\t", "\tG\t"))
    with pytest.raises(ScreeningError, match="taxonomic_level"):
        screen_taxa(bracken, screening.DEFAULT_CATALOG, 10)


def test_protein_evidence_preserves_hmm_na_and_validates_gene_join(tmp_path):
    _, gff, _, _ = fixture_files(tmp_path)
    raw = tmp_path / "raw.tsv"
    raw_report(raw)
    hits = read_protein_evidence(raw, gff, None)
    assert hits[0]["start"] == 1 and hits[0]["end"] == 12
    assert hits[0]["identity_percent"] == "NA"
    assert hits[0]["count"] == ""
    assert hits[0]["method"] == "HMM"
    raw_report(raw, "unknown")
    with pytest.raises(ScreeningError, match="not in predicted genes"):
        read_protein_evidence(raw, gff, None)


def test_screening_resume_requires_content_and_database_match(tmp_path, monkeypatch):
    proteins, gff, bracken, database = fixture_files(tmp_path)
    calls = []

    def run(command, **kwargs):
        if "--version" in command:
            return SimpleNamespace(stdout="4.2.7", stderr="", returncode=0)
        calls.append(command)
        raw_report(Path(command[command.index("--output") + 1]))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(screening.subprocess, "run", run)
    monkeypatch.setattr(screening, "run_command", run)
    first = run_screening(proteins, gff, bracken, tmp_path, database, total_input=10, count_unit="fragments")
    assert first["protein_types"] == {"AMR": 1}
    assert first["screened_genes"] == 2
    assert not first["reused"]
    again = run_screening(proteins, gff, bracken, tmp_path, database, total_input=10, count_unit="fragments", reuse=True)
    assert again["reused"] and len(calls) == 1
    (tmp_path / "screening/protein_evidence.tsv").write_text("corrupted")
    run_screening(proteins, gff, bracken, tmp_path, database, total_input=10, count_unit="fragments", reuse=True)
    assert len(calls) == 2
    (database / "AMRProt.fa").write_text("changed-v1\n")
    run_screening(proteins, gff, bracken, tmp_path, database, total_input=10, count_unit="fragments", reuse=True)
    assert len(calls) == 3
    assert json.loads((tmp_path / "screening/completion.json").read_text())["completed"]
    (tmp_path / "screening/completion.json").write_text("incomplete marker")
    run_screening(proteins, gff, bracken, tmp_path, database, total_input=10, count_unit="fragments", reuse=True)
    assert len(calls) == 4


def test_screening_is_explicit_and_requires_gene_prediction():
    config = load_config()
    config = replace(config, screening=replace(config.screening, enabled=True))
    names = [name for name, _ in build_default_pipeline(config)._stages]
    assert names[-2:] == ["Research Screening", "Reporting"]
    with pytest.raises(Exception, match="requires gene prediction"):
        build_default_pipeline(config, skip_annotation=True)


def test_protein_and_abundance_id_mismatches_are_not_hidden(tmp_path, monkeypatch):
    proteins, gff, bracken, database = fixture_files(tmp_path)
    monkeypatch.setattr(screening.subprocess, "run", lambda *a, **k: SimpleNamespace(stdout="4.2.7", stderr=""))
    proteins.write_text(">different\nMAAA\n")
    with pytest.raises(ScreeningError, match="FASTA and predicted GFF"):
        run_screening(proteins, gff, bracken, tmp_path, database, total_input=10, count_unit="reads")
    assert not (tmp_path / "screening/completion.json").exists()
    raw = tmp_path / "raw.tsv"
    raw_report(raw)
    abundance = tmp_path / "abundance.tsv"
    abundance.write_text("gene_id\tcount\tmean_depth\ttpm\ng1\t1\t2\t3\ng1\t2\t2\t3\ng2\t3\t2\t3\n")
    with pytest.raises(ScreeningError, match="Duplicate abundance"):
        read_protein_evidence(raw, gff, abundance)


def test_taxonomy_duplicate_ids_and_impossible_totals_fail(tmp_path):
    _, _, bracken, _ = fixture_files(tmp_path)
    bracken.write_text(bracken.read_text() + "Duplicate\t562\tS\t1\n")
    with pytest.raises(ScreeningError, match="Duplicate Bracken"):
        screen_taxa(bracken, screening.DEFAULT_CATALOG, 10)
    bracken.write_text("name\ttaxonomy_id\ttaxonomy_lvl\tnew_est_reads\nA\t562\tS\t8\nB\t999\tS\t8\n")
    with pytest.raises(ScreeningError, match="total exceeds"):
        screen_taxa(bracken, screening.DEFAULT_CATALOG, 10)


def test_empty_database_core_file_is_not_installed(tmp_path):
    _, _, _, database = fixture_files(tmp_path)
    (database / "AMR.LIB").write_text("")
    with pytest.raises(ScreeningError, match="incomplete"):
        screening.resolve_amr_database(database)


@pytest.mark.skipif(not os.environ.get("METAQUEST_TEST_AMR_DB") or not shutil.which("amrfinder"),
                    reason="set METAQUEST_TEST_AMR_DB to run the real AMRFinderPlus test")
def test_native_amrfinder_positive_and_negative_controls(tmp_path):
    from Bio import SeqIO

    database = screening.resolve_amr_database(Path(os.environ["METAQUEST_TEST_AMR_DB"]))
    reference = next(SeqIO.parse(database / "AMRProt.fa", "fasta"))
    proteins, gff, bracken, _ = fixture_files(tmp_path)
    proteins.write_text(f">g1\n{reference.seq}\n>g2\nM{'A' * 60}\n")
    gff.write_text(f"c1\tfixture\tCDS\t1\t{len(reference.seq) * 3}\t.\t+\t0\tID=g1;partial=00\n"
                   "c2\tfixture\tCDS\t1\t183\t.\t+\t0\tID=g2;partial=00\n")
    result = run_screening(proteins, gff, bracken, tmp_path, database,
                           total_input=10, count_unit="fragments", threads=2)
    with (tmp_path / result["protein_table"]).open() as handle:
        hits = list(csv.DictReader(handle, delimiter="\t"))
    assert hits and {row["gene_id"] for row in hits} == {"g1"}
    assert any("EXACT" in row["method"] or "ALLELE" in row["method"] for row in hits)


@pytest.mark.skipif(not os.environ.get("METAQUEST_TEST_AMR_DB") or not shutil.which("amrfinder"),
                    reason="set METAQUEST_TEST_AMR_DB to run the real AMRFinderPlus test")
def test_native_zero_hits_produce_completed_reports_with_empty_tables(tmp_path):
    from metaquest.pipeline.context import PipelineContext
    from metaquest.reporting.stable_reporter import generate_stable_reports

    proteins, gff, bracken, _ = fixture_files(tmp_path)
    proteins.write_text(f">g1\nM{'E' * 60}\n>g2\nM{'F' * 60}\n")
    gff.write_text("c1\tfixture\tCDS\t1\t183\t.\t+\t0\tID=g1;partial=00\n"
                   "c2\tfixture\tCDS\t20\t202\t.\t-\t0\tID=g2;partial=00\n")
    bracken.write_text("name\ttaxonomy_id\ttaxonomy_lvl\tnew_est_reads\nUncatalogued fixture\t999\tS\t10\n")
    result = run_screening(proteins, gff, bracken, tmp_path, Path(os.environ["METAQUEST_TEST_AMR_DB"]),
                           total_input=10, count_unit="reads", threads=2)
    assert result["taxon_candidates"] == result["protein_matches"] == 0
    for key in ("taxon_table", "protein_table"):
        with (tmp_path / result[key]).open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            assert reader.fieldnames and list(reader) == []
    # Only screening is exercised here; prediction/assembly must remain not run.
    context = PipelineContext(load_config(), [], tmp_path)
    context.screening = result
    generate_stable_reports(context)
    stage = json.loads((tmp_path / "reports/07_screening.json").read_text())
    assert stage["status"] == "completed"
    assert stage["metrics"]["protein_matches"] == 0
    assert json.loads((tmp_path / "reports/04_gene_prediction.json").read_text())["status"] == "not_run"
    assert (tmp_path / "plots/screening_proteins.png").is_file()
