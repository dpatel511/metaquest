"""Taxon-ID joins between observed Kraken assignments and Bracken estimates."""
import csv
from dataclasses import asdict
import hashlib
import html
import json
from pathlib import Path
import re

import pandas as pd

from metaquest.config import __version__


EVIDENCE_COLUMNS = (
    "taxonomy_id", "name", "rank", "count_unit", "kraken_direct_count",
    "kraken_clade_count", "bracken_assigned_count", "bracken_added_count",
    "bracken_estimated_count", "fraction_all_input", "fraction_bracken_profile",
)
RANK_NAMES = {"S": "species", "G": "genus", "F": "family"}


def _count(value):
    if not re.fullmatch(r"[0-9]+", str(value)):
        raise ValueError(f"Count must be a nonnegative integer: {value!r}")
    return int(value)


def _kraken_rows(path):
    rows = {}
    if path is not None and Path(path).is_file():
        # Kraken's documented six-column TSV has no library-specific parser.
        with Path(path).open(encoding="utf-8", newline="") as handle:
            for fields in csv.reader(handle, delimiter="\t"):
                if not fields:
                    continue
                if len(fields) != 6 or not re.fullmatch(r"[0-9]+", fields[4]):
                    raise ValueError("Malformed standard Kraken report")
                taxid = fields[4]
                if taxid in rows:
                    raise ValueError(f"Duplicate Kraken taxonomy ID: {taxid}")
                clade, direct = _count(fields[1]), _count(fields[2])
                if direct > clade:
                    raise ValueError(f"Kraken direct count exceeds clade count: {taxid}")
                rows[taxid] = (clade, direct)
        if not rows:
            raise ValueError("Empty Kraken report")
    return rows


def build_evidence(table, classification, unit, level):
    """Validate accounting and preserve missing support as null, never zero."""
    total = _count(classification.total_reads)
    classified = _count(classification.kraken_classified_reads)
    unclassified = _count(classification.unclassified_reads)
    assigned = _count(classification.species_assigned_reads)
    if classified + unclassified != total or not assigned <= classified <= total:
        raise ValueError("Taxonomy counts do not reconcile with their input denominator")
    kraken_path = getattr(classification, "kraken_report", None)
    kraken = _kraken_rows(kraken_path)
    if kraken and (kraken.get("1", (0, 0))[0] != classified
                   or kraken.get("0", (0, 0))[0] != unclassified):
        raise ValueError("Kraken source counts disagree with classification counts")
    if any(clade > classified for taxid, (clade, _) in kraken.items() if taxid != "0"):
        raise ValueError("Kraken clade count exceeds classified count")
    rows, seen = [], set()
    for _, item in table.iterrows():
        taxid = item.get("taxonomy_id")
        if taxid is not None:
            if not re.fullmatch(r"[1-9][0-9]*", str(taxid)):
                raise ValueError(f"Invalid Bracken taxonomy ID: {taxid!r}")
            if taxid in seen:
                raise ValueError(f"Duplicate Bracken taxonomy ID: {taxid}")
            seen.add(taxid)
        rank = item.get("taxonomy_lvl", level)
        if rank != level:
            raise ValueError("Bracken rank differs from recorded classification settings")
        estimate = _count(item["new_est_reads"])
        baseline = _count(item["kraken_assigned_reads"]) if "kraken_assigned_reads" in table else None
        added = _count(item["added_reads"]) if "added_reads" in table else None
        if baseline is not None and added is not None and baseline + added != estimate:
            raise ValueError(f"Bracken support counts do not reconcile: {taxid}")
        clade, direct = kraken.get(taxid, (None, None))
        rows.append(dict(zip(EVIDENCE_COLUMNS, (
            taxid, str(item["name"]), rank, unit, direct, clade, baseline, added, estimate,
            estimate / total if total else None, estimate / assigned if assigned else None,
        ))))
    if sum(row["bracken_estimated_count"] for row in rows) != assigned:
        raise ValueError("Bracken estimated counts disagree with classification counts")
    return pd.DataFrame(rows, columns=EVIDENCE_COLUMNS, dtype=object)


def _source(path):
    if path is None or not Path(path).is_file():
        return None
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path), "sha256": digest.hexdigest()}


def write_methods(ctx, evidence, unit, level):
    """Use recorded provenance; never substitute the currently installed tools."""
    path = Path(ctx.output_dir) / "analysis_metadata.json"
    metadata = json.loads(path.read_text()) if path.is_file() else {}
    provenance = metadata.get("tool_provenance", {})
    config = getattr(getattr(ctx, "config", None), "classification", None)
    settings = metadata.get("effective_config", {}).get("classification")
    methods = {
        "schema_version": 1, "metaquest_version": __version__, "count_unit": unit, "taxonomic_level": level,
        "effective_settings": settings if settings is not None else asdict(config) if config else {},
        "observed_read_length_bp": ctx.classification.observed_read_length,
        "bracken_model_read_length_bp": ctx.classification.bracken_read_length,
        "tool_versions": {name: provenance.get(name) for name in ("kraken2", "bracken")},
        "database": provenance.get("kraken_database"),
        "run_provenance": "analysis_metadata.json" if path.is_file() else None,
        "sources": {"kraken_report": _source(getattr(ctx.classification, "kraken_report", None)),
                    "bracken_report": _source(ctx.classification.bracken_file)},
        "taxa_without_kraken_row": int(evidence["kraken_clade_count"].isna().sum()),
        "denominators": {"fraction_all_input": ctx.classification.total_reads,
                         "fraction_bracken_profile": ctx.classification.species_assigned_reads},
        "interpretation": [
            "Kraken direct counts exclude descendants; clade counts include descendants.",
            "Bracken assigned counts are its starting clade counts, not Kraken direct counts.",
            "Bracken added counts are redistributed estimates, not independent observations.",
            "Classified counts outside the Bracken profile can include filtered taxa and higher ranks.",
            "Fractions are not confidence scores; missing support is unavailable, not zero.",
            "Database scope, related references and settings affect calls; no call does not establish absence.",
        ],
    }
    (Path(ctx.output_dir) / "taxonomy_methods.json").write_text(json.dumps(methods, indent=2) + "\n")
    return methods


def evidence_html(output, taxonomy, prefix=""):
    """Show a bounded evidence table; full data and methods remain downloadable."""
    table = pd.read_csv(Path(output) / taxonomy["evidence_table"], sep="\t", dtype=str, keep_default_na=False)
    columns = {"taxonomy_id": "Taxonomy ID", "name": "Taxon", "kraken_direct_count": "Kraken direct",
               "kraken_clade_count": "Kraken clade", "bracken_assigned_count": "Bracken assigned",
               "bracken_added_count": "Bracken added", "bracken_estimated_count": "Bracken estimate",
               "fraction_all_input": "% all input", "fraction_bracken_profile": "% Bracken profile"}
    shown = table.head(50)[list(columns)].copy()
    for name in ("fraction_all_input", "fraction_bracken_profile"):
        shown[name] = shown[name].map(lambda value: f"{float(value) * 100:.4f}" if value else "")
    rendered = shown.rename(columns=columns).replace("", "not available").to_html(
        index=False, escape=True, classes="taxonomy-evidence")
    methods = taxonomy["methods"]
    settings = methods["effective_settings"]
    database = methods["database"] or {}
    rows = {"Report software version": methods["metaquest_version"],
            "Kraken2 version": methods["tool_versions"]["kraken2"],
            "Bracken version": methods["tool_versions"]["bracken"],
            "Rank": methods["taxonomic_level"], "Minimum hit groups": settings.get("min_hit_groups"),
            "Bracken threshold": settings.get("bracken_threshold"),
            "Observed read length (bp)": methods["observed_read_length_bp"],
            "Bracken model length (bp)": methods["bracken_model_read_length_bp"],
            "Database release": database.get("manifest", {}).get("release") or database.get("configured_release")}
    details = "".join(f"<tr><th>{html.escape(key)}</th><td>{html.escape(str(value)) if value is not None else 'not available'}</td></tr>"
                      for key, value in rows.items())
    notes = "".join(f"<li>{html.escape(note)}</li>" for note in methods["interpretation"])
    return (f"<p>Counts are {html.escape(taxonomy['count_unit'])}. Showing up to 50 taxa. "
            "Wide tables scroll horizontally. "
            f"<a href='{prefix}taxonomy_evidence.tsv'>Full taxon evidence TSV</a> · "
            f"<a href='{prefix}taxonomy_methods.json'>Methods and source hashes JSON</a></p>"
            f"<div class='scroll'>{rendered}</div><h3>Methods and interpretation</h3>"
            f"<table>{details}</table><ul>{notes}</ul>")
