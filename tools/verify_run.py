#!/usr/bin/env python3
"""Check saved run accounting and offline report artifacts; no accuracy claims.

Gene coordinates use the pipeline's 1-based inclusive GFF3 convention.
Run this with the installed MetaQuest environment, including Biopython.
"""

import argparse
import csv
import json
import math
from html.parser import HTMLParser
from pathlib import Path

from Bio import SeqIO

from metaquest.core.abundance import read_gene_gff


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rows(path):
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


class ReportLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.images = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "a" and "href" in attrs:
            self.links.append(attrs["href"])
        if tag == "img":
            self.images.append(attrs.get("src", ""))


def verify(root):
    root = root.resolve()
    metadata = json.loads((root / "analysis_metadata.json").read_text())
    summary = json.loads((root / "analysis_summary.json").read_text())
    require(metadata.get("status") == "completed", "Run did not complete")
    taxonomy = summary["taxonomy"]
    total = taxonomy["total_input"]
    require(taxonomy["kraken_classified"] + taxonomy["unclassified_reads"] == total,
            "Taxonomy classified/unclassified counts do not reconcile")
    bracken = list(rows(root / "bracken_report.tsv"))
    require(len({item["taxonomy_id"] for item in bracken}) == len(bracken), "Duplicate taxonomy IDs")
    require(sum(int(item["new_est_reads"]) for item in bracken) == taxonomy["species_assigned"],
            "Bracken counts disagree with summary")
    require(0 <= taxonomy["species_assigned"] <= taxonomy["kraken_classified"] <= total,
            "Taxonomy counts exceed their input denominator")
    genes = set()
    if "gene_prediction" in summary:
        proteins = list(record.id for record in SeqIO.parse(root / "gene_prediction/genes.faa", "fasta"))
        genes = set(proteins)
        require(len(genes) == len(proteins), "Duplicate protein IDs")
        predicted = summary["annotation"]["predicted_genes"]
        require(len(genes) == predicted, "Protein count disagrees with predicted gene count")
        gff = read_gene_gff(root / "gene_prediction/genes.gff3")
        require(len(gff) == len(genes) and {gene.gene_id for gene in gff} == genes,
                "Protein and GFF gene IDs do not match")
    annotation = summary.get("annotation", {})
    if annotation.get("functional_status") == "completed":
        annotations = list(rows(root / annotation["functional_annotations"]))
        require(len(annotations) == len(genes) and {row["query"] for row in annotations} == genes,
                "Annotation and predicted gene IDs do not match")
        require(sum(row["annotation_status"] == "annotated" for row in annotations) == annotation["annotated_genes"],
                "Annotated gene count disagrees with summary")
    abundance = summary.get("abundance")
    if abundance:
        values = list(rows(root / abundance["gene_abundance"]))
        require(len(values) == len(genes) and {row["gene_id"] for row in values} == genes,
                "Abundance and predicted gene IDs do not match")
        for row in values:
            require(all(math.isfinite(float(row[key])) and float(row[key]) >= 0
                        for key in ("count", "mean_depth", "tpm")), "Invalid abundance values")
        unit = abundance["count_unit"]
        assigned = abundance[f"{unit}_in_genes"]
        require(0 <= assigned <= abundance[f"mapped_{unit}"] <= abundance[f"input_{unit}"],
                "Mapping and assignment counts exceed their denominator")
        require(math.isclose(sum(float(row["count"]) for row in values), assigned,
                             rel_tol=0, abs_tol=max(.001, len(values) * .00005)),
                "Gene count sum disagrees with assigned reads/fragments")
        require(math.isclose(sum(float(row["tpm"]) for row in values), 1e6 if assigned else 0,
                             rel_tol=0, abs_tol=max(.001, len(values) * .0000005)),
                "Gene TPM does not reconcile with its normalization denominator")
    screening = summary.get("screening")
    if screening:
        evidence = list(rows(root / screening["protein_table"]))
        require(len(evidence) == screening["protein_matches"], "Protein evidence row count disagrees")
        require({row["gene_id"] for row in evidence}.issubset(genes), "Unknown screening gene IDs")
        candidates = list(rows(root / screening["taxon_table"]))
        require(len(candidates) == screening["taxon_candidates"], "Taxon candidate count disagrees")
    reports = summary["stage_reports"]
    require(len(reports) == 8, "Expected eight stage status reports")
    html_paths = [root / "report.html"]
    for report in reports:
        for key in ("html", "json", "table"):
            require((root / report[key]).is_file(), f"Missing stage artifact: {report[key]}")
        payload = json.loads((root / report["json"]).read_text())
        require(payload["status"] == report["status"], "Stage status differs between report and summary")
        html_paths.append(root / report["html"])
    for path in html_paths:
        parser = ReportLinks()
        parser.feed(path.read_text(encoding="utf-8"))
        require(all(src.startswith("data:image/") for src in parser.images), "Report has nonembedded images")
        for href in parser.links:
            if href.startswith(("#", "http://", "https://")):
                continue
            require((path.parent / href).is_file(), f"Broken report link: {href}")
    for path in summary.get("figures", []):
        require((root / path).is_file(), f"Missing figure: {path}")
    return {"status": "passed", "metaquest_version": metadata["metaquest_version"],
            "input_count": total, "count_unit": taxonomy["count_unit"],
            "predicted_genes": len(genes), "annotated_genes": annotation.get("annotated_genes"),
            "protein_evidence_rows": screening["protein_matches"] if screening else None,
            "stage_reports": len(reports), "figure_files": len(summary.get("figures", [])),
            "interpretation": "Accounting and artifact checks; not a biological accuracy benchmark."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="Completed MetaQuest output directory")
    args = parser.parse_args()
    try:
        print(json.dumps(verify(args.results), indent=2))
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f"Acceptance check failed: {exc}\n")
