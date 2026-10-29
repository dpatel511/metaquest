"""Research screening: explicit taxon targets and curated AMRFinderPlus evidence.

Taxonomic candidates and protein matches are independent evidence. This module
does not assign a protein to a species or infer clinical risk from co-occurrence.
All genomic coordinates are 1-based inclusive, as in the Pyrodigal GFF3.
"""

from __future__ import annotations

import csv
import json
import subprocess
from collections import Counter
from pathlib import Path
from metaquest.io.processes import run_command

from metaquest.core.abundance import _sha256, read_gene_gff
from metaquest.exceptions import MetaQuestError


class ScreeningError(MetaQuestError):
    """Screening could not be executed or its evidence could not be validated."""


DEFAULT_CATALOG = Path(__file__).parents[1] / "data" / "human_associated_taxa.tsv"
AMR_FILES = ("version.txt", "AMRProt.fa", "AMR.LIB", "fam.tsv")
TAXON_COLUMNS = ("taxid", "name", "estimated_count", "fraction_all_input", "source", "note")
PROTEIN_COLUMNS = (
    "gene_id", "contig_id", "start", "end", "strand", "partial", "element_type",
    "subtype", "symbol", "description", "scope", "class", "subclass", "method",
    "identity_percent", "reference_coverage_percent", "reference_accession",
    "hmm_accession", "count", "mean_depth", "tpm",
)
LIMITATIONS = [
    "Taxonomic candidates are species-level classifications, not confirmation of infection or pathogenic strains.",
    "Protein evidence uses AMRFinderPlus curated AMR and selected virulence/stress families; coverage is organism dependent.",
    "Protein-only mode does not screen nucleotide mutations, expression or phenotype.",
    "Protein hits are linked to genes and contigs, not attributed to the taxonomic candidates by co-occurrence.",
    "No reported hit does not establish absence. Unassembled, unpredicted and database-unrepresented sequences are not assessed.",
    "The bundled bacterial target catalog is a transparent research subset, not a complete human-pathogen catalog.",
]


def resolve_amr_database(path: Path) -> Path:
    path = Path(path).expanduser().resolve()
    candidate = (path / "latest").resolve() if (path / "latest").is_dir() else path
    missing = [name for name in AMR_FILES
               if not (candidate / name).is_file() or (candidate / name).stat().st_size == 0]
    if missing:
        raise ScreeningError(
            f"AMRFinderPlus database is incomplete at {candidate}: {', '.join(missing)}",
            suggestions=["Run 'metaquest databases --database amr' using the same database root."],
        )
    return candidate


def _write_table(path: Path, columns: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def screen_taxa(bracken: Path, catalog: Path, total_input: int) -> tuple[list[dict], int]:
    with Path(catalog).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"taxid", "name", "source", "note"}.issubset(reader.fieldnames or ()):
            raise ScreeningError("Taxon catalog requires taxid, name, source and note columns")
        targets = {}
        for row in reader:
            taxid = row["taxid"]
            if not taxid.isdigit() or taxid in targets or not row["name"] or not row["source"]:
                raise ScreeningError("Taxon catalog has invalid/duplicate taxids or missing provenance")
            targets[taxid] = row
    if not targets:
        raise ScreeningError("Taxon catalog is empty")
    found = []
    with Path(bracken).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"taxonomy_id", "taxonomy_lvl", "new_est_reads"}.issubset(reader.fieldnames or ()):
            raise ScreeningError("Taxon screening requires species-level Bracken taxonomy IDs and counts")
        seen = set()
        assigned = 0
        for row in reader:
            if row["taxonomy_lvl"] != "S":
                raise ScreeningError("Taxon screening requires classification.taxonomic_level: S")
            if row["taxonomy_id"] in seen:
                raise ScreeningError("Duplicate Bracken taxonomy IDs prevent taxon evidence joins")
            seen.add(row["taxonomy_id"])
            target = targets.get(row["taxonomy_id"])
            count = float(row["new_est_reads"])
            if count < 0 or not count.is_integer() or count > total_input:
                raise ScreeningError("Bracken counts are invalid for the reported input denominator")
            assigned += int(count)
            if target and count > 0:
                found.append({
                    "taxid": row["taxonomy_id"], "name": target["name"],
                    "estimated_count": int(count),
                    "fraction_all_input": count / total_input if total_input else 0.0,
                    "source": target["source"], "note": target["note"],
                })
    if assigned > total_input:
        raise ScreeningError("Bracken count total exceeds the reported input denominator")
    return sorted(found, key=lambda row: (-row["estimated_count"], row["taxid"])), len(targets)


def read_protein_evidence(raw_table: Path, gff: Path, abundance: Path | None) -> list[dict]:
    genes = read_gene_gff(gff)
    coordinates = {gene.gene_id: gene for gene in genes}
    if len(coordinates) != len(genes):
        raise ScreeningError("Duplicate GFF gene IDs prevent protein evidence joins")
    quantities = {}
    if abundance:
        with Path(abundance).open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not {"gene_id", "count", "mean_depth", "tpm"}.issubset(reader.fieldnames or ()):
                raise ScreeningError("Gene abundance table lacks required quantification columns")
            for row in reader:
                if row["gene_id"] in quantities:
                    raise ScreeningError("Duplicate abundance gene IDs prevent protein evidence joins")
                quantities[row["gene_id"]] = row
        if set(quantities) != set(coordinates):
            raise ScreeningError("Abundance and predicted-gene IDs do not match")
    evidence = []
    with Path(raw_table).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        symbol = "Element symbol" if "Element symbol" in fields else "Gene symbol"
        description = "Element name" if "Element name" in fields else "Sequence name"
        element_type = "Type" if "Type" in fields else "Element type"
        required = {"Protein id", "Method", symbol, description, element_type, "Scope"}
        if not required.issubset(fields):
            raise ScreeningError(f"AMRFinderPlus output lacks required columns: {sorted(required - fields)}")
        for row in reader:
            gene_id = row["Protein id"]
            if gene_id not in coordinates:
                raise ScreeningError(f"AMRFinderPlus protein ID is not in predicted genes: {gene_id}")
            gene = coordinates[gene_id]
            quantity = quantities.get(gene_id, {})
            evidence.append({
                "gene_id": gene_id, "contig_id": gene.contig_id, "start": gene.start,
                "end": gene.end, "strand": gene.strand, "partial": gene.partial,
                "element_type": row[element_type], "subtype": row.get("Subtype", row.get("Element subtype", "")),
                "symbol": row[symbol], "description": row[description], "scope": row["Scope"],
                "class": row.get("Class", ""), "subclass": row.get("Subclass", ""),
                "method": row["Method"], "identity_percent": row.get("% Identity to reference", ""),
                "reference_coverage_percent": row.get("% Coverage of reference", ""),
                "reference_accession": row.get("Closest reference accession", ""),
                "hmm_accession": row.get("HMM accession", ""),
                **{key: quantity.get(key, "") for key in ("count", "mean_depth", "tpm")},
            })
    return evidence


def run_screening(
    proteins: Path, gff: Path, bracken: Path, output_dir: Path, database: Path,
    *, total_input: int, count_unit: str, catalog: Path | None = None,
    abundance: Path | None = None, threads: int = 4, reuse: bool = False,
) -> dict:
    catalog = Path(catalog or DEFAULT_CATALOG).resolve()
    database = resolve_amr_database(database)
    output = Path(output_dir).resolve() / "screening"
    output.mkdir(parents=True, exist_ok=True)
    version_run = subprocess.run(["amrfinder", "--version"], capture_output=True, text=True, check=True)
    version = version_run.stdout.strip() or version_run.stderr.strip()
    state = {
        "schema": 1, "tool_version": version,
        "inputs": {str(Path(path).resolve()): _sha256(path) for path in (proteins, gff, bracken, catalog)},
        "abundance_sha256": _sha256(abundance) if abundance else None,
        "database_files": {name: _sha256(database / name) for name in AMR_FILES},
        "total_input": total_input, "count_unit": count_unit,
    }
    completion = output / "completion.json"
    summary_path = output / "summary.json"
    raw = output / "amrfinder.tsv"
    taxon_path, protein_path = output / "taxon_candidates.tsv", output / "protein_evidence.tsv"
    if reuse and all(path.is_file() for path in (completion, summary_path, raw, taxon_path, protein_path)):
        try:
            previous = json.loads(completion.read_text())
        except (OSError, ValueError):
            previous = {}
        output_hashes = {path.name: _sha256(path) for path in (summary_path, raw, taxon_path, protein_path)}
        if previous.get("state") == state and previous.get("completed") is True and previous.get("outputs") == output_hashes:
            result = json.loads(summary_path.read_text())
            result["reused"] = True
            return result
    completion.unlink(missing_ok=True)
    candidates, target_count = screen_taxa(bracken, catalog, total_input)
    genes = read_gene_gff(gff)
    from Bio import SeqIO

    protein_ids = set()
    for protein in SeqIO.parse(proteins, "fasta"):
        if protein.id in protein_ids:
            raise ScreeningError("Duplicate protein FASTA IDs prevent predicted-gene joins")
        protein_ids.add(protein.id)
    if len({gene.gene_id for gene in genes}) != len(genes) or protein_ids != {gene.gene_id for gene in genes}:
        raise ScreeningError("Protein FASTA and predicted GFF gene IDs do not match")
    if genes:
        command = ["amrfinder", "--protein", str(Path(proteins).resolve()), "--plus",
                   "--database", str(database), "--threads", str(threads), "--output", str(raw)]
        with (output / "amrfinder.log").open("w", encoding="utf-8") as log:
            log.write("[MetaQuest] " + " ".join(command) + "\n")
            log.flush()
            process = run_command(command, stdout=log, stderr=log)
        if process.returncode:
            raise ScreeningError(f"AMRFinderPlus failed (exit {process.returncode}); see {output / 'amrfinder.log'}")
    else:
        raw.write_text("Protein id\tElement symbol\tElement name\tType\tScope\tMethod\n")
    evidence = read_protein_evidence(raw, gff, abundance)
    _write_table(taxon_path, TAXON_COLUMNS, candidates)
    _write_table(protein_path, PROTEIN_COLUMNS, evidence)
    summary = {
        "status": "completed", "research_use_only": True, "tool": "AMRFinderPlus",
        "tool_version": version, "database_release": (database / "version.txt").read_text().strip(),
        "database_sha256": state["database_files"], "catalog_sha256": _sha256(catalog),
        "catalog_path": str(catalog), "catalog_targets": target_count,
        "taxon_candidates": len(candidates), "protein_matches": len(evidence),
        "unique_matched_genes": len({row["gene_id"] for row in evidence}),
        "protein_types": dict(Counter(row["element_type"] for row in evidence)),
        "screened_genes": len(genes), "count_unit": count_unit, "total_input": total_input,
        "coordinates": "1-based inclusive", "protein_mode": "protein-only --plus; no organism-specific mutation model",
        "limitations": LIMITATIONS, "reused": False,
        "taxon_table": "screening/taxon_candidates.tsv", "protein_table": "screening/protein_evidence.tsv",
        "raw_table": "screening/amrfinder.tsv",
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    output_hashes = {path.name: _sha256(path) for path in (summary_path, raw, taxon_path, protein_path)}
    completion.write_text(json.dumps({"state": state, "completed": True, "outputs": output_hashes}, indent=2) + "\n")
    return summary
