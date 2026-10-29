"""Metagenomic gene prediction and eggNOG functional annotation."""

import csv
import hashlib
import json
import os
import subprocess
import tempfile
import threading
import time
from collections import Counter
from pathlib import Path
from metaquest.io.processes import run_command

from Bio import SeqIO

from ..exceptions import AnnotationError


EGGNOG_OUTPUT_COLUMNS = (
    "query",
    "annotation_status",
    "seed_ortholog",
    "evalue",
    "score",
    "max_annot_lvl",
    "cog_categories",
    "description",
    "preferred_name",
    "gos",
    "ecs",
    "kos",
    "pfams",
)


def run_gene_prediction(
    fasta_path: Path,
    output_dir: Path,
    *,
    min_contig_length: int = 200,
) -> tuple[Path, int]:
    """Predict genes from metagenomic contigs with bounded-memory Pyrodigal."""
    try:
        import pyrodigal
    except ImportError as exc:
        raise AnnotationError(
            "Pyrodigal is required for gene prediction; install pyrodigal>=3.7.1"
        ) from exc

    prediction_dir = Path(output_dir) / "gene_prediction"
    prediction_dir.mkdir(parents=True, exist_ok=True)
    proteins_path = prediction_dir / "genes.faa"
    genes_path = prediction_dir / "genes.fna"
    gff_path = prediction_dir / "genes.gff3"
    contig_map_path = prediction_dir / "contig_id_map.tsv"
    stable_contigs_path = prediction_dir / "contigs.stable.fasta"

    finder = pyrodigal.GeneFinder(meta=True)
    contigs_seen = 0
    contigs_processed = 0
    genes_predicted = 0
    contig_ids: Counter[str] = Counter()

    try:
        with proteins_path.open("w", encoding="utf-8") as proteins_out, \
             genes_path.open("w", encoding="utf-8") as genes_out, \
             gff_path.open("w", encoding="utf-8") as gff_out, \
             contig_map_path.open("w", encoding="utf-8") as map_out, \
             stable_contigs_path.open("w", encoding="utf-8") as contigs_out:
            gff_out.write("##gff-version 3\n")
            map_out.write("stable_contig_id\toriginal_contig_id\tsequence_sha256\n")
            for record in SeqIO.parse(Path(fasta_path), "fasta"):
                contigs_seen += 1
                if len(record.seq) < min_contig_length:
                    continue
                digest = hashlib.sha256(bytes(record.seq)).hexdigest()[:16]
                base_id = f"contig_{digest}"
                contig_ids[base_id] += 1
                sequence_id = (
                    base_id
                    if contig_ids[base_id] == 1
                    else f"{base_id}_{contig_ids[base_id]}"
                )
                map_out.write(f"{sequence_id}\t{record.id}\t{hashlib.sha256(bytes(record.seq)).hexdigest()}\n")
                contigs_out.write(f">{sequence_id}\n{record.seq}\n")
                predictions = finder.find_genes(bytes(record.seq))
                contigs_processed += 1
                genes_predicted += len(predictions)
                predictions.write_translations(
                    proteins_out,
                    sequence_id=sequence_id,
                    include_stop=False,
                    full_id=True,
                )
                predictions.write_genes(
                    genes_out,
                    sequence_id=sequence_id,
                    full_id=True,
                )
                predictions.write_gff(
                    gff_out,
                    sequence_id=sequence_id,
                    header=False,
                    include_translation_table=True,
                    full_id=True,
                )
    except Exception as exc:
        raise AnnotationError(f"Pyrodigal gene prediction failed: {exc}") from exc

    summary = {
        "tool": "Pyrodigal",
        "tool_version": getattr(pyrodigal, "__version__", "unknown"),
        "mode": "metagenomic",
        "minimum_contig_length": min_contig_length,
        "contigs_seen": contigs_seen,
        "contigs_processed": contigs_processed,
        "genes_predicted": genes_predicted,
    }
    (prediction_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    return prediction_dir, genes_predicted


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_fasta_sha256(path: Path) -> str:
    """Hash FASTA records independently of record order."""
    digest = hashlib.sha256()
    records = sorted((record.id, str(record.seq)) for record in SeqIO.parse(path, "fasta"))
    for identifier, sequence in records:
        digest.update(identifier.encode())
        digest.update(b"\0")
        digest.update(sequence.encode())
        digest.update(b"\0")
    return digest.hexdigest()


def _emapper_version(database_dir: Path) -> str:
    try:
        result = subprocess.run(
            ["emapper.py", "--data_dir", str(database_dir), "--version"],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError as exc:
        raise AnnotationError(
            "eggNOG-mapper is required for functional annotation; "
            "install eggnog-mapper=2.1.15"
        ) from exc
    output = f"{result.stdout}\n{result.stderr}".strip()
    if result.returncode:
        raise AnnotationError(f"Could not determine eggNOG-mapper version: {output}")
    return output


def _read_emapper_annotations(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    header: list[str] | None = None
    annotations: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("#query"):
                header = line.removeprefix("#").split("\t")
                continue
            if not line or line.startswith("#"):
                continue
            if header is None:
                raise AnnotationError("eggNOG annotations are missing the #query header")
            values = line.split("\t")
            values.extend([""] * (len(header) - len(values)))
            row = dict(zip(header, values))
            annotations[row["query"]] = row
    if header is None:
        raise AnnotationError("eggNOG annotations are missing the #query header")
    return header, annotations


def _terms(value: str) -> list[str]:
    if not value or value == "-":
        return []
    return [term.strip() for term in value.split(",") if term.strip() and term != "-"]


def _write_functional_outputs(
    proteins_path: Path,
    raw_annotations: Path,
    output_dir: Path,
) -> tuple[Path, Path, int, int]:
    """Left-join eggNOG hits to every predicted protein and aggregate terms."""
    _, annotations = _read_emapper_annotations(raw_annotations)
    protein_ids = [record.id for record in SeqIO.parse(proteins_path, "fasta")]
    table_path = output_dir / "functional_annotations.tsv"
    counts: dict[str, Counter[str]] = {
        "COG": Counter(),
        "GO": Counter(),
        "EC": Counter(),
        "KO": Counter(),
    }
    annotated_count = 0

    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EGGNOG_OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for protein_id in protein_ids:
            source = annotations.get(protein_id, {})
            annotated = bool(source)
            annotated_count += int(annotated)
            cog_value = source.get("COG_category", "")
            go_value = source.get("GOs", "")
            ec_value = source.get("EC", "")
            ko_value = source.get("KEGG_ko", "")
            writer.writerow(
                {
                    "query": protein_id,
                    "annotation_status": "annotated" if annotated else "unannotated",
                    "seed_ortholog": source.get("seed_ortholog", ""),
                    "evalue": source.get("evalue", ""),
                    "score": source.get("score", ""),
                    "max_annot_lvl": source.get("max_annot_lvl", ""),
                    "cog_categories": cog_value,
                    "description": source.get("Description", ""),
                    "preferred_name": source.get("Preferred_name", ""),
                    "gos": go_value,
                    "ecs": ec_value,
                    "kos": ko_value,
                    "pfams": source.get("PFAMs", ""),
                }
            )
            for category in "".join(_terms(cog_value)):
                if category.isalpha():
                    counts["COG"][category] += 1
            for namespace, value in (("GO", go_value), ("EC", ec_value), ("KO", ko_value)):
                counts[namespace].update(_terms(value))

    category_path = output_dir / "functional_category_summary.tsv"
    with category_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("namespace", "term", "gene_count"))
        for namespace in ("COG", "KO", "EC", "GO"):
            for term, count in sorted(counts[namespace].items()):
                writer.writerow((namespace, term, count))
    return table_path, category_path, annotated_count, len(protein_ids)


def run_functional_annotation(
    proteins_path: Path,
    output_dir: Path,
    database_dir: Path,
    *,
    threads: int = 8,
    evalue: float = 1e-6,
    diamond_block_size: float = 0.5,
    tax_scope: str = "auto",
    expected_version: str = "2.1.15",
    database_release: str = "5.0.2",
    reuse: bool = True,
) -> tuple[Path, Path, int, bool]:
    """Annotate every predicted protein with pinned eggNOG-mapper v2 data."""
    proteins_path = Path(proteins_path).resolve()
    output_dir = (Path(output_dir) / "functional_annotation").resolve()
    database_dir = Path(database_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if not proteins_path.is_file():
        raise AnnotationError(f"Predicted protein file not found: {proteins_path}")
    raw_annotations = output_dir / "metaquest.emapper.annotations"
    summary_path = output_dir / "summary.json"
    completion_path = output_dir / "completion.json"
    protein_sha256 = _canonical_fasta_sha256(proteins_path)
    version_output = _emapper_version(database_dir)
    if expected_version not in version_output:
        raise AnnotationError(
            f"Unsupported eggNOG-mapper version; expected {expected_version}, got {version_output}"
        )
    mapper_build = version_output.split(" / ", 1)[0].removeprefix("emapper-")

    expected_state = {
        "protein_content_sha256": protein_sha256,
        "eggnog_mapper_version": expected_version,
        "eggnog_mapper_build": mapper_build,
        "eggnog_database_release": database_release,
        "tax_scope": tax_scope,
        "evalue": evalue,
        "diamond_block_size": diamond_block_size,
    }
    if reuse and completion_path.is_file() and raw_annotations.is_file() and summary_path.is_file():
        existing = json.loads(completion_path.read_text(encoding="utf-8"))
        if all(existing.get(key) == value for key, value in expected_state.items()):
            table, categories, annotated, total = _write_functional_outputs(
                proteins_path, raw_annotations, output_dir
            )
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            summary.update(
                {
                    "total_genes": total,
                    "annotated_genes": annotated,
                    "unannotated_genes": total - annotated,
                    "reused": True,
                }
            )
            summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
            return table, categories, annotated, True

    completion_path.unlink(missing_ok=True)
    required_database_files = ("eggnog.db", "eggnog.taxa.db", "eggnog_proteins.dmnd")
    missing = [name for name in required_database_files if not (database_dir / name).is_file()]
    if missing:
        raise AnnotationError(
            f"eggNOG database is incomplete at {database_dir}: missing {', '.join(missing)}"
        )

    if not any(SeqIO.parse(proteins_path, "fasta")):
        raw_annotations.write_text(
            "#query\tseed_ortholog\tevalue\tscore\tmax_annot_lvl\tCOG_category\t"
            "Description\tPreferred_name\tGOs\tEC\tKEGG_ko\tPFAMs\n",
            encoding="utf-8",
        )
        table_path, category_path, annotated_count, total_count = _write_functional_outputs(
            proteins_path, raw_annotations, output_dir
        )
        summary_path.write_text(
            json.dumps(
                {
                    "tool": "eggNOG-mapper",
                    "tool_version": mapper_build,
                    "database_release": database_release,
                    "method": "diamond",
                    "tax_scope": tax_scope,
                    "evalue": evalue,
                    "diamond_block_size": diamond_block_size,
                    "threads": threads,
                    "total_genes": total_count,
                    "annotated_genes": annotated_count,
                    "unannotated_genes": 0,
                    "raw_annotations": raw_annotations.name,
                    "reused": False,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        completion_path.write_text(
            json.dumps({**expected_state, "completed": True}, indent=2) + "\n",
            encoding="utf-8",
        )
        return table_path, category_path, annotated_count, False

    command = [
        "emapper.py",
        "-i", str(proteins_path),
        "--itype", "proteins",
        "-m", "diamond",
        "--data_dir", str(database_dir),
        "--output_dir", str(output_dir),
        "--output", "metaquest",
        "--cpu", str(threads),
        "--tax_scope", tax_scope,
        "--evalue", str(evalue),
        "--block_size", str(diamond_block_size),
        "--override",
    ]
    log_path = output_dir / "eggnog_mapper.log"
    environment = os.environ.copy()
    environment["PYTHONUNBUFFERED"] = "1"
    with tempfile.TemporaryDirectory(prefix=".emapper-tmp-", dir=output_dir) as temp_dir:
        with log_path.open("w", encoding="utf-8") as log:
            stop_heartbeat = threading.Event()
            started = time.monotonic()

            def heartbeat():
                while not stop_heartbeat.wait(30):
                    elapsed = time.monotonic() - started
                    log.write(
                        f"[MetaQuest] eggNOG-mapper still running "
                        f"({elapsed / 60:.1f} minutes elapsed)\n"
                    )
                    log.flush()

            heartbeat_thread = threading.Thread(target=heartbeat, daemon=True)
            heartbeat_thread.start()
            try:
                result = run_command(
                    command,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                    cwd=temp_dir,
                    env=environment,
                )
            finally:
                stop_heartbeat.set()
                heartbeat_thread.join(timeout=1)
    if result.returncode:
        tail = "\n".join(log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-20:])
        raise AnnotationError(
            f"eggNOG-mapper failed with exit code {result.returncode}. "
            f"See {log_path}. Last output:\n{tail}"
        )
    if not raw_annotations.is_file():
        raise AnnotationError(f"eggNOG-mapper did not create {raw_annotations}")

    table_path, category_path, annotated_count, total_count = _write_functional_outputs(
        proteins_path, raw_annotations, output_dir
    )
    summary = {
        "tool": "eggNOG-mapper",
        "tool_version": mapper_build,
        "database_release": database_release,
        "method": "diamond",
        "tax_scope": tax_scope,
        "evalue": evalue,
        "diamond_block_size": diamond_block_size,
        "threads": threads,
        "total_genes": total_count,
        "annotated_genes": annotated_count,
        "unannotated_genes": total_count - annotated_count,
        "raw_annotations": raw_annotations.name,
        "reused": False,
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    completion_path.write_text(
        json.dumps({**expected_state, "completed": True}, indent=2) + "\n",
        encoding="utf-8",
    )
    return table_path, category_path, annotated_count, False
