"""Annotation stage — metagenomic gene prediction with Pyrodigal."""

from __future__ import annotations

import logging
import json
from pathlib import Path

from metaquest.exceptions import AnnotationError
from metaquest.pipeline.context import AnnotationResult, PipelineContext

logger = logging.getLogger(__name__)


def run_annotation_stage(ctx: PipelineContext) -> PipelineContext:
    """Predict protein-coding genes from assembled metagenomic contigs."""
    from metaquest.core.functional_analysis import run_gene_prediction

    if ctx.assembly is None:
        raise AnnotationError("Assembly stage must complete before gene prediction")

    config = ctx.config.annotation
    logger.info("Running Pyrodigal in metagenomic mode")
    existing_dir = ctx.output_dir / "gene_prediction"
    existing_summary = existing_dir / "summary.json"
    if (
        ctx.resume
        and existing_summary.is_file()
        and (existing_dir / "genes.faa").is_file()
        and (existing_dir / "genes.fna").is_file()
        and (existing_dir / "genes.gff3").is_file()
    ):
        prediction_dir = existing_dir
        gene_count = int(
            json.loads(existing_summary.read_text(encoding="utf-8"))["genes_predicted"]
        )
        from metaquest.io.output_formatter import get_formatter
        get_formatter().info("Reusing completed Pyrodigal gene prediction")
    else:
        try:
            prediction_dir, gene_count = run_gene_prediction(
                ctx.assembly.contigs_fasta,
                ctx.output_dir,
                min_contig_length=config.min_contig_length,
            )
        except Exception as exc:
            raise AnnotationError(
                f"Pyrodigal gene prediction failed: {exc}",
                cause=exc,
            ) from exc

    prediction_dir = Path(prediction_dir)
    ctx.annotation = AnnotationResult(
        gene_prediction_dir=prediction_dir,
        protein_file=prediction_dir / "genes.faa",
        gene_count=gene_count,
    )
    logger.info("Gene prediction complete: %d genes", gene_count)
    return ctx


def run_functional_annotation_stage(ctx: PipelineContext) -> PipelineContext:
    """Annotate predicted proteins with eggNOG-mapper and eggNOG v5."""
    from metaquest.core.functional_analysis import run_functional_annotation
    from metaquest.io.output_formatter import get_formatter

    if ctx.annotation is None or ctx.annotation.protein_file is None:
        raise AnnotationError("Gene prediction must complete before functional annotation")

    config = ctx.config.annotation
    logger.info("Running eggNOG-mapper functional annotation")
    try:
        table, categories, annotated_count, reused = run_functional_annotation(
            ctx.annotation.protein_file,
            ctx.output_dir,
            ctx.config.databases.functional_dir,
            threads=config.threads,
            evalue=config.evalue,
            diamond_block_size=config.diamond_block_size,
            tax_scope=config.tax_scope,
            expected_version=config.eggnog_version,
            database_release=config.eggnog_database_release,
            reuse=ctx.resume,
        )
    except AnnotationError:
        raise
    except Exception as exc:
        raise AnnotationError(f"eggNOG functional annotation failed: {exc}", cause=exc) from exc

    ctx.annotation.functional_annotations = Path(table)
    ctx.annotation.functional_category_summary = Path(categories)
    ctx.annotation.annotated_count = annotated_count
    ctx.annotation.functional_reused = reused
    if reused:
        get_formatter().info("Reusing completed eggNOG functional annotation")
    if annotated_count == 0:
        get_formatter().warning(
            "eggNOG-mapper completed but assigned no annotations; all genes remain reported"
        )
    logger.info(
        "Functional annotation complete: %d/%d genes annotated%s",
        annotated_count,
        ctx.annotation.gene_count,
        " (reused)" if reused else "",
    )
    return ctx
