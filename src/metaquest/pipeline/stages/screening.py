"""Curated bacterial protein and human-associated taxon research screening."""

from metaquest.core.screening import ScreeningError, run_screening


def run_screening_stage(ctx):
    if ctx.annotation is None or ctx.classification is None:
        raise ScreeningError("Screening requires taxonomic classification and predicted proteins")
    config = ctx.config.screening
    ctx.screening = run_screening(
        ctx.annotation.protein_file, ctx.annotation.gene_prediction_dir / "genes.gff3",
        ctx.classification.bracken_file, ctx.output_dir,
        config.database_dir or ctx.config.databases.base_dir / "amr",
        total_input=ctx.classification.total_reads,
        count_unit="reads" if ctx.read_mode == "single" else "fragments",
        catalog=config.taxa_catalog, threads=config.threads, reuse=ctx.resume,
        abundance=ctx.abundance.gene_abundance if ctx.abundance else None,
    )
    return ctx
