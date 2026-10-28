"""Descriptive reports for the validated MetaQuest pipeline surface."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def _native(value: Any) -> Any:
    """Convert pandas/numpy scalar values into JSON-compatible values."""
    if pd.isna(value):
        return None
    return value.item() if hasattr(value, "item") else value


def _taxonomy_summary(classification, *, paired: bool = False, top_limit: int = 20) -> tuple[dict[str, Any], str]:
    bracken_file = classification.bracken_file
    table = pd.read_csv(bracken_file, sep="\t")
    required = {"name", "new_est_reads"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(
            f"Bracken report is missing required columns: {', '.join(sorted(missing))}"
        )

    species_assigned = classification.species_assigned_reads
    total_reads = classification.total_reads
    denominator_label = "fragments" if paired else "reads"
    table = table.sort_values("new_est_reads", ascending=False)
    top = []
    for _, row in table.head(top_limit).iterrows():
        top.append(
            {
                "name": str(row["name"]),
                "estimated_reads": int(row["new_est_reads"]),
                "estimated_count": int(row["new_est_reads"]),
                "fraction_species_assigned_reads": (
                    int(row["new_est_reads"]) / species_assigned if species_assigned else 0.0
                ),
                "fraction_all_input_reads": (
                    int(row["new_est_reads"]) / total_reads if total_reads else 0.0
                ),
                "fraction_all_input_fragments": (
                    int(row["new_est_reads"]) / total_reads if total_reads else 0.0
                ),
                "fraction_all_input": (
                    int(row["new_est_reads"]) / total_reads if total_reads else 0.0
                ),
            }
        )

    summary = {
        "reported_taxa": int(len(table)),
        "count_unit": denominator_label,
        "total_input": total_reads,
        "total_input_reads": total_reads,
        "kraken_classified": classification.kraken_classified_reads,
        "kraken_classified_reads": classification.kraken_classified_reads,
        "unclassified_reads": classification.unclassified_reads,
        "classification_rate": (
            classification.kraken_classified_reads / total_reads if total_reads else 0.0
        ),
        "species_assigned_reads": species_assigned,
        "species_assigned": species_assigned,
        "classified_above_species": max(
            0, classification.kraken_classified_reads - species_assigned
        ),
        "observed_read_length": classification.observed_read_length,
        "bracken_read_length": classification.bracken_read_length,
        "top_taxa": top,
        "primary_denominator": f"all cleaned input {denominator_label}",
    }

    lines = [
        "METAQUEST TAXONOMIC PROFILE",
        "============================",
        "",
        f"Reported taxa: {summary['reported_taxa']}",
        f"Total cleaned input {denominator_label}: {summary['total_input_reads']}",
        f"Kraken-classified {denominator_label}: {summary['kraken_classified_reads']}",
        f"Unclassified {denominator_label}: {summary['unclassified_reads']}",
        f"Classification rate: {summary['classification_rate']:.2%}",
        f"Species-assigned {denominator_label}: {summary['species_assigned_reads']}",
        f"Observed read length: {summary['observed_read_length']} bp",
        f"Bracken model length: {summary['bracken_read_length']} bp",
        "",
        "Top taxa",
        "--------",
        f"Taxon\tEstimated {denominator_label}\tFraction species-assigned\tFraction all input",
    ]
    lines.extend(
        f"{item['name']}\t{item['estimated_reads']}\t"
        f"{item['fraction_species_assigned_reads']:.6f}\t"
        f"{item['fraction_all_input_reads']:.6f}"
        for item in top
    )
    lines.extend(
        [
            "",
            "Interpretation note",
            "-------------------",
            "These are computational taxonomic abundance estimates from Kraken2/Bracken.",
            "They are research outputs and do not establish pathogenicity or clinical risk.",
        ]
    )
    return summary, "\n".join(lines) + "\n"


def generate_stable_reports(ctx) -> None:
    """Write descriptive text and JSON outputs without custom risk scoring."""
    output_dir = Path(ctx.output_dir).resolve()
    summary: dict[str, Any] = {
        "analysis_type": "short_read_metagenomics",
        "research_use_only": True,
        "completed_stages": [*ctx.completed_stages, "Reporting"],
        "validation": getattr(ctx, "metadata", {}).get("validation", {"status": "not_recorded"}),
    }

    if ctx.classification:
        reporting_config = getattr(getattr(ctx, "config", None), "reporting", None)
        taxonomy, report = _taxonomy_summary(
            ctx.classification,
            paired=getattr(ctx, "read_mode", "single") != "single",
            top_limit=getattr(reporting_config, "top_taxa", 20),
        )
        summary["taxonomy"] = taxonomy
        (output_dir / "01_taxonomic_report.txt").write_text(report, encoding="utf-8")

    if ctx.assembly:
        summary["assembly"] = {
            "total_contigs": ctx.assembly.total_contigs,
            "total_bases": ctx.assembly.total_bases,
            "n50": ctx.assembly.n50,
            "max_length": ctx.assembly.max_length,
            "mean_length": ctx.assembly.mean_length,
            "l50": ctx.assembly.l50,
            "n90": ctx.assembly.n90,
            "l90": ctx.assembly.l90,
            "gc_percent": ctx.assembly.gc_percent,
            "contigs_ge_1000": ctx.assembly.contigs_ge_1000,
            "contigs_ge_5000": ctx.assembly.contigs_ge_5000,
            "contigs_ge_10000": ctx.assembly.contigs_ge_10000,
        }
        assembly_lines = [
            "METAQUEST METAGENOMIC ASSEMBLY",
            "===============================",
            "",
            f"Contigs: {ctx.assembly.total_contigs}",
            f"Total bases: {ctx.assembly.total_bases}",
            f"N50: {ctx.assembly.n50} bp",
            f"L50: {ctx.assembly.l50} contigs",
            f"N90: {ctx.assembly.n90} bp",
            f"L90: {ctx.assembly.l90} contigs",
            f"GC content: {ctx.assembly.gc_percent:.2f}%",
            f"Maximum contig length: {ctx.assembly.max_length} bp",
            f"Mean contig length: {ctx.assembly.mean_length:.2f} bp",
            f"Contigs >= 1,000 bp: {ctx.assembly.contigs_ge_1000}",
            f"Contigs >= 5,000 bp: {ctx.assembly.contigs_ge_5000}",
            f"Contigs >= 10,000 bp: {ctx.assembly.contigs_ge_10000}",
            "",
            "Interpretation note",
            "-------------------",
            "Assembly continuity does not by itself establish assembly correctness.",
        ]
        (output_dir / "02_assembly_report.txt").write_text(
            "\n".join(assembly_lines) + "\n", encoding="utf-8"
        )

    if ctx.annotation:
        prediction_path = getattr(ctx.annotation, "gene_prediction_dir", output_dir / "gene_prediction") / "summary.json"
        summary["gene_prediction"] = (
            json.loads(prediction_path.read_text()) if prediction_path.is_file()
            else {"genes_predicted": ctx.annotation.gene_count}
        )
        functional_completed = bool(ctx.annotation.functional_annotations)
        annotation_summary = {
            "predicted_genes": ctx.annotation.gene_count,
            "annotated_genes": ctx.annotation.annotated_count if functional_completed else None,
            "unannotated_genes": ctx.annotation.gene_count - ctx.annotation.annotated_count if functional_completed else None,
            "functional_status": (
                "completed" if ctx.annotation.functional_annotations else "skipped"
            ),
        }
        if ctx.annotation.functional_annotations:
            functional_dir = ctx.annotation.functional_annotations.parent
            eggnog_summary_path = functional_dir / "summary.json"
            eggnog_summary = json.loads(eggnog_summary_path.read_text(encoding="utf-8"))
            annotation_summary.update(
                {
                    "functional_annotations": str(
                        ctx.annotation.functional_annotations.resolve().relative_to(output_dir)
                    ),
                    "functional_category_summary": str(
                        ctx.annotation.functional_category_summary.resolve().relative_to(output_dir)
                    ),
                    "functional_reused": ctx.annotation.functional_reused,
                    "eggnog_mapper_version": eggnog_summary["tool_version"],
                    "eggnog_database_release": eggnog_summary["database_release"],
                    "tax_scope": eggnog_summary["tax_scope"],
                }
            )
            lines = [
                "METAQUEST FUNCTIONAL ANNOTATION",
                "===============================",
                "",
                f"Predicted genes: {ctx.annotation.gene_count}",
                f"eggNOG-annotated genes: {ctx.annotation.annotated_count}",
                f"Unannotated genes: {ctx.annotation.gene_count - ctx.annotation.annotated_count}",
                f"Annotation rate: {ctx.annotation.annotated_count / ctx.annotation.gene_count:.2%}" if ctx.annotation.gene_count else "Annotation rate: unavailable",
                f"eggNOG-mapper: {eggnog_summary['tool_version']}",
                f"eggNOG database: {eggnog_summary['database_release']}",
                f"Taxonomic scope: {eggnog_summary['tax_scope']}",
                "",
                "Interpretation note",
                "-------------------",
                "Annotations are computational orthology-based assignments for research use.",
                "They do not establish gene expression, phenotype, or clinical significance.",
            ]
            (output_dir / "03_functional_report.txt").write_text(
                "\n".join(lines) + "\n", encoding="utf-8"
            )
        summary["annotation"] = annotation_summary
    abundance = getattr(ctx, "abundance", None)
    if abundance:
        abundance_summary = {
            key: value for key, value in abundance.summary.items() if key != "functional"
        }
        abundance_summary["gene_abundance"] = str(
            abundance.gene_abundance.resolve().relative_to(output_dir)
        )
        if abundance.functional_abundance:
            abundance_summary["functional_abundance"] = str(
                abundance.functional_abundance.resolve().relative_to(output_dir)
            )
            abundance_summary["functional"] = abundance.summary.get("functional", {})
        summary["abundance"] = abundance_summary
        unit = abundance_summary["count_unit"]
        lines = [
            "METAQUEST GENE ABUNDANCE",
            "========================",
            "",
            f"Mapper: BBMap {abundance_summary['tool_version']}",
            f"Minimum alignment identity: {abundance_summary['min_identity']}",
            f"Multi-mapping policy: {abundance_summary['ambiguous']}",
            f"Input {unit}: {abundance_summary[f'input_{unit}']}",
            f"Mapped {unit}: {abundance_summary[f'mapped_{unit}']}",
            f"Mapping rate: {abundance_summary['mapping_rate']:.2%}",
            f"Reads excluded for unknown identity: {abundance_summary.get('reads_unknown_identity', 0)}",
            f"{unit.capitalize()} assigned to genes: {abundance_summary[f'{unit}_in_genes']}",
            f"Intergenic {unit}: {abundance_summary[f'intergenic_{unit}']}",
            f"Gene assignment rate: {abundance_summary['gene_assignment_rate']:.2%}",
            f"Genes with counts: {abundance_summary['genes_with_counts']} / {abundance_summary['total_genes']}",
            f"Normalization: {abundance_summary['normalization']}",
        ]
        functional = abundance_summary.get("functional")
        if functional:
            lines += [
                f"Term attribution: {functional['term_attribution']}",
                f"TPM in eggNOG-annotated genes: {functional['tpm_in_annotated_genes']:.1f} / 1,000,000",
            ]
        lines += [
            "",
            "Interpretation note",
            "-------------------",
            "Abundance describes reads that map to assembled, predicted genes only.",
            "Unassembled and intergenic reads are excluded from TPM denominators.",
            "TPM values are compositional; cross-sample comparisons need appropriate",
            "compositional or count-based statistics.",
        ]
        (output_dir / "04_abundance_report.txt").write_text(
            "\n".join(lines) + "\n", encoding="utf-8"
        )
    if getattr(ctx, "preprocessing", None):
        summary["preprocessing"] = ctx.preprocessing
    screening = getattr(ctx, "screening", None)
    if screening:
        summary["screening"] = screening
        lines = ["METAQUEST RESEARCH SCREENING", "=" * 28,
                 f"Taxonomic candidates: {screening['taxon_candidates']} / {screening['catalog_targets']} targets",
                 f"Screened predicted proteins: {screening['screened_genes']}",
                 f"Protein evidence rows: {screening['protein_matches']}",
                 f"AMRFinderPlus version: {screening['tool_version']}",
                 f"Database release: {screening['database_release']}",
                 "", "Interpretation limits", *screening["limitations"]]
        (output_dir / "05_screening_report.txt").write_text("\n".join(lines) + "\n")
    if hasattr(ctx, "config"):
        from metaquest.reporting.plots import generate_plots
        figures = generate_plots(ctx, summary.get("taxonomy"))
        summary["figures"] = [
            str(path.resolve().relative_to(output_dir)) for path in figures
        ]
    else:
        summary["figures"] = []
    from metaquest.reporting.steps import write_step_reports

    summary["stage_reports"] = write_step_reports(output_dir, summary)
    summary_path = output_dir / "analysis_summary.json"
    temporary_summary = summary_path.with_suffix(".json.tmp")
    temporary_summary.write_text(
        json.dumps(summary, indent=2, default=_native) + "\n",
        encoding="utf-8",
    )
    temporary_summary.replace(summary_path)
    _write_html_report(output_dir, summary)


def _write_html_report(output_dir: Path, summary: dict[str, Any]) -> None:
    """Write a readable, dependency-free report with embedded SVG figures."""
    import base64
    import html

    def metric_table(rows: list[tuple[str, Any]]) -> str:
        rendered = "".join(
            f"<tr><th>{html.escape(label)}</th><td>{html.escape(str(value))}</td></tr>"
            for label, value in rows
        )
        return f'<table class="metrics">{rendered}</table>'

    sections = []
    stage_links = "".join(
        f'<tr><td><a href="{html.escape(item["html"], quote=True)}">{html.escape(item["stage"])}</a></td>'
        f'<td>{html.escape(item["status"])}</td></tr>' for item in summary.get("stage_reports", [])
    )
    sections.append('<section><h2>Workflow evidence</h2><table><tr><th>Step report</th><th>Status</th></tr>' + stage_links + '</table></section>')
    validation = summary.get("validation", {})
    sections.append(f'<section><h2>Input validation</h2>{metric_table([("Status", validation.get("status", "not_recorded")), ("Strict checks", validation.get("strict", False))])}</section>')
    preprocessing = summary.get("preprocessing")
    if preprocessing:
        before = preprocessing.get("before", {})
        after = preprocessing.get("after", {})
        before_reads = before.get("reads", 0)
        rows = [
            ("Tool version", preprocessing.get("fastp_version", "unknown")),
            ("Input reads", f"{before_reads:,}"),
            ("Retained reads", f"{after.get('reads', 0):,}"),
            (
                "Read retention",
                f"{after.get('reads', 0) / before_reads:.2%}"
                if before_reads
                else "unavailable",
            ),
            ("Mean read length", f"{after.get('mean_length', 0):.1f} bp"),
            ("Q30 bases after filtering", f"{after.get('q30_rate', 0):.2%}"),
        ]
        sections.append(f"<section><h2>Preprocessing</h2>{metric_table(rows)}</section>")

    taxonomy = summary.get("taxonomy")
    if taxonomy:
        unit = taxonomy["count_unit"]
        rows = [
            (f"Cleaned input {unit}", f"{taxonomy['total_input']:,}"),
            (f"Kraken-classified {unit}", f"{taxonomy['kraken_classified']:,}"),
            ("Classification rate", f"{taxonomy['classification_rate']:.2%}"),
            (f"Species-assigned {unit}", f"{taxonomy['species_assigned']:,}"),
            (f"Classified above species {unit}", f"{taxonomy['classified_above_species']:,}"),
            (f"Unclassified {unit}", f"{taxonomy['unclassified_reads']:,}"),
        ]
        sections.append(f"<section><h2>Taxonomy</h2>{metric_table(rows)}</section>")

    assembly = summary.get("assembly")
    if assembly:
        rows = [
            ("Contigs", f"{assembly['total_contigs']:,}"),
            ("Assembly size", f"{assembly['total_bases']:,} bp"),
            ("N50 / L50", f"{assembly['n50']:,} bp / {assembly['l50']:,} contigs"),
            ("N90 / L90", f"{assembly['n90']:,} bp / {assembly['l90']:,} contigs"),
            ("GC content", f"{assembly['gc_percent']:.2f}%"),
            ("Contigs ≥1 kb", f"{assembly['contigs_ge_1000']:,}"),
            ("Contigs ≥5 kb", f"{assembly['contigs_ge_5000']:,}"),
            ("Contigs ≥10 kb", f"{assembly['contigs_ge_10000']:,}"),
        ]
        sections.append(f"<section><h2>Assembly</h2>{metric_table(rows)}</section>")

    if summary.get("gene_prediction"):
        prediction = summary["gene_prediction"]
        sections.append(f'<section><h2>Gene prediction</h2>{metric_table([(key.replace("_", " "), value) for key, value in prediction.items() if not isinstance(value, (dict, list))])}</section>')
    annotation = summary.get("annotation")
    if annotation:
        genes = annotation["predicted_genes"]
        annotated = annotation["annotated_genes"]
        rows = [
            ("Predicted genes", f"{genes:,}"),
            ("Annotated genes", f"{annotated:,}" if annotated is not None else "not assessed"),
            ("Annotation rate", f"{annotated / genes:.2%}" if genes and annotated is not None else "not assessed"),
            ("Functional annotation", annotation["functional_status"]),
            ("eggNOG-mapper", annotation.get("eggnog_mapper_version", "not run")),
            ("eggNOG database", annotation.get("eggnog_database_release", "not run")),
        ]
        sections.append(f"<section><h2>Functional annotation</h2>{metric_table(rows)}</section>")

    abundance = summary.get("abundance")
    if abundance:
        unit = abundance["count_unit"]
        rows = [
            ("Mapper", f"BBMap {abundance['tool_version']}"),
            (f"Input {unit}", f"{abundance[f'input_{unit}']:,}"),
            ("Mapping rate", f"{abundance['mapping_rate']:.2%}"),
            ("Reads excluded for unknown identity", f"{abundance.get('reads_unknown_identity', 0):,}"),
            ("Gene assignment rate", f"{abundance['gene_assignment_rate']:.2%}"),
            ("Genes with counts", f"{abundance['genes_with_counts']:,} / {abundance['total_genes']:,}"),
            ("Normalization", abundance["normalization"]),
        ]
        functional = abundance.get("functional")
        if functional:
            rows.append(("Term attribution", functional["term_attribution"]))
        sections.append(f"<section><h2>Gene abundance</h2>{metric_table(rows)}</section>")

    screening = summary.get("screening")
    if screening:
        rows = [("Taxonomic candidates", screening["taxon_candidates"]),
                ("Protein evidence rows", screening["protein_matches"]),
                ("Screened genes", screening["screened_genes"]),
                ("AMRFinderPlus", screening["tool_version"]),
                ("Database release", screening["database_release"])]
        notes = "".join(f"<li>{html.escape(note)}</li>" for note in screening["limitations"])
        sections.append(f'<section><h2>Research screening</h2>{metric_table(rows)}<ul>{notes}</ul></section>')

    captions = {
        "validation_quality": "Input quality rates in the validation sample; analyzed read counts and sampling flags are in the source TSV.",
        "qc_metrics": "Base-quality rates before and after preprocessing.",
        "qc_retention": "Individual reads retained and filtered by fastp.",
        "qc_filtering_reasons": "Individual reads removed by fastp, grouped by reason.",
        "taxonomy_top": "Dominant taxa and mutually exclusive classification remainder.",
        "assembly_contig_lengths": "Assembly contig lengths on a logarithmic scale.",
        "assembly_cumulative_length": "Cumulative assembly size with L50 and L90 markers.",
        "abundance_cog": "COG categories ranked by summed gene TPM.",
        "abundance_ko": "KEGG orthologs ranked by summed gene TPM.",
        "abundance_ec": "Enzyme Commission terms ranked by summed gene TPM.",
        "abundance_go": "Gene Ontology terms ranked by summed gene TPM.",
        "gene_prediction_lengths": "Predicted CDS lengths and partial-gene status from 1-based inclusive GFF coordinates.",
        "screening_taxa": "Exact species-taxid catalog matches; these classifications do not establish pathogenic strains or infection.",
        "screening_proteins": "Native curated AMRFinderPlus protein evidence classes; counts are evidence rows, not clinical risk scores.",
    }
    figures = []
    for relative in summary.get("figures", []):
        path = output_dir / relative
        if path.suffix.lower() not in (".svg", ".png"):
            continue
        if path.suffix.lower() == ".png" and path.with_suffix(".svg").is_file():
            continue
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        stem = path.stem
        caption = captions.get(
            stem,
            stem.replace("functional_", "Top ").replace("_", " ").title()
            + " assignments.",
        )
        figures.append(
            '<figure><img src="data:' + ("image/svg+xml" if path.suffix == ".svg" else "image/png") + ';base64,'
            f'{encoded}" alt="{html.escape(caption)}">'
            f"<figcaption>{html.escape(caption)}</figcaption></figure>"
        )

    title = "MetaQuest analysis report"
    from metaquest.reporting.pages import STYLE

    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title>
<style>
{STYLE}
body{{font-family:system-ui,sans-serif;max-width:1100px;margin:auto;padding:2rem;color:#17202a}}
section,figure{{border:1px solid #d5d8dc;border-radius:.5rem;padding:1rem;margin:1rem 0}}
.metrics{{border-collapse:collapse;width:100%}}.metrics th,.metrics td{{padding:.55rem;border-bottom:1px solid #e5e7e9;text-align:left}}
.metrics th{{width:45%;color:#34495e}}img{{display:block;max-width:100%;height:auto;margin:auto}}
figcaption{{margin-top:.75rem;color:#566573}}.note{{background:#fff4cc;padding:1rem;border-radius:.5rem}}
</style></head><body><h1>{title}</h1>
<p class="note">Research-use descriptive output. No statistical significance, pathogenicity, or clinical risk is inferred.</p>
{"".join(sections)}
<h2>Figures</h2>{"".join(figures) if figures else "<p>No figures were generated.</p>"}
</body></html>"""
    (output_dir / "report.html").write_text(document, encoding="utf-8")
