"""Deterministic publication-oriented plots and their source TSV files."""

from __future__ import annotations

from pathlib import Path
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


COG_LABELS = {
    "A": "RNA processing and modification",
    "C": "Energy production and conversion",
    "D": "Cell cycle and division",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "J": "Translation and ribosome biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "M": "Cell wall and membrane biogenesis",
    "N": "Cell motility",
    "O": "Protein turnover and chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolite biosynthesis",
    "S": "Function unknown",
    "T": "Signal transduction",
    "U": "Intracellular trafficking and secretion",
    "V": "Defense mechanisms",
}


def _save(fig, base: Path, formats: tuple[str, ...], dpi: int) -> list[Path]:
    paths = []
    for extension in formats:
        path = base.with_suffix(f".{extension}")
        fig.savefig(path, dpi=dpi if extension == "png" else None, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def _style(palette: str) -> None:
    sns.set_theme(style="whitegrid", palette=palette, context="paper")


def generate_plots(ctx, taxonomy: dict | None) -> list[Path]:
    config = ctx.config.reporting
    output = Path(ctx.output_dir).resolve() / "plots"
    data_dir = output / "plot_data"
    data_dir.mkdir(parents=True, exist_ok=True)
    _style(config.color_palette)
    sample = Path(ctx.output_dir).name
    generated: list[Path] = []

    validation_rows = []
    for report in getattr(ctx, "metadata", {}).get("validation", {}).get("reports", []):
        bases = report.get("total_bases", 0)
        if bases:
            for threshold in (20, 30):
                validation_rows.append({"file": Path(report["file_path"]).name,
                                        "metric": f"Q{threshold} bases",
                                        "percent": report[f"q{threshold}_bases"] / bases * 100,
                                        "reads_analyzed": report["total_sequences"],
                                        "sampling_limited": report.get("sampling_limited", False)})
    if validation_rows:
        table = pd.DataFrame(validation_rows)
        table.to_csv(data_dir / "validation_quality.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(7.5, 4.2))
        sns.barplot(data=table, x="file", y="percent", hue="metric", ax=ax)
        ax.set(ylim=(0, 100), xlabel="Input file", ylabel="Analyzed bases meeting threshold (%)",
               title=f"{sample}: input FASTQ quality in the validation sample")
        ax.tick_params(axis="x", labelrotation=15)
        ax.legend(title="Base quality", bbox_to_anchor=(1.02, 1), loc="upper left")
        generated += _save(fig, output / "validation_quality", config.plot_formats, config.plot_dpi)

    if ctx.preprocessing and ctx.preprocessing.get("before"):
        before, after = ctx.preprocessing["before"], ctx.preprocessing["after"]
        table = pd.DataFrame([
            {"stage": "Before fastp", "metric": "Q20 bases", "percent": before["q20_rate"] * 100},
            {"stage": "After fastp", "metric": "Q20 bases", "percent": after["q20_rate"] * 100},
            {"stage": "Before fastp", "metric": "Q30 bases", "percent": before["q30_rate"] * 100},
            {"stage": "After fastp", "metric": "Q30 bases", "percent": after["q30_rate"] * 100},
        ])
        table.to_csv(data_dir / "qc_metrics.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        sns.barplot(data=table, x="metric", y="percent", hue="stage", ax=ax)
        ax.set(ylim=(0, 100), xlabel="Quality threshold", ylabel="Bases meeting threshold (%)")
        ax.set_title(f"{sample}: read quality before and after preprocessing")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
        for container in ax.containers:
            ax.bar_label(container, fmt="%.1f%%", padding=2, fontsize=8)
        generated += _save(fig, output / "qc_metrics", config.plot_formats, config.plot_dpi)

        if "reads" in before and "reads" in after:
            retention = pd.DataFrame([
                {"status": "Retained", "reads": after["reads"]},
                {"status": "Filtered", "reads": before["reads"] - after["reads"]},
            ])
            retention.to_csv(data_dir / "qc_retention.tsv", sep="\t", index=False)
            fig, ax = plt.subplots(figsize=(6.4, 4.2))
            sns.barplot(data=retention, x="status", y="reads", hue="status", ax=ax)
            ax.set(xlabel="Read status", ylabel="Individual reads", title=f"{sample}: preprocessing retention")
            for container in ax.containers:
                ax.bar_label(container, fmt="%d", padding=2, fontsize=8)
            generated += _save(fig, output / "qc_retention", config.plot_formats, config.plot_dpi)

        filtering_rows = [
            {"reason": key.removesuffix("_reads").replace("_", " ").title(), "reads": value}
            for key, value in ctx.preprocessing.get("filtering_result", {}).items()
            if key != "passed_filter_reads" and value > 0
        ]
        if filtering_rows:
            filtering = pd.DataFrame(filtering_rows).sort_values("reads", ascending=False)
            filtering.to_csv(data_dir / "qc_filtering_reasons.tsv", sep="\t", index=False)
            fig, ax = plt.subplots(figsize=(6.4, 4.2))
            sns.barplot(data=filtering, y="reason", x="reads", color=sns.color_palette()[0], ax=ax)
            ax.set(xlabel="Filtered individual reads", ylabel="Filtering reason", title=f"{sample}: fastp filtering reasons")
            generated += _save(fig, output / "qc_filtering_reasons", config.plot_formats, config.plot_dpi)

    if taxonomy and taxonomy.get("top_taxa"):
        rows = taxonomy["top_taxa"][: config.top_taxa]
        unit = taxonomy.get("count_unit", "reads")
        visible_rows = [
            row for row in rows if row["fraction_all_input"] >= 0.001
        ]
        top_count = sum(row["estimated_count"] for row in visible_rows)
        total = taxonomy["total_input"]
        plot_rows = [
            {"taxon": row["name"], "percent_all_input": row["fraction_all_input"] * 100}
            for row in visible_rows
        ]
        other_species = max(0, taxonomy["species_assigned"] - top_count)
        classified_above = taxonomy["classified_above_species"]
        for label, count in (
            ("Other species-assigned", other_species),
            ("Classified above species", classified_above),
            ("Unclassified", taxonomy["unclassified_reads"]),
        ):
            if count:
                plot_rows.append({"taxon": label, "percent_all_input": count / total * 100})
        table = pd.DataFrame(plot_rows)
        table.to_csv(data_dir / "taxonomy_top.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(7.2, max(4.2, len(table) * 0.28)))
        sns.barplot(data=table, y="taxon", x="percent_all_input", color=sns.color_palette()[0], ax=ax)
        ax.set(xlabel=f"All cleaned input {unit} (%)", ylabel="Taxon or classification status")
        ax.set_title(f"{sample}: dominant taxa and classification remainder")
        generated += _save(fig, output / "taxonomy_top", config.plot_formats, config.plot_dpi)

    if ctx.assembly and ctx.assembly.lengths:
        lengths = pd.DataFrame({"contig_length_bp": ctx.assembly.lengths})
        lengths.to_csv(data_dir / "assembly_contig_lengths.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        smallest = lengths["contig_length_bp"].min()
        largest = max(smallest + 1, lengths["contig_length_bp"].max())
        bins = np.geomspace(smallest, largest, 40)
        sns.histplot(data=lengths, x="contig_length_bp", bins=bins, ax=ax)
        ax.set_xscale("log")
        ax.set(xlabel="Contig length (bp)", ylabel="Contig count", title=f"{sample}: assembly contig-length distribution")
        generated += _save(fig, output / "assembly_contig_lengths", config.plot_formats, config.plot_dpi)
        cumulative = pd.DataFrame({"contig_length_bp": sorted(ctx.assembly.lengths, reverse=True)})
        cumulative["cumulative_assembly_bp"] = cumulative["contig_length_bp"].cumsum()
        cumulative.insert(0, "contig_rank", range(1, len(cumulative) + 1))
        cumulative.to_csv(data_dir / "assembly_cumulative_length.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        sns.lineplot(data=cumulative, x="contig_rank", y="cumulative_assembly_bp", ax=ax)
        ax.axvline(ctx.assembly.l50, color="tab:orange", linestyle="--", label="L50")
        ax.axvline(ctx.assembly.l90, color="tab:green", linestyle="--", label="L90")
        ax.set(xlabel="Contig rank (longest first)", ylabel="Cumulative assembly length (Mbp)", title=f"{sample}: cumulative assembly length")
        ax.ticklabel_format(axis="y", style="plain")
        ticks = ax.get_yticks()
        ax.set_yticks(ticks, [f"{tick / 1_000_000:g}" for tick in ticks])
        ax.legend()
        generated += _save(fig, output / "assembly_cumulative_length", config.plot_formats, config.plot_dpi)

    if ctx.annotation and ctx.annotation.functional_category_summary:
        source = pd.read_csv(ctx.annotation.functional_category_summary, sep="\t")
        for namespace in ("COG", "KO", "EC", "GO"):
            table = source[source["namespace"] == namespace].nlargest(config.top_functional_terms, "gene_count")
            if table.empty:
                continue
            table = table.copy()
            if namespace == "COG":
                table["display_term"] = table["term"].map(
                    lambda term: f"{term} — {COG_LABELS.get(term, 'Other COG category')}"
                )
            else:
                table["display_term"] = table["term"]
            table.to_csv(data_dir / f"functional_{namespace.lower()}.tsv", sep="\t", index=False)
            fig, ax = plt.subplots(figsize=(7.2, max(4.2, len(table) * 0.28)))
            sns.barplot(data=table, y="display_term", x="gene_count", color=sns.color_palette()[0], ax=ax)
            ax.set(xlabel="Genes assigned to term", ylabel=f"{namespace} term", title=f"{sample}: top {namespace} assignments")
            generated += _save(fig, output / f"functional_{namespace.lower()}", config.plot_formats, config.plot_dpi)
    abundance = getattr(ctx, "abundance", None)
    if abundance and abundance.functional_abundance:
        source = pd.read_csv(abundance.functional_abundance, sep="\t")
        for namespace in ("COG", "KO", "EC", "GO"):
            table = source[source["namespace"] == namespace].nlargest(config.top_functional_terms, "tpm")
            if table.empty:
                continue
            table = table.copy()
            if namespace == "COG":
                table["display_term"] = table["term"].map(
                    lambda term: f"{term} — {COG_LABELS.get(term, 'Other COG category')}"
                )
            else:
                table["display_term"] = table["term"]
            stem = f"abundance_{namespace.lower()}"
            table.to_csv(data_dir / f"{stem}.tsv", sep="\t", index=False)
            fig, ax = plt.subplots(figsize=(7.2, max(4.2, len(table) * 0.28)))
            sns.barplot(data=table, y="display_term", x="tpm", color=sns.color_palette()[1], ax=ax)
            ax.set(xlabel="Summed gene TPM", ylabel=f"{namespace} term", title=f"{sample}: most abundant {namespace} terms")
            generated += _save(fig, output / stem, config.plot_formats, config.plot_dpi)
    if ctx.annotation:
        gff = getattr(ctx.annotation, "gene_prediction_dir", Path(ctx.output_dir) / "gene_prediction") / "genes.gff3"
        lengths, partials = Counter(), Counter()
        if gff.is_file():
            with gff.open() as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue
                    fields = line.rstrip().split("\t")
                    if len(fields) != 9 or fields[2] != "CDS":
                        continue
                    length = int(fields[4]) - int(fields[3]) + 1
                    lengths[min(length // 100, 100)] += 1
                    attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
                    partials["Complete" if attrs.get("partial") == "00" else "Partial or unknown"] += 1
        if lengths:
            table = pd.DataFrame([{"length_bin_start_bp": key * 100, "gene_count": value}
                                  for key, value in sorted(lengths.items())])
            table.to_csv(data_dir / "gene_prediction_lengths.tsv", sep="\t", index=False)
            fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
            axes[0].bar(table["length_bin_start_bp"], table["gene_count"], width=90, color="#177d8d")
            axes[0].set(xlabel="CDS length bin start (bp; ≥10 kb pooled)", ylabel="Predicted genes", title="Gene lengths")
            pd.DataFrame(list(partials.items()), columns=["status", "gene_count"]).to_csv(data_dir / "gene_prediction_completeness.tsv", sep="\t", index=False)
            axes[1].bar(list(partials), list(partials.values()), color=["#177d8d", "#d3a54e"])
            axes[1].set(ylabel="Predicted genes", title="Gene completeness")
            fig.suptitle(f"{sample}: gene prediction (GFF coordinates: 1-based inclusive)")
            generated += _save(fig, output / "gene_prediction_lengths", config.plot_formats, config.plot_dpi)
    screening = getattr(ctx, "screening", None)
    if screening:
        candidates = pd.read_csv(Path(ctx.output_dir) / screening["taxon_table"], sep="\t")
        candidates.to_csv(data_dir / "screening_taxa.tsv", sep="\t", index=False)
        if not candidates.empty:
            table = candidates.head(config.top_taxa)
            fig, ax = plt.subplots(figsize=(8, max(4, len(table) * .3)))
            sns.barplot(data=table, y="name", x="fraction_all_input", color="#177d8d", ax=ax)
            ax.set(xlabel=f"Fraction of all cleaned input {screening['count_unit']}", ylabel="Catalog taxon candidate",
                   title=f"{sample}: human-associated bacterial screening candidates")
            generated += _save(fig, output / "screening_taxa", config.plot_formats, config.plot_dpi)
        types = pd.DataFrame([{"evidence_type": key, "matches": screening["protein_types"].get(key, 0)}
                              for key in ("AMR", "VIRULENCE", "STRESS")])
        types.to_csv(data_dir / "screening_proteins.tsv", sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(6.4, 4.2))
        sns.barplot(data=types, x="evidence_type", y="matches", color="#177d8d", ax=ax)
        ax.set(xlabel="Curated evidence class", ylabel="Protein evidence rows",
               title=f"{sample}: AMRFinderPlus matches (not clinical risk)")
        ax.set_ylim(0, max(1, types["matches"].max() * 1.15))
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        generated += _save(fig, output / "screening_proteins", config.plot_formats, config.plot_dpi)
    return generated
