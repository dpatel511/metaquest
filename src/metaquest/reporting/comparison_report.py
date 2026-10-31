"""Descriptive cohort graphics and explicitly labeled group-test outputs."""

import html

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from metaquest.reporting.pages import figure_html, write_page
from metaquest.reporting.plots import _save, _style


def write_comparison_report(output, relative, labels, qc, distances, tests, summary):
    plots = output / "plots"
    plots.mkdir(exist_ok=True)
    _style("colorblind")
    images = []
    top = relative.mean().nlargest(10).index
    composition = relative.loc[:, top].copy()
    composition["Other features"] = (1 - composition.sum(axis=1)).clip(lower=0)
    composition.columns = [labels.get(key, key) for key in composition.columns]
    fig, ax = plt.subplots(figsize=(max(7.5, len(relative) * .55), 5))
    composition.plot.bar(stacked=True, ax=ax, colormap="tab20", width=.8)
    ax.set(xlabel="Sample", ylabel="Fraction of reported profile total", ylim=(0, 1),
           title=f"{summary['feature']}: composition across samples")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    composition.to_csv(plots / "composition.tsv", sep="\t", index_label="sample_id")
    _save(fig, plots / "composition", ("svg", "png"), 200)
    images.append(figure_html(plots / "composition.svg", summary["relative_denominator"]))
    fig, ax = plt.subplots(figsize=(max(5, len(relative) * .5), max(4, len(relative) * .45)))
    sns.heatmap(distances, vmin=0, vmax=1, cmap="Blues", square=True, ax=ax,
                cbar_kws={"label": "Bray-Curtis dissimilarity"})
    ax.set(title="Between-sample composition differences")
    _save(fig, plots / "bray_curtis", ("svg", "png"), 200)
    images.append(figure_html(plots / "bray_curtis.svg", "Descriptive distances; no group significance is inferred from this heatmap."))
    heatmap = relative.loc[:, relative.mean().nlargest(20).index].T
    fig, ax = plt.subplots(figsize=(max(6, len(relative) * .5), max(4, len(heatmap) * .28)))
    sns.heatmap(heatmap, cmap="viridis", ax=ax, yticklabels=[labels.get(key, key) for key in heatmap.index],
                cbar_kws={"label": "Relative profile abundance"})
    ax.set(title="Dominant features", xlabel="Sample", ylabel="Feature")
    heatmap.to_csv(plots / "feature_heatmap.tsv", sep="\t", index_label="feature_id")
    _save(fig, plots / "feature_heatmap", ("svg", "png"), 200)
    images.append(figure_html(plots / "feature_heatmap.svg", "The full exported matrix, not these top features, is used for statistical testing."))
    overview = pd.DataFrame([
        {"metric": "Samples", "value": summary["sample_count"]},
        {"metric": "Aligned features", "value": summary["feature_count"]},
        {"metric": "Input values", "value": summary["input_unit"]},
        {"metric": "Method", "value": summary["statistics"]["method"]},
    ])
    sections = [("Study overview", overview.to_html(index=False, escape=True)),
                ("Sample quality and coverage", '<div class="scroll">' + qc.to_html(index=False, escape=True, na_rep="not run") + '</div>'),
                ("Composition and distances", "".join(images))]
    if tests is not None:
        effect = tests.assign(absolute_effect=tests["median_clr_difference"].abs()).nlargest(20, "absolute_effect")
        fig, ax = plt.subplots(figsize=(8, max(4, len(effect) * .28)))
        colors = ["#177d8d" if value else "#b9c7cd" for value in effect["passes_fdr"]]
        ax.barh(range(len(effect)), effect["median_clr_difference"], color=colors)
        ax.set_yticks(range(len(effect)), effect["label"])
        ax.invert_yaxis()
        ax.axvline(0, color="#526779", linewidth=.8)
        ax.set(xlabel="Median CLR difference (comparison − reference)", ylabel="Feature",
               title=f"{summary['statistics']['comparison_group']} vs {summary['statistics']['reference_group']}")
        effect.drop(columns="absolute_effect").to_csv(plots / "group_effects.tsv", sep="\t", index=False)
        _save(fig, plots / "group_effects", ("svg", "png"), 200)
        settings = pd.DataFrame(list(summary["statistics"].items()), columns=["setting", "value"])
        sections.append(("Exploratory independent-group tests", settings.to_html(index=False, escape=True)
                         + figure_html(plots / "group_effects.svg", "Teal bars meet the selected FDR threshold; gray bars do not. Effects are log-ratios, not absolute fold changes.")
                         + '<div class="scroll">' + tests.head(30).to_html(index=False, escape=True) + '</div>'))
    notes = "<ul>" + "".join(f"<li>{html.escape(note)}</li>" for note in summary["limitations"]) + "</ul>"
    sections.append(("Methods, provenance and interpretation", notes + '<p>Full results: feature_values.tsv, relative_abundance.tsv, sample_qc.tsv, features.tsv, bray_curtis.tsv and comparison_summary.json.</p>'))
    write_page(output / "report.html", "MetaQuest comparative analysis", sections,
               note="Research comparison of recorded profiles. Formal results, when requested, are exploratory log-ratio tests for two independent groups; they do not establish causality or clinical risk.")
    text = ["METAQUEST COMPARATIVE ANALYSIS", "=" * 30,
            f"Samples: {len(relative)}", f"Feature: {summary['feature']}",
            f"Features: {relative.shape[1]}", f"Method: {summary['statistics']['method']}"]
    if tests is not None:
        text += [f"Features passing FDR: {int(tests['passes_fdr'].sum())} / {len(tests)}",
                 f"Correction: {summary['statistics']['fdr_method']}; alpha={summary['statistics']['alpha']}",
                 f"Contrast: {summary['statistics']['comparison_group']} - {summary['statistics']['reference_group']}",
                 f"Transformation: {summary['statistics']['transformation']}",
                 f"Pseudocount: {summary['statistics']['pseudocount']} (proportion of full profile total)",
                 f"Prevalence threshold: {summary['statistics']['min_prevalence']}; tested={summary['statistics']['features_tested']}; excluded={summary['statistics']['features_excluded']}",
                 f"Design: {summary['statistics']['design']}",
                 f"P-value policy: {summary['statistics']['pvalue_policy']}",
                 f"Small tied features: {summary['statistics']['small_tied_features']}"]
    text += ["", "Interpretation limits", *summary["limitations"]]
    (output / "comparison_report.txt").write_text("\n".join(text) + "\n")
