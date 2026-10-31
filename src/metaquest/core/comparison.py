"""Aligned study tables and explicit exploratory independent-group comparisons.

Full profiles are joined by stable taxonomy IDs or functional term IDs, never
by the top-N rows displayed in a report. Formal tests compare CLR distributions,
not absolute abundance, and do not adjust for covariates or repeated measures.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from metaquest.core.abundance import _sha256
from metaquest.config import __version__
from metaquest.exceptions import MetaQuestError


class ComparisonError(MetaQuestError):
    """Study inputs, provenance or requested statistics are incompatible."""


def _json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ComparisonError(f"Cannot read MetaQuest results at {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ComparisonError(f"Expected a JSON object at {path}")
    return value


def _table(path: Path, *, separator: str = "\t", dtype=None) -> pd.DataFrame:
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            columns = next(csv.reader(handle, delimiter=separator), [])
        if len(columns) != len(set(columns)):
            raise ComparisonError(f"Duplicate columns in table: {path}")
        return pd.read_csv(path, sep=separator, dtype=dtype, keep_default_na=False,
                           encoding="utf-8-sig")
    except (OSError, UnicodeError, csv.Error, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        raise ComparisonError(f"Cannot read comparison table {path}: {exc}") from exc


def load_samples(path: Path, group_column: str | None = None) -> pd.DataFrame:
    path = Path(path).resolve()
    samples = _table(path, separator="\t" if path.suffix.lower() == ".tsv" else ",", dtype=str)
    required = {"sample_id", "result_dir"} | ({group_column} if group_column else set())
    if not required.issubset(samples.columns):
        raise ComparisonError(f"Sample sheet requires columns: {sorted(required)}")
    for column in required:
        samples[column] = samples[column].str.strip()
        if (samples[column] == "").any():
            raise ComparisonError(f"Missing sample-sheet values in {column}")
    if "subject_id" in samples:
        samples["subject_id"] = samples["subject_id"].str.strip()
    if len(samples) < 2 or samples["sample_id"].duplicated().any():
        raise ComparisonError("Comparison requires at least two unique sample IDs")
    paths = [Path(value).expanduser() for value in samples["result_dir"]]
    samples["result_dir"] = [str((p if p.is_absolute() else path.parent / p).resolve()) for p in paths]
    if samples["result_dir"].duplicated().any():
        raise ComparisonError("Multiple samples point to the same result directory")
    return samples.set_index("sample_id", drop=False)


def _profile(root: Path, feature: str) -> tuple[pd.Series, dict, dict, Path]:
    metadata = _json(root / "analysis_metadata.json")
    summary = _json(root / "analysis_summary.json")
    if metadata.get("status") != "completed":
        raise ComparisonError(f"Run is not completed: {root}")
    provenance = metadata.get("tool_provenance", {})
    config = metadata.get("effective_config", {})
    if not isinstance(provenance, dict) or not isinstance(config, dict):
        raise ComparisonError(f"Analysis provenance is incomplete for comparison: {root}")
    preprocessing = config.get("preprocessing")
    if not isinstance(preprocessing, dict) or "enabled" not in preprocessing:
        raise ComparisonError(f"Preprocessing provenance is incomplete: {root}")
    if feature == "taxonomy":
        path = root / "bracken_report.tsv"
        frame = _table(path, dtype={"taxonomy_id": str})
        if not {"taxonomy_id", "taxonomy_lvl", "name", "new_est_reads"}.issubset(frame.columns):
            raise ComparisonError(f"Bracken table lacks taxonomy IDs, names or counts: {path}")
        frame["taxonomy_id"] = frame["taxonomy_id"].str.strip()
        if not frame["taxonomy_id"].str.fullmatch(r"[1-9][0-9]*").all():
            raise ComparisonError(f"Profile has invalid taxonomy IDs: {path}")
        classification = config.get("classification", {})
        if not isinstance(classification, dict) or not classification.get("taxonomic_level"):
            raise ComparisonError(f"Classification provenance is incomplete: {root}")
        if not (frame["taxonomy_lvl"] == classification["taxonomic_level"]).all():
            raise ComparisonError(f"Profile rank differs from classification settings: {path}")
        identifiers, values = frame["taxonomy_id"], frame["new_est_reads"]
        labels = dict(zip(identifiers, frame["name"]))
        database = provenance.get("kraken_database", {}).get("manifest", {})
        database_identity = {"release": database.get("release"), "archive_md5": database.get("archive_md5")}
        if not all(database_identity.values()):
            raise ComparisonError(f"Taxonomy database provenance is incomplete: {root}")
        identity = {
            "database": database_identity, "classifier": provenance.get("kraken2"),
            "estimator": provenance.get("bracken"), "classification": classification,
            "bracken_model_read_length": provenance.get("bracken_model_read_length"),
            "count_unit": summary.get("taxonomy", {}).get("count_unit"),
        }
    else:
        abundance = summary.get("abundance", {})
        relative = abundance.get("functional_abundance")
        if not relative:
            raise ComparisonError(f"Functional abundance was not completed: {root}")
        path = (root / relative).resolve()
        if not path.is_relative_to(root.resolve()):
            raise ComparisonError("Functional table path leaves its result directory")
        frame = _table(path, dtype={"term": str})
        if not {"namespace", "term", "tpm"}.issubset(frame.columns):
            raise ComparisonError(f"Functional table lacks namespace, term or TPM: {path}")
        frame = frame.loc[frame["namespace"] == feature].copy()
        frame["term"] = frame["term"].str.strip()
        identifiers, values = frame["term"], frame["tpm"]
        labels = {key: key for key in identifiers}
        database_info = provenance.get("functional_database")
        if not isinstance(database_info, dict):
            raise ComparisonError(f"Functional database provenance is incomplete: {root}")
        database = database_info.get("manifest", {})
        artifacts = database.get("artifacts", [])
        if not database.get("release") or not artifacts or any(not item.get("sha256") for item in artifacts):
            raise ComparisonError(f"Functional database provenance is incomplete: {root}")
        identity = {
            "database_release": database["release"],
            "database_artifacts": sorted(item["sha256"] for item in artifacts),
            "mapper": abundance.get("tool_version"), "annotation": config.get("annotation"),
            "annotator": provenance.get("eggnog_mapper"),
            "gene_predictor": provenance.get("pyrodigal"), "assembler": provenance.get("megahit"),
            "min_identity": abundance.get("min_identity"), "ambiguous": abundance.get("ambiguous"),
            "term_attribution": abundance.get("functional", {}).get("term_attribution"),
            "count_unit": abundance.get("count_unit"),
        }
    if identifiers.isna().any() or identifiers.duplicated().any() or (identifiers.astype(str) == "").any():
        raise ComparisonError(f"Profile has missing or duplicate feature IDs: {path}")
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(numeric).all() or (numeric < 0).any():
        raise ComparisonError(f"Profile must contain finite nonnegative values and a positive total: {path}")
    try:
        total = math.fsum(numeric)
    except OverflowError as exc:
        raise ComparisonError(f"Profile total is not finite: {path}") from exc
    if total <= 0:
        raise ComparisonError(f"Profile must have a positive total: {path}")
    if feature == "taxonomy":
        if not (numeric == np.floor(numeric)).all():
            raise ComparisonError(f"Bracken estimated counts must be integers: {path}")
        taxonomy = summary.get("taxonomy", {})
        assigned, classified, input_total = (taxonomy.get(key) for key in
                                            ("species_assigned", "kraken_classified", "total_input"))
        counts = (assigned, classified, input_total)
        if (any(not isinstance(value, (int, float)) or not math.isfinite(value) for value in counts)
                or not 0 < total == assigned <= classified <= input_total):
            raise ComparisonError(f"Taxonomy counts disagree with summary: {path}")
    identity["preprocessor"] = provenance.get("fastp") if preprocessing["enabled"] else "not_run"
    if identity.get("count_unit") not in ("reads", "fragments"):
        raise ComparisonError(f"Count-unit provenance is incomplete: {root}")
    if any(value is None for value in identity.values()):
        raise ComparisonError(f"Analysis provenance is incomplete for comparison: {root}")
    for section in ("classification", "annotation"):
        if isinstance(identity.get(section), dict):
            identity[section] = {key: value for key, value in identity[section].items() if key != "threads"}
    identity["preprocessing"] = {key: value for key, value in preprocessing.items() if key != "threads"}
    series = pd.Series(numeric, index=identifiers.astype(str).to_numpy())
    return series, labels, {"identity": identity, "metadata": metadata, "summary": summary}, path


def group_tests(
    relative: pd.DataFrame, samples: pd.DataFrame, *, group_column: str,
    reference_group: str, pseudocount: float, fdr_method: str, alpha: float,
    min_prevalence: float,
) -> tuple[pd.DataFrame, dict]:
    from scipy.stats import false_discovery_control, mannwhitneyu

    if not 0 < pseudocount < 1 or not np.isfinite(pseudocount):
        raise ComparisonError("Formal CLR comparisons require an explicit proportion-scale pseudocount in (0, 1)")
    if fdr_method not in ("bh", "by") or not 0 < alpha < 1 or not 0 <= min_prevalence <= 1:
        raise ComparisonError("Invalid FDR method, alpha or prevalence threshold")
    if (samples.index.has_duplicates or relative.index.has_duplicates
            or set(relative.index) != set(samples.index)):
        raise ComparisonError("Profile sample IDs must match the unique sample-sheet IDs")
    if relative.columns.has_duplicates or relative.columns.isna().any() or (relative.columns == "").any():
        raise ComparisonError("Comparison requires unique nonempty feature IDs")
    values = relative.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if (not np.isfinite(values).all() or (values < 0).any() or (values > 1).any()
            or not np.allclose(values.sum(axis=1), 1, rtol=0, atol=1e-8)):
        raise ComparisonError("Formal tests require finite nonnegative proportions summing to one per sample")
    relative = pd.DataFrame(values, index=relative.index, columns=relative.columns).loc[samples.index]
    if group_column in samples and (samples[group_column].isna().any()
                                    or samples[group_column].astype(str).str.strip().eq("").any()):
        raise ComparisonError("Formal comparison requires nonempty sample group labels")
    if group_column not in samples or reference_group not in set(samples[group_column]):
        raise ComparisonError("Formal comparison requires a valid group column and reference group")
    groups = list(dict.fromkeys(samples[group_column]))
    if len(groups) != 2:
        raise ComparisonError("Formal comparison currently supports exactly two independent groups")
    if "subject_id" in samples and samples["subject_id"].astype("string").str.strip().replace("", pd.NA).dropna().duplicated().any():
        raise ComparisonError("Repeated subject IDs are incompatible with an independent-samples test")
    other = next(group for group in groups if group != reference_group)
    first = samples.index[samples[group_column] == other]
    second = samples.index[samples[group_column] == reference_group]
    if min(len(first), len(second)) < 3:
        raise ComparisonError("Formal comparison requires at least three independent samples per group")
    prevalence = (relative > 0).mean(axis=0)
    kept = relative.columns[(prevalence >= min_prevalence) & (prevalence > 0)]
    if len(kept) < 2:
        raise ComparisonError("CLR comparison requires at least two features after prevalence filtering")
    transformed = np.log(relative.loc[:, kept] + pseudocount)
    transformed = transformed.subtract(transformed.mean(axis=1), axis=0)
    rows = []
    for feature in kept:
        x = transformed.loc[first, feature].to_numpy()
        y = transformed.loc[second, feature].to_numpy()
        # With one pooled value every rank is tied: U is centered and there is
        # no evidence of a group difference. Some SciPy versions return NaN here.
        if np.all(x == x[0]) and np.all(y == x[0]):
            statistic, p_value = len(first) * len(second) / 2, 1.0
        else:
            if min(len(x), len(y)) < 10 and np.unique(np.concatenate((x, y))).size < len(x) + len(y):
                raise ComparisonError(
                    f"Feature {feature}: nonconstant tied CLR values in a small group require a permutation-based test; "
                    "this method requires at least 10 samples per group when ties occur"
                )
            result = mannwhitneyu(x, y, alternative="two-sided", method="auto", nan_policy="raise")
            statistic, p_value = float(result.statistic), float(result.pvalue)
        if not np.isfinite(p_value) or not 0 <= p_value <= 1:
            raise ComparisonError(f"Undefined Mann-Whitney p-value for feature {feature}; inspect the study profiles")
        rows.append({
            "feature_id": feature, "comparison_group": other, "reference_group": reference_group,
            "n_comparison": len(first), "n_reference": len(second),
            "prevalence": float(prevalence[feature]),
            "mean_relative_comparison": relative.loc[first, feature].mean(),
            "mean_relative_reference": relative.loc[second, feature].mean(),
            "median_clr_difference": float(np.median(x) - np.median(y)),
            "rank_biserial_effect": float(2 * statistic / (len(first) * len(second)) - 1),
            "u_statistic": float(statistic), "p_value": p_value,
        })
    results = pd.DataFrame(rows)
    results["q_value"] = false_discovery_control(results["p_value"].to_numpy(), method=fdr_method)
    results["passes_fdr"] = results["q_value"] <= alpha
    return results.sort_values(["q_value", "feature_id"]), {
        "method": "two-sided Mann-Whitney U on CLR abundances (SciPy auto exact/asymptotic)",
        "transformation": "natural log(proportion + pseudocount), centered across retained features per sample",
        "pseudocount": pseudocount, "pseudocount_scale": "proportion of full profile total",
        "fdr_method": fdr_method, "alpha": alpha, "min_prevalence": min_prevalence,
        "features_tested": len(kept), "features_excluded": relative.shape[1] - len(kept),
        "reference_group": reference_group, "comparison_group": other,
        "design": "independent unpaired samples; no covariate adjustment",
        "null": "equal distributions of feature CLR values across the two groups",
        "constant_features": "A single pooled CLR value is reported with U=n1*n2/2 and p=1",
        "small_tied_features": 0,
        "pvalue_policy": "SciPy auto: exact for an untied group of size <=8; otherwise tie-corrected asymptotic with continuity correction. Nonconstant tied CLR samples are rejected when either group has fewer than 10 observations; no permutation-based analysis is performed.",
    }


def run_comparison(
    sample_sheet: Path, output_dir: Path, *, feature: str = "taxonomy",
    group_column: str | None = None, test: str | None = None,
    reference_group: str | None = None, pseudocount: float | None = None,
    fdr_method: str = "by", alpha: float = 0.05, min_prevalence: float = 0.1,
) -> Path:
    if feature not in ("taxonomy", "COG", "KO", "EC", "GO"):
        raise ComparisonError("Feature must be taxonomy, COG, KO, EC or GO")
    samples = load_samples(sample_sheet, group_column)
    profiles, labels, qc, fingerprints = {}, {}, [], []
    identity = None
    seen_inputs = set()
    for sample_id, sample in samples.iterrows():
        series, names, details, source = _profile(Path(sample["result_dir"]), feature)
        if identity is None:
            identity = details["identity"]
        elif identity != details["identity"]:
            raise ComparisonError(f"Incompatible database, tool versions or analysis settings for sample {sample_id}")
        metadata, summary = details["metadata"], details["summary"]
        input_hashes = [item.get("sha256") or "" for item in metadata.get("input_files", [])]
        if test and (not input_hashes or any(not isinstance(value, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", value) for value in input_hashes)):
            raise ComparisonError(f"Complete SHA-256 input checksums required for formal comparison: {sample_id}")
        canonical_hashes = {value.lower() for value in input_hashes}
        if test and (len(canonical_hashes) != len(input_hashes) or canonical_hashes & seen_inputs):
            raise ComparisonError("Duplicate input data cannot be treated as independent biological replicates")
        seen_inputs.update(canonical_hashes)
        profiles[sample_id] = series
        if any(key in labels and labels[key] != value for key, value in names.items()):
            raise ComparisonError(f"Conflicting feature labels for sample {sample_id}")
        labels.update(names)
        taxonomy, abundance = summary.get("taxonomy", {}), summary.get("abundance", {})
        qc.append({
            "sample_id": sample_id, "group": sample[group_column] if group_column else "not provided",
            "features_reported": len(series), "profile_total": float(series.sum()),
            "total_input": taxonomy.get("total_input"), "classification_rate": taxonomy.get("classification_rate"),
            "assembly_n50": summary.get("assembly", {}).get("n50"),
            "mapping_rate": abundance.get("mapping_rate"), "gene_assignment_rate": abundance.get("gene_assignment_rate"),
        })
        fingerprints.append({"sample_id": sample_id, "source": str(source), "sha256": _sha256(source),
                             "metadata_sha256": _sha256(Path(sample["result_dir"]) / "analysis_metadata.json")})
    matrix = pd.DataFrame.from_dict(profiles, orient="index").fillna(0).sort_index(axis=1)
    matrix = matrix.loc[samples.index]
    relative = matrix.div(matrix.sum(axis=1), axis=0)
    statistics, method = None, {"method": "descriptive only"}
    if test:
        if test != "clr-mannwhitney" or not group_column or reference_group is None or pseudocount is None:
            raise ComparisonError("Formal comparisons require --test clr-mannwhitney, --group-column, --reference-group and --pseudocount")
        statistics, method = group_tests(relative, samples, group_column=group_column,
                                        reference_group=reference_group, pseudocount=pseudocount,
                                        fdr_method=fdr_method, alpha=alpha, min_prevalence=min_prevalence)
    output = Path(output_dir).resolve()
    if any(output == Path(path) or output in Path(path).parents for path in samples["result_dir"]):
        raise ComparisonError("Comparison output must not replace a sample result directory or its parent")
    if output.exists() and (not output.is_dir() or any(
            entry.name != "metaquest.log" or not entry.is_file() for entry in output.iterdir())):
        raise ComparisonError("Comparison output already contains results; choose a new output directory")
    output.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output / "feature_values.tsv", sep="\t", index_label="sample_id")
    relative.to_csv(output / "relative_abundance.tsv", sep="\t", index_label="sample_id")
    pd.DataFrame({"feature_id": matrix.columns, "label": [labels[key] for key in matrix.columns]}).to_csv(
        output / "features.tsv", sep="\t", index=False)
    qc_frame = pd.DataFrame(qc)
    qc_frame.to_csv(output / "sample_qc.tsv", sep="\t", index=False)
    samples.to_csv(output / "sample_sheet.tsv", sep="\t", index=False)
    if statistics is not None:
        statistics["label"] = statistics["feature_id"].map(labels)
        statistics.to_csv(output / "group_tests.tsv", sep="\t", index=False)
    from scipy import __version__ as scipy_version
    from scipy.spatial.distance import pdist, squareform

    distances = pd.DataFrame(squareform(pdist(relative.to_numpy(), metric="braycurtis")),
                             index=relative.index, columns=relative.index)
    distances.to_csv(output / "bray_curtis.tsv", sep="\t", index_label="sample_id")
    summary = {
        "analysis_type": "comparative_metagenomics", "research_use_only": True,
        "metaquest_version": __version__,
        "sample_count": len(samples), "feature_count": len(matrix.columns), "feature": feature,
        "input_unit": "Bracken estimated counts" if feature == "taxonomy" else "summed gene TPM",
        "relative_denominator": "sum of all reported profile features per sample; unclassified/unassembled data excluded",
        "group_column": group_column, "statistics": method, "scipy_version": scipy_version,
        "provenance": identity, "source_files": fingerprints,
        "sample_sheet_sha256": _sha256(Path(sample_sheet)),
        "profile_signature": hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest(),
        "limitations": [
            "Formal tests are exploratory distribution tests on log-ratios, not absolute abundance or causal effects.",
            "Samples must be independent biological replicates; technical replicates, pairing and covariates are not modeled.",
            "Zero replacement and prevalence filtering change the CLR reference; sensitivity analyses are required.",
            "BY controls FDR under arbitrary test dependence; BH requires independence or appropriate positive dependence.",
            "Separate assemblies recover different gene catalogs; functional differences may reflect assembly recovery.",
            "Full term attribution can count a gene in multiple terms; term totals are not disjoint DNA copy counts.",
        ],
    }
    from metaquest.reporting.comparison_report import write_comparison_report

    write_comparison_report(output, relative, labels, qc_frame, distances, statistics, summary)
    (output / "comparison_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return output / "comparison_summary.json"
