"""Run the frozen small synthetic checkpoint-4 studies; no reads/databases are used."""
import argparse
import csv
import hashlib
import inspect
import json
import math
from pathlib import Path
import platform
import time

import numpy as np
import pandas as pd
import scipy

from metaquest.core.comparison import group_tests


def wilson(successes, count):
    z = 1.959963984540054
    rate = successes / count
    denominator = 1 + z * z / count
    center = (rate + z * z / (2 * count)) / denominator
    radius = z * math.sqrt(rate * (1 - rate) / count + z * z / (4 * count**2)) / denominator
    return [max(0, center - radius), min(1, center + radius)]


def mean_bounds(values):
    values = np.asarray(values, dtype=float)
    radius = math.sqrt(math.log(20) / (2 * len(values)))
    return {"mean": float(values.mean()), "mc_standard_error": float(values.std(ddof=1) / math.sqrt(len(values))),
            "lower_one_sided_95": max(0, float(values.mean()) - radius),
            "upper_one_sided_95": min(1, float(values.mean()) + radius)}


def generate(protocol, scenario, replicate):
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([scenario["seed"], replicate])))
    count = scenario["n_reference"] + scenario["n_comparison"]
    width = protocol["features"]
    if scenario["name"] == "sparse_null":
        weights = np.exp(rng.uniform(-2, 2, (count, width)))
        dropout = rng.random((count, width)) < .65
        dropout[:, :2] = False
        weights[dropout] = 0
        weights /= weights.sum(axis=1, keepdims=True)
        values = np.array([rng.multinomial(1000, row) for row in weights], dtype=float) / 1000
    else:
        coordinates = rng.uniform(-1, 1, (count, width))
        for feature, shift in scenario.get("changed_features", {}).items():
            coordinates[scenario["n_reference"]:, int(feature[1:])] += shift
        weights = np.exp(coordinates - coordinates.max(axis=1, keepdims=True))
        weights /= weights.sum(axis=1, keepdims=True)
        values = (1 + width * protocol["parameters"]["pseudocount"]) * weights - protocol["parameters"]["pseudocount"]
        if (values <= 0).any():
            raise ValueError("Frozen inverse-closure design must remain positive")
    samples = pd.DataFrame({"group": ["control"] * scenario["n_reference"] + ["case"] * scenario["n_comparison"]},
                           index=[f"s{i:03d}" for i in range(count)])
    return pd.DataFrame(values, index=samples.index, columns=[f"f{i:02d}" for i in range(width)]), samples


def evaluate(values, samples, protocol, scenario, method, *, pseudocount=None, min_prevalence=None):
    settings = protocol["parameters"]
    result, info = group_tests(values, samples, group_column="group", reference_group="control",
                               pseudocount=settings["pseudocount"] if pseudocount is None else pseudocount,
                               min_prevalence=settings["min_prevalence"] if min_prevalence is None else min_prevalence,
                               alpha=settings["alpha"], fdr_method=method)
    if not result["p_value"].between(0, 1, inclusive="right").all():
        raise ValueError("Unexpected invalid/zero p-value for these small bounded studies")
    indexed = result.set_index("feature_id")
    detected = set(indexed.index[indexed["passes_fdr"]])
    truth = scenario.get("changed_features", {})
    return {"features_tested": info["features_tested"], "discoveries": len(detected),
            "false_discoveries": len(detected - truth.keys()),
            "fdp": len(detected - truth.keys()) / max(len(detected), 1),
            "target_recall": len(detected & truth.keys()) / len(truth) if truth else None,
            "both_target_directions_correct": all(indexed.loc[key, "median_clr_difference"] * value > 0
                                                   for key, value in truth.items()) if truth else None,
            "minimum_p": float(result["p_value"].min()), "minimum_q": float(result["q_value"].min()),
            "small_tied_features": info["small_tied_features"]}


def summarize(records, protocol):
    limits = protocol["acceptance_limits"]
    summary = []
    for scenario in protocol["scenarios"]:
        for method in protocol["parameters"]["fdr_methods"]:
            selected = [row for row in records if row["scenario"] == scenario["name"] and
                        row["method"] == method and row["setting"] == "primary"]
            if len(selected) != protocol["replicates_per_scenario"]:
                raise ValueError("Missing or duplicated primary study records")
            row = {"scenario": scenario["name"], "method": method, "replicates": len(selected)}
            if "changed_features" not in scenario:
                events = sum(item["discoveries"] > 0 for item in selected)
                interval = wilson(events, len(selected))
                row.update(global_null_discovery_studies=events, empirical_fdr=events / len(selected),
                           wilson_95=interval, passed=interval[1] <= limits["null_wilson_upper_max"])
            else:
                fdp = mean_bounds([item["fdp"] for item in selected])
                recall = mean_bounds([item["target_recall"] for item in selected])
                direction = wilson(sum(item["both_target_directions_correct"] for item in selected), len(selected))
                row.update(fdp=fdp, target_recall=recall, direction_wilson_95=direction,
                           passed=(fdp["upper_one_sided_95"] <= limits["effect_fdp_hoeffding_upper_max"] and
                                   recall["lower_one_sided_95"] >= limits["effect_recall_hoeffding_lower_min"] and
                                   direction[0] >= limits["effect_direction_wilson_lower_min"]))
            summary.append(row)
    return summary


def plots(summary, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.7))
    methods = [("by", "#197c8c"), ("bh", "#d47c24")]
    names = ["dense_null", "unbalanced_null", "sparse_null"]
    for index, (method, color) in enumerate(methods):
        for position, name in enumerate(names):
            row = next(item for item in summary if item["scenario"] == name and item["method"] == method)
            value = row["empirical_fdr"]
            axes[0].errorbar(position + (index - .5) * .16, value,
                            yerr=[[max(0, value - row["wilson_95"][0])], [row["wilson_95"][1] - value]],
                            fmt="o", color=color, capsize=3, label=method.upper() if position == 0 else None)
        row = next(item for item in summary if item["scenario"] == "balanced_clr_effect" and item["method"] == method)
        for axis, key, bound in ((axes[1], "fdp", "upper_one_sided_95"), (axes[2], "target_recall", "lower_one_sided_95")):
            metric = row[key]
            distance = abs(metric[bound] - metric["mean"])
            error = [[0], [distance]] if key == "fdp" else [[distance], [0]]
            axis.errorbar(index, metric["mean"], yerr=error, fmt="o", color=color, capsize=4)
    axes[0].set(xticks=range(3), xticklabels=["Dense", "Unequal n", "Sparse"], ylim=(0, .1),
                ylabel="Studies with any discovery", title="Global null: 95% Wilson intervals")
    axes[0].axhline(.05, color="#666666", linestyle="--", linewidth=1)
    axes[0].legend(frameon=False)
    for axis in axes[1:]:
        axis.set_xticks([0, 1], ["BY", "BH"])
    axes[1].set(ylim=(0, .12), ylabel="Mean false-discovery proportion", title="Known CLR effects: upper 95% bound")
    axes[1].axhline(.05, color="#666666", linestyle="--", linewidth=1)
    axes[2].set(ylim=(.7, 1.02), ylabel="Mean target recall", title="Known CLR effects: lower 95% bound")
    axes[2].axhline(.8, color="#666666", linestyle="--", linewidth=1)
    fig.suptitle("Checkpoint 4: frozen synthetic studies", fontsize=13)
    fig.text(.5, .015, "500 independent studies/scenario; effect bounds are one-sided Hoeffding bounds. No biological-accuracy claim.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .05, 1, .94))
    for extension in ("png", "svg", "pdf"):
        fig.savefig(output / f"comparison-audit.{extension}", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    if protocol["protocol_version"] != 1 or protocol["replicates_per_scenario"] < 2:
        raise ValueError("Unsupported simulation protocol")
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    records = []
    for scenario in protocol["scenarios"]:
        print(f'Running {scenario["name"]}: {protocol["replicates_per_scenario"]} studies', flush=True)
        for replicate in range(protocol["replicates_per_scenario"]):
            values, samples = generate(protocol, scenario, replicate)
            digest = hashlib.sha256(np.asarray(values, dtype="<f8").tobytes()).hexdigest()
            settings = [("primary", {})]
            if scenario["name"] == "balanced_clr_effect":
                settings += [("pseudocount_1e-4", {"pseudocount": 1e-4}), ("prevalence_0.5", {"min_prevalence": .5})]
            if scenario["name"] == "sparse_null" and replicate < 50:
                settings += [("pseudocount_1e-4", {"pseudocount": 1e-4})]
            for setting, options in settings:
                for method in protocol["parameters"]["fdr_methods"]:
                    result = evaluate(values, samples, protocol, scenario, method, **options)
                    # A changed pseudocount changes the hypotheses; sensitivity truth is descriptive.
                    if setting != "primary":
                        result["false_discoveries"] = result["fdp"] = None
                    records.append({"scenario": scenario["name"], "replicate": replicate, "method": method,
                                    "setting": setting, "input_sha256": digest, **result})
            if (replicate + 1) % 100 == 0:
                print(f'  {replicate + 1} completed', flush=True)
    with (args.output / "study-metrics.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=records[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    summary = summarize(records, protocol)
    implementation = Path(inspect.getsourcefile(group_tests))
    result = {"protocol_sha256": hashlib.sha256(args.protocol.read_bytes()).hexdigest(),
              "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
              "implementation": str(implementation), "implementation_sha256": hashlib.sha256(implementation.read_bytes()).hexdigest(),
              "elapsed_seconds": time.monotonic() - started, "criteria_passed": all(row["passed"] for row in summary),
              "primary": summary, "sensitivity_rows": sum(row["setting"] != "primary" for row in records),
              "limits": protocol["limits"]}
    (args.output / "summary.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    plots(summary, args.output)
    print(json.dumps(result, indent=2), flush=True)
    return 0 if result["criteria_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
