"""Exact two-point global-null calculations, including explicitly unsupported designs."""
import argparse
import hashlib
import inspect
import json
import math
from pathlib import Path

import pandas as pd

from metaquest.core.comparison import ComparisonError, group_tests


def calculate(protocol):
    settings = protocol["settings"]
    rows = []
    state_count = blocked_states = 0
    for nx, ny in protocol["n_comparison_reference"]:
        states = []
        for kx in range(nx + 1):
            for ky in range(ny + 1):
                values = pd.DataFrame({"f0": [.9] * kx + [.1] * (nx - kx) + [.9] * ky + [.1] * (ny - ky)})
                values["f1"] = 1 - values["f0"]
                samples = pd.DataFrame({"group": ["case"] * nx + ["control"] * ny})
                outcomes = {}
                for method in settings["methods"]:
                    try:
                        result, _ = group_tests(values, samples, group_column="group", reference_group="control",
                                                pseudocount=settings["pseudocount"], fdr_method=method,
                                                alpha=settings["alpha"], min_prevalence=settings["min_prevalence"])
                        outcomes[method] = bool(result["passes_fdr"].any())
                    except ComparisonError as exc:
                        if "permutation-based" not in str(exc) or min(nx, ny) >= 10 or kx + ky in (0, nx + ny):
                            raise
                        outcomes[method] = None
                        blocked_states += 1
                states.append((kx, ky, outcomes))
                state_count += 1
        for probability in protocol["population_probabilities"]:
            for method in settings["methods"]:
                weights = [math.comb(nx, kx) * math.comb(ny, ky) * probability**(kx + ky)
                           * (1 - probability)**(nx + ny - kx - ky) for kx, ky, _ in states]
                if abs(math.fsum(weights) - 1) > 1e-12:
                    raise ValueError("Binomial probability mass does not sum to one")
                unsupported = math.fsum(weight for weight, (_, _, outcome) in zip(weights, states) if outcome[method] is None)
                fdr = None if unsupported else math.fsum(weight for weight, (_, _, outcome) in zip(weights, states) if outcome[method])
                rows.append({"n_comparison": nx, "n_reference": ny, "probability": probability, "method": method,
                             "unsupported_probability_mass": unsupported, "exact_global_null_fdr": fdr,
                             "passes_nominal_alpha": None if fdr is None else fdr <= settings["alpha"] + 1e-12})
    return {"count_pair_states": state_count, "blocked_state_method_evaluations": blocked_states,
            "population_evaluations": rows,
            "fully_supported_calibration_passed": all(row["passes_nominal_alpha"] for row in rows if row["passes_nominal_alpha"] is not None)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Choose a new output JSON path")
    result = calculate(json.loads(args.protocol.read_text()))
    result["protocol_sha256"] = hashlib.sha256(args.protocol.read_bytes()).hexdigest()
    source = Path(inspect.getsourcefile(group_tests))
    result["implementation"] = str(source)
    result["implementation_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    result["interpretation"] = "Unsupported probability mass is not counted as zero FDR. Exact population calculations for a two-state family are not universal calibration."
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "population_evaluations"}, indent=2))
    return 0 if result["fully_supported_calibration_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
