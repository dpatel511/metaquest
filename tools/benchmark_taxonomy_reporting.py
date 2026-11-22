"""Repeat report generation on identical cached taxonomy outputs (no classifier run).

Run once with the frozen baseline source on PYTHONPATH and once with the candidate
installed wheel. Output must be a fresh directory. Timers do not establish accuracy.
"""
import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import resource
import shutil
import sys
import time

import yaml


def benchmark(source, output):
    import metaquest
    from metaquest.config import __version__
    from metaquest.core.preprocessing import _summary_metrics
    from metaquest.pipeline.context import ClassificationResult, PipelineContext
    from metaquest.reporting.stable_reporter import generate_stable_reports
    from metaquest.settings import load_config

    if output.exists():
        raise ValueError("Reporting benchmark requires a fresh output directory")
    metadata = json.loads((source / "analysis_metadata.json").read_text())
    if metadata.get("status") != "completed" or not metadata["workflow"]["taxonomy_only"]:
        raise ValueError("Supply a completed taxonomy-only run")
    if importlib.metadata.version("metaquest-bio") != __version__:
        raise ValueError("Installed package metadata and imported source version differ")
    output.mkdir(parents=True)
    config_path = output / "effective_config.yaml"
    config_path.write_text(yaml.safe_dump(metadata["effective_config"]))
    config = load_config(config_path=config_path)
    counts = {}
    with (source / "kraken_report.txt").open() as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if fields[3] in ("U", "R"):
                counts[fields[3]] = int(fields[1])
    with (source / "bracken_report.tsv").open() as handle:
        taxa = list(csv.DictReader(handle, delimiter="\t"))
    model = json.loads((source / "analysis_summary.json").read_text())["taxonomy"]
    timings = []
    for repeat in range(3):
        directory = output / f"repeat-{repeat + 1}"
        directory.mkdir()
        for name in ("analysis_metadata.json", "kraken_report.txt", "bracken_report.tsv"):
            shutil.copyfile(source / name, directory / name)
        (directory / "preprocessing").mkdir()
        shutil.copyfile(source / "preprocessing/fastp.json", directory / "preprocessing/fastp.json")
        ctx = PipelineContext(config, [], directory, read_mode=metadata["workflow"]["read_mode"], skip_annotation=True)
        ctx.preprocessing = _summary_metrics(json.loads((source / "preprocessing/fastp.json").read_text()))
        ctx.metadata["validation"] = metadata["validation"]
        ctx.classification = ClassificationResult(directory / "kraken_report.txt", directory / "bracken_report.tsv",
            species_count=len(taxa), species_assigned_reads=sum(int(row["new_est_reads"]) for row in taxa),
            total_reads=counts.get("U", 0) + counts.get("R", 0), kraken_classified_reads=counts.get("R", 0),
            unclassified_reads=counts.get("U", 0), observed_read_length=model["observed_read_length"],
            bracken_read_length=model["bracken_read_length"])
        wall, cpu = time.perf_counter_ns(), time.process_time_ns()
        generate_stable_reports(ctx)
        timings.append({"elapsed_seconds": (time.perf_counter_ns() - wall) / 1e9,
                        "cpu_seconds": (time.process_time_ns() - cpu) / 1e9,
                        "process_lifetime_max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    result = {"version": __version__, "module": metaquest.__file__, "python": sys.executable,
              "timings": timings, "source_sha256": {name: hashlib.sha256((source / name).read_bytes()).hexdigest()
                for name in ("kraken_report.txt", "bracken_report.tsv", "preprocessing/fastp.json")},
              "scope": "Three reporting repetitions on one cached run; shared-host timings, no speedup or accuracy claim"}
    (output / "measurements.json").write_text(json.dumps(result, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(benchmark(args.run.resolve(), args.output.resolve()), indent=2))


if __name__ == "__main__":
    main()
