"""Score a completed paired Zymo run independently of MetaQuest's report parsers.

No reads or databases are downloaded. Nominal DNA abundance is not read truth.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_counts(path):
    counts = {}
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require({"taxonomy_id", "taxonomy_lvl", "name", "new_est_reads"} <= set(reader.fieldnames or []),
                "Missing Bracken columns")
        for row in reader:
            taxid, count = int(row["taxonomy_id"]), int(row["new_est_reads"])
            require(taxid > 0 and count >= 0, "Invalid Bracken taxid/count")
            require(row["taxonomy_lvl"] == "S", "Benchmark requires species-level Bracken output")
            require(taxid not in counts, "Duplicate Bracken taxid")
            counts[taxid] = {"name": row["name"], "count": count}
    return counts


def kraken_counts(path):
    values = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            require(len(fields) >= 6, "Malformed Kraken report")
            rank = fields[3].strip()
            if rank in ("U", "R"):
                require(rank not in values, "Duplicate Kraken root/unclassified row")
                values[rank] = int(fields[1])
                require(values[rank] >= 0, "Negative Kraken count")
    require("R" in values, "Missing Kraken root row")
    return values["R"], values.get("U", 0)


def score(counts, total, members, threshold):
    require(total > 0, "Empty classification denominator")
    require(math.isfinite(threshold) and 0 < threshold <= 1, "Invalid detection threshold")
    require(members and any(member["bacterial"] for member in members), "Missing bacterial truth")
    taxids = [taxid for member in members for taxid in member["taxids"]]
    require(len(taxids) == len(set(taxids)), "Overlapping truth taxids")
    require(all(member["taxids"] for member in members), "Empty truth taxids")
    require(sum(row["count"] for row in counts.values()) <= total, "Bracken counts exceed input")
    known = set(taxids)
    expected = []
    for member in members:
        count = sum(counts.get(taxid, {"count": 0})["count"] for taxid in member["taxids"])
        expected.append({**member, "count": count, "fraction_of_input": count / total,
                         "detected": count / total >= threshold})
    unexpected = [
        {"taxid": taxid, **row, "fraction_of_input": row["count"] / total,
         "detected": row["count"] / total >= threshold}
        for taxid, row in sorted(counts.items()) if taxid not in known
    ]
    tp = sum(row["detected"] for row in expected)
    fp = sum(row["detected"] for row in unexpected)
    bacterial = [row for row in expected if row["bacterial"]]
    return {
        "threshold": threshold,
        "community_true_positives": tp,
        "unexpected_reported_species": fp,
        "community_false_negatives": len(expected) - tp,
        "community_precision": tp / (tp + fp) if tp + fp else None,
        "community_recall": tp / len(expected),
        "bacterial_recall": sum(row["detected"] for row in bacterial) / len(bacterial),
        "unexpected_fraction_of_input_all_abundances": sum(row["count"] for row in unexpected) / total,
        "expected": expected,
        "unexpected": unexpected,
    }


def evaluate(directory, protocol_path):
    protocol = json.loads(protocol_path.read_text())
    require(protocol["protocol_version"] == 1, "Unsupported benchmark protocol")
    metadata_path = directory / "analysis_metadata.json"
    summary_path = directory / "analysis_summary.json"
    metadata = json.loads(metadata_path.read_text())
    summary = json.loads(summary_path.read_text())
    require(metadata["status"] == "completed", "Run is not completed")
    require(metadata["workflow"]["read_mode"] == "paired", "Benchmark requires paired fragments")
    require(not metadata["validation"].get("bypassed", False), "Input validation was bypassed")
    require(metadata["validation"]["status"] == "passed", "Input validation did not pass")
    reports = metadata["validation"]["reports"]
    require(len(reports) == 2, "Expected two validated mates")
    for report, source in zip(reports, protocol["inputs"]):
        require(report["md5_checksum"] == source["md5"], "Input MD5 does not match frozen dataset")
        analyzed = report["total_sequences"]
        require(0 < analyzed <= protocol["raw_read_pairs"], "Invalid analyzed read count")
        require(report.get("sampling_limited", False) or analyzed == protocol["raw_read_pairs"],
                "Unsampled input read count differs")
    database = metadata["tool_provenance"]["kraken_database"]["manifest"]
    require(database["release"] == protocol["taxonomy_database"]["release"], "Database release differs")
    require(protocol["taxonomy_database"]["collection"] in database["name"], "Database collection differs")
    bracken_path, kraken_path = directory / "bracken_report.tsv", directory / "kraken_report.txt"
    counts = read_counts(bracken_path)
    classified, unclassified = kraken_counts(kraken_path)
    total = classified + unclassified
    require(total > 0, "Empty Kraken input")
    require(sum(row["count"] for row in counts.values()) <= classified, "Bracken counts exceed classified fragments")
    require(summary["taxonomy"]["total_input"] == total, "Summary/input counts disagree")
    require(summary["taxonomy"]["kraken_classified"] == classified, "Classified counts disagree")
    require(summary["taxonomy"]["unclassified_reads"] == unclassified, "Unclassified counts disagree")
    require(summary["taxonomy"]["species_assigned"] == sum(row["count"] for row in counts.values()),
            "Summary/Bracken counts disagree")
    fastp_path = directory / "preprocessing" / "fastp.json"
    fastp = json.loads(fastp_path.read_text())
    require(fastp["summary"]["before_filtering"]["total_reads"] == protocol["raw_read_pairs"] * 2,
            "fastp raw read count does not match frozen dataset")
    require(fastp["summary"]["after_filtering"]["total_reads"] == total * 2,
            "Kraken fragment denominator does not match fastp surviving mates")
    metrics = [score(counts, total, protocol["expected_members"], threshold)
               for threshold in [protocol["primary_threshold"], *protocol["secondary_thresholds"]]]
    primary = metrics[0]
    passed = (primary["bacterial_recall"] >= protocol["minimum_bacterial_recall"]
              and primary["community_precision"] is not None
              and primary["community_precision"] >= protocol["minimum_community_precision"])
    unassigned = 1 - sum(row["count"] for row in counts.values()) / total
    nominal_l1 = (sum(abs(row["fraction_of_input"] - row["dna_fraction"])
                      for row in primary["expected"])
                  + primary["unexpected_fraction_of_input_all_abundances"] + unassigned)
    return {
        "dataset": protocol["dataset"], "role": protocol["role"],
        "protocol_sha256": sha256(protocol_path),
        "artifacts_sha256": {path.name: sha256(path) for path in
                             (metadata_path, summary_path, bracken_path, kraken_path, fastp_path)},
        "classification_input_fragments": total,
        "kraken_classified_fragments": classified, "unclassified_fragments": unclassified,
        "unclassified_fraction": unclassified / total,
        "bracken_assigned_fragments": sum(row["count"] for row in counts.values()),
        "nominal_dna_l1_distance_not_read_accuracy": nominal_l1,
        "bracken_unassigned_fraction": unassigned,
        "metrics": metrics, "bacterial_baseline_gate_passed": passed,
        "whole_checkpoint_passed": False,
        "limits": [protocol["abundance_truth"], protocol["interpretation"],
                   "The two yeasts are outside the Standard-8 reference scope.",
                   "The Cryptococcus hybrid is one grouped community member, not two independent species.",
                   "No reference-defined negative universe: specificity is not estimated.",
                   "Unexpected calls are relative to nominal composition; contamination or reference ambiguity can contribute.",
                   "Single development sample: no biological replicate uncertainty estimate."],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--protocol", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = evaluate(args.run, args.protocol)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        print(f"Benchmark rejected: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["bacterial_baseline_gate_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
