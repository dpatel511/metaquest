"""Summarize one QUAST assembly against the checksum-pinned Zymo isolate drafts.

QUAST 5.3.0 counts aligned non-N bases from 1-based inclusive coordinates. This
tool aggregates its counted-base output; it does not convert genomic intervals.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import sys

from Bio import SeqIO


MEMBERS = {"bs": "B. spizizenii [paper: B. subtilis]", "ef": "E. faecalis", "ec": "E. coli",
           "lf": "L. fermentum", "lm": "L. monocytogenes", "pa": "P. aeruginosa",
           "se": "S. enterica", "sa": "S. aureus", "sc": "S. cerevisiae",
           "cn": "Cryptococcus hybrid"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fasta_records(path):
    records = {}
    normalization = str.maketrans("MKRYWSVBHD", "NNNNNNNNNN")
    for record in SeqIO.parse(path, "fasta"):
        require(record.id not in records, "Duplicate FASTA reference ID")
        # Match QUAST's documented uppercase / IUPAC-to-N normalization.
        sequence = str(record.seq).upper().translate(normalization)
        require(set(sequence) <= set("ACGTN"), "Invalid reference sequence alphabet")
        records[record.id] = (len(sequence), len(sequence) - sequence.count("N"),
                              hashlib.sha256(sequence.encode()).hexdigest())
    require(records, "Empty reference FASTA")
    return records


def aggregate(directory, reference):
    expected = fasta_records(reference)
    corrected = directory / "quast_corrected_input" / reference.name
    require(corrected.is_file(), "Missing retained reference; run QUAST with --debug")
    require(fasta_records(corrected) == expected, "QUAST reference content differs from supplied FASTA")
    with (directory / "report.tsv").open(encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    require(rows and rows[0][0] == "Assembly" and len(rows[0]) == 2,
            "Exactly one evaluated assembly required; maxima across assemblies cannot be combined")
    report = {row[0]: row[1] for row in rows if len(row) == 2}
    values = {key: {"member": label, "reference_non_N_bases": 0, "covered_non_N_bases": 0,
                    "reference_fragments": 0} for key, label in MEMBERS.items()}
    pattern = re.compile(r"^\s+(\S+) \(total length: (\d+) bp, total length without N's: (\d+) bp, maximal covered length: (\d+) bp\)$")
    seen = set()
    with (directory / "genome_stats/genome_info.txt").open(encoding="utf-8") as handle:
        for line in handle:
            match = pattern.fullmatch(line.rstrip("\r\n"))
            if not match:
                continue
            name, total, non_n, covered = match.groups()
            prefix = name.split("_", 1)[0]
            require(prefix in values, f"Unknown reference member: {name}")
            require(name not in seen, "Duplicate reference fragment")
            seen.add(name)
            require(0 <= int(covered) <= int(non_n) <= int(total), "Invalid reference coverage counts")
            require(name in expected and (int(total), int(non_n)) == expected[name][:2],
                    "QUAST reference statistics differ from supplied FASTA")
            values[prefix]["reference_non_N_bases"] += int(non_n)
            values[prefix]["covered_non_N_bases"] += int(covered)
            values[prefix]["reference_fragments"] += 1
    require(seen == expected.keys(), "Missing reference fragments in QUAST coverage output")
    for row in values.values():
        require(row["reference_non_N_bases"] > 0, f'Missing reference member: {row["member"]}')
        row["coverage_percent"] = row["covered_non_N_bases"] / row["reference_non_N_bases"] * 100
    total = sum(row["reference_non_N_bases"] for row in values.values())
    covered = sum(row["covered_non_N_bases"] for row in values.values())
    reported = report["Genome fraction (%)"]
    precision = len(reported.partition(".")[2])
    require(abs(covered / total * 100 - float(reported)) <= .5 * 10**(-precision) + 1e-8,
            "Per-member coverage disagrees with QUAST's aggregate fraction")
    return {"quast_report": report, "per_member": values,
            "total_reference_non_N_bases": total, "total_covered_non_N_bases": covered,
            "whole_checkpoint_passed": False, "reference_content_verified": True,
            "interpretation": "Agreement with filtered isolate draft assemblies; not finished-genome error truth. One assembly only; per-fragment maxima equal that assembly's coverage."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("quast", type=Path)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    args = parser.parse_args()
    try:
        protocol = json.loads(args.protocol.read_text())
        digest = hashlib.md5()
        with args.reference.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        require(digest.hexdigest() == protocol["reference"]["md5"], "Reference MD5 differs from protocol")
        result = aggregate(args.quast, args.reference)
        result["reference_md5"] = digest.hexdigest()
        result["protocol_sha256"] = hashlib.sha256(args.protocol.read_bytes()).hexdigest()
    except (ValueError, KeyError, TypeError, OSError) as exc:
        print(f"Reference assessment rejected: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
