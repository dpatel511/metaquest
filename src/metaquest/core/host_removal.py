"""Stream Bowtie2 SAM through vetted Samtools conversion/counting tools.

SAM text positions are 1-based; filtering uses flags, not genomic intervals.
"""
from contextlib import ExitStack
import gzip
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import re
import subprocess
import tempfile

from Bio.SeqIO.QualityIO import FastqGeneralIterator

from metaquest.exceptions import ConfigError, MetaQuestError
from metaquest.config import __version__
from metaquest.io.processes import managed_process


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def index_files(prefix):
    prefix = Path(prefix).expanduser().resolve()
    groups = [[Path(f"{prefix}.{part}.{kind}") for part in ("1", "2", "3", "4", "rev.1", "rev.2")]
              for kind in ("bt2", "bt2l")]
    present = [group for group in groups if any(p.exists() for p in group)]
    if len(present) != 1 or not all(p.is_file() and p.stat().st_size for p in present[0]):
        raise ConfigError(f"Host index must have one complete, nonempty .bt2 or .bt2l set: {prefix}. Use its prefix, not a FASTA or individual index file.")
    return present[0]


def index_identity(prefix):
    files = index_files(prefix)
    prefix = Path(prefix).expanduser().resolve()
    return {"prefix": str(prefix), "sha256": {str(p)[len(str(prefix)) + 1:]: _sha256(p) for p in files}}


def tool_versions():
    versions = {}
    for name, pattern in (("bowtie2", r"\bversion\s+(\S+)"), ("samtools", r"^samtools\s+(\S+)")):
        try:
            result = subprocess.run([name, "--version"], capture_output=True, text=True, timeout=15)
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise ConfigError(f"Host removal requires executable {name}") from exc
        match = re.search(pattern, result.stdout + result.stderr)
        if result.returncode or match is None:
            raise ConfigError(f"Cannot determine host-filter {name} version")
        versions[name] = match.group(1)
        release = re.match(r"^(\d+)\.(\d+)(?:\.|$)", versions[name])
        bounds = ((2, 5), (3, 0)) if name == "bowtie2" else ((1, 18), (2, 0))
        if release is None or not bounds[0] <= tuple(map(int, release.groups())) < bounds[1]:
            raise ConfigError(f"Unsupported host-filter {name} version: {versions[name]}")
    return versions


def count_fastq(paths, paired):
    if len(paths) != (2 if paired else 1):
        raise MetaQuestError("Host-filter FASTQ count requires one single file or two paired files")
    with ExitStack() as stack:
        handles = [stack.enter_context((gzip.open if str(path).endswith(".gz") else open)(path, "rt")) for path in paths]
        count = 0
        for records in zip_longest(*(FastqGeneralIterator(handle) for handle in handles)):
            if any(row is None for row in records):
                raise MetaQuestError("Host-filter paired FASTQ counts differ")
            if paired:
                names = [row[0].split()[0].removesuffix("/1").removesuffix("/2") for row in records]
                if names[0] != names[1]:
                    raise MetaQuestError("Host-filter paired FASTQ identifiers differ")
            count += 1
    return count


def run_host_removal(reads, output_dir, index, *, paired, threads, reuse=False):
    """Retain only primary unmapped reads, or pairs with both mates unmapped."""
    reads = [Path(path).resolve() for path in reads]
    if not all(path.is_file() for path in reads):
        raise ConfigError("Host-filter inputs must be regular FASTQ files")
    if threads < 1:
        raise ConfigError("Host-filter threads must be >= 1")
    reference = index_identity(index)
    root = Path(output_dir).resolve()
    output = root / "preprocessing" / "host_removal"
    if output.is_symlink() or not output.resolve().is_relative_to(root):
        raise ConfigError("Host-filter output path leaves the analysis directory")
    output.mkdir(parents=True, exist_ok=True)
    names = [f"retained_R{number}.fastq.gz" for number in range(1, 3 if paired else 2)]
    targets = [output / name for name in names]
    if any(path.exists() and not path.is_file() for path in targets):
        raise ConfigError("Host-filter output must be a regular file")
    state = {"inputs": [{"path": str(p), "sha256": _sha256(p)} for p in reads],
             "index": reference, "tool_versions": tool_versions(), "paired": paired,
             "threads": threads, "preset": "sensitive", "seed": 0,
             "pair_policy": "remove_pair_if_either_maps"}
    completion = output / "summary.json"
    if reuse and completion.is_file() and all(p.is_file() for p in targets):
        try:
            cached = json.loads(completion.read_text())
        except ValueError:
            cached = {}
        if (cached.get("state") == state and cached.get("status") == "completed"
                and cached.get("output_sha256") == {p.name: _sha256(p) for p in targets}
                and cached.get("evidence_sha256")
                and all((output / name).resolve().is_relative_to(output)
                        and (output / name).is_file() and _sha256(output / name) == digest
                        for name, digest in cached["evidence_sha256"].items())):
            return targets, {**cached, "reused": True}
    before = count_fastq(reads, paired)
    attempt = Path(tempfile.mkdtemp(prefix="attempt-", dir=output))
    if completion.is_file():
        (attempt / "previous-summary.json").write_bytes(completion.read_bytes())
        completion.unlink()
    staged = [attempt / name for name in names]
    mapper = ["bowtie2", "-x", reference["prefix"], "--sensitive", "-k", "1", "--mm", "--reorder", "--seed", "0", "-p", str(threads)]
    mapper += ["-1", str(reads[0]), "-2", str(reads[1])] if paired else ["-U", str(reads[0])]
    converter = ["samtools", "fastq", "-@", "0", "-c", "6", "-F", "2304", "-f", "12" if paired else "4"]
    converter += (["-N", "-1", str(staged[0]), "-2", str(staged[1]), "-0", str(attempt / "other.fastq"), "-s", str(attempt / "singletons.fastq")]
                  if paired else ["-n", "-0", str(staged[0])])
    converter.append("-")
    counter = ["samtools", "flagstat", "-@", "0", "-O", "json", "-"]
    try:
        with ExitStack() as stack:
            logs = [stack.enter_context((attempt / name).open("wb")) for name in ("bowtie2.log", "fastq.log", "flagstat.log", "flagstat.json")]
            align = stack.enter_context(managed_process(mapper, stdout=subprocess.PIPE, stderr=logs[0]))
            fastq = stack.enter_context(managed_process(converter, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=logs[1], bufsize=0))
            stats = stack.enter_context(managed_process(counter, stdin=subprocess.PIPE, stdout=logs[3], stderr=logs[2], bufsize=0))
            try:
                for block in iter(lambda: align.stdout.read(65536), b""):
                    for pipe in (fastq.stdin, stats.stdin):
                        pending = memoryview(block)
                        while pending:
                            pending = pending[pipe.write(pending):]
            finally:
                fastq.stdin.close()
                stats.stdin.close()
            codes = {"bowtie2": align.wait(), "fastq": fastq.wait(), "flagstat": stats.wait()}
        if any(codes.values()):
            raise MetaQuestError(f"Host-filter commands failed: {codes}; see {attempt}")
    except OSError as exc:
        raise MetaQuestError(f"Host-filter streaming failed; see {attempt}: {exc}") from exc
    metrics = json.loads((attempt / "flagstat.json").read_text())
    primary = sum(part["primary"] for part in metrics.values())
    mapped = sum(part["primary mapped"] for part in metrics.values())
    singletons = sum(part["singletons"] for part in metrics.values()) if paired else 0
    units = 2 if paired else 1
    retained_reads = primary - mapped - singletons
    retained = count_fastq(staged, paired)
    paired_records = sum(part["paired in sequencing"] for part in metrics.values())
    if (primary != before * units or paired_records != (primary if paired else 0)
            or retained_reads < 0 or retained_reads != retained * units
            or any((attempt / name).is_file() and (attempt / name).stat().st_size for name in ("other.fastq", "singletons.fastq"))):
        raise MetaQuestError(f"Host-filter FASTQ count disagrees with primary mapping counts; check duplicate identifiers or malformed alignments. See {attempt}")
    for source, target in zip(staged, targets):
        source.replace(target)
    summary = {"schema_version": 1, "metaquest_version": __version__, "status": "completed", "state": state,
               "count_unit": "fragments" if paired else "reads", "input": before,
               "retained": retained, "removed": before - retained,
               "retained_fraction": retained / before if before else None,
               "output_files": [str(p.relative_to(root)) for p in targets],
               "output_sha256": {p.name: _sha256(p) for p in targets},
               "commands": [mapper, converter, counter], "attempt_dir": str(attempt.relative_to(root)),
               "evidence_sha256": {str(path.relative_to(output)): _sha256(path) for path in
                                   (attempt / name for name in ("bowtie2.log", "fastq.log", "flagstat.log", "flagstat.json"))},
               "methods_file": "preprocessing/host_removal/summary.json", "reused": False,
               "interpretation": ["Matching is defined by the supplied index and Bowtie2 settings, without a MAPQ cutoff.",
                                  "A pair is removed when either mate aligns; retained paired files stay synchronized.",
                                  "Shared microbial/reference sequence can be removed; retained reads do not establish host absence.",
                                  "Index hashes identify the supplied index files, not an original FASTA or biological origin.",
                                  "SAM conversion normalizes mate suffixes and omits FASTQ header comments."]}
    temporary = attempt / "summary.json"
    temporary.write_text(json.dumps(summary, indent=2) + "\n")
    temporary.replace(completion)
    return targets, summary
