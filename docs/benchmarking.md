# Biological benchmark protocol

Checkpoint 3 is **in progress**. Runtime tests, a development-data baseline and
held-out biological validation are separate kinds of evidence.
The [recorded first baseline](../benchmarks/zymo-even/RESULTS.md) retains the failed
primary taxonomy gate and unresolved structural/timing discrepancies.

## First baseline: ERR2984773

Use the complete paired Illumina Even mock from
[Nicholls et al. (2019)](https://doi.org/10.1093/gigascience/giz043), study
`PRJEB29504`, sample `SAMEA5157164`. The reads contain 8,785,731 pairs and total
1,956,674,499 compressed bytes. This sample was already used during MetaQuest
development; it is **not held out**. The authors release their data under
CC BY 4.0. Attribute the paper when reusing them.

The frozen [protocol](../benchmarks/zymo-even/protocol.json) contains input MD5s,
the reference URL/checksum, truth identities, thresholds and interpretation limits.
Record its SHA-256 before processing; retain that snapshot with every result.
Never edit thresholds after seeing results and present the run as prospective.
Any subsequent protocol revision needs a new recorded identity and rationale.

The current baseline tests the default Standard-8 collection (2026-06-26).
[Its documented contents](https://benlangmead.github.io/aws-indexes/k2) exclude
fungi. Keep both yeasts in the truth table and report their detection separately;
failure to detect them exposes database scope, not evidence of sample absence.

Taxonomy identities require care:

- The paper's *Bacillus subtilis* ATCC 6633 is now
  [*B. spizizenii*](https://www.atcc.org/products/6633), NCBI taxid 96241.
  Use that current species identity. Do not silently merge other *Bacillus* species.
- *Lactobacillus fermentum* is now *Limosilactobacillus fermentum*, taxid 1613.
- The paper describes a *Cryptococcus neoformans × C. deneoformans* hybrid.
  Group taxids 5207/40410 as **one community member** for whole-community reporting.
  These metrics are therefore member-level, not a strict ten-species score.

## Scoring rules frozen before this baseline

Primary detection threshold: estimated member count divided by **all surviving
paired fragments**, at least 0.001 (0.1%). Secondary descriptive thresholds are
0.0001 and 0.01. Bracken's rounded fraction column and classified-only fractions
are not the denominator. Sum the two hybrid parents before testing its threshold.

The proposed project gate requires all eight bacterial members to be detected
and community precision of at least 0.90. Precision is detected expected members
divided by detected expected members plus unexpected species calls. Also retain
whole-community recall, every missing member, every unexpected species (including
below-threshold calls), their total assigned mass and unclassified fractions.
No calls gives undefined precision (`null`), not perfect precision. There is no
reference-defined negative universe, so specificity is not estimated.
Unexpected calls are relative to the declared mixture; laboratory contamination
and reference ambiguity can contribute and must not automatically be called
pipeline implementation bugs.

Nominal proportions are DNA fractions: eight bacteria at 0.12 each, two yeast
members at 0.02 each. They are **not exact truth for read/fragment fractions**.
The scorer's nominal-DNA L1 distance includes unexpected and Bracken-unassigned
mass as additional categories. It is descriptive, not an abundance-accuracy gate.
A read-origin truth dataset is still needed to establish abundance error.
One mock sample does not provide biological replication or population uncertainty.

## Run and preserve evidence

Use a repository checkout for these development tools. Obtain the full reads
with the download and MD5 commands in [acceptance checks](acceptance.md), without
the 10,000-pair subsetting step. Check sufficient disk before starting. The first
local baseline uses four threads, a 2 GB BBMap Java heap and memory-mapped Kraken;
reserve about 15 GB working space. Assembly uses MEGAHIT's default memory policy.
The observed resource record, rather than this planning estimate, is authoritative.
No new taxonomy or functional database download is needed if the recorded release
is installed already. Keep MEGAHIT scratch on a POSIX filesystem under WSL.

```bash
sha256sum benchmarks/zymo-even/protocol.json
metaquest --version
/usr/bin/time -v -o pipeline-resources.txt \
  metaquest --no-color run --config benchmarks/zymo-even/config.yaml \
  --db-dir "$DB_ROOT" --low-memory --skip-functional \
  --paired "$READS/ERR2984773_1.fastq.gz" "$READS/ERR2984773_2.fastq.gz" \
  --output "$OUT"
python tools/verify_run.py "$OUT" > accounting.json
python tools/benchmark_taxonomy.py "$OUT" \
  --protocol benchmarks/zymo-even/protocol.json > taxonomy-benchmark.json
```

Set `DB_ROOT`, `READS` and a fresh `OUT` to your local paths first. The scorer exits
0 for a passing bacterial baseline gate, 1 for a failed gate **while retaining the
metrics**, and 2 for invalid artifacts/provenance. Keep its exit code and stderr.
It checks input MD5s, full raw counts from fastp, post-filter paired denominators,
database collection/release and agreement between source tables and the summary.
Validation's sampled QC count is not mistaken for the full read count.

Record the source revision and uncommitted patch, imported module path, package
and native versions, effective settings, database manifests, input SHA-256s,
commands, wall time, maximum RSS reported by GNU time and output size. On Linux,
the child RSS maximum is the largest child value, not aggregate simultaneous
memory use; it also does not measure the machine's total filesystem-cache use.
Do not commit reads, native databases or large generated outputs.

For integration checking, independently run Kraken2 and Bracken on the exact
cleaned mates, with the same database, hit-group threshold, species rank,
read-length model, read threshold and thread count recorded in metadata. Compare
taxid/count maps and Kraken accounting with MetaQuest; this checks orchestration
and parsing, not the independent accuracy of Kraken/Bracken. Likewise, compare
direct gene calls on identical contigs using the recorded Pyrodigal release.
MEGAHIT record ordering/headers can vary: compare sequence content, not FASTA bytes.

For this four-thread, 151 bp-input protocol, the closest installed Bracken model
is 150 bp. Confirm the recorded model in `analysis_metadata.json` before repeating:

```bash
mkdir -p "$DIRECT"
kraken2 --db "$DB_ROOT/taxonomy" --threads 4 --minimum-hit-groups 2 \
  --memory-mapping --report "$DIRECT/kraken_report.txt" --output /dev/null \
  --paired "$OUT/preprocessing/cleaned_R1.fastq.gz" "$OUT/preprocessing/cleaned_R2.fastq.gz"
bracken -d "$DB_ROOT/taxonomy" -i "$DIRECT/kraken_report.txt" \
  -o "$DIRECT/bracken_report.tsv" -w "$DIRECT/bracken_report.txt" -r 150 -l S -t 10
cmp "$OUT/kraken_report.txt" "$DIRECT/kraken_report.txt"
cmp "$OUT/bracken_report.tsv" "$DIRECT/bracken_report.tsv"
```

Set `DIRECT` to a fresh directory. A failed byte comparison warrants a taxid/count
comparison to distinguish ordering changes from scientific discrepancies.

The authors provide a checksum-pinned, 71,433,493-byte (68.1 MiB)
[isolate draft reference](https://lomanlab.github.io/mockcommunity/).
Use QUAST with an explicit minimum contig length and fragmented-reference mode
to report reference coverage, alignment identity, unaligned sequence and apparent
misassemblies. Draft-reference disagreements can originate in either assembly;
they are not automatically MetaQuest errors. MetaQuest GFF3 and derived evidence
use **1-based inclusive** coordinates. If reading PAF/BED alignments for additional
checks, those intervals are **0-based half-open** and require explicit conversion.

After downloading and checking the reference MD5 in the protocol, this baseline
uses QUAST 5.3.0 (an evaluation dependency, not a MetaQuest runtime dependency):

```bash
quast.py "$OUT/megahit_assembly/final.contigs.fa" -r reference.fasta \
  -o "$QUAST_OUT" -t 4 --min-contig 200 --fragmented --no-html --no-icarus --no-plots --debug
python tools/benchmark_reference.py "$QUAST_OUT" --reference reference.fasta \
  --protocol benchmarks/zymo-even/protocol.json > reference-assessment.json
```

Use a fresh `QUAST_OUT` directory and retain `report.tsv`, `quast.log`, alignments
and GNU time output. Report per-member reference coverage when interpreting the
combined draft reference; the low-coverage yeasts affect its aggregate fraction.
The summarizer weights coverage by non-N reference bases, checks that the member
totals reproduce QUAST's rounded aggregate, and rejects multi-assembly maxima.
`--debug` retains QUAST's normalized reference; its sequence identities and base
counts must match the supplied checksum-pinned FASTA, allowing QUAST's explicit
uppercase/IUPAC-to-N normalization. This guards against using the wrong output folder.
Its successful exit confirms a valid assessment, not biological correctness.

## Remaining checkpoint-3 gates

- An independent held-out community, chosen before tuning, with declared truth,
  licensing, database coverage and acceptance limits.
- Read-compatible abundance truth and an error metric; DNA proportions alone
  cannot establish it.
- Reference-supported assembly and coding-sequence assessment. Direct Pyrodigal
  agreement is an integration check, not coding-gene truth. Its prokaryotic model
  does not establish yeast gene prediction accuracy.
- Functional annotation assessment against independently curated reference
  labels. Annotation rate alone is insufficient. This first core baseline skips
  functional annotation and screening; retain existing native regression controls.
- Measured resource baselines and budgets before future regression comparisons.
- Results, failures and residual limits published together in a focused PR.

Passing the bacterial gate does not close checkpoint 3 or justify a stable release.

One candidate for the next independent benchmark is
[MBARC-26](https://doi.org/10.1038/sdata.2016.81), which provides 23 bacterial and
three archaeal members with finished reference genomes. Its input checksums,
reference versions, abundance units and resource budget still need verification
before selecting a held-out protocol. No MBARC reads have been downloaded or scored.
