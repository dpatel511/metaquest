# Optional host-reference filtering

This feature belongs to **2.0.0b5 (unreleased)**. It is disabled by default.
MetaQuest uses Bowtie2 to align reads against a supplied local reference index
and Samtools to retain unmapped reads. It does not download a host genome or
infer a reference from taxonomic results.

## Setup and execution

Install the checkout's core Conda environment, including Bowtie2 `>=2.5,<3` and
Samtools `>=1.18,<2`, then install the beta-5 checkout/wheel. Keep screening in
its separate environment. Build an index from a reference you have selected
and reviewed; record its source, assembly and any masking separately:

```bash
bowtie2-build /data/references/selected-host.fa /data/references/selected-host
metaquest check --db-dir /data/metaquest-db --host-index /data/references/selected-host
metaquest run --db-dir /data/metaquest-db \
  --host-index /data/references/selected-host \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --taxonomy-only --output results/sample-filtered
```

Pass the **prefix**, not the FASTA or an individual `.bt2` file. Exactly one
complete, nonempty six-file `.bt2` or `.bt2l` index is required. Paths with spaces
are supported. Genome size, indexing storage and mapping cost depend on your
reference; no human-sized reference has been benchmarked here.

The equivalent configuration is:

```yaml
preprocessing:
  enabled: true
  host_index: /data/references/selected-host
  threads: 2
```

`--host-index` overrides YAML. Omit the flag and use `host_index: null` to
disable filtering. Setting fastp `enabled: false` leaves host filtering enabled
when an index is configured. Single, paired and interleaved input modes work;
interleaved reads are split into synchronized mate files.

## Policy and output contract

Filtering follows fastp and precedes every downstream read-consuming stage.
Bowtie2 runs with `--sensitive -k 1 --mm --reorder --seed 0` and the configured
preprocessing thread count. There is no MAPQ cutoff. The fixed policy removes
the **whole pair if either mate aligns**, including discordant or mixed
alignments. It retains only primary records with both mates unmapped; single-end
mode retains primary unmapped reads. Using Bowtie2 `--un-conc` alone would not
implement this policy, because failure to align concordantly does not mean both
mates are unmapped. SAM text positions use **1-based** coordinates; filtering
uses SAM flags rather than genomic interval arithmetic.

The filter streams SAM through Samtools `fastq` and `flagstat`; it does not
materialize a SAM/BAM file. Samtools conversion normalizes mate suffixes and
omits FASTQ header comments. Sequences and qualities are retained. Output:

- `preprocessing/host_removal/retained_R1.fastq.gz` and, for pairs, `retained_R2.fastq.gz`;
- `preprocessing/host_removal/summary.json`, schema 1;
- `attempt-*/` command logs and primary mapping statistics, including failed attempts.

The JSON records `state.inputs` (absolute paths/full hashes), `state.index`
(prefix/full hashes of all six shards), tool versions, paired mode, threads,
preset, seed and pair policy. Exact command vectors, relative output/evidence
paths and full hashes are retained. Index hashes identify the supplied binary
files; they do not establish the original FASTA's source or biological identity.
Hashing is streaming but adds reference-sized I/O, including during strict resume.

`input`, `retained` and `removed` count **fragments** in paired/interleaved mode,
and **reads** in single-end mode. `input = retained + removed`; `retained_fraction`
is 0–1, or null for zero input. Input here means after fastp, or original reads
when fastp is disabled. fastp's counts remain individual reads. Downstream taxonomy
and abundance use the retained read files; taxonomy's all-input denominator is
the retained fragment/read count. Reports retain eight stages, with host
evidence in preprocessing and links to the full methods JSON.

Conversion output counts are checked independently against primary mapping
counts. Samtools can collapse consecutive records sharing a QNAME; disagreement
fails before completion rather than silently discarding fragments. If every
fragment is removed, the filter records zero retained and analysis stops with
an error before taxonomy. Completion is atomic; partial attempts are preserved.

## Resume, comparison and verification

Strict resume requires unchanged index bytes, input hashes, effective settings
and native versions. Output and native evidence hashes must match before reuse.
Changed retained files/evidence trigger recomputation. Use a new output directory
when the reference or settings change; do not edit recorded provenance.
When filtering is disabled, legacy provenance without `host_index` remains
compatible. Enabled/disabled cohort profiles cannot be mixed. Comparisons require
the same index hashes, native versions and mapping policy/settings; local index
paths and worker counts may differ. These checks do not introduce batch adjustment.

```bash
python -m pytest -q tests/test_host_removal.py
python tools/benchmark_host_removal.py --output /tmp/host-controls-fresh
python tools/verify_run.py results/sample-filtered
```

The resource protocol is versioned in
[`benchmarks/host-removal/protocol.json`](../benchmarks/host-removal/protocol.json).
Native tests require Bowtie2, Bowtie2-build, Samtools and the core tools; pip-only
test environments skip those controls explicitly. Generated controls test record
handling, not held-out human-read removal accuracy. The resource runner requires
GNU `/usr/bin/time`, an installed non-editable beta-5 wheel and a fresh output
directory; it verifies that imports come from the selected environment.

## Interpretation limits and component references

Mapping depends on the selected reference, masking and alignment settings.
Reads from unrepresented or divergent host sequence may remain. Shared
microbial/reference sequence can be removed; a microbial label does not protect
a matching read. Retained reads do not establish host absence, privacy protection,
infection status or clinical suitability. No sensitivity/specificity claim or
biological improvement follows from the deterministic controls.

The policy follows the documented meanings of
[Bowtie2's alignment/output options](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
and [Samtools FASTQ conversion](https://www.htslib.org/doc/samtools-fasta.html).
[Hostile](https://pmc.ncbi.nlm.nih.gov/articles/PMC10749771/) provides related
published precedent for Bowtie2/Samtools host filtering and discusses reference
masking; its performance and accuracy results do not transfer to MetaQuest.
