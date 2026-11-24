# Host-removal integration and resource controls

Measured 2026-11-23 against unreleased `2.0.0b5`, based on published beta 4
commit `ebd860161218941219687ce56cad9113009216fa`.
These are deterministic software controls and a small accounting smoke run.
They do not measure held-out host-removal accuracy or close the failed
Zymo biological benchmark gate.

## Protocol and environment

[`protocol.json`](protocol.json) was written before measurement. The reproducible
runner is [`tools/benchmark_host_removal.py`](../../tools/benchmark_host_removal.py).
Filtering uses an installed wheel outside the source tree: Python 3.11.6,
Bowtie2 2.5.4 and Samtools 1.18, on Ubuntu 24.04/WSL2 Linux
6.18.58.1, Intel Core i5-13500H (16 visible logical CPUs), with two mapping
threads. Both implementations use the same inputs, index, native versions,
preset, seed, primary flags, compression level and decoded-record verifier.
Each measurement starts a fresh Python worker. No performance measurements
overlapped other acceptance tests.

The fixed seed generates a 50,000 bp reference, an unrelated 50,000 bp sequence
and 10,000 pairs of 150 bp reads. Four groups of 2,500 pairs match the reference
in both mates, R1 only, R2 only, or neither. Synthetic slicing uses 0-based
Python offsets. Filtering uses flags from 1-based SAM records, not intervals.
The expected retained set is the neither-matching group: **2,500 fragments**.
All six runs matched this set, sequence and quality content exactly.
This round count follows construction of the controls; it is not an accuracy
estimate or independent biological evidence.

Input SHA-256:

| File | SHA-256 |
| --- | --- |
| host.fa | `cadb2926de3e1595b4add1cc638ef6d01b0a59797d4197ebd816352f530c60d6` |
| R1.fastq | `889df0b45dd35583955ad931e30f9990a13d585b2c03cd484d6515ec708bda10` |
| R2.fastq | `c49bdb955165e66ba059b8f493d1b5dbc20638eac21fb8991901c014c5853b2f` |
| truth.json | `58e66ca9ffff3237e3ef24fbcbb007760c4df6d5836f5b5558e280b70774734f` |

## Resource measurements

The direct native worker pipes Bowtie2 to Samtools FASTQ conversion. MetaQuest
also hashes inputs/index/outputs, tees SAM into independent flag statistics,
checks FASTQ counts and writes provenance. Both workers include the same
decoded-record verification and Python startup. Index construction is excluded
from filtering: GNU time recorded 1.62 s elapsed, 0.34 s CPU and 95,120 KiB
maximum RSS for that separate step.

| Order | Worker | Monotonic elapsed (s) | Process-tree CPU (s) | GNU maximum RSS (KiB) |
| --- | --- | ---: | ---: | ---: |
| 1 | Direct native | 2.007122 | 3.329023 | 47,960 |
| 2 | MetaQuest | 2.717831 | 3.433886 | 49,088 |
| 3 | MetaQuest | 2.834921 | 3.574303 | 49,500 |
| 4 | Direct native | 2.025011 | 3.258945 | 48,068 |
| 5 | Direct native | 1.925896 | 3.274101 | 48,244 |
| 6 | MetaQuest | 1.949571 | 2.821680 | 49,352 |

Median elapsed was **2.717831 s** for MetaQuest and **2.007122 s** for the native
worker. The predefined limit `max(2 × native median, native median + 1 s)` was
**4.014243 s**, so this small resource-regression gate passed. CPU and memory
are reported without independent pass thresholds. GNU maximum RSS is its
reported maximum, not a sum of simultaneous process memory. Raw GNU elapsed
timers and monotonic measurements are both retained locally; rounded GNU times
differ by approximately 0.03–0.04 s from parent timing. No speedup is claimed.
Filesystem/cache noise and these tiny references limit extrapolation.

## Cached public-data accounting smoke

The same native dependency environment ran the exact public beta-4 wheel and
the beta-5 candidate on 1,000 previously hash-verified ENA `ERR2984773` pairs
plus 100 generated reference-matching pairs. No download was needed. The cached
Standard-8 (2026-06-26) Kraken/Bracken database and settings were unchanged.

- Disabled filtering: Kraken report, Bracken profile, taxonomy evidence and
  taxonomy plot-data TSV were byte-identical between beta 4 and beta 5.
- Enabled filtering: fastp retained 1,050 fragments; the filter removed 100 and
  retained 950. Every retained decoded FASTQ record matched direct Bowtie2/
  Samtools execution on the same cleaned inputs. Taxonomy used 950 as its
  input denominator; accounting and report-link checks passed.
- Unchanged resume preserved the completion JSON and scientific tables. A
  deliberately altered native log invalidated reuse; recomputation retained the
  altered attempt and produced unchanged scientific tables. A missing index
  failed preflight.

This smoke covers taxonomy accounting, not full annotation, abundance accuracy,
human-genome memory requirements, privacy protection or clinical validation.
Generated controls and repeated smoke subsets are not held-out biological data.
Shared microbial/reference sequence can be removed, and divergent host sequence
can remain. No new database was selected or classifier default tuned here.

Additional installed-wheel checks covered single reads (1,090 after fastp,
100 removed, 990 retained) and interleaved pairs with fastp disabled (1,100
fragments in, 100 removed, 1,000 retained). Both passed accounting and preserved
unchanged completion JSON on resume.

A separate cached 10,000-pair ENA subset plus 100 generated reference-matching
pairs exercised taxonomy, MEGAHIT, Pyrodigal and BBMap gene abundance. fastp
retained 9,669 fragments, filtering removed 100 and retained 9,569. The run
produced 100 contigs and 106 predicted genes and passed the output/ID/accounting
checker. Only 3.53% of fragments mapped to the shallow assembly, so abundance
describes that assembled fraction; these results do not validate abundance
accuracy. Functional annotation and research screening were skipped in this
additional pipeline smoke. The full software suite separately exercised
controlled eggNOG subprocess fixtures and real AMRFinderPlus positive/negative
and zero-hit controls; it did not run a new native eggNOG database search.

## Software checks

The feature source and independently installed wheel each passed **407 tests,
zero skips**, including all 41 host-filter cases and the existing native
screening control using AMRFinderPlus data 2026-10-29.1. A fresh pip-only wheel
environment and extracted source distribution each passed **389 tests, 18 native
controls skipped**. Fresh pip dependencies included NumPy 2.4.6, pandas 3.1.1,
Matplotlib 3.11.2 and Biopython 1.88. The Conda environment emitted 276 existing
Matplotlib/Pyparsing deprecation warnings; the fresh pip runs emitted none.
Warnings were retained rather than suppressed.

Ruff, strict Twine metadata checks, direct-dependency license checks, `pip check`,
43 packaged runtime/config/catalog file comparisons and independent sdist rebuild
checks passed. Regressions cover single/paired/interleaved modes, one-mate
matches, all/zero retained, duplicate-name loss, unsupported versions, changed
reference/output/evidence, cohort compatibility, safe symlink/FIFO output
handling and SIGINT/SIGTERM child-process cleanup. Interleaved resume reuses only
validated split files, and preprocessing HTML keeps counts ahead of detailed
downloadable provenance.

Offline Chromium 156 and Firefox 157 checks passed all 54 desktop, 375 px narrow
and print views across the nine report pages, with no page overflow, missing
images or broken local links. Representative narrow preprocessing and rendered
PDF host sections were visually inspected for count readability, wrapping and
complete interpretation text. Physical printing was not tested.

Raw commands, hashes, environment solves, failed controls, JUnit records, native
logs and browser captures are retained in the ignored local
`test_results/host-removal/` audit directory. The protocol and runner are shipped
so the generated controls can be reproduced independently.
