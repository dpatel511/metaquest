# ERR2984773 development baseline

Status: full core run completed; **the frozen primary taxonomy gate failed**.
Checkpoint 3 remains open. This is development data, not held-out validation.

The [protocol](protocol.json) was saved before the full-data run with SHA-256
`30c859006c71a974a99baefba3644161afa1fc8fd695fd439851548f27b07ea0`.
The pipeline source is merged revision
`48c2ceedcd8d0c5db815ad1b93cfd65acb36d085`, version `2.0.0b1`; this branch adds
development tools/documentation without changing the pipeline engine.

## Software checks

| Check | Observed result |
| --- | --- |
| New scorer regressions | 39 passed |
| Final source suite with native controls, clean Python dependencies | 212 passed; no skips or warnings |
| Earlier prototype suites using inherited developer dependencies | 201 source / 201 wheel tests passed; 276 Matplotlib/Pyparsing deprecation warnings each |
| Inherited-environment dependency check | Failed: eggNOG-mapper metadata pins older Biopython, psutil and XlsxWriter |
| Final fresh isolated installed-wheel suite with native controls | 212 passed; no skips or warnings |
| Fresh wheel dependencies, direct dependency licenses, strict metadata | Passed |
| Packaged source/config/catalog consistency | All 40 files matched the checkout; AGENTS.md excluded |
| Ruff | Passed |

The clean environment's wheel imports were verified outside the source tree.
The candidate wheel is local and unreleased; it does not replace the public b1.
The failed inherited-environment check remains recorded; it is not evidence of
a clean Conda environment. Functional annotation is excluded from this core run.

## Native benchmark

Both full cached input files matched their ENA sizes and MD5s. The author-provided
isolate draft reference matched MD5 `e7b60972c37869e0e91513a06bdf9d33`.
The run used all 8,785,731 raw pairs, four configured processing threads,
memory-mapped Kraken, a 2 GB BBMap Java heap, and the [recorded config](config.yaml).
Functional annotation and screening were not requested for this core baseline.

| Measurement | Observed |
| --- | ---: |
| Post-fastp paired fragments | 8,504,800 |
| Kraken-classified fragments | 7,573,632 |
| Unclassified fragments | 931,168 (10.9487%) |
| Bracken species-assigned fragments | 7,569,048 |
| Raw Bracken reported taxa | 192 |
| Assembly | 19,526 contigs; 40,468,399 bp; N50 89,640 bp |
| Candidate CDS / protein sequences | 49,573 |
| Fragments mapped to assembly | 8,068,261 (94.8671%) |
| Fragments assigned to genes | 7,154,552 (84.1237%) |
| Accounting/artifact check | Passed: eight stage reports, 16 figure files |

Mapping and assignment are fragment-level accounting, not abundance-accuracy
estimates. Candidate CDS counts do not establish coding-gene truth.

## Taxonomy result: retain the failure

At the frozen primary threshold of **0.1% of all surviving fragments**, all eight
bacterial members were detected, with five unexpected species calls. Nominal
community precision is **8/13 = 61.54%**, below the required 90%. Whole-community
recall is **8/10 = 80%**; neither yeast was reported. Standard-8 excludes fungal
libraries, so these missing calls do not establish biological absence.

| Unexpected primary-threshold call | Taxid | Estimated fragments | % of all surviving fragments |
| --- | ---: | ---: | ---: |
| Shigella flexneri | 623 | 19,340 | 0.2274 |
| Bacillus subtilis | 1423 | 30,587 | 0.3596 |
| Enterococcus sp. AZ037 | 2774924 | 9,066 | 0.1066 |
| Enterococcus sp. DIV1316a | 2774959 | 56,234 | 0.6612 |
| Enterococcus sp. DIV1202 | 2774964 | 34,007 | 0.3999 |

All unexpected taxa, including those below the primary threshold, account for
**2.6982% of input fragments**. This is a different metric from detection precision;
61.54% is not per-read classification accuracy. Unexpected calls are relative to
the declared mixture; reference ambiguity, algorithm behaviour and contamination
need adjudication before assigning a cause.

The two predeclared secondary thresholds are retained as descriptive results:

| Threshold | Bacterial recall | Community recall | Community precision | Unexpected calls |
| --- | ---: | ---: | ---: | ---: |
| **0.1% (primary)** | 8/8 | 8/10 | 8/13 (61.54%) | 5 |
| 0.01% | 8/8 | 8/10 | 8/34 (23.53%) | 26 |
| 1% | 8/8 | 8/10 | 8/8 (100%) | 0 |

The secondary result does not replace the failed primary gate. No thresholds or
pipeline defaults were tuned after viewing this run. Nominal-DNA L1 distance is
0.399442, **descriptive only** because DNA fractions are not exact fragment truth.

## Independent component and reference checks

Direct Kraken2 and Bracken runs on the identical cleaned mates reproduced both
raw tables **byte-for-byte**, including the unexpected calls. Taxid-keyed rows
also matched. Direct Pyrodigal calls matched all 49,573 protein sequences and
CDS coordinate tuples. GFF3/Pyrodigal coordinates are **1-based inclusive**.
These establish integration consistency for the checked components, not the
independent biological correctness of their algorithms or databases.

QUAST 5.3.0 used the pinned isolate drafts, minimum contig length 200 bp,
minimum alignment 65 bp, minimum identity 95%, one ambiguity assignment,
fragmented-reference mode and retained debug artifacts. The normalized reference
sequence identities and per-fragment base counts were independently reconciled.

| Reference member | Non-N reference bases covered (%) |
| --- | ---: |
| B. spizizenii (paper: B. subtilis ATCC 6633) | 99.81 |
| E. faecalis | 99.95 |
| E. coli | 99.80 |
| L. fermentum | 99.71 |
| L. monocytogenes | 98.58 |
| P. aeruginosa | 99.84 |
| S. enterica | 99.67 |
| S. aureus | 99.90 |
| S. cerevisiae | 50.00 |
| Cryptococcus hybrid | 12.97 |

The aggregate is **56.785%**: 39,830,420 covered / 70,142,935 non-N reference bases.
QUAST also reports **1,197 apparent extensive misassemblies**, 78 local events,
97,364 unaligned bases, 80.45 mismatches and 6.58 indel events per 100 kbp.
The 1,124 extensively flagged contig IDs reconcile with the native report:
810 have best alignments only to yeast references, 107 only to bacterial
references and 207 to both groups. This association is not an error adjudication.
The draft/hybrid reference and genuine assembly errors both need consideration;
high coverage alone does not establish correctness. Retain these disagreements.

## Resources and timing limitations

The WSL environment exposed about 7.6 GiB RAM and 16 logical CPUs. Data and final
outputs were on Windows-mounted storage; MEGAHIT scratch used Linux storage.

| Resource | Observed |
| --- | --- |
| GNU time elapsed | 1:10:32 (4,232 seconds) |
| CLI completion text | 1h 15m |
| GNU time maximum RSS | 6,750,492 KiB (~6.44 GiB); not summed concurrent process-tree memory |
| Retained pipeline output bytes | 3,461,900,211 |
| Pipeline / accounting / taxonomy scorer exits | 0 / 0 / 1 |

The elapsed timers **disagree**. The CLI uses a monotonic clock; the environment
and timer sources need a controlled repeat before establishing a performance
budget. Do not silently select one timing as validated. GNU time's zero Linux
block-I/O counters also do not represent the files read/written through the
Windows-mounted filesystem. Some software verification overlapped initial QC;
this is not an isolated performance experiment.

The initial 45-minute observer was replaced with a 90-minute total observation
window while the same native CLI continued uninterrupted. Its metadata contains
one completed attempt, with no forced cache/provenance edits. This observer
handoff and both timing outputs remain in the local evidence.

Local evidence is retained under ignored `test_results/benchmark-check/`, including
the frozen protocol/input/source snapshot, commands, full logs and failed checks.
Large reads, reference sequences and generated results are excluded from Git.
PNG/SVG/PDF benchmark plots and source TSVs were generated; the PNG renderings
were visually checked. Native report HTML passed structural checks; browser
layout review remains checkpoint 5.

## Limits

This previously used mock is development data, not held-out validation. Standard-8
excludes fungal libraries. DNA proportions are not exact read/fragment truth.
The isolate reference is a filtered draft assembly, not finished-genome truth.
Direct component agreement is an integration check. Prokaryotic gene prediction
does not establish yeast gene accuracy. Independent abundance/functional truth,
held-out data and resource regression budgets remain open.

Next: adjudicate unexpected taxonomy calls and structural reference disagreements,
freeze an independent held-out protocol before evaluating it, and repeat timing
under controlled conditions. The current evidence does not justify closing
checkpoint 3 or declaring stable biological accuracy.
