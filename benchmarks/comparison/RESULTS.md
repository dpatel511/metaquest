# Checkpoint 4: comparative-analysis audit

Status: local validation criteria met after correcting input/output contracts and
restricting unsupported small tied studies. The safeguards are included in
**2.0.0b3** after [PR #34](https://github.com/dpatel511/metaquest/pull/34) merged and
[main CI](https://github.com/dpatel511/metaquest/actions/runs/41244598939) passed.
Original `2.0.0b2` artifacts are unchanged.

## Scope and frozen protocols

The source baseline is `cd1050fc8a6f9c1fe7766a42e3dacd3337e03bce`. Local checks
used Python 3.11.9, NumPy 2.4.6, pandas 3.1.0 and SciPy 1.17.1. The additional
supported Python versions passed the merged audit CI matrix (3.10/3.11/3.12).
The [primary protocol](protocol.json) was frozen before simulation with SHA-256
`44b236dd6e3bbd97c0cab31d1f59e0e61c3351cfc4ccfc939881e2adb6253b5a`.
The supplemental [exact tied-null protocol](tied-null-protocol.json), added after
noticing the sparse CLR studies had no transformed ties, has SHA-256
`5905dfa525aa219dfd383cc6b7802180fc6ca5c857032b71213a7995e6a2cef5`. Neither was edited
after scoring. The supplemental resource-planning text says 474 states; the actual
complete count is 374. This counting typo does not change the population grid or criteria.

All studies are synthetic calculation models. They are not biological replicates
or independent mock-community validation. Details and reproduction commands are in
[comparison validation](../../docs/comparison-validation.md).

## Defects reproduced and fixed

- Repeated subjects hidden by whitespace and even one shared raw-file hash across samples.
- Incomplete preprocessing/Bracken-model provenance, malformed profile IDs/quantities,
  duplicate headers, rank/count-summary disagreement and conflicting feature labels.
- Missing sample/feature alignment checks at the statistical API boundary.
- Existing reports overwritten through the library API. The CLI diagnostic log is
  allowed and preserved; fresh/forced CLI output and backups have dedicated tests.
- Text methods omitted actual transformation, pseudocount, filtering and contrast.
- The exact diagnostic now reads declared protocol settings; a regression verifies
  that changing alpha changes its calculations. Frozen-protocol findings are unchanged.
- Small nonconstant tied CLR data were admitted to an unreliable asymptotic test.
  They now fail with a permutation-test explanation when either group has <10 samples.
  Constant pooled features retain p=1; untied small studies and larger tied studies remain supported.

The initial contract reproductions retained 30 failing cases. A malformed checksum
case, the initially introduced CLI/log guard regression and the tied-calibration
failure were subsequently reproduced and corrected. Failed checks remain in local evidence.

## Primary Monte Carlo results

500 independent studies/scenario, 24 features, alpha 0.05, pseudocount 1e-6,
overall prevalence 0.1. Both corrections evaluate the same studies. The final run
reproduced all **6,100 per-study metric records byte-for-byte** from the initial run
(4,000 primary and 2,100 sensitivity evaluations). Repeating a run does not add
independent studies or narrow the original confidence intervals.

| Global-null scenario | Correction | Studies with any discovery | FDR estimate | 95% Wilson interval |
| --- | --- | ---: | ---: | --- |
| dense_null | BY | 1/500 | 0.20% | 0.04–1.12% |
| dense_null | BH | 17/500 | 3.40% | 2.13–5.38% |
| unbalanced_null | BY | 9/500 | 1.80% | 0.95–3.39% |
| unbalanced_null | BH | 20/500 | 4.00% | 2.60–6.10% |
| sparse_null | BY | 0/500 | 0.00% | 0.00–0.76% |
| sparse_null | BH | 14/500 | 2.80% | 1.68–4.64% |

| Known-effect scenario | Mean FDP | Upper one-sided 95% bound | Mean target recall | Lower one-sided 95% bound |
| --- | ---: | ---: | ---: | ---: |
| BY | 0.60% | 6.07% | 97.70% | 92.23% |
| BH | 3.29% | 8.76% | 99.20% | 93.73% |

Both target effect directions were correct in 500/500 studies for each correction;
the Wilson lower bound is 99.24%, not proof of universally perfect direction recovery.
BY sparse-null 0/500 discoveries still has a 0.76% upper interval bound.
All proposed screening limits passed. The nominal target remains 5%; the declared
Monte Carlo margins are not proof of universal 5% FDR control. BH dense/unequal-null
intervals and the conservative effect-FDP bounds extend above 5%.

Sensitivity settings left the dense-effect target recall and discovery counts unchanged.
That population is entirely positive, so its prevalence-threshold check cannot establish
robustness for sparse data. The first 50 sparse-null studies at pseudocount 1e-4 had
0 BY and 1 BH discovery-bearing studies; this is descriptive, not a tuned gate.

## Tied-data failure: retain the public-beta result

The main sparse model had zeros but **no transformed ties**. Exact binomial-weighted
calculations over a two-point null population exposed the missing coverage. The public
wheel independently reproduced every one of the 108 population/correction results,
including these four failures:

| Groups | High-state probability | Correction | Exact global-null FDR at alpha 5% |
| --- | ---: | --- | ---: |
| 3 / 30 | 0.05 | BH | 8.2074% |
| 3 / 30 | 0.1 | BH | 6.7900% |
| 3 / 30 | 0.9 | BH | 6.7900% |
| 3 / 30 | 0.95 | BH | 8.2074% |

These failures arise in the existing small-sample asymptotic p-value policy. The
candidate rejects 486 nonconstant small/tied count-pair/correction states rather
than silently dropping features or returning zero discoveries. Unsupported population
probability mass is recorded as unsupported and its FDR is **null**, not 0.
All 18 fully supported population/correction calculations (10 / 10 samples) meet
the nominal threshold in this particular two-point family. This is not universal calibration.

## Software and distribution verification

| Check | Result |
| --- | --- |
| source full suite | 280 passed, 2 skipped, 0 failures, 0 errors |
| wheel full suite | 280 passed, 2 skipped, 0 failures, 0 errors |
| sdist full suite | 280 passed, 2 skipped, 0 failures, 0 errors |
| Skipped controls | Two unchanged native AMRFinderPlus controls: local screening DB unavailable |
| Lint, strict wheel/source metadata, clean dependencies and licenses | Passed |
| Packaged runtime/config/catalog | 40 files match the checkout; AGENTS.md excluded |
| Actual installed-wheel CLI | Valid comparison, repeated-input failure streams, backups and small-tie errors verified |
| Real b2 taxonomy output | Accepted unchanged: 5 features, 623 assigned fragments |
| Public b2 numerical parity | 80 admitted study/correction evaluations; 1,898 feature rows agree within 1e-12 |
| Protocol and final simulation source | Frozen protocol and all 40 source/config/catalog hashes unchanged |

During the local audit, the same b2 version was used only for an ignored candidate
wheel; it was not a replacement upload of public b2. The local skips above remain
recorded. The subsequent merged native CI ran all 282 tests without skips, covering
the AMRFinderPlus controls. Beta 3 carries the fixes in a new version. Source
distributions include tools, protocols and test support.

## Resources and remaining limits

Final recorded resource output is retained in `test_results/comparison-check/`.
GNU time and the Python monotonic timer disagree again; retain both rather than
choosing one as authoritative. These are observed laptop runs, not a controlled
performance baseline. Software verification overlapped part of the repeated run.
The initial run recorded GNU time 5:47.92, harness elapsed 373.68 seconds, and
maximum RSS 191,516 KiB (~187 MiB). The final run recorded GNU time 5:39.71,
harness elapsed 363.26 seconds and maximum RSS 191,128 KiB (~187 MiB). The repeated
run remains the same 2,000 studies, with no additional confidence claimed.

No larger databases or sequence downloads were used. Biological abundance truth,
batch/covariate adjustment, pairing, independent biology, browser rendering and
user acceptance remain outside this local statistical audit. Checkpoint 3 remains
open and deferred; the proposed next local work is checkpoint 5 report/data review.
