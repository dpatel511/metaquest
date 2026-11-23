# Beta 4 release acceptance

Target: `2.0.0b4`, research beta. This record accompanies release preparation;
it does not assert publication or completion of stable-release checkpoints.

## Changes included

- Strict YAML configuration handling rejects unknown keys before analysis.
- The complete environment selects corrected public eggNOG-mapper build
  `2.1.15=pyhdfd78af_1`; version probing and annotation execution do not borrow
  the caller's Git hash.
- Taxonomy evidence schema 1 separates Kraken direct/clade observations from
  Bracken assigned/redistributed estimates, records denominators and source
  hashes, and exposes methods in stage and consolidated reports.

Classifier algorithms and defaults are unchanged. Host removal is a separate
planned feature and is not included in this release.

## Verification

The implementation merged in [PR #40](https://github.com/dpatel511/metaquest/pull/40).
All six checks passed on merge revision
`336a3e5a68b3404e75262ea8ab717e5459a37a87` in
[main CI](https://github.com/dpatel511/metaquest/actions/runs/45060683494).
The native integration artifact retains JUnit results and its screening
database manifest. Subsequent release-preparation revisions need their own CI.

The [taxonomy evidence results](../benchmarks/taxonomy-evidence/RESULTS.md)
retain the 200,000-pair integration comparison and repeated reporting timings.
Kraken and Bracken tables were byte-identical across the reporting change;
the reporting median met the predeclared resource budget. This does not
establish a speedup or improved biological accuracy.

Additional report review passed 54 offline views across Chromium 156.0.8078.4
and Firefox 157.0: desktop 1440 px, narrow 375 px and print media 794 px.
Representative Firefox screenshots and Chromium A4 PDF taxonomy pages were
visually inspected. Evidence columns fit; methods continued across page breaks.
Native print dialogs and physical printers were not tested.

The full installed-wheel smoke exposed a remaining execution-path provenance
bug: eggNOG's raw headers borrowed the output repository's Git hash even though
the isolated version probe was correct. [PR #41](https://github.com/dpatel511/metaquest/pull/41)
fixes this at revision `249c4479f1f9c7c7bf60b97c9dd8e9589504f19e`.
Both new subprocess regressions failed before the fix and passed afterward.
All six [fix CI jobs](https://github.com/dpatel511/metaquest/actions/runs/45063112565)
passed; native CI passed 366 tests without skips.

Local checks used Python 3.11.6 on Linux/WSL, the public-build core Conda
environment and a separately resolved AMRFinderPlus 4.2.7 environment with
database `2026-10-29.1`. The isolated wheel imported from its own site-packages
outside the checkout but shared Conda dependency packages; this was not a fresh
installation of every dependency.

| Check | Observed result |
| --- | --- |
| Fixed source, installed wheel, extracted sdist | Each passed 366 tests, no skips; 276 existing plotting deprecation warnings per suite |
| Ruff, direct dependency licenses, `pip check`, strict distribution metadata | Passed |
| Fresh full eight-stage smoke before the execution fix | 10,000 raw pairs; 9,569 retained pairs; 100 contigs; 106 genes; 95 eggNOG annotations; two curated protein-evidence rows; eight reports and 36 figure files; accounting checker passed |
| Fixed wheel resumes that full run | Passed; checked scientific tables retained their original SHA-256 hashes |
| Fresh annotation through the fixed wheel | Three recorded positive-control proteins; native log and annotation headers both report `2.1.15`; exported annotation/category tables byte-identical to recorded controls |
| Fixed wheel: fresh single-end taxonomy | Passed accounting; 9,886 retained reads |
| Fixed wheel: fresh interleaved taxonomy | Passed accounting; 9,569 retained paired fragments |
| Fixed wheel: skip functional annotation and abundance, retain screening | Passed accounting; 106 predicted genes; annotation explicitly not run |
| Existing output, missing database, mismatched paired IDs | Failed as intended with nonzero exits |
| Native positive/negative screening and completed zero-hit reports | Passed in the fixed full suites; no database-dependent skips |
| Complete smoke reports in Chromium/Firefox | Another 54 offline views passed; representative narrow screenshots and A4 PDF taxonomy/screening pages visually inspected |

The smoke uses the first 10,000 records per mate of cached ENA `ERR2984773`.
Both complete raw files matched the MD5 checksums in the
[acceptance procedure](acceptance.md). Paired mates were checked explicitly;
interleaved data placed each R1 record immediately before its R2 partner.
The full run's 3.53% mapping rate limits abundance interpretation to the small
assembled fraction. Its original raw eggNOG headers remain preserved as failing
evidence; the fixed fresh native annotation check covers three controls rather
than a second fresh eight-stage run. The pre-fix wheel also failed both new
regressions; those failures were retained alongside the passing fixed suites.

GNU time recorded 12:24.93 elapsed and 2,476,512 KiB maximum RSS for the full
smoke, while the CLI reported 12m56s. Retain both timers. Concurrent verification
work and the previously recorded WSL clock discrepancy preclude a speedup or
causal timing claim from this run. This is runtime/accounting acceptance,
not a biological accuracy benchmark or an independent comparison.

Ignored evidence is retained under `test_results/release-b4/` (commands, exit
codes, input/artifact hashes, JUnit, native logs, manifests and browser artifacts)
and `test_results/taxonomy-evidence/`; it is not bundled into distributions.

## Residual limits

- The frozen full-input mock benchmark detected eight expected bacteria and
  five unexpected species: 8/13 (61.54%) nominal detection precision, below the
  90% project gate. The reused 200,000-pair subset scored 8/9 (88.89%); neither
  is held-out validation. Standard-8 excludes fungi.
- Assembly/reference disagreements and conflicting elapsed measurements remain
  unresolved. Reporting and smoke acceptance do not close checkpoint 3.
- Screening is a 16-species bacterial target subset plus selected curated
  protein families. Co-occurrence does not establish species attribution,
  pathogenicity or resistance phenotype; zero hits do not establish absence.
- Formal comparisons support two independent groups with CLR normalization
  and recorded filtering, pseudocount and BY/BH correction. They do not adjust
  for pairing, repeated subjects or batch effects.
- Independent installation/user acceptance and broader biological validation
  remain outstanding. This is a beta release, not stable 2.0.0.
