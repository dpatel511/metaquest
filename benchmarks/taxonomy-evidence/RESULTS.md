# Taxonomy evidence validation — 2026-11-22

Candidate source version: **2.0.0b4 (unreleased)**, based on main
`757e387d6a107f161ca4afc872d81280b2ef80fa`. The [protocol](protocol.json) fixes
input hashes, native settings, reporting repetitions and the overhead budget.
Evidence is retained in ignored `test_results/taxonomy-evidence/`; it is not
bundled into distributions or substituted for the original Zymo benchmark.

## Native comparison

Both revisions processed the same cached **200,000 paired fragments** before
filtering (44,574,567 compressed bytes across the two FASTQs), using four threads,
memory-mapped Kraken2 and Standard-8 release 2026-06-26. The 151 bp observed read
length selected the same installed 150 bp Bracken model, with the mismatch
warning retained. No data/database download or default/threshold tuning occurred.
Native banners reported Kraken 2.1.3, Bracken 3.0.1 and fastp 1.4.0. Conda package
identities were Kraken2 `2.1.3-pl5321hdcf5f25_3`, Bracken `3.1-h9948957_0` and
fastp `1.4.0-h43da1c4_0`; the Bracken package and its banner differ, so both are
retained rather than silently treating them as interchangeable version strings.

| Observation | Baseline and candidate |
| --- | ---: |
| Surviving paired fragments | 193,504 |
| Kraken-classified fragments | 172,108 |
| Unclassified fragments | 21,396 |
| Bracken-estimated fragments in the reported species profile | 171,754 |
| Classified fragments outside that profile | 354 |
| Reported taxa before display thresholds | 23 |
| Kraken report / Bracken table comparison | Both byte-identical |
| Independent accounting, source-hash and offline-link checks | Passed |

Additional native Bracken component controls reused the same Kraken report at
genus/family ranks. The evidence joins and count reconciliation passed for all
14 genus rows (171,988 estimated fragments) and 10 family rows (172,034 estimated
fragments), respectively. These were rank-output controls, not additional
independent classification or biological benchmark runs.

The new evidence makes differences such as **Bacillus subtilis** (taxid 1423)
visible: 30 Kraken direct assignments, 31 clade assignments, 483 Bracken-added
fragments and a 514-fragment estimate. These counts have different meanings;
redistribution is not independent observation or a confidence score.

At the existing 0.1%-of-all-cleaned-fragments threshold, this reused subset
detected all eight expected bacteria and one unexpected species call, giving
nominal precision **8/9 = 88.89%**, below the inherited 90% criterion. Neither
yeast was detected; this database excludes fungi. These are descriptive
nominal-composition results, not independent accuracy estimates or evidence of
absence. No cause is assigned to the unexpected call. The reporting feature
does not change this result.

This subset does not replace the frozen full-input benchmark: its **8/13 =
61.54% precision failure**, missing fungal coverage and other limitations remain
recorded in [the original results](../zymo-even/RESULTS.md).

## Reporting resource measurements

`tools/benchmark_taxonomy_reporting.py` generated reports three times per revision
from the **same completed baseline outputs, effective settings and native source
hashes**. The candidate reports retain the original classifier/tool provenance
and separately record their beta-4 report-generator version. Elapsed and CPU
timers are recorded independently. RSS is Linux KiB and is a process-lifetime
high-water mark, not incremental memory per repetition.

| Revision | Elapsed seconds, repetitions 1 / 2 / 3 | CPU seconds, repetitions 1 / 2 / 3 | Median elapsed |
| --- | --- | --- | ---: |
| Main / beta 3 | 3.970 / 2.111 / 2.309 | 2.673 / 1.338 / 1.426 | 2.309 s |
| Beta-4 candidate | 2.956 / 1.727 / 2.051 | 2.104 / 1.044 / 1.170 | 2.051 s |

The predeclared budget was `max(2 × baseline median, baseline median + 1 s)`:
**4.619 seconds**, and the candidate remained within it. Process-lifetime maximum
RSS reached 217,176 KiB for the baseline and 217,832 KiB for the candidate.
These short, sequential shared-host measurements do **not** establish a speedup.

Native runs were executed once each, so their resource logs are context, not a
performance-regression benchmark. GNU time recorded 104.10 s / 86.31 s and peak
RSS 6,769,716 / 6,768,044 KiB for baseline/candidate. The CLI separately printed
1m 47s / 1m 30s. Both timer records are retained; the previously noted WSL
elapsed-timer discrepancy and WSL clock-jump possibility remain unresolved. No speedup or causal attribution
is made from these conflicting elapsed measurements.

Classifier/preprocessing source content is unchanged in canonical Git text.
The archive/worktree source files differ in CRLF/LF line endings; raw hashes and
the explicit canonical comparison are retained. Biological input/output hashes
were never normalized or edited.

## Software and rendering checks

- The initial feature regressions failed before implementation; the final suite
  includes ID/denominator fixtures, malformed inputs, missing support, legacy
  data, zero counts, non-species ranks and independently detected tampering.
- The full source suite passed **362 tests**, with two AMRFinderPlus database
  controls skipped and 276 existing plotting dependency deprecation warnings.
- The installed-wheel suite passed **362 tests**, with the same two skips and
  warnings. Imported package paths and version metadata were checked outside the
  source tree in an isolated environment sharing the verified public Conda
  dependency packages. Distribution metadata and `pip check` passed.
- The extracted-sdist suite also passed **362 tests**, with the same two skips
  and 276 warnings. Source, wheel and sdist checks exercised the same cases; these
  are not three independent biological validations.
- The initial unsupported `python -m metaquest` benchmark launch and the initial
  raw-source hash comparison failure were retained, investigated and corrected;
  the actual console launcher and explicit line-ending accounting were used.
- Actual stage and consolidated HTML were inspected in the in-app browser at
  desktop and 390×844 viewport settings. A narrow-layout table compression bug
  was fixed and rechecked; horizontal scrolling and methods navigation worked
  without body overflow. Four JPEG screenshots are retained locally.

The first source test invocation lacked the native tool PATH and skipped BBMap;
the final source/wheel runs used explicit executable paths and exercised that
control. Screening database controls, native print pagination, a second browser
engine and held-out biological validation remain outstanding. No tag, release,
package publication or stable-release claim is made by this local preparation.

## Reproduce the report measurements

With the original native run and frozen baseline available locally:

```bash
# Run from outside the source tree; use the same Conda dependencies for both.
PYTHONPATH=/path/to/frozen-main/src /path/to/baseline-python \
  /path/to/current/tools/benchmark_taxonomy_reporting.py /path/to/native-run /tmp/report-before
/path/to/beta4-wheel-python /path/to/current/tools/benchmark_taxonomy_reporting.py \
  /path/to/native-run /tmp/report-after
python tools/verify_run.py /tmp/report-after/repeat-1
```

Use fresh output directories. Native comparison commands and GNU time logs are
retained with the ignored evidence. The complete raw ENA inputs remain governed
by the original [biological benchmark protocol](../../docs/benchmarking.md).
