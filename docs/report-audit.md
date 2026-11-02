# Checkpoint 5 local report audit — 2026-11-02

Branch: `fix/report-contract-audit`, based on published b3. No version bump,
publication, new biological dataset or database download is part of this audit.

## Findings and changes

1. Taxonomy top-feature summaries lost taxids; identical names could become the
   same seaborn category, averaging distinct counts into one bar. Preserve the
   source ID and disambiguate repeated names and reserved remainder labels.
2. Comparison composition exports used display names as columns, losing identity
   when names repeated. Add `plots/composition_by_id.tsv`; retain the existing
   display export and distinguish repeated labels in composition, heatmap and
   group-effect plots. Full primary matrices and statistical calculations remain
   unchanged.
3. Actual narrow-screen inspection found status words split across lines and
   figures reduced until labels were difficult to read. Keep workflow statuses
   on one line and allow full-size figures to scroll horizontally within their
   containers. Desktop/print figures continue to scale to their containers.

Two identifier regressions failed before their fixes, then passed. A third
regression covers skipped/not-run/completed-zero states, typed nulls and escaped
HTML. Tests use deterministic numerical values, not image hashes. See
[report data contracts](report-data-contracts.md) for units, aliases, schemas,
null handling, ID joins and additive migration expectations.

## Evidence

- Cached b3 taxonomy smoke: acceptance checker passed; 950 cleaned fragments,
  eight stage reports and ten figure files. This is accounting/runtime evidence,
  not a biological accuracy estimate. Original outputs were preserved.
- Baseline browser probes: 60 views across the cached b3 smoke and formal
  comparison report. Updated report probes: 168 views covering consolidated and
  per-stage smoke reports, duplicate-taxid/name data, comparison and edge-state
  fixtures. A separate 54-view run covers completed zero-hit screening.
- Engines: Chromium 156.0.8078.4 and Firefox 157.0, Playwright 1.64.1. Viewports:
  1440 px desktop, 375 px narrow, 794 px print-media; 950 px viewport height.
  Browsers ran headlessly with networking disabled. No page-level horizontal
  overflow, broken embedded images, external assets or script errors were found.
  Final reusable-tool rerun: all 168 views passed, including local companion-link checks. Figures intentionally scroll within narrow-screen containers.
- Representative narrow, print-media and comparison screenshots were visually
  inspected. Chromium A4 PDFs were generated; Firefox print-media screenshots
  were retained. Native print dialogs and physical-printer pagination were not
  exercised. The older cached full-workflow HTML was also rendered separately;
  that does not substitute for a new full-workflow biological run.
- Final source and isolated installed-wheel suites each passed 287 tests with
  two AMRFinderPlus database-dependent skips. The BBMap native control ran.
  Wheel tests imported from their own site-packages outside the checkout;
  dependency packages were shared with the existing Conda environment, so this
  is not a new clean-dependency installation claim.
- Ruff and wheel/source distribution metadata checks passed. Existing plotting
  dependencies emitted 276 Matplotlib/pyparsing deprecation warnings per suite;
  warnings were recorded, not suppressed or repaired by unrelated dependency
  changes.

Pre-commit verification repeated the source/wheel and 168-view browser checks.
The existing native Conda environment's `pip check` reports eggNOG-mapper 2.1.15
pins for Biopython 1.76, psutil 5.7.0 and xlsxwriter 1.4.3 conflicting with the
installed 1.88, 7.2.2 and 3.2.9 versions. Those native-environment constraints
were not changed by this report patch. A separate fresh wheel environment with
no system-site packages passes `pip check` and all 287 tests (two native database-dependent skips), with no warnings reported. See the additional clean-wheel test
evidence in `test_results/checkpoint5/clean-wheel-tests.log`.

Ignored evidence: `test_results/checkpoint5/` contains per-browser screenshots,
Chromium PDFs, JSON probe records, test logs/JUnit XML, builds and original/fixed
fixture outputs. These are local artifacts, not bundled biological benchmarks.

## Repeating browser probes

In a development environment, install the optional test driver and browser
runtimes (these are not MetaQuest runtime requirements):

```bash
python -m pip install playwright
python -m playwright install chromium firefox
python tools/verify_report_browser.py /tmp/report-evidence RESULTS_DIR [MORE_RESULTS_DIRS ...]
```

The command records versions, layouts, missing local companion links, embedded
images, external assets and script errors, and exits nonzero on probe failures.
Inspect its screenshots separately; passing geometry probes alone missed the
original narrow-screen readability issues.

This audit covers the stated fixtures and cached reports. It does not close
checkpoint 3, independent user acceptance, all possible cohort sizes or the
native printer/pagination review. Do not describe these as universal biological
or rendering guarantees.
