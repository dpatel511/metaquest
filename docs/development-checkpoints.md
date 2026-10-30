# Development checkpoints after 2.0.0b1

Keep the current beta scope fixed while improving reliability, evidence and
usability. A passing software test suite does not establish biological accuracy.
This plan defines evidence required before proposing a stable 2.0.0 release;
it does not declare the remaining checkpoints complete.

## Basis and attribution

Torsten Seemann's 2013 paper,
[Ten recommendations for creating usable bioinformatics command line software](https://pmc.ncbi.nlm.nih.gov/articles/PMC4076505/),
addresses CLI usability and explicitly sets algorithmic correctness aside.
Its recommendations become checkpoint 1 below. The other checkpoints are
MetaQuest-specific acceptance criteria, not a validation standard authored or
endorsed by Seemann.

For broader guidance on standard formats, provenance, demonstration data and
user feedback, see List, Ebert and Albrecht's
[Ten Simple Rules for Developing Usable Software in Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005265).

## Recorded baseline

- Released package: `metaquest-bio==2.0.0b1`.
- Prior release evidence: all six required CI jobs passed; the downloaded
  production wheel passed 95 tests outside the source checkout. See
  [release CI](https://github.com/dpatel511/metaquest/actions/runs/40753444237)
  and [publishing workflow](https://github.com/dpatel511/metaquest/actions/runs/40753931047).
- Native smoke runs covered the main workflow, optional modes and resume;
  report structure, data joins and selected PNG plots were inspected.
- Those runs were not truth-based biological benchmarks. Full browser layout
  review and independent user acceptance remain outstanding.
- On 2026-10-29, 23 subprocess CLI probes against the cached, previously
  hash-verified production wheel produced 20 passing cases and three failures.
  All expected exit codes passed. Missing input, quiet missing input and no
  database selection wrote their errors to stdout instead of stderr.
  The isolated wheel shared the existing environment's dependency packages;
  this was not a new clean-install test of all dependencies.
- Local audit evidence is in ignored `test_results/seemann-audit/cli-audit.json`
  with separate stdout/stderr captures and the probe script. It is not included
  in the distribution. The source formatter also uses stdout for progress.

## 1. Command-line contract: first fix checkpoint

Seemann's recommendations, paraphrased and mapped to current evidence:

| Recommendation | Current evidence / required check |
| --- | --- |
| Useful response without arguments | Usage text observed; retain a subprocess regression. |
| Standard help options | Root and six subcommands passed both `-h` and `--help`. |
| Readable version option | Both `-v` and `--version` passed. |
| Diagnostics on stderr | Fails three error probes; progress also needs routing. |
| Nonzero exit on failure | Six probed failures returned the expected nonzero codes; expand coverage. |
| Validate arguments | Missing/unknown arguments tested; numerical boundaries need audit. |
| Avoid fixed installation paths | Configurable paths exist; relocation still needs testing. |
| Keep helper names off PATH | Packaging defines one MetaQuest console command; check installed files. |
| Check dependencies early | Preflight exists; verify missing tools fail before expensive work. |
| Provide launchers for Java tools | No direct MetaQuest JAR distribution; inspect wrapped-tool failures. |

Local implementation on `fix/cli-contract` (2026-10-29): terminal diagnostics,
validation errors, debug output and spinner activity now use stderr. Fatal
validation messages survive quiet mode; invalid validation thresholds and thread
counts fail before output creation. This is a local candidate, not an update to
the already published beta.

Verification added 44 CLI/formatter cases. The source suite passed 136 tests
with three native tests skipped and 276 deprecation warnings in the older
developer environment. A clean installed-wheel environment with current
dependencies and existing native fixtures passed all 139 tests, with no skips
or warnings. Lint, dependency consistency, strict wheel metadata and direct
dependency license checks passed. Evidence is retained locally in ignored
`test_results/seemann-audit/`. CI includes these new tests through its existing
source and isolated-wheel jobs; remote CI has not been run for this branch.

Interrupt/child-process recovery remains to be exercised with checkpoint 2;
these CLI checks alone do not close every recovery or biological-validation gate.

Pass criteria:

- Help/version output stays on stdout; diagnostics, banners, progress, warnings
  and errors go to stderr. File-producing commands must not pollute stdout.
- `--quiet` suppresses progress but preserves actionable errors on stderr;
  `--no-color` and redirected streams contain no ANSI control sequences.
- Add subprocess regressions to CI for these observable behaviours, including
  invalid values and paths containing spaces. Test the installed wheel too.
- Missing input/config/tool/database and malformed FASTQ fail clearly without
  misleading success text or unnecessary expensive work. Interrupts have a
  defined nonzero exit status and leave recoverable outputs.

## 2. Failure recovery and reproducibility

Implementation and observed recovery behaviour are documented in
[download and run recovery](recovery.md). The audit retains the failed initial
MEGAHIT byte comparison: all sequence content matched, while ordering/headers
differed. SIGKILL/power-loss cleanup and concurrent writers are explicitly outside
the handled-signal guarantees.

Local verification on 2026-10-29 added 34 recovery/process cases. The complete
source and independently installed-wheel suites each passed 173 tests with no
skips or warnings, including the existing native controls. Lint, strict wheel
metadata, dependency consistency and license checks passed; 40 packaged source,
config and catalog files matched the checkout. Remote CI and publication remain
pending for this checkpoint branch. Evidence is retained under ignored
`test_results/recovery-check/`.

Verification: use disposable output/database directories to exercise interrupted
downloads, truncated archives, insufficient space and interrupted stages.
Compare a resumed run to an uninterrupted run using the same input hashes,
database releases, effective settings and native tool versions.

Pass criteria:

- Failed installation cannot replace an existing valid database; partial files
  and integrity failures are clearly distinguished from installed data.
- Changed input/config/database invalidates affected cached work. An interrupted
  stage is never mistaken for a completed stage.
- Scientific tables agree within declared numerical tolerances. Compare stable
  fields rather than timestamps, paths or image-file bytes. Record any native
  tool nondeterminism rather than demanding unjustified byte equality.
- Provenance records input checksums, commands/settings, package/tool versions,
  database identities and attempts. The documented two-environment installation
  works from a new directory without relying on the developer's shell setup.

## 3. Biological benchmarking

The first development-data baseline follows the frozen
[benchmark protocol](benchmarking.md). Keep its bacterial detection gate separate
from whole-community scope, held-out validation and the remaining accuracy gates.
The [full-data baseline results](../benchmarks/zymo-even/RESULTS.md) retain a failed
primary taxonomy gate, structural reference disagreements and a timing discrepancy;
checkpoint 3 remains open despite passing software and accounting checks.

Verification: choose public mock-community data with documented composition and
compatible ground truth. Record accessions, checksums, truth provenance and
licenses. Freeze thresholds and acceptance limits before scoring held-out data.
Keep the small runtime fixture separate from an adequately covered benchmark.

Pass criteria:

- Taxonomy: report species detection precision/recall and false positives at
  declared abundance thresholds; include unclassified reads and abundance error
  against compatible truth. Do not equate cell, DNA-mass and read proportions.
- Assembly/gene/functional stages: use available reference truth and independent
  checks to assess errors, completeness and assignments. Contig count, N50 and
  annotation rate alone cannot establish correctness. Where truth is unavailable,
  label consistency checks as such.
- Compare MetaQuest's parsed outputs with direct runs of its component tools
  using identical data/settings; distinguish integration correctness from the
  accuracy of the underlying algorithms and databases.
- Keep existing curated-screening regression controls and report their scope.
  Zero reported matches cannot establish biological absence or clinical safety.
- Record wall time, peak memory, disk use, hardware, threads and tool/database
  versions. Set resource regression budgets from a measured baseline.
- Publish expected, observed and discrepant results, including failures. Do not
  tune parameters on the final evaluation data and call it independent validation.

## 4. Comparative-analysis validity

Verification: retain deterministic calculation fixtures and add independent
null/effect simulations for the supported two-independent-group design.

Pass criteria:

- Check sample/feature ID overlap and row counts explicitly. Duplicate samples,
  repeated subjects, invalid quantities and incompatible provenance must fail.
- State normalization, prevalence filtering, pseudocount, CLR definition,
  reference direction and multiple-testing method in outputs.
- Use simulation with recorded seeds and uncertainty to assess null behaviour
  and effect recovery; define tolerances before judging results. Synthetic
  calculation fixtures are not biological replicates.
- Document independence/exchangeability assumptions and the lack of covariate
  adjustment. A batch-confounded comparison needs a different study design/model;
  the current exploratory test must not silently claim to adjust for it.

## 5. Report and data-contract review

Verification: run the [acceptance checker](acceptance.md), then inspect actual
HTML in Firefox and Chromium at desktop/narrow widths and print preview. Retain
screenshots and browser versions. Browser rendering is still an open checkpoint.

Pass criteria:

- Eight stage reports and the consolidated report agree with source tables and
  plot-data TSVs on counts, denominators, units and status.
- Empty results, skipped stages and stages not requested are distinguishable;
  labels/tables remain readable and report images/links work offline.
- Document stable TSV/JSON column names, units, null semantics and identifier
  relationships. GFF3 and derived genomic evidence use **1-based inclusive**
  coordinates; specify conversion explicitly wherever another format is used.
- Keep representative expected outputs under regression checks without brittle
  timestamp/image-byte comparisons. Define compatibility/migration expectations
  before changing an output consumed by users or `compare`/`--resume`.

## 6. Independent acceptance and release decision

Proposed project criterion: at least two users outside the development setup
complete installation and the documented example on fresh supported Linux/WSL
environments, recording package/tool versions, commands, time and blockers.
This is a usability checkpoint, not a statistically powered validation study.

Pass criteria:

- Installation, database setup, example run and report interpretation require no
  undocumented manual repair. Repeat the affected path after each blocker fix.
- Each bug gets a minimal reproduction, regression test, focused PR and changelog
  entry. Each PR records what passed and what remains untested.
- Keep the Python matrix, license, isolated-wheel and native integration checks
  mandatory. Large benchmark runs can be manual release checks with retained
  artifacts; avoid forcing multi-GB downloads into every small PR.
- Use another beta for fixes while these checkpoints remain open. Consider
  stable 2.0.0 only after recording evidence for every gate and disclosing residual
  limitations. A passing test count alone is insufficient.
- Reuse production Trusted Publishing. TestPyPI remains an optional rehearsal;
  the [release guide](releasing.md) describes the established publishing route.

Immediate order: CLI contract, recovery, biological benchmarks, comparison
validation, browser/data-contract review, then independent acceptance. Some
verification can run in parallel, but no checkpoint is passed by planning it.
