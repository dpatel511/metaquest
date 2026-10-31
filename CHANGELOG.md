# Changelog

## Unreleased

## [2.0.0b3] — 2026-10-31

Corrective research beta for comparative-analysis safeguards. In a declared
small tied global-null model, beta 2 reached 8.21% FDR under BH at nominal 5%.
Beta 3 rejects nonconstant tied CLR features when either group has fewer than
10 samples; constant pooled features retain p=1 and untied studies retain the
three-sample input minimum. No permutation or batch-adjusted model is added.

Existing result files are preserved. Comparison inputs now require complete,
compatible model/preprocessing evidence; reprocess missing legacy evidence
rather than editing metadata. The synthetic audit does not close biological,
browser-review or independent-acceptance checkpoints.

- Preserve frozen benchmark protocol checksums across Windows/Linux checkouts
  using explicit LF attributes and byte-identity regressions.
- Audit independent-group comparisons: reject shared input files and repeated
  subjects hidden by whitespace, incompatible Bracken models/preprocessing,
  malformed profiles and nonempty comparison outputs; reconcile taxonomic
  counts and preserve explicit transformation/design settings in text reports.
- Add frozen null/effect simulations and independent CLR, U-test, FDR, ID-join
  and uncertainty checks. Retain small-sample/tie, compositional and design limits;
  public 2.0.0b2 artifacts are unchanged. Reject nonconstant tied CLR studies
  below 10 samples per group after an exact null calculation exposed small-sample
  asymptotic miscalibration; preserve the failed public-beta results.

## [2.0.0b2] — 2026-10-30

Reliability beta with CLI and recovery fixes plus a reproducible development
benchmark. Biological benchmarking remains incomplete; this release does not
claim improved classification accuracy. See the
[recorded results](benchmarks/zymo-even/RESULTS.md) for the failed primary taxonomy
gate, unresolved structural disagreements and timing discrepancy. Larger-database
testing is deferred because of available hardware.

Older outputs without the new database identities cannot be reused with strict
`--resume`. Preserve the original results and rerun into a new directory.

- Include tests, benchmark support, documentation and environment definitions in
  the source distribution; run its extracted test suite in CI.

- Start checkpoint 3 with a frozen full-mock benchmark protocol and independent
  taxonomy scorer, explicit paired-fragment denominators, database coverage and
  development-versus-held-out limits; add calculation and provenance regressions.
- Add a checksum/content-verified QUAST coverage summarizer and retain the first
  baseline's failed taxonomy gate, apparent structural disagreements and timer
  discrepancy alongside passing integration checks.

- Preserve database installations across failed/interrupted promotion, recover
  durable replacement backups, validate manifests and discard corrupt functional
  archives so retries can download them again.
- Bound HTTP operations and verify range bodies/completed partial downloads.
- Record interrupted attempts, check database/tool identity and essential output
  hashes before resume, and invalidate downstream caches after damaged stages.
- Terminate native process groups on SIGINT/SIGTERM and return exit 130/143;
  document recovery guarantees and the observed native MEGAHIT ordering variation.

- Route terminal diagnostics, validation errors, debug output and animation to
  stderr; retain help/version on stdout and fatal validation details under `--quiet`.
- Reject invalid validation thresholds and annotation thread counts before output
  creation; add subprocess CLI regressions exercised against source and wheels.
- Record the post-beta development checkpoints and clarify that TestPyPI is optional.

## [2.0.0b1] — 2026-10-28

Research beta prepared for publication. Screening and formal comparisons remain
limited to their documented scope; runtime checks do not establish biological accuracy.

- Mark the package and documentation as a beta, correct the issue-report link,
  and record MetaQuest's version in comparative summaries.
- Test isolated wheel installations and native AMRFinderPlus controls in CI;
  gate publishing on successful CI for the exact release commit and verify the
  exact TestPyPI package after upload.
- Add reproducible acceptance commands and a checker for result accounting,
  gene-ID joins, per-stage artifacts and embedded report links/images.
- Resolve core and screening native tools in separate environments to avoid
  eggNOG/DIAMOND and AMRFinderPlus curl/zlib conflicts. Handle the all-tied CLR
  case explicitly and reject other undefined statistical p-values.

- Add optional `--screen-pathogens` research screening with exact species-taxid
  targets and native AMRFinderPlus AMR/selected virulence/stress evidence;
  preserve gene/contig coordinates, method, provenance and interpretation limits.
- Install tool-compatible AMRFinderPlus databases through `metaquest databases
  --database amr`, with staged validation and recorded release/file hashes.
- Generate HTML, JSON and TSV reports for every workflow stage, gene-prediction
  and screening plots, and explicit skipped/not-run states in the consolidated report.
- Restore `compare` using complete ID-aligned taxonomy or functional-abundance
  profiles, descriptive plots, and explicit exploratory CLR/Mann–Whitney tests
  for two independent groups with BY/BH multiple-testing correction.
- Fix database setup elapsed-time reporting and diagnose database subdirectories
  supplied in place of the database root.
- Keep MEGAHIT scratch on the system temporary filesystem so compressed-read
  named pipes work when results are on Windows-mounted WSL drives.
- Honor configured annotation threads unless explicitly overridden on the CLI.
- Retain the complete MD5 in validation evidence rather than a truncated prefix.

- Count only aligned bases for gene depth and midpoint assignment, exclude
  alignments with unverifiable identity, and hash read contents before reusing
  gene abundance results.
- Add a Gene Abundance stage that maps quality-controlled reads to the
  stable-ID contigs with BBMap, assigns fragments to Pyrodigal CDS features,
  and writes per-gene counts, mean depth, RPK, and TPM
  (`gene_abundance/gene_abundance.tsv`).
- Aggregate gene TPM into COG, KO, EC, and GO functional abundance
  (`gene_abundance/functional_abundance.tsv`) with configurable full or split
  attribution for multi-term genes.
- Report mapping, gene-assignment, and annotated-TPM denominators in
  `analysis_summary.json`, `04_abundance_report.txt`, and the HTML report.
- Add `--skip-abundance` and the `abundance` configuration section.
- Continue interrupted database downloads across multiple validated HTTP range
  responses instead of rejecting each completed chunk as an unexpected size.
- Direct database help and missing-data diagnostics to MetaQuest's managed
  `databases` commands.
- Suggest the correct `eggnog-mapper` Conda package when `emapper.py` is missing.

All notable changes to MetaQuest are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## Version-history notice

The formal public Git history currently runs from `v1.0.0` through `v1.3.0`.
The former 3.x, 4.x, and 5.x numbers below were untagged prototype milestones,
not public semantic-version releases. They are retained for development
provenance and are explicitly labelled as prototypes.

MetaQuest 2.0 resumes semantic-version continuity from the public 1.x series.
Pre-release identifiers are used until the refactored workflow and scientific
claims are validated.

---

## [2.0.0a1] — 2026-08-22

### Release highlights

- Add fastp preprocessing before classification and assembly.
- Add resumable start metadata, full validation records, and stage timings.
- Add stable-ID contig FASTA aligned with the gene-prediction GFF3.
- Add N50/L50, N90/L90, GC, and contig-threshold assembly metrics.
- Add Matplotlib/Seaborn figures, plotted-data TSVs, and offline HTML reporting.
- Correct mean Phred quality display and paired-end fragment terminology.
- Infer read length and select the closest installed Bracken model.
- Report classified and unclassified denominators explicitly.
- Add safe `--resume` and recoverable `--force` output modes.
- Validate complete FASTQ structure and synchronized paired reads.
- Treat quality findings as warnings unless `--strict-validation` is used.
- Stabilize gene identifiers and functional-cache keys across record ordering.
- Keep eggNOG temporary files inside the selected result directory.
- Add assembly reporting and expanded reproducibility metadata.
- Remove configuration options not implemented by the stable alpha runtime.

### Changed

- Reset the stable runtime to a research-use-only, validation-first surface.
- Introduced `run` and `databases` as the primary command names while retaining
  `analyze` and `setup-db` as compatibility aliases.
- Excluded custom pathogen scoring, ML prediction, HMM, ESM, pathogenicity
  islands, Bayesian integration, and clinical-risk reporting from the default
  pipeline pending independent validation.
- Made taxonomy-only analysis skip assembly and annotation.
- Replaced active risk-oriented reports with descriptive text and JSON output.
- Made Bracken failure explicit instead of substituting a Kraken report with a
  different schema.
- Removed global `pkill` and automatic `tbl2asn` process termination.
- Made dependency checks conditional on the selected workflow.
- Replaced the provisional annotation path with Pyrodigal gene prediction and
  pinned eggNOG-mapper 2.1.15 functional annotation against eggNOG 5.0.2.
- Added complete per-gene annotation tables, explicit unannotated genes, COG,
  KO, EC, and GO summaries, and checksum-based restart-safe reuse.
- Added `--skip-functional` and `--taxonomy-only`; retained
  `--skip-annotation` as a deprecated taxonomy-only alias.

### Added

- Stable-runtime boundary tests.
- Explicit experimental-feature quarantine documentation.
- Clear research-use-only and non-clinical interpretation statements.
- Linux CI for Python 3.10-3.12, a complete Bioconda environment check, and
  wheel/source-distribution verification.
- TestPyPI and PyPI Trusted Publishing workflows with protected GitHub
  environments and release-tag/version validation.
- GPL-3.0-or-later licensing, third-party provenance, and a direct Python
  dependency license gate.

### Removed

- Unsupported perfect sensitivity/specificity and production-readiness claims
  from the primary project documentation.
- ML model artifacts from default package data.

---

## [Prototype 5.0.0 — untagged] — 2026-03-03

### Added
- **`reporting/validation_engine.py`** — Multi-evidence pathogen validation engine. Vectorized gene-to-species linkage (pandas merge, 5–10× faster). Bracken confidence scoring: assigns HIGH/MODERATE/LOW confidence to species based on ratio of direct Kraken hits vs. inferred reads. Reduces risk scores by 30% for low-confidence dominant taxa.
- **`io/data_loaders.py`** — Centralized loaders for annotation files, ML predictions, and protein sequences. Streaming FASTA loader pulls only sequences with database hits (~99% memory reduction).
- **`io/text_parsers.py`** — Shared text parsing utilities (`extract_organism_name`, Prokka GFF parsers) extracted from multiple modules.
- **`core/assembly.py`** — Dedicated assembly module wrapping MEGAHIT/SPAdes with contig filtering, coverage extraction, and quality metrics.
- **`config/pathogen_config.json`** — Pathogen configuration externalized from `config.py` into a versioned config directory.

### Changed
- **`core/comparative_analysis.py`**
  - Rarefaction deprecated: warns with McMurdie & Holmes (2014) citation, falls back to TSS
  - CLR validation added: rejects data with >95% sparsity or negative values
  - Power analysis corrected using Noether (1987) nonparametric formula
  - Reproducibility: `random_seed` parameter, saves `analysis_metadata.json`
  - Memory: explicit `del` + cleanup of intermediate DataFrames
  - OutputFormatter integration replaces bare `print()` calls
  - Raw counts preserved separately for alpha diversity (which requires counts, not normalized values)
- **`reporting/risk_scoring.py`**
  - Scientific model: z-score baseline comparison replaced with absolute functional density scoring (sample-agnostic)
  - `DENSITY_THRESHOLDS`: virulence 2%, AMR 3%, transposase 5% of proteome
  - Fixed divide-by-zero when `total_cds` is `None`, `0`, or non-numeric
  - `_analyze_pathogen_db_hits()` vectorized: pandas masks replace `iterrows()` (10–50× faster)
- **`reporting/taxonomy_reporter.py`**
  - Classification rate now derived from Kraken2 `unclassified` row (Bracken estimates always showed ~100%)
  - BSL lookup reads from `pathogen_config.json` (was hardcoded dict)
  - Confidence column (HIGH/MODERATE/LOW) added using `ValidationEngine` data
  - Non-bacterial reads section added with accurate fractions
- **`reporting/functional_reporter.py`**
  - Annotation counts clarified: DIAMOND reports multiple hits per protein
  - Section 2.1 now shows "Annotated CDS" (unique) vs "Total annotations"
  - Added `_generate_methodology_note()` explaining DIAMOND multi-match behavior
- **`reporting/pathogen_risk_reporter.py`**
  - Section 5 interpretation made sample-agnostic (valid for clinical, environmental, food samples)
- **`reporting/base_reporter.py`**
  - `VERSION` constant corrected: `"5.0.0"` (was `"5.0.1"`, embedded in all report footers)
- **`core/pathogen_analysis.py`**
  - 7-factor fragment-aware confidence scoring replaces simple identity/coverage filter
  - Factors: fragment completeness (sigmoid), identity×coverage, organism abundance (log₁₀), contig depth, critical motifs, Prokka cross-validation, signal peptide penalty
  - Tiers: HIGH ≥70, MEDIUM 45–69, LOW <45
  - Outputs: `pathogen_detections_validated.json` (full) + `pathogen_results.tsv` (legacy)
- **`cli.py`** — Version display reads from `config.__version__` dynamically; `--normalization` help updated to note rarefaction deprecation
- **`config.py`** — `__version__ = "5.0.0"`, `__release_date__ = "2026-03-03"`

### Fixed
- `setup.py`: `check_external_tools()` moved inside `if __name__ == "__main__"` guard (was running at import time)
- `setup.py`: `package_data` path corrected to `config/pathogen_config.json`
- `README.md`: All Table of Contents links fixed (were pointing to `google.com/search` URLs)
- `ml/feature_extractor.py`: Removed `✅` symbols from `print()` statements being emitted to stdout

### Removed
- Development comment markers (`✅`) cleaned from all source files

---

## [Prototype 4.1.0 — untagged] — 2025-11-01

### Added
- **Comparative Analysis v2.0**: Full compositional data normalization (TSS, CLR, rarefaction)
- **Statistical Power Assessment**: Cohort size analysis with detectable effect size estimates
- **Multi-group Alpha Diversity**: Kruskal-Wallis with automated pairwise post-hoc tests
- **Effect Size Metrics**: Cliff's Delta (differential abundance), epsilon-squared (alpha diversity), R² (PERMANOVA)
- **ML Biomarker Discovery**: SelectKBest feature selection + GridSearchCV hyperparameter tuning
- **Prevalence Filtering**: Filters rare species before differential abundance testing

### Changed
- **Pathogen DB v4.1**: Rebuilt to filter housekeeping genes (*rpoB*, *gyrA*, etc.)
  - Sensitivity: 33.3% → 100% | Specificity: 33.3% → 100%
  - Database reduced from 70,752 to 21,152 sequences (70% smaller, 2.6× faster)
- PERMANOVA results now include R² effect size
- ANOSIM results include group separability interpretation

---

## [Prototype 4.0.0 — untagged] — 2025-10-01

### Added
- Enhanced COG + SwissProt dual-database annotation
- ML model v2.1 with updated training data
- Advanced gene prediction controls (contig filtering, tbl2asn timeout)
- Professional logging system (`OutputFormatter`) with `--debug` mode

---

## [Prototype 3.6.0 — untagged] — 2025-10-01

### Added
- SPAdes metagenomic assembly support
- Modular visualization system with specialized plotters (`taxonomic_visualizer`, `functional_visualizer`, `pathogenic_visualizer`)

---

## [Prototype 3.5.0 — untagged] — 2025-08-01

### Added
- Enhanced statistical testing for comparative analysis
- Advanced ML biomarker discovery pipeline
- Comprehensive comparative visualizations

---

## [Prototype 3.3.0 — untagged] — 2025-08-01

### Added
- Comparative analysis pipeline (initial)
- Beta diversity analysis and PCoA visualization
- Alpha diversity visualization
