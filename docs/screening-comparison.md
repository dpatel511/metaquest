# Research screening and comparative analysis

These features are included in `2.0.0b1`. Before beta publication, install this
checkout with `python -m pip install -e .`; the older PyPI alpha does not contain them.

## Run the complete workflow with screening

Use the same database root for installation, checking and analysis:

```bash
metaquest databases --db-dir databases --database taxonomy
metaquest databases --db-dir databases --database functional
metaquest databases --db-dir databases --database amr
metaquest check --db-dir databases --screen-pathogens
metaquest run --db-dir databases \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --screen-pathogens --output results/sample
```

The repository Conda environment includes `ncbi-amrfinderplus`. Existing
environments can install it with
`conda install -c conda-forge -c bioconda 'ncbi-amrfinderplus>=4.2.7,<5'`.
Screening is optional by default. It needs species-level classification (`S`)
and gene prediction, so it cannot be combined with `--taxonomy-only`.
It can run with `--skip-functional` or `--skip-abundance`; unavailable abundance
is left blank rather than reported as zero.

For an existing AMRFinderPlus database or a custom taxon target table:

```yaml
screening:
  enabled: true
  database_dir: /data/amrfinder/2026-08-07.1
  taxa_catalog: /data/study/species_targets.tsv
  threads: 4
```

Custom catalogs require unique numeric `taxid`, `name`, `source` and `note`
columns in a tab-separated file. Bracken IDs are matched exactly at species
rank; genus/name similarity does not create a hit. The bundled table is a
16-species human-associated bacterial research subset informed by the
[WHO bacterial priority pathogen list](https://www.who.int/publications/i/item/9789240093461).
It is neither a complete pathogen catalog nor a classifier for resistant strains.
WHO's resistance-specific priority categories are not assigned to species hits.

`screening/taxon_candidates.tsv` reports estimated species counts and fractions
of all cleaned input reads/fragments. `screening/protein_evidence.tsv` preserves
native AMRFinderPlus types, methods, reference identity/coverage and gene IDs.
Gene coordinates are **1-based inclusive**, matching the Pyrodigal GFF3. Native
HMM values such as `NA` remain `NA`; they are not fabricated percent identities.
Gene abundance joins require complete, unique ID agreement.

AMRFinderPlus uses protein-only `--plus` mode and its curated thresholds. This
covers AMR plus selected virulence/stress families with organism-dependent
coverage. It does not assess nucleotide mutations, expression, phenotype, all
pathogens, or all virulence mechanisms. In a mixed metagenome a protein hit is
linked to its contig, not attributed to a species by co-occurrence. No hit does
not establish absence. See the upstream [output interpretation](https://github.com/ncbi/amr/wiki/Interpreting-results).

## Inspect the evidence

- `report.html`: consolidated report with links to every stage.
- `reports/00_validation.html` through `reports/07_screening.html`: stage reports
  with JSON and TSV counterparts, explicit status, metrics and relevant figures.
- `plots/`: PNG/SVG figures (or configured formats); `plots/plot_data/`: plotted data.
- `analysis_metadata.json`: inputs/checksums, settings, tool/database versions and timings.
- `screening/summary.json`, `amrfinder.tsv`, `amrfinder.log`: screening provenance
  and native results. The species and protein evidence tables are independent.

Skipped, not run, and zero detected hits are distinct states. If input validation
was bypassed, its report says so. Sequence evidence is for research interpretation,
not a clinical diagnosis or a calculated pathogen-risk score.

## Descriptive sample comparison

Create `samples.tsv`; paths are relative to the sheet. Each directory must be a
completed MetaQuest run with compatible recorded database/tool/settings provenance.

```text
sample_id	result_dir	group	subject_id
C01	results/C01	control	person01
C02	results/C02	control	person02
C03	results/C03	control	person03
T01	results/T01	case	person04
T02	results/T02	case	person05
T03	results/T03	case	person06
```

Use real tabs in a TSV or a comma-separated CSV. `group` and `subject_id` are
optional for descriptive comparisons. Leading zeros in sample IDs are preserved.

```bash
metaquest compare --sample-sheet samples.tsv --feature taxonomy \
  --group-column group --output comparisons/taxonomy
metaquest compare --sample-sheet samples.tsv --feature KO \
  --group-column group --output comparisons/KO
```

Features are aligned by taxonomy ID or functional term ID using full source
profiles, never a report's top-N display. Functional inputs are summed gene TPM
from the BBMap abundance stage. Namespace choices are `COG`, `KO`, `EC`, `GO`.
Raw values, relative profiles, feature labels, sample QC, Bray–Curtis distances,
offline HTML and source-data plots are exported. Relative profiles sum to one
within reported features; unclassified/unassembled material is excluded.
They are compositional, not absolute organism loads or gene copy numbers.

## Exploratory tests for two independent groups

```bash
metaquest compare --sample-sheet samples.tsv --feature taxonomy \
  --group-column group --test clr-mannwhitney \
  --reference-group control --pseudocount 0.000001 \
  --fdr-method by --alpha 0.05 --min-prevalence 0.1 \
  --output comparisons/taxonomy-groups
```

The pseudocount is an explicit user choice on the relative-proportion scale,
not a default biological detection limit. The method filters by overall
prevalence without using group labels, retains at least two features, adds the
pseudocount, then computes the centered log-ratio (natural log minus the mean
log across retained features for each sample). Each feature receives a
two-sided Mann–Whitney U test on these CLR values. Its null concerns equal
distributions; a median difference alone does not define that null.

`group_tests.tsv` records U, raw p, FDR-adjusted q, median CLR difference
(comparison minus reference), rank-biserial effect, prevalence and sample counts.
The default Benjamini–Yekutieli (`by`) correction allows arbitrary dependence
among tests; `bh` assumes independence or appropriate positive dependence.
The correction includes every tested feature, not just displayed features.
The report records the retained feature count, settings and SciPy version.

The initial design requires exactly two independent groups with at least three
samples per group. Three is an input minimum, not a power recommendation;
small groups have very limited power after correction. Duplicate input hashes
and repeated nonempty `subject_id` values are rejected for formal tests.
Users must still verify biological independence, batch balance, confounding and
replicate identity: different files do not prove different biological samples.
Paired, longitudinal, covariate-adjusted and multi-group designs are unsupported.

Treat these as exploratory log-ratio distribution tests. Check sensitivity to
pseudocount and prevalence, and inspect QC and effect sizes before interpreting
q values. Separate sample assemblies recover different gene catalogs, so apparent
functional changes can reflect assembly recovery. Full multi-term attribution
can count a gene in several terms. These tests do not establish absolute fold
changes, causality or clinical risk. See the primary
[Mann–Whitney documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html)
and [FDR correction documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html).
