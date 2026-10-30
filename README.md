# MetaQuest

MetaQuest is a research-use command-line pipeline for short-read metagenomic
FASTQ analysis. The `2.0.0b2` research beta provides a
maintained workflow for profiling, assembly, annotation, abundance and reporting.
It does not make clinical or pathogen-risk claims.

Screening and restored comparison are included in the beta. Screening covers
16 human-associated bacterial species targets and curated, organism-dependent
AMRFinderPlus protein families. Formal two-group comparisons are exploratory.
The older `2.0.0a1` package does not contain these additions. Beta 2 adds CLI and
recovery fixes and records the first full development benchmark.

## Current capabilities

- full FASTQ structural validation and sampled quality warnings
- fastp adapter and quality preprocessing with retained-read QC
- synchronized paired-read identifier and count validation
- Kraken2 taxonomic classification
- Bracken abundance re-estimation
- optional taxonomy-only execution
- MEGAHIT metagenomic assembly
- Pyrodigal metagenomic gene prediction
- eggNOG-mapper orthology-based functional annotation
- per-gene annotations and aggregated COG, KO, EC, and GO counts
- read-mapping gene abundance (BBMap, TPM) with TPM-weighted COG, KO, EC, and
  GO functional abundance
- descriptive text, JSON, offline HTML, and publication figure reporting
- explicit classified/unclassified denominators and reproducibility metadata
- versioned taxonomy and eggNOG database installation
- optional species-target and AMRFinderPlus protein evidence screening
- per-stage HTML, JSON and TSV reports with explicit skipped/not-run states
- aligned multi-sample profiles and exploratory comparisons of two independent groups

MetaQuest uses fastp for read trimming. It does not currently remove host reads,
diagnose infection, establish resistance phenotypes, or recommend treatment.
Screening reports curated sequence evidence with database provenance; its
biological sensitivity and specificity require independent validation.

## Pipeline

~~~text
FASTQ
  ├── validation
  ├── fastp preprocessing
  ├── Kraken2 → Bracken
  ├── MEGAHIT
  ├── Pyrodigal
  ├── eggNOG-mapper → COG / KO / EC / GO
  ├── BBMap read mapping → gene TPM → functional abundance
  ├── optional taxon targets + AMRFinderPlus protein evidence
  └── per-stage + consolidated reporting → HTML + figures + plotted-data TSV
~~~

The maintained pipeline currently runs directly in Python. Migration to
Snakemake and optional Bowtie2 host filtering remain planned.

## Quick start

Install the exact Python distribution from PyPI:

~~~bash
python -m pip install --upgrade --pre 'metaquest-bio==2.0.0b2'
metaquest --version
~~~

The distribution is named `metaquest-bio`; the installed Python package and
command remain `metaquest`. The PyPI distribution does not bundle Kraken2,
Bracken, MEGAHIT, DIAMOND, or eggNOG-mapper. For the complete runtime, clone
the repository and create its Conda environment before installing the package:

~~~bash
git clone --branch v2.0.0b2 https://github.com/dpatel511/metaquest.git
cd metaquest
conda env create -f environment/environment.yml
conda activate metaquest
python -m pip install --no-deps 'metaquest-bio==2.0.0b2'
~~~

For optional screening, also install `environment/screening.yml` separately and
append its executables to PATH using the [installation guide](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/installation.md).
It isolates AMRFinderPlus's native libraries from eggNOG's older DIAMOND stack.

Inspect and install the taxonomy database:

~~~bash
metaquest databases --list
metaquest databases --database taxonomy
metaquest databases --database functional
metaquest check
~~~

By default, MetaQuest installs reference data in `./databases`. Use
`--db-dir` or `METAQUEST_DB_DIR` only when you intentionally want another
location.

Run taxonomic profiling:

~~~bash
metaquest run \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --output results/sample
~~~

Use `--skip-functional` to skip eggNOG annotation while retaining gene prediction
and gene abundance, `--skip-abundance` to skip
read mapping and gene abundance, or `--taxonomy-only` to skip assembly, gene
prediction, functional annotation, and abundance. The former
`--skip-annotation` spelling remains as a deprecated alias for
`--taxonomy-only`.

Each fresh run requires a new or empty output directory. Use `--resume` only
to reuse a matching MetaQuest run, or `--force` to move an existing output
directory to a timestamped backup before starting fresh.

## Commands

~~~text
metaquest check       # verify tools and databases
metaquest validate    # validate FASTQ input
metaquest run         # run the analysis pipeline
metaquest compare     # align profiles; optionally test two independent groups
metaquest databases   # inspect or install reference data
metaquest init-config # create a configuration file
~~~

The previous `analyze` and `setup-db` names remain available as aliases for
backward compatibility.

Global output controls may appear before or after the command:

~~~bash
metaquest run --verbose --single reads.fastq.gz --output results/
metaquest databases --no-color --list
metaquest --quiet check
~~~

For systems that cannot load the Kraken2 database fully into RAM, enable
Kraken2 memory mapping with `--low-memory`:

~~~bash
metaquest run --low-memory --single reads.fastq.gz --output results/
~~~

Currently this flag only adds Kraken2's `--memory-mapping` option. It does not
change resource settings for Bracken, MEGAHIT, or gene prediction.

## Database storage

Reference databases are not stored in Git. From the repository root, install
them directly into the default `databases/` directory:

~~~bash
metaquest databases --database taxonomy
metaquest databases --database functional
metaquest check
~~~

The database path resolution order is:

1. `--db-dir`
2. `METAQUEST_DB_DIR`
3. `databases.base_dir` in YAML
4. `./databases`

## Documentation

- [Installation](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/installation.md)
- [Database management](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/databases.md)
- [Usage](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/usage.md)
- [Gene prediction and functional annotation](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/annotation.md)
- [Research screening and comparative analysis](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/screening-comparison.md)
- [Release process](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/releasing.md)
- [Reproducible acceptance checks](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/acceptance.md)
- [Biological benchmark protocol](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/benchmarking.md)
- [Development checkpoints](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/development-checkpoints.md)
- [Changelog](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/CHANGELOG.md)

## Validation status

The complete ERR2984773 development mock reproduced direct Kraken2, Bracken and
Pyrodigal outputs. At the frozen 0.1% detection threshold, all eight expected
bacterial members were detected alongside five unexpected species calls:
8/13 (61.54%) detection precision, below the proposed 90% project gate. This is
not per-read classification accuracy. The Standard-8 database excludes fungi.

The benchmark also retains unresolved assembly/reference disagreements and
conflicting elapsed times. It skipped functional annotation and screening.
See the [full benchmark results](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/benchmarks/zymo-even/RESULTS.md) and
[protocol](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/benchmarking.md) before interpreting these results.

MetaQuest remains beta software. Biological benchmarking, comparative-method
validation, browser/data-contract review and independent user acceptance remain
open in the [development checkpoints](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/development-checkpoints.md). Larger
database comparisons are deferred; passing software tests does not close these
gates. Older outputs lacking current database identities require a new run
instead of strict `--resume`; see [recovery](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/docs/recovery.md).

## License

MetaQuest is distributed under the GNU General Public License v3.0 or later.
Copyright (c) 2026 Dev Patel. External programs and reference databases are not
included in the Python distribution and remain subject to their own terms. See
[Third-party licenses](https://github.com/dpatel511/metaquest/blob/v2.0.0b2/THIRD_PARTY_LICENSES.md) for the maintained dependency
and provenance inventory.
