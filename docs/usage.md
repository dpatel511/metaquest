# Command-line usage

Run `metaquest COMMAND --help` for the authoritative option list.

## Global options

Global options may appear before or after the command:

~~~bash
metaquest --version
metaquest run --config metaquest.yaml ...
metaquest run --verbose ...
metaquest run --debug ...
metaquest --quiet check
metaquest databases --no-color --list
~~~

| Option | Meaning |
|---|---|
| `--config FILE` | YAML configuration override |
| `--verbose` | additional progress context |
| `--debug` | commands and diagnostic output |
| `--quiet` | errors only |
| `--no-color` | disable ANSI colors and animation |
| `--low-memory` | use Kraken2 database memory mapping |

Help and version output use stdout. Terminal status, summaries, progress,
warnings and errors use stderr; redirect them with `2>metaquest.stderr.log`.
`--quiet` suppresses progress while preserving fatal errors. Colour and animation
follow stderr's terminal status, and `--no-color` disables both explicitly.

## Check

~~~bash
metaquest check --db-dir /data/metaquest-db
~~~

`check` validates required executables, Python packages, and database files.
The full check validates the installed eggNOG core database.

## Validate FASTQ

~~~bash
metaquest validate --single reads.fastq.gz

metaquest validate \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz

metaquest validate --interleaved sample.fastq.gz
~~~

Optional thresholds:

~~~text
--min-quality SCORE
--min-sequences NUM
--overrep-threshold FRACTION
--strict-validation
~~~

`--min-quality` must be nonnegative, `--min-sequences` must be at least 1,
and `--overrep-threshold` must be a finite fraction in `[0, 1]`.
An explicit `--annotation-threads` value on `run` must be at least 1.

FASTQ structure is checked across the complete input. Quality statistics are
sampled. Low mean quality, low Q30 content, adapters, and overrepresented
sequences are warnings by default; `--strict-validation` makes them fatal.
Paired files must have equal counts and synchronized read identifiers.
Validation itself does not trim reads. `run` subsequently uses fastp, and its
cleaned reads feed both Kraken2 and MEGAHIT.

## Run

### Taxonomy-only

~~~bash
metaquest run \
  --db-dir /data/metaquest-db \
  --paired sample_R1.fastq.gz sample_R2.fastq.gz \
  --taxonomy-only \
  --output results/sample
~~~

This runs validation, fastp preprocessing, Kraken2, Bracken, and reporting.

### Assembly, gene prediction, and functional annotation

~~~bash
metaquest run \
  --db-dir /data/metaquest-db \
  --single reads.fastq.gz \
  --output results/sample
~~~

This additionally runs MEGAHIT, Pyrodigal, eggNOG-mapper, and BBMap gene
abundance. Use `--skip-functional` to skip eggNOG annotation while retaining
gene prediction and gene abundance, and `--skip-abundance` to skip read mapping.

Useful options:

| Option | Meaning |
|---|---|
| `--output DIR` | result directory |
| `--db-dir DIR` | database root override |
| `--skip-validation` | bypass input validation |
| `--taxonomy-only` | taxonomy-only execution |
| `--skip-annotation` | deprecated alias for `--taxonomy-only` |
| `--skip-functional` | skip eggNOG annotation; retain gene prediction and abundance |
| `--screen-pathogens` | species-target and curated AMRFinderPlus protein evidence |
| `--skip-abundance` | skip read mapping and gene/functional abundance |
| `--annotation-threads N` | eggNOG-mapper worker threads |
| `--plot-formats svg png pdf` | override configured figure formats |
| `--resume` | reuse stages only when inputs and effective configuration match |
| `--force` | preserve the existing output as a timestamped backup and start fresh |

Without `--resume` or `--force`, a nonempty output directory is rejected.
MEGAHIT scratch uses the system temporary filesystem (or `TMPDIR`). Under WSL,
keep `TMPDIR` on a Linux filesystem with space for assembly intermediates;
Windows-mounted result directories do not support its compressed-read FIFOs.

## Compare

Provide a CSV/TSV sample sheet with unique `sample_id` and `result_dir` columns.
Paths are relative to the sheet. The descriptive mode aligns full profiles by
feature ID and exports composition, QC and Bray–Curtis distances:

~~~bash
metaquest compare --sample-sheet samples.tsv --output comparisons/taxonomy
~~~

Formal testing initially supports two independent groups. Supply a group column,
reference label and an explicit proportion-scale pseudocount. See
[screening and comparison](screening-comparison.md) for assumptions and examples.

## Setup databases

~~~bash
metaquest databases --db-dir /data/metaquest-db --list
metaquest databases --db-dir /data/metaquest-db --database taxonomy
~~~

Use `--force` only when deliberately replacing an existing invalid or older
installation.

See [Database management](databases.md).

## Initialize configuration

~~~bash
metaquest init-config --output metaquest.yaml
~~~

Database root resolution:

~~~text
--db-dir
  → METAQUEST_DB_DIR
    → databases.base_dir in YAML
      → ./databases
~~~

## Current output layout

~~~text
results/sample/
├── preprocessing/
│   ├── cleaned_R1.fastq.gz
│   ├── fastp.json
│   └── fastp.html
├── kraken_classified.txt
├── kraken_report.txt
├── bracken_report.tsv
├── megahit_assembly/
├── gene_prediction/
│   ├── genes.faa
│   ├── genes.fna
│   ├── genes.gff3
│   ├── contigs.stable.fasta
│   ├── contig_id_map.tsv
│   └── summary.json
├── functional_annotation/
│   ├── functional_annotations.tsv
│   ├── functional_category_summary.tsv
│   ├── eggnog_mapper.log
│   ├── summary.json
│   └── completion.json
├── gene_abundance/
│   ├── gene_abundance.tsv
│   └── functional_abundance.tsv
├── screening/                   # when --screen-pathogens is enabled
│   ├── taxon_candidates.tsv
│   ├── protein_evidence.tsv
│   ├── amrfinder.tsv
│   └── summary.json
├── reports/                     # one HTML, JSON and TSV per stage
├── 01_taxonomic_report.txt
├── 02_assembly_report.txt
├── 03_functional_report.txt
├── analysis_summary.json
├── analysis_metadata.json
├── report.html
├── plots/
│   └── plot_data/*.tsv
└── metaquest.log
~~~

Taxonomy-only runs omit assembly and annotation outputs.

Taxonomy JSON reports both the fraction of species-assigned reads and the
fraction of all input reads. It also records Kraken-classified reads,
unclassified reads, classification rate, observed read length, and the
Bracken model length selected from the installed database.

## Exit behavior

- `0`: command completed successfully
- `1`: validation, dependency, database, or pipeline failure
- `130`: interrupted with Ctrl+C

Errors include an actionable message. Detailed external-tool output is retained
in the run log when logging is enabled.
