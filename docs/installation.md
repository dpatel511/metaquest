# Installation

This guide covers the MetaQuest `2.0.0b1` research beta.

## Supported environment

- Linux or WSL2
- x86-64
- Python 3.10, 3.11, or 3.12
- conda or mamba
- at least 16 GB RAM recommended for the current Standard-8 profile
- at least 15 GB free during taxonomy database installation

MetaQuest has not yet been validated on native Windows or macOS.

## Install from PyPI

The distribution is named `metaquest-bio`; it installs the `metaquest` Python
package and CLI command:

~~~bash
python -m pip install --pre 'metaquest-bio==2.0.0b1'
metaquest --version
~~~

Use that command after beta publication; until then install this checkout below.
PyPI installs the Python layer only. The external bioinformatics tools listed
below must be installed separately. The recommended complete environment is the
repository's Conda environment.

## Install the complete Conda environment

From the repository root:

~~~bash
conda env create -f environment/environment.yml
conda activate metaquest
python -m pip install -e .
~~~

For a regular non-editable installation:

~~~bash
python -m pip install .
~~~

Confirm the CLI:

~~~bash
metaquest --version
metaquest --help
~~~

## Runtime components

The minimal environment contains:

| Component | Role |
|---|---|
| fastp | adapter and quality preprocessing |
| Kraken2 | read-level taxonomic classification |
| Bracken | abundance re-estimation |
| MEGAHIT | metagenomic assembly |
| Pyrodigal | metagenomic gene prediction |
| eggNOG-mapper and DIAMOND | orthology-based functional annotation |
| BBMap | gene abundance mapping and interleaved FASTQ splitting |
| AMRFinderPlus | optional curated protein evidence screening |

Pyrodigal is a Python dependency and does not require a separate Prodigal or
Prokka executable.

## Configure databases

Choose a database root outside Git:

~~~bash
metaquest databases --db-dir /data/metaquest-db --list
metaquest databases --db-dir /data/metaquest-db --database taxonomy
metaquest databases --db-dir /data/metaquest-db --database functional
~~~

Then verify the installed tools and databases:

~~~bash
metaquest check --db-dir /data/metaquest-db
~~~

For repository-local convenience with physical storage on another disk:

~~~bash
mkdir -p /mnt/e/metaquest/databases
ln -s /mnt/e/metaquest/databases databases
metaquest databases --database taxonomy
~~~

See [Database management](databases.md) for manifests, checksums, and recovery.

## Configuration

Generate a YAML file:

~~~bash
metaquest init-config --output metaquest.yaml
~~~

Global options may appear before or after the command:

~~~bash
metaquest run --config metaquest.yaml \
  --single reads.fastq.gz \
  --taxonomy-only \
  --output results/
~~~

## Development dependencies

Install the test extras:

~~~bash
python -m pip install -e '.[test]'
pytest -q
~~~

## Troubleshooting

### Database not found

Confirm the selected root:

~~~bash
metaquest databases --db-dir /data/metaquest-db --list
ls -lh /data/metaquest-db/taxonomy/hash.k2d
~~~

Pass the same root to analysis or export it:

~~~bash
export METAQUEST_DB_DIR=/data/metaquest-db
~~~

### Interrupted database download

Run the same setup command again. Partial archives are retained under
`.downloads/` and resumed when the server supports HTTP range requests.

### Invalid existing database

MetaQuest will not silently overwrite an invalid directory. Inspect or move it,
or deliberately replace it:

~~~bash
metaquest databases \
  --db-dir /data/metaquest-db \
  --database taxonomy \
  --force
~~~

### Functional annotation fails preflight

The maintained functional stage requires eggNOG-mapper 2.1.15, DIAMOND, and all
three eggNOG 5.0.2 database files. Run `metaquest check` for exact missing
components.
