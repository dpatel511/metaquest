# Installation

This guide covers the MetaQuest `2.0.0b4` research beta.

The source checkout also contains unreleased beta-5 host filtering. To test that
feature, use the checkout's `environment/environment.yml` and install the
checkout with `python -m pip install --no-deps .`; the published beta-4 wheel
does not contain it. Bowtie2 `>=2.5,<3` and Samtools `>=1.18,<2` are added to the
core Conda environment. They are required at runtime only when a host index is
configured. See [host reference setup](host-removal.md). Existing beta-4
environments need those tools installed before testing the new wheel; selecting
a different Python executable alone does not change which native tools PATH uses.

Check the [release page](https://github.com/dpatel511/metaquest/releases) for
publication status before using the version-pinned index/tag commands below.
For an unpublished checkout, build a wheel with `python -m build`, install that
wheel with `--no-deps` inside the complete environment, and verify
`metaquest --version` and `pip check`.

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
python -m pip install --upgrade --pre 'metaquest-bio==2.0.0b4'
metaquest --version
~~~

PyPI installs the Python layer only. The external bioinformatics tools listed
below must be installed separately. The recommended complete environment is the
repository's Conda environment.

## Install the complete Conda environment

Use the environment files matching the package release:

~~~bash
git clone --branch v2.0.0b4 https://github.com/dpatel511/metaquest.git
cd metaquest
conda env create -f environment/environment.yml
conda activate metaquest
conda install --strict-channel-priority -c conda-forge -c bioconda 'eggnog-mapper=2.1.15=pyhdfd78af_1'
python -m pip install --no-deps 'metaquest-bio==2.0.0b4'
python -m pip check
~~~

`--no-deps` preserves the Python packages selected by Conda. For source development,
use `python -m pip install --no-deps -e .` from the intended checkout; use
`python -m pip install --no-deps .` for a non-editable source installation.

The maintained environment pins eggNOG-mapper's corrected public build
`2.1.15=pyhdfd78af_1`. The explicit installation above also selects that build when
using an older release's environment file. Build 0 has conflicting Python
dependency metadata; upgrade the package rather than editing installed metadata.
The [environment audit](environment-validation.md) records the original failure,
upstream correction and fresh public-build verification. `pip check` must pass
after installing MetaQuest; Conda resolving alone does not establish consistency.

Optional screening uses a separate native-tool environment. eggNOG-mapper's
DIAMOND `<2.1` requires zlib `<1.3`, while current AMRFinderPlus requires a newer
curl/zlib stack; installing both into one environment does not solve reliably.
Keep the core environment active and append the screening executables:

~~~bash
conda env create -f environment/screening.yml
conda activate metaquest
export PATH="$PATH:$(conda run -n metaquest-screening bash -c 'printf "%s/bin" "$CONDA_PREFIX"')"
amrfinder --version
metaquest databases --db-dir /data/metaquest-db --database amr
metaquest check --db-dir /data/metaquest-db --screen-pathogens
~~~

Repeat that PATH export in new shells where screening is needed. Appending keeps
the core environment's Python first. AMRFinderPlus retains its own native library
and helper-program paths. Do not downgrade the screening environment's shared
libraries to force it into the core environment.

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

## Upgrading from earlier betas

Beta 2 fixes CLI stream handling and interrupted download/run recovery; it does
not change the underlying profiling algorithms. Preserve earlier outputs. Strict
`--resume` rejects outputs lacking the new database identities; rerun into a fresh
output directory rather than modifying metadata. See [recovery](recovery.md).

Beta 3 adds comparison safeguards. Nonconstant tied CLR features are rejected
when either group has fewer than 10 samples; untied studies still require at
least three. All-pooled-constant features retain p=1. The current method does
not substitute a permutation, paired or batch-adjusted test.

Comparisons require compatible selected Bracken models, preprocessing/tool
versions and count/rank evidence. Formal studies also require complete raw-file
SHA-256 values and independent subjects/files. Existing b2 results with complete
provenance can be reused as comparison inputs; incomplete legacy outputs need
reprocessing. Preserve original results rather than editing their metadata.
See [comparison validation](comparison-validation.md) for the measured b2 failure,
new boundary and interpretation limits.

The [development benchmark](benchmarking.md) retains a failed taxonomy gate and
unresolved reference/timing discrepancies. Biological validation is incomplete.

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
