# Database management

MetaQuest keeps reference data outside the Python package and records the exact
database release used by an analysis.

## Available profiles

| Key | Database | Release | Download | Installed | Status |
|---|---|---:|---:|---:|---|
| `taxonomy` | Kraken2 Standard-8 with Bracken | 2026-06-26 | ~5.5 GB | ~7.5–8.1 GB | available |
| `functional` | eggNOG-mapper core data | 5.0.2 | ~11.3 GB | ~53 GB | available |
| `amr` | NCBI AMRFinderPlus | tool-compatible release | varies | varies | available |
| `virulence` | VFDB Set A candidate | — | — | — | planned |

The taxonomy archive comes from the official Kraken2 index catalog:
<https://benlangmead.github.io/aws-indexes/k2>.

Planned entries are visible so users can understand the intended layout, but
MetaQuest refuses to install them until their versions and validation contracts
are finalized. The functional profile installs eggNOG annotation, taxonomy, and
DIAMOND data only; it excludes Pfam, MMseqs, HMM, and novel-family data.

## Inspect status

~~~bash
metaquest databases --db-dir /data/metaquest-db --list
~~~

Example:

~~~text
taxonomy: installed | 2026-06-26 | 5.5 GB download / 7.5 GB installed
functional: available | 5.0.2 | 11.3 GB download / 53 GB installed
amr: available | tool-compatible release | upstream release size
virulence: planned | pending | not finalized
~~~

## Install taxonomy

~~~bash
metaquest databases \
  --db-dir /data/metaquest-db \
  --database taxonomy
~~~

The installer:

1. resolves and creates the database root;
2. checks free disk space;
3. downloads the upstream checksum list;
4. downloads or resumes the archive;
5. selects the checksum matching the archive filename;
6. verifies the complete archive;
7. rejects unsafe archive paths and links;
8. extracts into a temporary directory;
9. verifies required Kraken2 index files;
10. writes `metaquest-db.json`;
11. atomically publishes `taxonomy/`;
12. removes the completed download cache.

An existing valid installation is reused. An invalid directory is never
silently overwritten; replacement requires `--force`.

## Install functional data

~~~bash
metaquest databases \
  --db-dir /data/metaquest-db \
  --database functional
~~~

The functional installer keeps resumable partial files under
`.downloads/functional/`, checks the official compressed byte sizes, records a
local SHA-256 for each artifact, and validates both SQLite files plus the
DIAMOND database before publishing the installation. Partial files are retained
after a failure and removed only after a successful install.

This checkpoint makes eggNOG 5.0.2 installation reproducible. The full pipeline
uses these files immediately after Pyrodigal gene prediction.

## Install screening data

Install the isolated screening tool environment, then use MetaQuest's database
command while the core environment is active:

~~~bash
conda env create -f environment/screening.yml
conda activate metaquest
export PATH="$PATH:$(conda run -n metaquest-screening bash -c 'printf "%s/bin" "$CONDA_PREFIX"')"
metaquest databases --db-dir /data/metaquest-db --database amr
metaquest check --db-dir /data/metaquest-db --screen-pathogens
~~~

MetaQuest invokes NCBI's native updater to download and index a compatible
release. It validates the staged database with AMRFinderPlus, records its actual
version and core-file SHA-256 hashes, then publishes `ROOT/amr/` atomically.
Its `latest` link selects the installed version. Reuse that installation across
all samples in a study. `--force` updates it deliberately; compatible releases
can change upstream, so a new installation is not a historical version pin.
The updater needs network access and at least 2 GiB of working space.

## Directory layout

~~~text
/data/metaquest-db/
├── taxonomy/
    ├── hash.k2d
    ├── opts.k2d
    ├── taxo.k2d
    ├── database50mers.kmer_distrib
    ├── database75mers.kmer_distrib
    ├── database100mers.kmer_distrib
    ├── database150mers.kmer_distrib
    ├── database200mers.kmer_distrib
    ├── database250mers.kmer_distrib
    ├── database300mers.kmer_distrib
│   └── metaquest-db.json
└── functional/
    ├── eggnog.db
    ├── eggnog.taxa.db
    ├── eggnog_proteins.dmnd
    └── metaquest-db.json
~~~

## Manifest

`metaquest-db.json` records:

- MetaQuest database key and descriptive name;
- upstream release;
- source URL;
- download and installed-size estimates;
- installation timestamp;
- upstream artifact URL and compressed byte size;
- verified taxonomy archive MD5 or local functional-artifact SHA-256.

MD5 is used because it is the checksum published by the upstream index catalog.
It protects transfer integrity; it is not treated as a cryptographic signature.

## Path resolution

MetaQuest resolves the root in this order:

1. command-level `--db-dir`;
2. `METAQUEST_DB_DIR`;
3. `databases.base_dir` in the selected YAML file;
4. `./databases`.

The taxonomy index defaults to `ROOT/taxonomy`.
Pass the root to `--db-dir`, for example `databases/`, rather than
`databases/taxonomy/`. Passing the installed taxonomy/functional subdirectory
now produces a targeted error instead of appending a second subdirectory.

## External disk with a repository-local path

~~~bash
mkdir -p /mnt/e/metaquest/databases
ln -s /mnt/e/metaquest/databases databases
metaquest databases --database taxonomy
~~~

The repository ignores `/databases`, including this symlink.

## Updating a release

Taxonomy and functional database releases are pinned in source. AMRFinderPlus
uses its native compatible-release updater and records the installed release.
Changing a release requires:

- updating the catalog entry and expected size;
- adding or updating checksum/parser tests;
- verifying Kraken2 and all bundled Bracken read lengths;
- recording the change in the changelog;
- rerunning benchmark datasets before making the new release the default.
