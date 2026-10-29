# Reproducible beta acceptance checks

Use an installed beta wheel and the repository's complete Conda environment.
For screening, also install the isolated screening environment and append its
bin directory as described in [installation](installation.md).
Run `metaquest --version` and save its output. Databases and native tools must
be the same across samples in a comparison. This checks runtime behavior and
accounting; it does not estimate sensitivity, specificity or clinical utility.

## Small public-data smoke run

The ENA run `ERR2984773` was used for the existing Zymo mock-community smoke
checks. Its paired downloads total about 1.96 GB. The following checksums are
from the saved ENA manifest; verify them before subsetting:

```bash
mkdir -p acceptance-data
curl --fail --location --retry 5 \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR298/003/ERR2984773/ERR2984773_1.fastq.gz \
  --output acceptance-data/ERR2984773_1.fastq.gz
curl --fail --location --retry 5 \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR298/003/ERR2984773/ERR2984773_2.fastq.gz \
  --output acceptance-data/ERR2984773_2.fastq.gz
md5sum -c <<'CHECKSUMS'
f5794d7ab28de74d803e9dad609a150a  acceptance-data/ERR2984773_1.fastq.gz
78bf00772b0cd8168ce37c16df345f7d  acceptance-data/ERR2984773_2.fastq.gz
CHECKSUMS
```

Create the first 10,000 read pairs using Biopython's streaming FASTQ iterator:

```bash
python - <<'PY'
import gzip
from itertools import islice
from pathlib import Path
from Bio.SeqIO.QualityIO import FastqGeneralIterator

directory = Path("acceptance-data")
for mate in (1, 2):
    with gzip.open(directory / f"ERR2984773_{mate}.fastq.gz", "rt") as source:
        with gzip.open(directory / f"smoke_R{mate}.fastq.gz", "wt") as output:
            count = 0
            for title, sequence, quality in islice(FastqGeneralIterator(source), 10000):
                output.write(f"@{title}\n{sequence}\n+\n{quality}\n")
                count += 1
            if count != 10000:
                raise SystemExit(f"Expected 10,000 reads in mate {mate}; got {count}")
PY
```

MetaQuest validates that the selected read identifiers and counts are paired.
Install/check all database profiles and run the full workflow:

```bash
metaquest databases --db-dir databases --database taxonomy
metaquest databases --db-dir databases --database functional
metaquest databases --db-dir databases --database amr
metaquest check --db-dir databases --screen-pathogens
/usr/bin/time -v -o acceptance-data/resources.txt \
  metaquest --no-color run --db-dir databases --low-memory \
  --paired acceptance-data/smoke_R1.fastq.gz acceptance-data/smoke_R2.fastq.gz \
  --screen-pathogens --output results/beta-smoke
python tools/verify_run.py results/beta-smoke > acceptance-data/accounting.json
```

If databases are installed elsewhere, use their common root throughout. The
reference data require roughly 11.3 GB compressed functional download and
53 GB installed functional storage, plus taxonomy, screening and working space.
MEGAHIT scratch must be on a filesystem supporting POSIX named pipes; under
WSL leave `TMPDIR` on Linux storage. `/usr/bin/time -v` records elapsed time and
the process-tree resource measurements reported by GNU time; compare resources
only with tool versions, input hashes, settings and hardware recorded.

The earlier 10,000-pair test used a subset of an existing subsample, rather than
this first-read prefix. It retained 9,569 pairs, assembled 100 contigs, predicted
106 genes and annotated 95. These are an observed runtime example, not expected
biological truth or fixed acceptance counts for a different subset. Its 3.53%
mapping rate demonstrates why a shallow smoke run cannot validate abundance accuracy.
The full reference search took about ten minutes on that machine despite the
small query. Larger, truth-labeled mock data are needed for biological validation.

## Resume and optional modes

Repeat the exact command with `--resume`. Input hashes, effective settings and
workflow flags must match. Check `analysis_metadata.json` attempt history and
the logs for cache reuse. Then rerun `tools/verify_run.py`.
See [download and run recovery](recovery.md) for integrity checks, legacy-output
compatibility, signal exit codes and hard-termination limitations.

Use fresh directories for additional checks:

```bash
metaquest run --db-dir databases --low-memory --single acceptance-data/smoke_R1.fastq.gz \
  --taxonomy-only --output results/beta-single
python tools/verify_run.py results/beta-single

metaquest run --db-dir databases --low-memory \
  --paired acceptance-data/smoke_R1.fastq.gz acceptance-data/smoke_R2.fastq.gz \
  --skip-functional --skip-abundance --screen-pathogens --output results/beta-skip
python tools/verify_run.py results/beta-skip
```

For interleaved inputs, place each R1 record immediately before its matching R2
record and use `--interleaved FILE` with another fresh output. Never reuse the
same biological sample as independent replicates in a formal group test.
Use the [comparison sample-sheet example](screening-comparison.md#descriptive-sample-comparison)
for real independent samples; synthetic tests validate calculations only.

## Inspect reports and failures

The acceptance checker verifies completed status, classified/unclassified and
Bracken denominators, unique/overlapping gene IDs, annotation counts, finite
nonnegative gene quantities, TPM normalization, screening joins, eight stage
reports, local report links and embedded image sources. Coordinates follow
1-based inclusive GFF3. It exits nonzero on an inconsistent artifact.

Open `report.html` and all eight `reports/*.html` pages in a browser. Check narrow
screens, printing, readable wide tables and the distinction between not run,
skipped and screened with no reported hits. PNG/SVG figures and their source TSVs
must agree. Structural checks do not replace browser rendering or manual review.

A missing database, mismatched paired IDs or a nonempty output without `--resume`
or `--force` must fail clearly. Screening with zero evidence rows must still
produce valid headers and an explicit completed state; it does not establish absence.
Retain `analysis_metadata.json`, the acceptance JSON, resource log, native logs
and report screenshots with the beta release evidence. Run `python -m pytest -q`
for the failure/zero-hit regressions, and set `METAQUEST_TEST_AMR_DB` to the
installed `databases/amr` directory to enable native positive/negative controls.
