# Download and run recovery

Recovery checks protect existing data and recorded provenance. They do not
establish biological accuracy. Use one writer per output/database directory.

## Database downloads and replacement

Repeat the same managed setup command after a connection failure:

```bash
metaquest databases --db-dir /data/metaquest --database taxonomy
metaquest databases --db-dir /data/metaquest --database functional
```

Partial downloads remain under the database root's `.downloads` directory.
Requests have a 30-second network-operation timeout and a bounded retry count;
this is not a deadline for downloading an entire multi-GB file. Byte-range
responses must agree with their declared offsets and total size. A server that
ignores a resume range is rejected without overwriting the existing prefix.
A complete partial file is accepted on HTTP 416 only when the reported total
matches its local size and any configured expected size.

Taxonomy archives are checked against the upstream MD5. Functional archives use
pinned compressed sizes, archive integrity checks, SQLite/native-tool validation,
and locally recorded SHA-256 values. The locally calculated SHA-256 values are
provenance, not independent upstream authenticity checks. A corrupt functional
archive rejected during extraction is discarded so a retry can download it again;
an extraction filesystem error retains downloaded archives for a later retry.

Replacement is staged and validated before promotion. `--force` permits replacing
an existing installation; it does not delete that installation before validation.
The old installation is temporarily retained at `.taxonomy.previous`,
`.functional.previous` or `.amr.previous` under the common database root.
These are reserved installer backup locations. Ordinary exceptions and handled
signals restore the old installation if promotion has not completed. On a later
setup invocation, a backup left with no target is restored before installation
checks. If both exist, a valid target permits backup cleanup; an invalid target
causes an error retaining the backup for inspection.

Free-space preflight cannot guarantee that another process will not fill the
disk later. Keep the existing installation and downloaded files when resolving
a space error; do not edit manifests to force an installation to appear valid.

## Pipeline resume

Repeat the original `metaquest run` command with its original input paths,
configuration, workflow options and output directory, adding `--resume`.
For example, if the original command was:

```bash
metaquest run --paired R1.fastq.gz R2.fastq.gz \
  --db-dir /data/metaquest --output results/sample
```

resume with:

```bash
metaquest run --paired R1.fastq.gz R2.fastq.gz \
  --db-dir /data/metaquest --output results/sample --resume
```

Input hashes, effective configuration and workflow flags must match. Database
identity records manifest hashes and index paths, sizes and modification times;
normal database replacements invalidate resume. This avoids rehashing tens of GB
on every metadata save, but cannot detect deliberate index modifications that
preserve both size and modification time. Fastp, Kraken2, Bracken and MEGAHIT
version changes are rejected when their previous versions were available.

Each completed stage records hashes for its essential reusable outputs.
Only completed stages with matching artifacts may reuse work; a failed stage,
missing file or changed hash invalidates that stage and downstream reuse.
Interleaved split reads are included when preprocessing is disabled. Functional
annotation also honours the runner's reuse decision and clears stale completion
markers before a new attempt. Reports are regenerated.

Older outputs without the new database identities are rejected for strict resume.
Retain those outputs and reprocess in a new directory. Never insert missing
fingerprints or edit configuration/provenance to bypass the check.

Ctrl+C exits 130; SIGTERM exits 143. Both unwind native process groups and record
an interrupted attempt. Native descendants are terminated, including fixtures
that ignore SIGTERM. Metadata replacement uses a temporary file, preserving the
last complete JSON if promotion fails. Attempt history retains earlier failures;
a previous unfinished attempt is marked interrupted on the next resume.

SIGKILL and power loss cannot execute Python cleanup. Backup-layout recovery and
unfinished-attempt handling cover their persisted state, not live-process cleanup.
Ensure surviving tools have stopped before retrying; overlapping writers are not
supported. This audit does not simulate a physical filesystem/power failure.

## Checkpoint 2 verification

The regression suite uses small disposable installations, a localhost HTTP server
and owned process trees. It covers disconnects, invalid ranges, completed partials,
network timeouts, corrupt archives/manifests, low/full disks, failed promotion,
backup recovery, metadata atomicity, interrupted attempts, cache invalidation and
SIGINT/SIGTERM cleanup. Run it in the documented tool environments:

```bash
python -m pytest -q tests/test_recovery.py tests/test_process_recovery.py
```

On 2026-10-29, all 34 recovery/process tests passed. The full source and installed
wheel suites each passed 173 tests, including existing native fixtures, without
skips or warnings. Lint, strict wheel metadata, dependency and license checks also
passed. The changes were merged in
[PR #31](https://github.com/dpatel511/metaquest/pull/31); the merged main commit
passed [remote CI](https://github.com/dpatel511/metaquest/actions/runs/40891232601).
These recovery guarantees apply to beta 2; prior public b1 artifacts are unchanged.

The native control used the existing 10,000-pair mock-community smoke fixture,
omitting functional search, abundance and screening. An interruption before gene
prediction exited 130; resume reused fastp, taxonomy and assembly and completed.
Both runs retained 9,569 fragments and produced 100 contigs and 106 genes. Their
taxonomy tables matched byte-for-byte, their contig sequence multisets and gene
sequence maps matched, and their GFF3 records matched using **1-based inclusive**
coordinates. Independent MEGAHIT runs differed in raw FASTA ordering/headers;
the initial bytewise comparison failed and was retained as an audit finding.
A repeated resume left the cached scientific files byte-identical.

Both outputs passed the acceptance checker, including eight stage reports and
16 figure files. These are accounting/structure checks; browser rendering remains
a separate checkpoint. Local logs and scripts are retained under ignored
`test_results/recovery-check/`.
