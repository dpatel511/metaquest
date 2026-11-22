# eggNOG environment audit — 2026-11-22

The original eggNOG-mapper 2.1.15 Conda build resolved and installed with conflicting
Python dependency declarations. The correction was accepted in
[Bioconda PR #72011](https://github.com/bioconda/bioconda-recipes/pull/72011) and
published as `2.1.15=pyhdfd78af_1`. MetaQuest's maintained environment now pins that
build and requires `python -m pip check` in Conda CI. The original failure and
local correction below are retained as historical evidence. A passing test suite
alone does not establish dependency consistency or biological annotation accuracy.

## Reproduced conflict

The fresh Linux environment was resolved from `environment/environment.yml` on
top of MetaQuest revision `ab1bfa2`. Its 264 Conda packages included Python
3.11.6, eggNOG-mapper 2.1.15, Biopython 1.88, psutil 7.2.2 and XlsxWriter 3.2.9.
`python -m pip check` exited 1 with three conflicts.

| Dependency | eggNOG v2.1.15 Python metadata | Bioconda runtime requirement |
| --- | --- | --- |
| Biopython | `==1.76` | `>=1.76` |
| psutil | `==5.7.0` | `>=5.7.0` |
| XlsxWriter | `==1.4.3` | `>=1.4.3` |

These exact pins are present in the
[upstream v2.1.15 source](https://github.com/eggnogdb/eggnog-mapper/blob/v2.1.15/setup.cfg),
while the [original Bioconda recipe](https://github.com/bioconda/bioconda-recipes/blob/2ed02b79ee0f469484041e2433a1be79ea8752ba/recipes/eggnog-mapper/meta.yaml)
declares minimum versions. MetaQuest requires Biopython >=1.81, so the exact
upstream pin cannot be satisfied in its core Python environment. Do not downgrade
MetaQuest dependencies, delete package metadata or suppress this failure.
Upgrading to eggNOG-mapper's current development branch would also change its
annotation engine and database format; that is outside this reliability fix.

## Local package correction

A candidate Conda build of eggNOG-mapper 2.1.15 changes its three Python metadata
requirements to the same minimums already declared by Bioconda. The recipe build
number is incremented and `python -m pip check` is added to its tests. The source
archive remains pinned to SHA-256
`e2fa4402c15cfdcb85d9055450b974a100a2f74697a2c5b739a287c5182fa23b`.

In a second fresh environment, all 264 Conda package identities matched the first
environment except the candidate eggNOG package. Its 56 Python source files
matched the original package byte-for-byte. Dependency consistency then passed
with `pip check` exit 0. This is a source-package correction; installed metadata
was not edited, and no annotation algorithm or database was upgraded.

The initial package was built with `--no-test`; the submission candidate was then
rebuilt with all recipe tests enabled and passed, as did official Bioconda lint.
`python -m pip check` runs in `run_test.sh` with pip as a test-only dependency.
The earlier description of `download_eggnog_data.py -sfy` as a database download
test was incorrect: `-s` simulates downloads. That simulation passed in the full
recipe test. Local checks were Linux-only; upstream CI and review are tracked in
[Bioconda PR #72011](https://github.com/bioconda/bioconda-recipes/pull/72011).
The mismatch was reported in
[eggNOG-mapper issue #625](https://github.com/eggnogdb/eggnog-mapper/issues/625).

## Public build verification

On 2026-11-22, a fresh Linux x86-64 environment under WSL2 was created from the
maintained environment file using public Conda channels. It selected
`eggnog-mapper=2.1.15=pyhdfd78af_1`, Python 3.11.6, Biopython 1.88, psutil 7.2.2
and XlsxWriter 3.2.9. `python -m pip check` passed before and after installing a
locally built MetaQuest wheel from this working tree, based on revision `ab1bfa2`.
The package version remains `2.0.0b3`; this local wheel is not a new PyPI release.

The source, installed-wheel and extracted-sdist suites each passed 334 tests,
with two native AMRFinderPlus database controls skipped and 276 existing plotting
dependency deprecation warnings. All four version-provenance regressions were
collected by name and passed. The wheel was tested outside the source tree in the
fresh Conda environment, with its imported module path and package metadata
verified; it did not use the developer environment or a source editable install.
Lint, direct dependency licenses and strict distribution metadata checks passed.
The separate screening environment solved independently; its runtime/database
controls remain for complete-tool CI.

The public-build native control uses the installed wheel and the same three
cached positive protein IDs described below. It checks unique query joins,
annotation fields, the recorded eggNOG version and byte-identical annotation and
category tables against the earlier local candidate. This remains a scoped
integration check, not independent evidence of biological accuracy.

## Version provenance

eggNOG 2.x first probes Git in its installation directory, but its fallback
`git rev-parse` inherits the caller's directory. In a MetaQuest checkout this
reported `emapper-2.1.15-ab1bfa2`, falsely borrowing MetaQuest's revision.

MetaQuest now probes the version from a temporary directory, removes inherited
Git repository overrides, and prevents discovery of parent repositories even
when temporary storage is inside a checkout. A real eggNOG source checkout may
still report its own Git identity. Regressions exercise an unrelated caller
repository, repository-local temporary storage and inherited Git overrides.

Older outputs may contain a caller-derived eggNOG build suffix. Preserve those
outputs; do not rewrite their provenance to force cache reuse or comparison.
New outputs can have a different recorded tool identity and require reprocessing
when strict compatibility checks reject older evidence.

## Verification scope and next gate

The native control uses three previously annotated proteins from the cached
ERR2984773 development run and the existing eggNOG 5.0.2 database. It checks input
IDs, annotation status, seed orthologs, COG/KO/EC/GO/PFAM fields and output joins.
These are integration/consistency checks, not independent functional truth,
sensitivity estimates or performance benchmarks. No databases or sequencing data
were downloaded for this control.

Ignored evidence is retained under `test_results/eggnog-environment/`: resolved
package specifications, original failures, recipe/patch, build logs, package
hashes, dependency audits, native controls and test records. Those local candidate
artifacts remain distinct from the subsequently published Bioconda build. The
developer's existing environment was not modified.

Fresh public-build evidence is retained separately under
`test_results/eggnog-public-build/`, including the public package URL and checksum,
resolved environment specification, dependency checks, package tests and native
control. Keep the core and screening environments separate. Independent user
acceptance and held-out biological validation remain separate checkpoints.
