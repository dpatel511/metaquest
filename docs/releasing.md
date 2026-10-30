# Release process

MetaQuest publishes Linux-supported Python distributions through GitHub Actions
and PyPI Trusted Publishing. No long-lived PyPI API token is stored in GitHub.

## One-time repository configuration

1. Reuse the existing production `pypi` environment and Trusted Publisher.
   If an optional TestPyPI rehearsal is wanted, also configure `testpypi`.
2. Require a maintainer approval for publishing environments. Restrict `pypi` to
   protected tags matching `v*` if the repository plan exposes that control.
3. Protect `main` and require these CI checks before merging:
   `Python 3.10`, `Python 3.11`, `Python 3.12`,
   `Direct dependency licenses`, `Build distributions`, and
   `Conda integration (Python 3.11)`.
4. For an optional rehearsal, configure a TestPyPI trusted publisher for:
   - owner: `dpatel511`
   - repository: `metaquest`
   - workflow: `testpypi.yml`
   - environment: `testpypi`
5. The production PyPI Trusted Publisher uses the same owner and repository,
   workflow `release.yml` and environment `pypi`. Verify this existing configuration
   before changing it; a separate TestPyPI account is not required for production.

## Optional TestPyPI rehearsal

This is an additional packaging check, not a prerequisite for an established
production Trusted Publishing route. It requires its own TestPyPI account-side
publisher setup. If omitted, verify the built wheel in isolation before release
and download/install the exact public PyPI artifacts after publication.

First run CI on the exact commit to be uploaded and wait for all six required
checks to pass. Both publishing workflows refuse to build for publication if
the latest push/manual CI run for that commit has not succeeded. PR CI checks a
merge revision, so dispatch CI on a feature branch before its TestPyPI rehearsal;
after merging, wait for the `main` push CI run before tagging that commit.

Run the **Publish to TestPyPI** workflow manually. Approve the `testpypi`
environment deployment, then inspect both the project page and uploaded wheel
and source distribution. Versions on package indexes are immutable, so bump the
pre-release serial before repeating an upload.

The workflow installs runtime dependencies from PyPI, then the exact MetaQuest
version from TestPyPI with `--no-deps`. This keeps the package source unambiguous.
Its final job checks the installed version, packaged species catalog and CLI help
from outside the checkout. A green upload alone is not a successful rehearsal;
the installation-verification job must also pass.

In a fresh complete Conda environment, repeat the documented small pipeline run
using that uploaded package instead of `pip install -e .`:

```bash
python -m pip install --no-deps --index-url https://test.pypi.org/simple/ 'metaquest-bio==2.0.0b2'
metaquest --version
```

Follow [acceptance checks](acceptance.md), retain the resulting evidence and
inspect HTML in a browser. Reusing existing valid reference databases is fine;
record their manifests and actual tool versions.

## Beta release gate

The research-beta scope is frozen around FASTQ validation/preprocessing, taxonomy,
assembly, gene prediction, annotation, abundance, reporting and the documented
optional research screening/comparison modes. During beta prioritize fixes,
installation reliability, documentation and independent validation.
Track post-beta evidence against the [development checkpoints](development-checkpoints.md).

Before publishing:

- Update version, release date, Beta classifier, installation commands and changelog together.
- Require the Python 3.10/3.11/3.12, license, distribution and complete Conda checks.
  The distribution job tests a wheel in an isolated environment outside the source
  tree. Conda CI installs screening data with MetaQuest and sets
  `METAQUEST_TEST_AMR_DB`, so the native positive/negative test runs.
  The core and screening tool environments resolve separately to avoid their
  incompatible DIAMOND/zlib and AMRFinderPlus/curl constraints.
- Inspect the native-integration artifact for its JUnit results and database
  manifest. The AMRFinderPlus cache key includes the actual tool version; increment
  its cache revision deliberately when refreshing the selected compatible data.
- Finish isolated-wheel verification and acceptance checks for read modes,
  skip/resume, missing-data errors and zero-hit reporting. TestPyPI can provide
  an additional optional packaging rehearsal.
- Review report rendering on desktop/narrow screens and in print. Automated
  link/image checks and inspected PNGs do not prove browser layout correctness.
- Describe screening as a 16-species bacterial target subset plus selected curated
  protein families; describe formal comparisons as exploratory two-group CLR tests.
  Runtime validation is not a sensitivity/specificity or clinical-validation claim.
- Record dependency deprecation warnings and verify with fresh supported packages.
  Do not suppress failures or broad classes of warnings to make checks green.

Older outputs lacking complete provenance may be unsuitable for `compare` or
strict `--resume`. Reprocess into a new directory rather than editing provenance
to force compatibility. Keep original output directories intact.

## Beta 2 scope and residual limits

`2.0.0b2` carries the merged CLI contract and recovery fixes plus the frozen
ERR2984773 development benchmark. Checkpoint 3 remains open: its primary taxonomy
gate failed, structural reference disagreements are unresolved, and elapsed
measurements disagree. The benchmark omitted functional annotation and screening.
Held-out biological validation, comparative-method simulations, browser review
and independent user acceptance remain outstanding. Larger-database testing is
deferred because of available hardware. Retain these limits in release notes;
a reliability beta does not require claiming that stable-release gates passed.

## Production release

1. Update `src/metaquest/config.py` to the canonical PEP 440 version and add the
   corresponding changelog entry.
2. Merge the release preparation PR after all required checks pass.
3. Create a signed or annotated tag with the same canonical version prefixed by
   `v`, for example `v2.0.0b2`.
4. Create and publish a GitHub Release from that tag. Mark beta releases as
   prereleases; the upload workflow checks that this flag matches the package version.
5. Complete any configured `pypi` environment approval. The current production
   environment publishes automatically after its required jobs pass.

Publishing the GitHub Release triggers `release.yml`. The workflow verifies
that the tag and package versions match, builds and checks the wheel and source
distribution in one job, and gives only the final upload job an OIDC token.
Draft releases do not publish to PyPI.

After publication, download the exact wheel and source distribution from PyPI,
verify their SHA-256 hashes against index metadata and the publishing artifacts,
and install the wheel outside the checkout in a fresh environment. Check
`metaquest --version`, root/subcommand help, packaged catalog resources, CLI failure
streams and the isolated-wheel tests. Run `pip check`; record any native controls
that require separate tools/data. Build the source distribution independently to
confirm it contains the files needed to reproduce its packaged runtime.
