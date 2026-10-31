# Comparative-analysis validation

Checkpoint 4 evaluates the existing exploratory CLR/Mann–Whitney method for
**two independent groups**. This is an unreleased audit after `2.0.0b2`;
public PyPI artifacts are unchanged. It does not introduce a new statistical
model or establish biological differential-abundance accuracy.

## Frozen experiment

The [protocol](../benchmarks/comparison/protocol.json) was saved before simulation
with SHA-256 `44b236dd6e3bbd97c0cab31d1f59e0e61c3351cfc4ccfc939881e2adb6253b5a`.
It specifies 500 independent synthetic studies in each of four scenarios:

| Scenario | Samples per group | Purpose |
| --- | --- | --- |
| Dense global null | 12 / 12 | Equal CLR distributions with every feature present |
| Unequal-size global null | 6 / 18 | Exchangeable populations with unequal group sizes |
| Sparse global null | 6 / 6 | Multinomial counts with independent dropout and the same population in both groups |
| Balanced CLR effects | 24 / 24 | Two known shifts, +1 and -1, with 22 exact CLR distribution nulls |

Each study has 24 features. BY and BH use the same generated study; their results
are paired evaluations, not independent repetitions. NumPy PCG64 and independent
`SeedSequence([scenario_seed, replicate_index])` streams make the inputs reproducible.
There are 4,000 primary method evaluations and 2,100 predeclared sensitivity
evaluations; input hashes and all per-study metrics are retained.

For continuous scenarios, independent bounded latent log coordinates generate
proportions by inverse additive-pseudocount closure:
`p = (1 + D * epsilon) * softmax(z) - epsilon`.
This deliberately makes `CLR(p + epsilon) = z - mean(z)` exact, with positive
proportions under the chosen bounds. Equal-and-opposite shifts leave the other
22 CLR-coordinate distributions unchanged. This is a controlled calculation
model; it does not model all sequencing or biological variability. The sparse
null separately tests zeros and finite-depth counts under exchangeable groups.

Settings are frozen at alpha 0.05, proportion-scale pseudocount 1e-6 and overall
prevalence 0.1. Sensitivity checks use pseudocount 1e-4 and prevalence 0.5 on the
same effect studies, plus pseudocount 1e-4 for the first 50 sparse-null studies.
They are descriptive: changing the pseudocount/reference changes the hypotheses.
All dense-effect features are present, so changing their prevalence threshold
alone is expected to leave their retained set unchanged.

## Metrics and acceptance limits

Under a global null, FDR equals the probability of any discovery in a study.
Report that event count, its rate and a two-sided 95% Wilson interval over the
500 independent studies. Do not treat correlated features within a study as
independent replicates.

For the effect scenario, report the mean study-level FDP (false discoveries /
max(discoveries, 1)), target recall and their Monte Carlo standard errors.
The explicit Wilson calculation is regression-checked against SciPy binomial
intervals. One-sided 95% Hoeffding bounds use the bounded [0,1] study outcomes and require
no independence among features. Correct direction requires both target median
CLR effects to have their declared signs; its uncertainty uses a Wilson interval.
Zero observed FDP still has a nonzero upper uncertainty bound.

The proposed project screening limits were fixed before scoring:

- Global-null Wilson upper bound <= 0.08 for each scenario/method.
- Effect FDP Hoeffding upper bound <= 0.10.
- Target-recall Hoeffding lower bound >= 0.80.
- Correct-direction Wilson lower bound >= 0.95.

These margins accommodate finite Monte Carlo uncertainty. They are **not** a
claim that FDR of 8–10% is the nominal target, nor proof of universal 5% control.
Always report the observed estimates and nominal alpha alongside the screening
result. A failed criterion remains a failure; do not tune seeds, effects or
thresholds to erase it. See [recorded results](../benchmarks/comparison/RESULTS.md).

Reproduce from the intended source checkout in an environment with test/runtime
dependencies installed:

```bash
sha256sum benchmarks/comparison/protocol.json
PYTHONPATH=src python -m pytest -q tests/test_comparison.py \
  tests/test_comparison_statistics.py tests/test_comparison_simulation.py
PYTHONPATH=src /usr/bin/time -v -o comparison-resources.txt \
  python tools/benchmark_comparison.py \
  --protocol benchmarks/comparison/protocol.json --output comparison-validation
```

Use a new output directory. The runner writes per-study TSV, summary JSON and
standalone PNG/SVG/PDF plots. It exits 1 when a scientific screening criterion
fails while retaining the results. A successful harness run is not a biological
benchmark pass or a stable-release decision.

## Independent calculation and integration checks

The regression suite checks exact small-sample U p-values by exhaustive label
permutations, the documented U=17 example, and larger tied samples against an
independent continuity/tie-corrected normal calculation. It checks known BH/BY
step-up answers across every tested feature, explicit two-feature CLR log-ratios,
prevalence reference changes, contrast reversal and row/column reorder invariance.
Profile integration checks retain leading-zero sample IDs, verify the full union
of feature IDs and zero-filled missing features, and test row-scaling invariance.
These calculation fixtures are not biological replicates.

Input regressions cover repeated subjects after whitespace normalization,
shared/duplicated raw-file hashes, incomplete provenance, mismatched Bracken
models or fastp versions, duplicate headers/IDs, invalid quantities/ranks,
summary/count disagreement, conflicting feature labels and nonempty outputs.
Text, JSON and HTML methods retain the actual contrast, transformation,
pseudocount, prevalence, correction and supported design. Fresh/forced CLI output
is tested separately so its own diagnostic log cannot invalidate a fresh run.

## Interpretation and compatibility

The CLR reference is the retained feature set. The pseudocount is added on the
full reported-profile proportion scale before centering; retained features are
not renormalized before that addition. Unclassified/unassembled material is
excluded from the reported-profile denominator. Effects concern relative
log-ratio distributions, not organism loads, expression or absolute fold changes.
Prevalence and pseudocount sensitivity can change both effect sizes and nulls.

Mann–Whitney tests equality of CLR distributions. Interpreting it solely as a
median test requires additional distributional assumptions. The existing SciPy
`auto` policy uses an exact untied test if either group has <=8 observations;
otherwise it uses a tie-corrected normal approximation with continuity correction.
A supplemental exact two-point null calculation exposed miscalibration of that
approximation with 3 versus 30 tied samples: the public beta reaches 8.21% FDR
at nominal 5% under BH in one declared population. The candidate therefore
**rejects nonconstant tied CLR values if either group has fewer than 10 samples**,
with an actionable permutation-test message. No features are silently discarded.
All-pooled-constant features still receive p=1; untied small studies retain the
exact test. Larger tied studies retain the tie-corrected asymptotic calculation.
The zero small/tied count in successful summaries means the validity guard passed.
Three samples per group is an input minimum for untied data, not a power recommendation.
A permutation method is not introduced in this audit.

The separately frozen [tied-null protocol](../benchmarks/comparison/tied-null-protocol.json)
keeps the original Monte Carlo protocol unchanged. It enumerates independent
binomial count pairs for two-point profiles, across six group-size designs and
nine population probabilities. Its planning text overcounted the state total;
there are 374 actual count-pair states, not 474. All 108 population/correction
calculations are retained. Reproduce with the intended installed package (run
outside its source tree) or `PYTHONPATH=src` for the audit candidate:

```bash
python /path/to/metaquest/tools/benchmark_comparison_ties.py \
  --protocol /path/to/metaquest/benchmarks/comparison/tied-null-protocol.json \
  --output tied-null.json
```

The unchanged public `2.0.0b2` wheel reproduces four failed BH cases. The candidate
blocks the unsupported nonconstant small/tied states and reports their probability
mass as unsupported, **not as zero FDR**. Fully supported cases are reported
separately. The main sparse Monte Carlo model contained zeros but no transformed
ties; this is why its successful results did not expose the problem. The original
simulations are repeated against the final guard using identical seeds; repeated
runs are not counted as additional independent studies.

BY allows arbitrary test dependence when the individual p-values are valid.
BH requires independence or appropriate positive dependence; CLR dependence is
not automatically guaranteed to meet that condition. Simulations assess only
the declared families and do not remove these assumptions. Batch effects,
covariates, pairing and repeated measures remain unsupported. Different input
hashes and unique subject labels cannot prove biological independence; users
must verify study design and metadata themselves.

Current comparisons require complete preprocessing provenance and compatible
fastp versions when preprocessing ran. Taxonomy comparisons additionally require
the selected Bracken model length and rank, integer estimated counts and agreement
with assigned/classified/input summary counts. Missing legacy evidence is an
error; keep the original outputs and reprocess rather than inventing provenance.
Source tables with conflicting labels for the same stable ID require inspection.
Comparison output directories must be new or empty before the CLI starts; its
reserved `metaquest.log` is allowed by the library guard and is preserved. No
prior result tables or reports are overwritten. `--force` backs up prior CLI
output before preparing a fresh directory.

Primary sources: [SciPy Mann–Whitney documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html),
[SciPy FDR documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html),
and [Gloor et al. (2017)](https://doi.org/10.3389/fmicb.2017.02224).
The numerical limits above are MetaQuest project criteria, not criteria from
those sources. Browser layout review and independent biological/user acceptance
remain separate checkpoints.
