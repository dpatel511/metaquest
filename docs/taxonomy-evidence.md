# Taxonomy evidence — beta 4 preparation

MetaQuest `2.0.0b4` adds evidence schema **1**. This source version is prepared but
unpublished. Kraken2 and Bracken algorithms, default settings and original output
files are unchanged. Existing JSON count fields remain readable; the new evidence
files and canonical rank-neutral fields are additive.

## Full taxon evidence

`taxonomy_evidence.tsv` contains every row in the reported Bracken profile,
including rows omitted from the top-taxa display. Rows join the original Kraken
report by NCBI taxonomic ID, never by display name. Duplicate IDs and invalid,
nonintegral or inconsistent counts fail reporting instead of being silently
aggregated, truncated or clamped. This table does not list every Kraken-only
taxon or constitute a reference-defined negative universe.

| Column | Meaning |
| --- | --- |
| `taxonomy_id`, `name`, `rank` | Bracken's ID, display name and requested rank; IDs remain strings. |
| `count_unit` | `fragments` for paired/interleaved input; `reads` for single input. |
| `kraken_direct_count` | Kraken assignments directly to this ID, excluding descendants. |
| `kraken_clade_count` | Kraken assignments to this ID and descendants. Do not sum across nested clades. |
| `bracken_assigned_count` | Bracken's starting count (`kraken_assigned_reads`), including clade support; this is not the Kraken direct column. |
| `bracken_added_count` | Bracken's redistributed count (`added_reads`), not independent read observations. |
| `bracken_estimated_count` | Final estimate (`new_est_reads`); assigned + added when both source columns exist. |
| `fraction_all_input` | Estimate / all cleaned classification input, including unclassified input. |
| `fraction_bracken_profile` | Estimate / sum of exported integer Bracken estimates at the requested rank. |

Fractions use 0–1, while HTML displays percentages. The profile fraction is
recomputed from exported integer counts, so rounding can differ from Bracken's
`fraction_total_reads`. Blank support counts mean unavailable, never observed
zero. A missing Kraken row or a legacy report without support columns remains
unavailable. Zero denominators produce null canonical fractions, not confidence
estimates or claims of absence. Literal names such as `NA` retain their identity.

## Methods and compatibility

`taxonomy_methods.json` records `schema_version: 1`, the report generator's
`metaquest_version`, count units and denominators,
requested rank, recorded classification settings, observed/model read lengths,
tool/database provenance and SHA-256 hashes of the two native source tables.
It reads tool/database identity from saved run provenance rather than substituting
currently installed tools. Unrecorded provenance is null. The final run metadata
is linked without hashing it before the pipeline's final metadata update.

`analysis_summary.json` and stage JSON include `evidence_schema_version: 1`,
`evidence_table`, `methods_file`, `methods`, `taxonomic_level`, `rank_assigned`,
`classified_not_in_bracken`, `classified_fraction` and `unclassified_fraction`.
The latter two are null for empty input. Legacy `classification_rate` retains
its previous zero-input behavior, and `species_assigned*` / `classified_above_species`
remain compatibility aliases. Their names do not prove a species rank or explain
why classified input was excluded from the final Bracken profile. That remainder
can include threshold-filtered taxa as well as higher-rank assignments.

Consolidated and taxonomy-stage HTML show up to 50 evidence rows with links to
the complete TSV and methods JSON. Generic stage TSV/JSON retain all metrics and
legacy aliases. Top-taxon limits and plot display thresholds do not filter the
full evidence table. The plot remainder labels describe other Bracken-assigned
taxa, classified input outside the Bracken profile, and unclassified input.

These are research estimates, not confidence scores, cell counts, pathogenicity
findings or clinical predictions. Related references, database scope and settings
can affect calls. No reported match does not establish biological absence.

## Verification and benchmark scope

Regressions cover ID-based joins with identical display names, read/fragment
units, all-input/profile denominators, top-list truncation, missing support,
legacy reports, empty input, genus/family ranks, escaped labels, invalid counts,
duplicate IDs and independent detection of altered evidence/source hashes.
`tools/verify_run.py` cross-checks schema-1 evidence against native tables using
an independent CSV-based check. Existing runs without the new namespace retain
the previous acceptance path.

The scoped [benchmark protocol](../benchmarks/taxonomy-evidence/protocol.json)
uses the cached 200,000-pair ERR2984773 development subset, the existing
Standard-8 database, four threads and memory mapping. Native classifier outputs
are compared against frozen main, then reporting is measured three times per
revision on identical source outputs/settings. Results and limitations belong in
the [benchmark record](../benchmarks/taxonomy-evidence/RESULTS.md).

This is not the frozen 8,785,731-pair biological benchmark or a held-out sample.
The original full-input precision gate failure remains recorded in
[Zymo results](../benchmarks/zymo-even/RESULTS.md). No classifier defaults or
acceptance thresholds are tuned against it, and this reporting change does not
claim to improve sensitivity, specificity, abundance accuracy or checkpoint 3.

Definitions follow the [Kraken2 report specification](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown)
and [Bracken's output calculations](https://github.com/jenniferlu717/Bracken/blob/master/src/est_abundance.py).
