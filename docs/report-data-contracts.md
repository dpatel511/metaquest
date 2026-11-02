# Report and data contracts

These contracts describe research outputs, not biological validation. Prefer
stable identifiers to display labels. Read `analysis_metadata.json` for version,
effective settings, input hashes, tool/database identity and run status.

## Pipeline outputs

- `analysis_summary.json`: stage-level evidence. Absent stage keys mean no
  evidence was recorded; they do not imply a negative finding. Eight
  `stage_reports` entries link HTML, JSON and TSV evidence. `completed_stages`
  describes the workflow, not proof of biological accuracy.
- `reports/NN_stage.json`: `stage`, `status`, `metrics`. Statuses such as
  `completed`, `skipped`, `not_run`, `not_recorded` and validation `passed`
  retain their recorded meaning. Completed screening with zero matches is
  distinct from skipped/unrequested screening. JSON null means unavailable or
  unassessed, never zero.
- Companion stage TSV: `metric`, `value`; nested paths use ` / ` and lists use
  one-based positions. Null displays as `not available`. This is a human-readable
  flattened view; use JSON for typed values. HTML shows up to 150 metrics and
  screening tables up to 100 rows; companion files retain all rows.
- Taxonomy: `count_unit` is `fragments` for paired/interleaved reads and `reads`
  for single-end data. `total_input` is the cleaned classification input.
  `kraken_classified + unclassified_reads = total_input`; `species_assigned`
  is the Bracken estimated-count sum. Fractions are 0–1, not percentages.
  `top_taxa` is a display subset, not the full profile. New `taxonomy_id`
  retains source IDs as strings (null for legacy inputs lacking that column).
  Compatibility aliases containing `reads` or `fragments` retain their old
  names; interpret their counts using `count_unit`.
- `plots/plot_data/taxonomy_top.tsv`: `taxon`, additive `taxonomy_id`,
  `percent_all_input` (0–100). Remainder rows have no taxonomy ID. Repeated
  display names include IDs to prevent categorical aggregation; never use the
  label as a join key. Other species, above-species and unclassified remainders
  make the displayed percentages sum to 100 for nonempty inputs.
- fastp totals and retention refer to individual reads, including in paired
  mode. Q20/Q30 are fractions of bases in JSON; plot TSVs use percentages.
  Validation sampling indicators and analyzed counts accompany sample estimates.
- Assembly lengths and N50/N90 are bp; L50/L90 are contig counts. Continuity
  alone does not establish assembly correctness.
- Predicted protein FASTA IDs, GFF3 `ID`, annotation `query` and abundance
  `gene_id` refer to the same genes. GFF3 and derived evidence use **1-based
  inclusive** intervals: length = end - start + 1. BED conversion requires
  subtracting one from start, leaving end unchanged.
- Annotation counts are null if functional annotation was skipped. Completed
  annotation may have zero annotated genes. Mapping/assignment rates use the
  recorded read/fragment unit. Gene TPM sums to one million for assigned data,
  zero with no assignments; functional term attribution can count a gene in
  multiple terms, so term totals need not sum to one million.
- Screening protein counts are evidence rows, potentially multiple per gene.
  IDs must occur in the predicted genes. Co-occurrence does not assign a protein
  to a detected species or establish pathogenicity, phenotype or absence.

## Comparison outputs

`sample_id` remains a string, including leading zeros. Matrix rows follow the
sample sheet. `feature_values.tsv` contains full aligned feature IDs and counts
or summed gene TPM; missing features are zero-filled by ID, never label.
`features.tsv` maps `feature_id` to display `label`. `relative_abundance.tsv`
divides each row by its full reported profile total, excluding unclassified or
unassembled material. `bray_curtis.tsv` is a symmetric sample-by-sample matrix
in [0,1] with a zero diagonal. QC fractions use 0–1; blank QC values mean
unavailable, shown as `not run` in HTML.

`group_tests.tsv` is present only for formal tests; its effect is comparison
minus reference median CLR, not an absolute fold change. `comparison_summary.json`
records the actual transformation, pseudocount, prevalence filtering, test,
contrast, FDR correction and alpha. Read its design assumptions and limitations.

`plots/composition.tsv` retains its display-label format; repeated labels are
now disambiguated with feature IDs. New `plots/composition_by_id.tsv` is the
machine-readable plotted top-feature matrix plus `Other features` remainder.
Use full `relative_abundance.tsv` for inference. Heatmap source data uses feature
IDs, and group-effect source rows preserve formal-test identifiers.

## Compatibility and offline reports

These changes are additive: existing primary matrices and provenance contracts
are unchanged. Older summaries may lack the new taxonomy IDs and older reports
lack composition-by-ID exports. Do not fabricate missing provenance or IDs to
force compatibility. Keep original runs intact and regenerate into a new
output directory when needed. Strict `compare`/`--resume` checks still apply.

All report images are embedded and styles are inline. Companion links require
the output directory structure; copying HTML alone preserves images but not its
linked source evidence. Numerical/display regressions compare values and IDs,
not timestamps or image bytes. Browser evidence is separate from biological
accuracy and must state the browser versions, viewport and print method.
