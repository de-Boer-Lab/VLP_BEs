# results/crispresso_tables/

Curated CRISPResso2 output tables consumed by
[`crispresso_analysis.R`](../../workflow/scripts/README.md), covering all 13
benchmarked base editors × 3 targets × 3 replicates at every VLP dose.

- **`editing_summary.tsv`** — one row per sample × target, with aligned read counts
  and insertion/deletion/substitution counts (see `Snakefile_crispresso`'s
  `aggregate_results` rule for exactly how it's built).
- **`modifications/*.txt`** — one file per base editor × target × replicate; rows are
  per-position `Insertions`/`Deletions`/`Substitutions` counts across all doses (dose
  is encoded in each row's `Folder` value, not the filename).
- **`substitutions/*_subs.txt`** — same layout, rows are per-position `A`/`T`/`C`/`G`
  call frequencies (the base-composition matrix substitutions are computed from).
- **`plots/`** — empty until you run the R script; PDF figures are written here.
