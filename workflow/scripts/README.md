# workflow/scripts/

Helper scripts used by the Snakefiles, plus the standalone R plotting script.

- **`count_fastqs.py`** — counts reads in gzipped FASTQs; called by `Snakefile_demux`.
- **`combine_read_counts.py`** — merges per-sample replicate counts with sample
  metadata; called by `Snakefile_demux`.
- **`combine_cispresso_summary.py`** — not called by any current Snakefile rule; its
  expected input/output layout predates `Snakefile_crispresso`'s current structure.
  Kept for reference only.
- **`crispresso_analysis.R`** — standalone plotting script (run manually, not via
  Snakemake).

