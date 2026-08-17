# workflow/scripts/

Helper scripts called by the Snakefiles' `shell:` directives, plus the standalone R
plotting script. All Python scripts are CLI tools (`argparse`) — run `<script>.py -h`
for the full flag list.

## `count_fastqs.py`
Called by `Snakefile_demux`'s `count_replicate_fastqs` rule. Counts reads in a set of
gzipped FASTQs (line count / 4) and writes a `replicate` / `read_count` TSV.

```bash
python workflow/scripts/count_fastqs.py \
  --inputs results/replicates/SAMPLE/SAMPLE_rep1.fastq.gz ... \
  --output results/replicates/SAMPLE/SAMPLE.replicate_counts.tsv
```

## `combine_read_counts.py`
Called by `Snakefile_demux`'s `combine_read_count_summaries` rule. Merges every
sample's per-replicate count TSV (output of `count_fastqs.py`) with
`config/samples.tsv` metadata (`base_editor`, `dose`) into one combined table.

```bash
python workflow/scripts/combine_read_counts.py \
  --sample-sheet config/samples.tsv \
  --replicate-counts results/replicates/*/*.replicate_counts.tsv \
  --output results/replicates/summaries/combined_read_counts.tsv
```

## `combine_cispresso_summary.py`
**Not called by any current Snakefile rule.** It expects a config TSV with
`sample`/`target`/`replicate`/`fastq`/`amplicon_seq`/`guide_seq`/`quant_window_center`/
`quant_window_size` columns and checks for CRISPResso output at
`results/crispresso/{sample}/{target}/{replicate}/...` — a per-target directory layout
that doesn't match the sample × replicate layout `Snakefile_crispresso` actually
produces (`results/crispresso/{sample}/rep{N}/...`). Kept for reference; the equivalent,
currently-used logic lives inline in `Snakefile_crispresso`'s `aggregate_results` rule.

## `crispresso_analysis.R`
Standalone plotting script (not invoked by Snakemake). Loads the CRISPResso summary
tables and produces every figure used for base-editor comparison. Requires R with
`tidyverse`, `tools`, `ggh4x`, `patchwork`, `ggtext`, and `RColorBrewer`.

**Expected working directory / folder layout:** the script sets

```r
data_dir  <- "./Data"
plots_dir <- "./plots"
```

relative to wherever R's working directory is when you `source()` it — it expects
`./Data/modifications/`, `./Data/substitutions/`, and `./Data/editing_summary.tsv` to
exist, and writes PDFs to `./plots/`. This repo's tracked data lives directly under
[`results/crispresso_tables/`](../../results/crispresso_tables/README.md) (no `Data/`
subfolder). To run it as-is, either:

- `cd results/crispresso_tables && ln -s . Data` (one-time symlink) and make sure
  `plots/` exists there, then `Rscript ../../workflow/scripts/crispresso_analysis.R`; or
- edit `data_dir`/`plots_dir` at the top of the script to point at
  `results/crispresso_tables` and `results/crispresso_tables/plots` directly.

Outputs (all PDF, written to `plots_dir`):

| File | Content |
|---|---|
| `substitutions_total_compressed.pdf` | Per-position substitution % by dosage, adenine vs. cytosine editor panels |
| `D5_adenine_cytosine_heatmap.pdf` | Heatmap of D5 (10 µL) editing by position, adenine vs. cytosine editors |
| `abe_cbe_composition.pdf` | Substitution base composition (A/T/C/G) at D5 |
| `substitutions_spacers_lineplot_compressed_D5.pdf` | D5 substitution frequency across the full spacer, all editors/targets |
| `substitutions_percentage_dosage.pdf` | % substitutions vs. VLP dose, per editor/target |
| `deletions_spacers_lineplot_compressed_D5.pdf` | D5 deletion frequency across the spacer |
| `substitutions_vs_dels_scatter.pdf` | % substitutions vs. % deletions scatter |
| `del_substitution_ratio_D5.pdf` | Deletion:substitution ratio at D5 |

The script also builds `sub_composition_summary` (canonical C→T vs. non-canonical CBE
editing outcomes) in-memory but does not currently write it to a file.
