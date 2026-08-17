# results/crispresso_tables/

Curated CRISPResso2 output tables, tracked in git so that
[`crispresso_analysis.R`](../../workflow/scripts/README.md) can regenerate all figures
without re-running the sequencing pipeline. Covers all 13 base editors × 3 targets ×
3 replicates at every VLP dose.

## `editing_summary.tsv`

One row per sample × target (i.e. per `CRISPRessoPooled` run). Columns:

| column | description |
|---|---|
| `Sample` | Sample name (e.g. `ABE8e-SpRY_D0`) |
| `Base_Editor` | Base editor name |
| `Dose` | VLP dose code, `D0`–`D5` |
| `Replicate` | `rep1` / `rep2` / `rep3` |
| `Name` | Genomic target — `HEK3`, `PDCD1`, or `B2M` |
| `Unmodified%` / `Modified%` | % of aligned reads without/with any edit |
| `Reads_total` / `Reads_aligned` | Total reads in / reads aligned to the amplicon |
| `Unmodified` / `Modified` / `Discarded` | Read counts |
| `Insertions` / `Deletions` / `Substitutions` | Read counts with each outcome (not mutually exclusive) |
| `Only Insertions` … `Insertions Deletions and Substitutions` | Read counts by exclusive outcome combination |

## `modifications/*.txt` and `substitutions/*.txt`

One file per sample × target × replicate, named `{base_editor}_{dose}_rep{N}_{target}.txt`
(e.g. `ABE8e_B2M_rep1.txt`) — this repo's copies are named without the dose segment
where all doses were pooled per file; check the `Folder` column inside each file for
the exact source run. Both are per-nucleotide-position CRISPResso quantification
tables, wide format: a `Folder`/`Modification`(or `Nucleotide`) ID column followed by
one column per position in the sgRNA spacer.

- **`modifications/`** — rows are `Insertions`, `Deletions`, `Substitutions` (and
  `Only_*` variants); values are per-position edit counts.
- **`substitutions/`** (`*_subs.txt`) — rows are `A`/`T`/`C`/`G`/`-` (the nucleotide
  observed at that position); values are per-position frequencies of that base being
  called, i.e. the raw base-composition matrix substitutions are computed from.

`crispresso_analysis.R` parses the `Folder` column to recover `base_editor`, `dosage`,
`replicate`, and `target`, reshapes both tables to long format, and combines them with
`editing_summary.tsv` to build every figure (see
[`../../workflow/scripts/README.md`](../../workflow/scripts/README.md) for the full
output list and how to invoke it).

## `plots/`

Empty in the repo — this is the output directory `crispresso_analysis.R` writes its
PDF figures to. Populate it by running the R script.
