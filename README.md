# VLP Base Editor Benchmarking

Snakemake pipeline benchmarking 12 base editors delivered via engineered virus-like
particles (VLPs) across three genomic targets (**HEK3**, **B2M**, **PDCD1**), six VLP
doses (**D0(0 uL)**–**D5(10 uL)**), and three replicates. Raw amplicon sequencing reads are merged,
QC'd, demultiplexed by replicate barcode, quantified with CRISPResso2, and summarized
into editing/substitution/deletion tables and figures.

## Pipeline overview

The pipeline runs as three separate Snakefiles, in order, followed by an R plotting
script:

| Step | File | What it does |
|---|---|---|
| 1 | `workflow/Snakefile_Fastqc` | Raw-read FastQC, merges paired-end reads with `pear`, MultiQC report |
| 2 | `workflow/Snakefile_demux` | Splits each merged per-sample FASTQ into per-replicate FASTQs via 5′-anchored barcode matching (`cutadapt`) |
| 3 | `workflow/Snakefile_crispresso` | Runs `CRISPRessoPooled` per sample × replicate, then `CRISPRessoAggregate` per base editor, and builds a combined editing summary table |
| 4 | `workflow/scripts/crispresso_analysis.R` | Loads the CRISPResso output tables and produces all summary figures (PDF) |

See each subfolder's own README for details:

- [`config/`](config/README.md) — sample sheet, barcode map, amplicon/target definitions, pipeline parameters
- [`envs/`](envs/README.md) — conda environments
- [`workflow/`](workflow/README.md) — Snakefiles and helper scripts
- [`data/`](data/README.md) — raw sequencing data
- [`results/`](results/README.md) — pipeline outputs

## Setup

Two conda environments are used (see [`envs/README.md`](envs/README.md) for exact
contents):

```bash
conda env create -f envs/be_vlp_crispresso.yaml   
```

## Running the pipeline

Run each Snakefile from the repository root, activating the matching environment:

```bash
# 1. QC + read merging
conda activate be_vlp_crispresso
snakemake -s workflow/Snakefile_Fastqc --cores 8

# 2. Replicate demultiplexing
conda activate be_vlp_crispresso
snakemake -s workflow/Snakefile_demux --cores 8

# 3. CRISPResso quantification
conda activate be_vlp_crispresso
snakemake -s workflow/Snakefile_crispresso --cores 8
```

Add `-n` to any command to do a dry run first. Paths and parameters (raw data
location, barcode error rate, CRISPResso quantification window, etc.) are set in
[`config/config.yaml`](config/config.yaml).

Once `results/crispresso_tables/` is populated, generate the figures with:

```r
# from within results/crispresso_tables/ — see workflow/scripts/README.md
# for the exact folder layout this script expects
source("../../workflow/scripts/crispresso_analysis.R")
```

## Data availability

Raw FASTQ files are not included in this repository (see
[`data/README.md`](data/README.md)). The tracked deliverable is
`results/crispresso_tables/` — the per-sample CRISPResso modification/substitution
tables and the combined `editing_summary.tsv` — which is sufficient to reproduce all
figures via the R script without re-running the upstream sequencing pipeline.
