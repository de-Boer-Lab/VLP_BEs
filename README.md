# VLP Base Editor Collection Benchmarking

Data analysis pipeline for benchamrking base editors delivered via engineered virus-like
particles (VLPs) across three genomic targets (**HEK3**, **B2M**, **PDCD1**), six VLP
doses (**D0** = 0 µL – **D5** = 10 µL), and three replicates. Raw amplicon sequencing
reads are downloaded from SRA, merged, QC'd, demultiplexed by replicate barcode,
quantified with CRISPResso2, and summarized into editing/substitution/deletion tables
and figures.

## Pipeline overview

The pipeline runs as four separate Snakefiles, in order, followed by an R plotting
script:

| Step | File | What it does |
|---|---|---|
| 0 | `workflow/Snakefile_sra_download` | Downloads raw FASTQs from SRA into `data/seq_data/` using `config/sra_accession_map.tsv` |
| 1 | `workflow/Snakefile_Fastqc` | Raw-read FastQC, merges paired-end reads with `pear`, MultiQC report |
| 2 | `workflow/Snakefile_demux` | Splits each merged per-sample FASTQ into per-replicate FASTQs via 5′-anchored barcode matching (`cutadapt`) |
| 3 | `workflow/Snakefile_crispresso` | Runs `CRISPRessoPooled` per sample × replicate, then `CRISPRessoAggregate` per base editor, and builds a combined editing summary table |
| 4 | `workflow/scripts/crispresso_analysis.R` | Loads the CRISPResso output tables and produces all summary figures (PDF) |

See each subfolder's own README for details:

- [`config/`](config/README.md) — sample sheet, SRA accession map, barcode map, amplicon/target definitions, pipeline parameters
- [`workflow/`](workflow/README.md) — Snakefiles and helper scripts
- [`data/`](data/README.md) — raw sequencing data
- [`results/`](results/README.md) — pipeline outputs

## Setup

Two conda environments are used, plus `sra-tools` for the download step:

```bash
conda env create -f envs/be_analysis.yml
```

## Running the pipeline

Run each Snakefile from the repository root, activating the matching environment:

```bash
conda activate be_analysis

# 0. Download raw reads from SRA
snakemake -s workflow/Snakefile_sra_download --cores 8

# 1. QC + read merging
snakemake -s workflow/Snakefile_Fastqc --cores 8

# 2. Replicate demultiplexing
snakemake -s workflow/Snakefile_demux --cores 8

# 3. CRISPResso quantification
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

`results/crispresso_tables/` contains the curated set of CRISPResso output
tables, so the figures can be regenerated directly without rerunning
steps 0–3.
