# workflow/

Snakemake rules for the pipeline, split into three separate Snakefiles that are run in
order (each Snakefile has its own `rule all` and is invoked independently — they are
not `include:`d into one master file). All commands below are run from the repository
root.

## Run order

### 0. `Snakefile_sra_download` — fetch raw reads from SRA
Environment: `be_analysis`

```bash
snakemake -s workflow/Snakefile_sra_download --cores 8
```

Reads `config/sra_accession_map.tsv` (sample → SRR accession), then for each sample:
`prefetch`s the run and `fasterq-dump`s + gzips it into
`data/seq_data/{sample}/{sample}_R1.fastq.gz` / `_R2.fastq.gz`. Skip this step if you
already have raw FASTQs in that layout (see [`../data/README.md`](../data/README.md)).

### 1. `Snakefile_Fastqc` — QC and read merging
Environment: `be_preprocessing` ([`../envs/preprocessing.yaml`](../envs/preprocessing.yaml))

```bash
snakemake -s workflow/Snakefile_Fastqc --cores 8
```

For each sample in `config/samples.tsv`:
- `fastqc_raw` — FastQC on raw R1/R2 (`results/qc/raw/`)
- `merge_reads` — merges R1/R2 with `pear` into `results/merged/{sample}.merged.fastq.gz`
- `fastqc_merged` — FastQC on the merged reads (`results/qc/merged/`)
- `multiqc` — aggregates all FastQC reports into `results/qc/multiqc_report.html`

### 2. `Snakefile_demux` — replicate demultiplexing
Environment: `be_preprocessing`

```bash
snakemake -s workflow/Snakefile_demux --cores 8
```

Requires `results/merged/{sample}.merged.fastq.gz` from step 1. For each sample:
- `demux_replicates` — uses `cutadapt` to split the merged FASTQ by the 5′-anchored
  barcodes in `config/barcodes.tsv`, writing one FASTQ per replicate
  (`results/replicates/{sample}/{sample}_{rep}.fastq.gz`) plus an `_unassigned.fastq.gz`
- `count_replicate_fastqs` — read counts per replicate/unassigned file, via
  `workflow/scripts/count_fastqs.py`
- `combine_read_count_summaries` — merges all per-sample counts into
  `results/replicates/summaries/combined_read_counts.tsv`, via
  `workflow/scripts/combine_read_counts.py`

### 3. `Snakefile_crispresso` — editing quantification
Environment: `be_crispresso` ([`../envs/crispresso.yaml`](../envs/crispresso.yaml))

```bash
snakemake -s workflow/Snakefile_crispresso --cores 8
```

Requires the per-replicate FASTQs from step 2. For each sample × replicate:
- `crispresso_pooled` — runs `CRISPRessoPooled --base_editor_output` against
  `config/amplicons.txt`, writing to `results/crispresso/{sample}/rep{N}/`
- `aggregate_results` — collects every sample's `SAMPLES_QUANTIFICATION_SUMMARY.txt`
  into `results/crispresso/summary/editing_summary.tsv`
- `crispresso_aggregate` — for each base editor × replicate, symlinks all dose folders
  for that editor into a staging dir and runs `CRISPRessoAggregate`, writing to
  `results/crispresso_aggregated/{base_editor}/rep{N}/`

### 4. Plotting (not a Snakemake rule)
`workflow/scripts/crispresso_analysis.R` consumes the tracked, curated tables in
[`../results/crispresso_tables/`](../results/crispresso_tables/README.md) (a hand-picked
copy/subset of the `Snakefile_crispresso` output) and produces all summary figures.
Run manually — see [`workflow/scripts/README.md`](scripts/README.md) for the exact
folder layout it expects.

## scripts/
Helper Python and R scripts called by the rules above (or run standalone for
plotting). See [`workflow/scripts/README.md`](scripts/README.md).
