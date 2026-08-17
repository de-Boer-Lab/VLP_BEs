# data/

Raw sequencing input for the pipeline. `seq_data/` is populated by running
[`workflow/Snakefile_sra_download`](../workflow/README.md), which pulls each sample's
FASTQ pair from SRA using [`../config/sra_accession_map.tsv`](../config/sra_accession_map.tsv).
You can also point `raw_base` in [`../config/config.yaml`](../config/config.yaml) at
your own pre-downloaded data instead.

## `seq_data/`

Root directory referenced by `config.yaml`'s `raw_base: "./data/seq_data"`, read by
both `Snakefile_sra_download` and `Snakefile_Fastqc` via shared `get_r1`/`get_r2`
helpers:

```python
def get_r1(sample):
    return f"{RAW_BASE}/{sample}/{sample}_R1.fastq.gz"
def get_r2(sample):
    return f"{RAW_BASE}/{sample}/{sample}_R2.fastq.gz"
```

Expected layout — one subfolder per sample (matching the `sample` column in
[`../config/samples.tsv`](../config/samples.tsv)), each containing a gzipped
paired-end FASTQ pair:

```
data/seq_data/
├── ABE8e_D0/
│   ├── ABE8e_D0_R1.fastq.gz
│   └── ABE8e_D0_R2.fastq.gz
├── ABE8e_D1/
│   ├── ABE8e_D1_R1.fastq.gz
│   └── ABE8e_D1_R2.fastq.gz
...
```

76 sample folders total — 13 base editors × 6 doses (D0–D5), excpet AYBE and
SPACE-DeltaCBE (D1–D5, no D0), each with 3 replicates and three target loci
pooled together in one FASTQ pair — replicates are split out downstream by
`Snakefile_demux` using the barcodes in
[`../config/barcodes.tsv`](../config/barcodes.tsv), not by folder structure here.

`Snakefile_Fastqc`'s `merge_reads` rule reads directly from this folder and writes
merged output to `results/merged/`; nothing under `data/` is ever modified in place.
