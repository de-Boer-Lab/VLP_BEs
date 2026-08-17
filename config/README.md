# config/

Configuration inputs for all three Snakefiles (`Snakefile_Fastqc`, `Snakefile_demux`,
`Snakefile_crispresso`). Nothing in this folder is generated — these are hand-maintained
inputs that define samples, barcodes, amplicons, and pipeline parameters.

## Files

### `config.yaml`
Central config loaded by every Snakefile via `configfile: "config/config.yaml"`.

- `raw_base` — root directory of raw per-sample FASTQ pairs (see [`../data/README.md`](../data/README.md))
- `mapq` — mapping quality threshold (currently unused downstream but reserved)
- `barcode_error_rate` — max mismatch rate `cutadapt` allows when matching replicate barcodes (default `0.17`)
- `replicates_dir` — where demuxed per-replicate FASTQs live (output of `Snakefile_demux`)
- `amplicon_file` — path to `amplicons.txt`, passed to `CRISPRessoPooled`
- `crispresso_outdir` / `crispresso_aggregate_outdir` — output roots for `CRISPRessoPooled` and `CRISPRessoAggregate`
- `replicates` — list of replicate indices expected in every sample folder (`1, 2, 3`)
- `crispresso_params` — passed straight through to `CRISPRessoPooled`: `min_reads_to_use_region`, `min_frequency_alleles_around_cut`, `plot_window_size`, `quantification_window_center`, `quantification_window_size`
- `crispresso_threads` — threads per `CRISPRessoPooled` job

### `samples.tsv`
Sample sheet, one row per sample (tab-separated). Columns:

| column | description |
|---|---|
| `sample` | Sample name, matches the raw FASTQ folder/file prefix, e.g. `ABE8e-SpRY_D0` |
| `base_editor` | Base editor name, e.g. `ABE8e-SpRY` |
| `dose` | VLP dose code `D0`–`D5` (D0 = untreated control; see `crispresso_analysis.R` for the D0–D5 → µL mapping) |

13 base editors × 6 doses = 78 samples.

### `barcodes.tsv`
Maps replicate names to the 5′ barcode sequence used to demultiplex each sample's
merged FASTQ into per-replicate FASTQs (`rep1`, `rep2`, `rep3`).

### `amplicons.txt`
Tab-separated amplicon definitions passed to `CRISPRessoPooled --amplicons_file`.
Columns: `Amplicon_Name`, `Amplicon_Sequence`, `sgRNA`. Defines the three targeted loci
— `HEK3`, `PDCD1`, `B2M` — and their respective sgRNA spacer sequence.

### `targets.fa`
The same three amplicon sequences in FASTA format (reference for primer/amplicon
design, not consumed directly by the Snakefiles). Note the `PDCD1` header is spelled
`PCDC1` here — a pre-existing typo, kept as-is since it isn't parsed by any script.
