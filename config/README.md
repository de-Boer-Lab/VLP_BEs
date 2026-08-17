# config/

Hand-maintained inputs shared by all four Snakefiles: sample sheet, SRA accessions,
barcodes, amplicons, and pipeline parameters.

## Files

### `config.yaml`
Central config, loaded via `configfile: "./config/config.yaml"`.

- `raw_base` — root directory of raw per-sample FASTQ pairs, `./data/seq_data`
- `mapq` — mapping quality threshold (currently unused downstream but reserved)
- `barcode_error_rate` — max mismatch rate `cutadapt` allows when matching replicate barcodes (default `0.17`)
- `replicates_dir` — where demuxed per-replicate FASTQs live (output of `Snakefile_demux`)
- `amplicon_file` — path to `amplicons.txt`, passed to `CRISPRessoPooled`
- `crispresso_outdir` / `crispresso_aggregate_outdir` — output roots for `CRISPRessoPooled` and `CRISPRessoAggregate`
- `replicates` — list of replicate indices expected in every sample folder (`1, 2, 3`)
- `crispresso_params` — passed straight through to `CRISPRessoPooled`: `min_reads_to_use_region`, `min_frequency_alleles_around_cut`, `plot_window_size`, `quantification_window_center`, `quantification_window_size`
- `crispresso_threads` — threads per `CRISPRessoPooled` job
- `sra_map` — optional override for the accession map path used by `Snakefile_sra_download` (defaults to `./config/sra_accession_map.tsv`)

### `samples.tsv`
Sample sheet, one row per sample (tab-separated): `sample`, `base_editor`, `dose`
(`D0`–`D5`). (11 base editors × 6 doses) + (2 base editors × 5 doses) = 76 rows.

### `sra_accession_map.tsv`
Maps `sample` → SRA run accession (`srr`) + `biosample`, used by
`Snakefile_sra_download` to fetch raw reads.

### `barcodes.tsv`
Maps replicate names to the 5′ barcode sequence used to demultiplex each sample's
merged FASTQ into per-replicate FASTQs (`rep1`, `rep2`, `rep3`).

### `amplicons.txt` / `targets.fa`
Amplicon + sgRNA definitions for the three targeted loci (`HEK3`, `PDCD1`, `B2M`).
`amplicons.txt` (tab-separated: `Amplicon_Name`, `Amplicon_Sequence`, `sgRNA`) is what
`CRISPRessoPooled --amplicons_file` actually reads; `targets.fa` is the same sequences
as FASTA, for reference.
