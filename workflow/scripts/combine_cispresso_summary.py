#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd


def expected_json(sample, target, replicate):
    return Path(f"results/crispresso/{sample}/{target}/{replicate}/CRISPResso_on_{replicate}/CRISPResso2_info.json")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.config, sep="\t")

    rows = []
    for _, row in df.iterrows():
        out_json = expected_json(row["sample"], row["target"], row["replicate"])
        rows.append({
            "sample": row["sample"],
            "target": row["target"],
            "replicate": row["replicate"],
            "fastq": row["fastq"],
            "amplicon_seq": row["amplicon_seq"],
            "guide_seq": row["guide_seq"],
            "quant_window_center": row["quant_window_center"],
            "quant_window_size": row["quant_window_size"],
            "crispresso_info_json": str(out_json),
            "exists": out_json.exists()
        })

    out_df = pd.DataFrame(rows)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()