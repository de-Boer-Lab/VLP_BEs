#!/usr/bin/env python3

import argparse
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-sheet", required=True)
    parser.add_argument("--replicate-counts", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    samples_df = pd.read_csv(args.sample_sheet, sep="\t")
    sample_meta = samples_df.set_index("sample")[["base_editor", "dose"]].to_dict("index")

    rep_rows = []
    for f in args.replicate_counts:
        # Path structure: results/replicates/{sample}/{sample}.replicate_counts.tsv
        sample = f.split("/")[2]

        df = pd.read_csv(f, sep="\t")
        for _, row in df.iterrows():
            rep_rows.append({
                "sample": sample,
                "base_editor": sample_meta[sample]["base_editor"],
                "dose": sample_meta[sample]["dose"],
                "replicate": row["replicate"],
                "read_count": int(row["read_count"]),
            })

    final_df = pd.DataFrame(rep_rows)
    final_df = final_df[
        ["sample", "base_editor", "dose", "replicate", "read_count"]
    ].sort_values(["sample", "replicate"])

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    final_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()