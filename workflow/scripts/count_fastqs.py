#!/usr/bin/env python3

import argparse
import gzip
import os
import pandas as pd


def count_reads_in_fastq(path):
    n = 0
    with gzip.open(path, "rt") as fh:
        for _ in fh:
            n += 1
    if n % 4 != 0:
        raise ValueError(f"{path} does not look like a valid FASTQ: line count {n} not divisible by 4")
    return n // 4


def label_from_path(path):
    base = os.path.basename(path)
    if base.endswith(".fastq.gz"):
        return base[:-9]
    return base


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = []
    for fq in args.inputs:
        rows.append({
            "replicate": label_from_path(fq),
            "read_count": count_reads_in_fastq(fq)
        })

    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()