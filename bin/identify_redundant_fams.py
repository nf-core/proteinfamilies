#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import sys
import pandas as pd
import argparse


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--mapping",
        required=True,
        metavar="FILE",
        type=str,
        help="CSV metadata mapping input.",
    )
    parser.add_argument(
        "-d",
        "--domtbl",
        required=True,
        metavar="FILE",
        type=str,
        help="TSV hmmsearch domtbl out results for filtering.",
    )
    parser.add_argument(
        "-l",
        "--length_threshold",
        required=True,
        metavar="FLOAT",
        type=float,
        help="Minimum length percentage threshold of annotated domain (env) against query to keep.",
    )
    parser.add_argument(
        "-o",
        "--out_file",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output file with redundant family ids.",
    )
    return parser.parse_args(args)


def remove_self_hits(domtbl_df, representative_to_family):
    domtbl_df["target name"] = domtbl_df["target name"].map(representative_to_family)
    domtbl_df = domtbl_df[domtbl_df["target name"] != domtbl_df["query name"]]

    return domtbl_df


def filter_by_length(domtbl_df, length_threshold):
    domtbl_df = domtbl_df[
        (domtbl_df["env to"] - domtbl_df["env from"] + 1) / domtbl_df["qlen"]
        >= length_threshold
    ]

    return domtbl_df


def remove_redundant_fams(mapping, domtbl, length_threshold, out_file):
    mapping_df = pd.read_csv(
        mapping, comment="#", usecols=["Family Id", "Size", "Representative Id"]
    )
    domtbl_df = pd.read_csv(
        domtbl, sep=r"\s+", comment="#", header=None, usecols=[0, 3, 5, 19, 20]
    ).rename(
        columns={
            0: "target name",
            3: "query name",
            5: "qlen",
            19: "env from",
            20: "env to",
        }
    )

    representative_to_family = dict(
        zip(mapping_df["Representative Id"], mapping_df["Family Id"])
    )
    family_to_size = dict(zip(mapping_df["Family Id"], mapping_df["Size"]))

    domtbl_df = remove_self_hits(domtbl_df, representative_to_family)
    domtbl_df = filter_by_length(domtbl_df, length_threshold)
    domtbl_df = domtbl_df.drop(columns=["qlen", "env from", "env to"])
    domtbl_df["query size"] = domtbl_df["query name"].map(family_to_size)
    domtbl_df["target size"] = domtbl_df["target name"].map(family_to_size)

    redundant_fam_names = set()
    for _, row in domtbl_df.iterrows():
        query = row["query name"]
        target = row["target name"]
        query_size = int(row["query size"])
        target_size = int(row["target size"])

        if query_size < target_size:
            redundant_fam_names.add(query)
        elif target_size < query_size:
            redundant_fam_names.add(target)
        else: # sizes equal, keep alphabetically first as non-redundant to avoid triangular deletions
            redundant_fam_names.add(max(query, target))

    with open(out_file, "w") as f:
        for name in sorted(redundant_fam_names):
            f.write(name + "\n")


def main(args=None):
    args = parse_args(args)

    remove_redundant_fams(
        args.mapping,
        args.domtbl,
        args.length_threshold,
        args.out_file,
    )


if __name__ == "__main__":
    sys.exit(main())
