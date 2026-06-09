#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.
"""
Identifies redundant and similar protein families from an all-vs-all hmmsearch of family
representatives. Redundant families (fully contained by another) are flagged for removal;
similar families (partial overlap) are paired for potential merging.
"""

import sys
import pandas as pd
import argparse
from typing import Sequence


def parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
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
        "--skip_family_redundancy_removal",
        action="store_true",
        help="If set, skip filtering of similar families above redundancy threshold.",
    )
    parser.add_argument(
        "-r",
        "--redundancy_length_threshold",
        default=1.0,
        metavar="FLOAT",
        type=float,
        help="Minimum length percentage threshold of annotated domain (env) against query to keep.",
    )
    parser.add_argument(
        "-s",
        "--similarity_length_threshold",
        required=True,
        default=0.9,
        metavar="FLOAT",
        type=float,
        help="Secondary similarity length threshold; families above this but below main threshold will be reported.",
    )
    parser.add_argument(
        "--redundant_ids_file",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output file with redundant family ids.",
    )
    parser.add_argument(
        "--similar_ids_file",
        required=True,
        metavar="FILE",
        type=str,
        help="Output file for IDs of families appearing in similarities.csv.",
    )
    parser.add_argument(
        "--pairwise_similarities_file",
        default="similarities.csv",
        metavar="FILE",
        type=str,
        help="Output file for similar but non-redundant families.",
    )
    return parser.parse_args(args)


def remove_self_hits(
    domtbl_df: pd.DataFrame, representative_to_family: dict[str, str]
) -> pd.DataFrame:
    """
    Map target sequence IDs to their family IDs, then drop rows where a family's HMM
    found its own representative — these self-hits would otherwise be flagged as redundant.

    Args:
        domtbl_df (pandas.DataFrame): Parsed hmmsearch domain hits.
        representative_to_family (dict[str, str]): Mapping from representative ID to
            family ID.

    Returns:
        pandas.DataFrame: Filtered hit table with self-hits removed.
    """
    domtbl_df["target name"] = domtbl_df["target name"].map(representative_to_family)
    domtbl_df = domtbl_df[domtbl_df["target name"] != domtbl_df["query name"]]

    return domtbl_df


def filter_and_label_similar(
    domtbl_df: pd.DataFrame,
    redundancy_length_threshold: float,
    similarity_length_threshold: float,
    skip_family_redundancy_removal: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Return two DataFrames:
    - redundant candidates (>= redundancy_length_threshold)
    - similar but non-redundant (>= similarity_length_threshold and < redundancy_length_threshold)

    Args:
        domtbl_df (pandas.DataFrame): Parsed hmmsearch domain hits.
        redundancy_length_threshold (float): Coverage threshold for declaring redundancy.
        similarity_length_threshold (float): Coverage threshold for declaring similarity.
        skip_family_redundancy_removal (bool): Whether to suppress redundant-family output.

    Returns:
        tuple[pandas.DataFrame, pandas.DataFrame]: Redundant candidates and similar
            non-redundant hits.
    """
    domtbl_df["similarity_score"] = (
        (domtbl_df["env to"] - domtbl_df["env from"] + 1) / domtbl_df["qlen"]
    )

    if skip_family_redundancy_removal:
        redundant_df = pd.DataFrame(columns=domtbl_df.columns) # empty df
        similar_df = domtbl_df[
            domtbl_df["similarity_score"] >= similarity_length_threshold
        ]
    else:
        redundant_df = domtbl_df[
            domtbl_df["similarity_score"] >= redundancy_length_threshold
        ]
        similar_df = domtbl_df[
            (domtbl_df["similarity_score"] >= similarity_length_threshold)
            & (domtbl_df["similarity_score"] < redundancy_length_threshold)
        ]

    return redundant_df, similar_df


def process_redundant(
    redundant_df: pd.DataFrame,
    family_to_size: dict[str, int],
    redundant_ids_file: str,
    skip_family_redundancy_removal: bool,
) -> set[str]:
    """
    Determine which family to discard for each redundant pair: the smaller one is marked
    redundant. For equal-sized pairs, the alphabetically later name is chosen — this
    ensures deterministic, collision-free deduplication without marking both directions.
    Returns the set of redundant family names.

    Args:
        redundant_df (pandas.DataFrame): Redundant candidate hits.
        family_to_size (dict[str, int]): Mapping from family ID to family size.
        redundant_ids_file (str): Output path for redundant family IDs.
        skip_family_redundancy_removal (bool): Whether to suppress redundant-family output.

    Returns:
        set[str]: Family IDs selected for removal.
    """
    redundant_fam_names = set()

    if skip_family_redundancy_removal:
        return redundant_fam_names

    redundant_df = redundant_df.drop(columns=["qlen", "env from", "env to"])
    redundant_df["query size"] = redundant_df["query name"].map(family_to_size)
    redundant_df["target size"] = redundant_df["target name"].map(family_to_size)

    for _, row in redundant_df.iterrows():
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

    with open(redundant_ids_file, "w") as f:
        for name in sorted(redundant_fam_names):
            f.write(name + "\n")

    return redundant_fam_names


def process_similar(
    similar_df: pd.DataFrame,
    redundant_fam_names: set[str],
    pairwise_similarities_file: str,
    similar_ids_file: str,
) -> None:
    """
    Write pairwise similarity pairs (excluding already-redundant families) to CSV, and
    the set of all family IDs that appear in at least one similarity pair to a text file.
    Returns early without writing if no similar pairs remain.

    Args:
        similar_df (pandas.DataFrame): Similar but non-redundant hit pairs.
        redundant_fam_names (set[str]): Families already marked redundant.
        pairwise_similarities_file (str): Output CSV for surviving pairwise similarities.
        similar_ids_file (str): Output path for family IDs that appear in similarities.

    Returns:
        None: Outputs are written for their side effects only.
    """
    if similar_df.empty:
        return

    # remove any similarity rows involving families already marked redundant
    similar_df = similar_df[
        ~similar_df["query name"].isin(redundant_fam_names)
        & ~similar_df["target name"].isin(redundant_fam_names)
    ]

    if similar_df.empty:
        return

    similar_out = similar_df[["query name", "target name", "similarity_score"]].copy()
    similar_out = similar_out.rename(
        columns={"query name": "family_1", "target name": "family_2"}
    )
    similar_out.to_csv(pairwise_similarities_file, index=False)

    # collect all unique family IDs that appear at least once in similarities
    similar_ids = set(similar_out["family_1"]) | set(similar_out["family_2"])
    with open(similar_ids_file, "w") as f:
        for fam in sorted(similar_ids):
            f.write(f"{fam}\n")


def process_family_similarity(
    mapping: str,
    domtbl: str,
    redundancy_length_threshold: float,
    similarity_length_threshold: float,
    redundant_ids_file: str,
    similar_ids_file: str,
    pairwise_similarities_file: str,
    skip_family_redundancy_removal: bool,
) -> None:
    """
    Classify family relationships as redundant or similar from hmmsearch hits.

    Args:
        mapping (str): Metadata CSV containing family IDs, sizes, and representatives.
        domtbl (str): Whitespace-delimited hmmsearch domain table.
        redundancy_length_threshold (float): Coverage threshold for declaring redundancy.
        similarity_length_threshold (float): Coverage threshold for declaring similarity.
        redundant_ids_file (str): Output path for redundant family IDs.
        similar_ids_file (str): Output path for family IDs in similarity relationships.
        pairwise_similarities_file (str): Output CSV for similar non-redundant pairs.
        skip_family_redundancy_removal (bool): Whether to suppress redundant-family output.
    """
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

    redundant_df, similar_df = filter_and_label_similar(domtbl_df, redundancy_length_threshold, similarity_length_threshold, skip_family_redundancy_removal)

    redundant_fam_names = process_redundant(redundant_df, family_to_size, redundant_ids_file, skip_family_redundancy_removal)

    process_similar(similar_df, redundant_fam_names, pairwise_similarities_file, similar_ids_file)


def main(args: Sequence[str] | None = None) -> None:
    args = parse_args(args)

    # Pre-create both output files as empty so Nextflow sees them even if neither
    # redundant nor similar families are found and the functions return early.
    open(args.redundant_ids_file, "w").close()
    open(args.similar_ids_file, "w").close()

    process_family_similarity(
        args.mapping,
        args.domtbl,
        args.redundancy_length_threshold,
        args.similarity_length_threshold,
        args.redundant_ids_file,
        args.similar_ids_file,
        args.pairwise_similarities_file,
        args.skip_family_redundancy_removal
    )


if __name__ == "__main__":
    sys.exit(main())
