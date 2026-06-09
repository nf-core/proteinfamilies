#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.
"""
Groups similar protein families into pools using graph connected-component analysis.
Each pool is a maximal set of transitively similar families that will be merged into
a single super-family.
"""

import sys
import argparse
import pandas as pd
import networkx as nx
from typing import Sequence


def parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Group connected families based on similarity relationships using NetworkX."
    )
    parser.add_argument(
        "-i",
        "--input_csv",
        required=True,
        metavar="FILE",
        type=str,
        help="Input CSV file with columns: family_1,family_2,similarity_score.",
    )
    parser.add_argument(
        "-o",
        "--out_file",
        required=True,
        metavar="FILE",
        type=str,
        help="Output file listing family pools (connected components).",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        default=0.0,
        metavar="FLOAT",
        type=float,
        help="Minimum similarity score to consider families connected (default: 0.0).",
    )
    return parser.parse_args(args)


def build_pools(similarities_csv: str, threshold: float) -> list[list[str]]:
    """
    Build a graph of pairwise family similarities and extract connected components as pools.

    Uses transitive closure: if A is similar to B and B is similar to C, all three are
    pooled even if A and C have no direct similarity relationship. The undirected graph
    treats similarity as symmetric regardless of which direction the hmmsearch hit landed.

    Args:
        similarities_csv (str): CSV file containing ``family_1``, ``family_2``, and
            ``similarity_score`` columns.
        threshold (float): Minimum similarity score required to keep an edge.

    Returns:
        list[list[str]]: Connected components sorted internally and by first family ID.
    """
    # Load CSV
    df = pd.read_csv(similarities_csv)

    # Filter edges by threshold if specified
    df = df[df["similarity_score"] >= threshold]

    # Build undirected graph
    G = nx.from_pandas_edgelist(df, "family_1", "family_2")

    # Extract connected components (each is a pool)
    pools = [sorted(list(component)) for component in nx.connected_components(G)]

    # Sort pools by their first element for consistent output
    pools = sorted(pools, key=lambda x: x[0])

    return pools


def write_pools(pools: list[list[str]], out_file: str) -> None:
    """
    Write connected family pools to a newline-delimited text file.

    Args:
        pools (list[list[str]]): Connected components to write.
        out_file (str): Destination path for the pool listing.
    """
    with open(out_file, "w") as f:
        for pool in pools:
            f.write(",".join(pool) + "\n")


def main(args: Sequence[str] | None = None) -> None:
    args = parse_args(args)

    pools = build_pools(args.input_csv, args.threshold)
    write_pools(pools, args.out_file)

    print(f"[INFO] Found {len(pools)} connected family pools.")
    print(f"[INFO] Results written to: {args.out_file}")


if __name__ == "__main__":
    sys.exit(main())
