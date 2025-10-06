#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import sys
import argparse
import pandas as pd
import networkx as nx


def parse_args(args=None):
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


def build_pools(similarities_csv, threshold):
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


def write_pools(pools, out_file):
    with open(out_file, "w") as f:
        for pool in pools:
            f.write(",".join(pool) + "\n")


def main(args=None):
    args = parse_args(args)

    pools = build_pools(args.input_csv, args.threshold)
    write_pools(pools, args.out_file)

    print(f"[INFO] Found {len(pools)} connected family pools.")
    print(f"[INFO] Results written to: {args.out_file}")


if __name__ == "__main__":
    sys.exit(main())
