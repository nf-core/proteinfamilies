#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import sys
import os
import argparse
from collections import defaultdict
import csv
import pyfastx


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c",
        "--clustering",
        required=True,
        metavar="FILE",
        type=str,
        help="TSV clustering file input.",
    )
    parser.add_argument(
        "-s",
        "--sequences",
        required=True,
        metavar="FILE",
        type=str,
        help="Initial sequences FASTA file.",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        required=True,
        metavar="INT",
        type=int,
        help="Minimum cluster size to keep.",
    )
    parser.add_argument(
        "-o",
        "--out_folder",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Name of the output folder to be created.",
    )
    return parser.parse_args(args)


def collect_clusters(clustering_file, threshold):
    # Collect clusters with a size threshold, storing in a defaultdict
    clusters = defaultdict(list)

    with open(clustering_file) as f:
        csv_reader = csv.reader(f, delimiter="\t")
        for row in csv_reader:
            rep, member = row
            clusters[rep].append(member)

    # Filter clusters by threshold
    return {
        rep: members for rep, members in clusters.items() if len(members) >= threshold
    }


def main(args=None):
    args = parse_args(args)

    # Create output directory if it doesn't exist
    os.makedirs(args.out_folder, exist_ok=True)

    # Collect clusters that meet the threshold
    clusters = collect_clusters(args.clustering, int(args.threshold))

    # Use pyfastx for indexed FASTA access (no need to load sequences into memory)
    fasta = pyfastx.Fasta(args.sequences)

    # Write output files
    for chunk_num, (rep, members) in enumerate(clusters.items(), start=1):
        output_file = os.path.join(args.out_folder, f"{chunk_num}.fasta")
        with open(output_file, "w") as fasta_out:
            for member in members:
                if member in fasta:
                    fasta_out.write(f">{member}\n{fasta[member]}\n")


if __name__ == "__main__":
    sys.exit(main())
