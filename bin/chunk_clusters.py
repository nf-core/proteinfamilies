#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.
"""
Splits MMSeqs2 cluster output into per-cluster FASTA files. Clusters smaller than
--threshold are discarded; surviving clusters are written as numbered FASTA chunks
used as inputs to per-family alignment.
"""

import sys
import os
import argparse
from collections import defaultdict
import csv
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from Bio.SeqRecord import SeqRecord
from typing import Sequence


def parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--clustering", required=True, metavar="FILE", help="TSV clustering file input."
    )
    parser.add_argument(
        "-s", "--sequences", required=True, metavar="FILE", help="Initial sequences FASTA file."
    )
    parser.add_argument(
        "-t", "--threshold", required=True, metavar="INT", type=int, help="Minimum cluster size to keep."
    )
    parser.add_argument(
        "-p", "--threads", metavar="INT", type=int, default=4, help="Number of threads for parallel fasta writing (default: 4)."
    )
    parser.add_argument(
        "-o", "--out_folder", required=True, metavar="FOLDER", help="Name of the output folder to be created."
    )
    return parser.parse_args(args)


def collect_clusters(clustering_file: str, threshold: int) -> dict[str, list[str]]:
    """
    Read an MMSeqs2 TSV (col 0 = representative, col 1 = member) and return clusters at or
    above the size threshold. The representative is counted as a member of its own cluster,
    so a singleton cluster has size 1.

    Args:
        clustering_file (str): MMSeqs2 clustering TSV file.
        threshold (int): Minimum cluster size to retain.

    Returns:
        dict[str, list[str]]: Surviving clusters keyed by representative ID.
    """
    clusters = defaultdict(list)
    with open(clustering_file) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            rep, member = row
            clusters[rep].append(member)
    return {
        rep: members for rep, members in clusters.items() if len(members) >= threshold
    }


def load_sequences(fasta_file: str, needed_ids: set[str]) -> dict[str, SeqRecord]:
    """
    Load only the sequences that belong to surviving clusters, avoiding a full FASTA scan
    of sequences that will be discarded anyway.

    Args:
        fasta_file (str): Input FASTA containing all clustered sequences.
        needed_ids (set[str]): Sequence IDs that belong to retained clusters.

    Returns:
        dict[str, Bio.SeqRecord.SeqRecord]: Loaded sequences keyed by record ID.
    """
    sequences = {}
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in needed_ids:
                sequences[record.id] = record
    return sequences


def write_cluster(
    prefix: str,
    chunk_num: int,
    members: list[str],
    sequences: dict[str, SeqRecord],
    out_folder: str,
) -> str:
    """
    Write cluster members to '{prefix}_{chunk_num}.faa'. The sequential chunk_num (not
    the representative ID) names the file; downstream Nextflow tokenizes on '_' to
    extract this number as the chunk identifier.

    Args:
        prefix (str): Sample prefix derived from the clustering filename.
        chunk_num (int): Sequential cluster number used in the output filename.
        members (list[str]): Sequence IDs to include in the cluster FASTA.
        sequences (dict[str, Bio.SeqRecord.SeqRecord]): Loaded sequences keyed by ID.
        out_folder (str): Destination directory for per-cluster FASTA files.

    Returns:
        str: Path to the written cluster FASTA file.
    """
    output_file = os.path.join(out_folder, f"{prefix}_{chunk_num}.faa")
    with open(output_file, "w") as out_handle:
        for member in members:
            if member in sequences:
                SeqIO.write(sequences[member], out_handle, "fasta")
    return output_file


def main(args: Sequence[str] | None = None) -> None:
    args = parse_args(args)
    os.makedirs(args.out_folder, exist_ok=True)

    # Prefix comes from the clustering filename stem (e.g., 'sample_1' from 'sample_1.tsv'),
    # which becomes the sample name embedded in each output chunk filename.
    prefix = os.path.splitext(os.path.basename(args.clustering))[0]

    # Step 1: Parse clusters
    clusters = collect_clusters(args.clustering, args.threshold)
    print(f"Clusters filtered.")

    # Step 2: Collect sequence IDs only from clusters above threshold
    needed_ids = set()
    for members in clusters.values():
        needed_ids.update(members)
    print("Required sequence IDs collected.")

    # Step 3: Load only sequences of clusters above threshold
    sequences = load_sequences(args.sequences, needed_ids)
    print("Required sequences loaded.")

    # Step 4: Parallel output writing
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for chunk_num, (_, members) in enumerate(clusters.items(), 1):
            futures.append(executor.submit(write_cluster, prefix, chunk_num, members, sequences, args.out_folder))

    print(f"Done. {len(clusters)} clusters written to {args.out_folder}")


if __name__ == "__main__":
    sys.exit(main())
