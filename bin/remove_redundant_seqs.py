#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import sys
import gzip
import argparse
import csv
from Bio import SeqIO


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
        "-o",
        "--out_fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output fasta file with family representative sequences.",
    )
    return parser.parse_args(args)


def extract_rep_sequences(clustering, sequences, out_fasta):
    # Read the clustering file and extract unique values from column 1
    unique_representatives = set()
    with open(clustering, "r") as tsv_file:
        reader = csv.reader(tsv_file, delimiter="\t")
        for row in reader:
            if row:  # Ensure the row is not empty
                unique_representatives.add(row[0])

    # Read the sequences file and filter for representatives
    matching_records = []
    with (
        gzip.open(sequences, "rt")
        if sequences.endswith(".gz")
        else open(sequences, "r")
    ) as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in unique_representatives:
                matching_records.append(record)

    # Write the matching sequences to the output fasta file
    with open(out_fasta, "w") as output_file:
        SeqIO.write(matching_records, output_file, "fasta")


def main(args=None):
    args = parse_args(args)
    extract_rep_sequences(args.clustering, args.sequences, args.out_fasta)


if __name__ == "__main__":
    sys.exit(main())
