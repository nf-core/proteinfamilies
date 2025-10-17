#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import sys
import argparse
import gzip
import re
from Bio import SeqIO


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--domtbl",
        required=True,
        metavar="FILE",
        type=str,
        help="Domain summary annotations result from hmmsearch.",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Input fasta file containing all protein names with their amino acid sequence for mapping.",
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
        "--out_fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output fasta file with the fasta converted sequences (no gaps).",
    )
    return parser.parse_args(args)


def filter_sequences(domtbl, length_threshold):
    filtered_sequences = []

    with gzip.open(domtbl, "rt", encoding="utf-8") as file:
        for line in file:
            if line.startswith("#"):
                continue  # Skip comments

            columns = line.split()
            try:
                qlen = float(columns[5])
                env_from = int(columns[19])
                env_to = int(columns[20])
                env_length = env_to - env_from + 1

                if env_length >= length_threshold * qlen:
                    sequence_name = columns[0]

                    filtered_sequences.append(f"{sequence_name}/{env_from}-{env_to}")
            except (IndexError, ValueError):
                continue  # Skip malformed lines

    return filtered_sequences


def validate_and_parse_hit_name(hit):
    """
    Validates and parses a hit string.
    The hit must contain a string, at least one '/', and a valid range (integer-integer) after the last '/'.

    Args:
        hit (str): The hit string to validate and parse.

    Returns:
        tuple: (sequence_name, env_from, env_to) if the hit is valid.

    Raises:
        ValueError: If the hit is invalid.
    """
    # Define the regex pattern
    pattern = r"^(.*)/(\d+)-(\d+)$"

    # Match the pattern
    match = re.match(pattern, hit)
    if not match:
        raise ValueError(f"Skipping hit with invalid format: {hit}.")

    # Extract components
    sequence_name = match.group(1)  # Everything before the last '/'
    env_from = int(match.group(2))  # First integer in the range
    env_to = int(match.group(3))    # Second integer in the range

    return sequence_name, env_from, env_to


def extract_fasta_subset(filtered_sequences, fasta, out_fasta):
    if filtered_sequences:
        open_func = gzip.open if fasta.endswith(".gz") else open
        with open_func(fasta, "rt") as in_fasta:
            fasta_dict = {record.id: str(record.seq) for record in SeqIO.parse(in_fasta, "fasta")}

        with gzip.open(out_fasta, "wt") as out_file:
            for filtered_sequence in filtered_sequences:
                try:
                    sequence_name, env_from, env_to = validate_and_parse_hit_name(filtered_sequence)

                    # Get the original sequence
                    original_record = fasta_dict[sequence_name]

                    # Extract the specific range (adjust indices for 0-based indexing)
                    extracted_seq = original_record[env_from-1:env_to]

                    # Determine the new sequence ID
                    if len(extracted_seq) == len(original_record):
                        new_id = sequence_name  # Omit range if full-length
                    else:
                        new_id = f"{sequence_name}/{env_from}-{env_to}"

                    out_file.write(f">{new_id}\n{extracted_seq}\n")
                except KeyError:
                    print(f"Sequence {sequence_name} not found in the input FASTA.", file=sys.stderr)
                except ValueError as e:
                    print(e, file=sys.stderr)
    else:
        print("No filtered sequences remained to write out. Skipping out_fasta file creation.")


def filter_recruited(domtbl, fasta, length_threshold, out_fasta):
    filtered_sequences = filter_sequences(domtbl, length_threshold)
    extract_fasta_subset(filtered_sequences, fasta, out_fasta)


def main(args=None):
    args = parse_args(args)
    filter_recruited(args.domtbl, args.fasta, args.length_threshold, args.out_fasta)


if __name__ == "__main__":
    sys.exit(main())
