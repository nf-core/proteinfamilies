#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import sys
import argparse
import os
import gzip
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Input fasta file.",
    )
    parser.add_argument(
        "-d",
        "--domtbl",
        required=True,
        metavar="FILE",
        type=str,
        help="Domain summary annotations result from hmmsearch.",
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
        "-H",
        "--hits",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Name of the output folder with hit fasta files (one file per family, where the filename is the family id).",
    )
    parser.add_argument(
        "-n",
        "--non_hits",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output fasta file with the non hit sequences.",
    )
    return parser.parse_args(args)


def filter_sequences(domtbl, length_threshold):
    results = {}

    # Open the domtbl file (supporting gzip)
    open_func = gzip.open if domtbl.endswith(".gz") else open
    with open_func(domtbl, "rt") as file:
        for line in file:
            if line.startswith("#"):
                continue # Skip comments

            columns = line.split()
            try:
                qlen = float(columns[5])
                env_from = int(columns[19])
                env_to = int(columns[20])
                env_length = env_to - env_from + 1

                if env_length >= length_threshold * qlen:
                    sequence_name = columns[0]
                    query_name = columns[3]

                    if query_name not in results:
                        results[query_name] = set()
                    results[query_name].add(f"{sequence_name}/{env_from}-{env_to}")
            except (IndexError, ValueError):
                continue  # Skip malformed lines

    return results


# Open the file with gzip if it's gzipped, otherwise open normally
def open_fasta(file_path):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")  # Open gzipped file in text mode
    return open(file_path, "rt")  # Open plain text file


# Parse the file
def parse_fasta(file_path):
    with open_fasta(file_path) as file:
        return {record.id: record for record in SeqIO.parse(file, "fasta")}


def write_non_hit_sequences(filtered_sequences, sequences, non_hits):
    # Determine the non-hit sequences
    hit_sequence_names = {hit.split("/")[0] for hits in filtered_sequences.values() for hit in hits}
    non_hit_records = [
        record for name, record in sequences.items()
        if name not in hit_sequence_names
    ]

    # Write the non-hit sequences to a gzipped file
    with gzip.open(non_hits, "wt") as non_hits_file:  # 'wt' mode for text writing
        SeqIO.write(non_hit_records, non_hits_file, "fasta")
    print(f"Written {len(non_hit_records)} non-hit sequences to {non_hits}")


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


def write_family_fastas(results, sequences, output_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    for family, hits in results.items():
        family_records = []

        for hit in hits:
            try:
                sequence_name, env_from, env_to = validate_and_parse_hit_name(hit)

                # Get the original sequence
                original_record = sequences[sequence_name]

                # Extract the specific range (adjust indices for 0-based indexing)
                extracted_seq = original_record.seq[env_from-1:env_to]

                # Determine the new sequence ID
                if len(extracted_seq) == len(original_record.seq):
                    new_id = sequence_name  # Omit range if full-length
                else:
                    new_id = f"{sequence_name}/{env_from}-{env_to}"

                # Create a new SeqRecord for the extracted range
                new_record = SeqRecord(
                    Seq(extracted_seq),
                    id=new_id,
                    description=family
                )
                family_records.append(new_record)
            except KeyError:
                print(f"Sequence {sequence_name} not found in the input FASTA.", file=sys.stderr)
            except ValueError as e:
                print(e, file=sys.stderr)

        # Write the extracted sequences to a FASTA file for the family
        if family_records:
            family_fasta_path = os.path.join(output_dir, f"{family}.fasta")
            SeqIO.write(family_records, family_fasta_path, "fasta")
            print(f"Written {len(family_records)} sequences to {family_fasta_path}")


def filter_recruited(fasta, domtbl, length_threshold, hits, non_hits):
    filtered_sequences = filter_sequences(domtbl, length_threshold)
    sequences = parse_fasta(fasta)
    write_non_hit_sequences(filtered_sequences, sequences, non_hits)
    write_family_fastas(filtered_sequences, sequences, hits)


def main(args=None):
    args = parse_args(args)
    filter_recruited(
        args.fasta, args.domtbl, args.length_threshold, args.hits, args.non_hits
    )


if __name__ == "__main__":
    sys.exit(main())
