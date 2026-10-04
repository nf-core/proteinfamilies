#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.
"""
Assigns searched sequences to existing families using hmmsearch domain hits. Sequences whose
domain envelope covers at least --length_threshold of the query HMM are written, cut to that
envelope, into per-family FASTA files.
"""

import sys
import argparse
import os
import gzip
import re
from io import TextIOWrapper
from typing import Sequence
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Searched fasta file, that hits are cut from.",
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
    return parser.parse_args(args)


def filter_sequences(domtbl: str, length_threshold: float) -> dict[str, set[str]]:
    """
    Parse an hmmsearch domain table and return hits passing the length threshold.

    The envelope region (cols 19–20) must cover at least length_threshold × qlen (col 5).
    Returns a dict of query_name (family ID) → set of "seqname/env_from-env_to" strings.

    Args:
        domtbl (str): Path to the hmmsearch domain table, optionally gzipped.
        length_threshold (float): Minimum envelope coverage ratio relative to query length.

    Returns:
        dict[str, set[str]]: Passing hits grouped by family ID.
    """
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
def open_fasta(file_path: str) -> TextIOWrapper:
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    return open(file_path, "rt")


def parse_fasta(file_path: str) -> dict[str, SeqRecord]:
    """
    Parse a FASTA file into a record dictionary keyed by sequence ID.

    Args:
        file_path (str): Path to a plain-text or gzipped FASTA file.

    Returns:
        dict[str, Bio.SeqRecord.SeqRecord]: Parsed sequences indexed by record ID.
    """
    with open_fasta(file_path) as file:
        return {record.id: record for record in SeqIO.parse(file, "fasta")}


def validate_and_parse_hit_name(hit: str) -> tuple[str, int, int]:
    """
    Validates and parses a hit string.
    The hit must contain a string, at least one '/', and a valid range (integer-integer) after the last '/'.

    Args:
        hit (str): The hit string to validate and parse.

    Returns:
        tuple[str, int, int]: Parsed sequence name and envelope coordinates.

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


def slice_name(name: str, start: int, end: int) -> str:
    """
    Name a slice `name/start-end`. A name that already holds a `/s-e` range is itself a slice,
    so the new range is given in its parent's coordinates: `seq/10-200` [3, 180] -> `seq/12-189`.
    """
    match = re.match(r"^(.*)/(\d+)-(\d+)$", name)
    if match:
        offset = int(match.group(2)) - 1
        return f"{match.group(1)}/{start + offset}-{end + offset}"
    return f"{name}/{start}-{end}"


def write_family_fastas(
    results: dict[str, set[str]],
    sequences: dict[str, SeqRecord],
    output_dir: str,
) -> None:
    """
    Write per-family FASTA files with sequences cropped to their hit envelope.

    Coordinates are 1-based (HMMER convention), converted to 0-based for slicing.
    If the extracted range spans the full sequence, the ID omits the /from-to suffix.
    One file is written per family, named <family_id>.fasta.

    Args:
        results (dict[str, set[str]]): Passing hits grouped by family ID.
        sequences (dict[str, Bio.SeqRecord.SeqRecord]): Input sequences keyed by ID.
        output_dir (str): Directory where per-family FASTA files are written.
    """
    os.makedirs(output_dir, exist_ok=True)

    for family, hits in results.items():
        family_records = []

        for hit in sorted(hits):  # sets iterate in a per-run order
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
                    new_id = slice_name(sequence_name, env_from, env_to)

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


def split_family_hits(fasta: str, domtbl: str, length_threshold: float, hits: str) -> None:
    """
    Write the hits passing the length threshold into per-family FASTA files.

    Args:
        fasta (str): Searched FASTA that hits are cut from.
        domtbl (str): Path to the hmmsearch domain table.
        length_threshold (float): Minimum envelope coverage ratio relative to query length.
        hits (str): Output directory for per-family hit FASTA files.
    """
    write_family_fastas(filter_sequences(domtbl, length_threshold), parse_fasta(fasta), hits)


def main(args: Sequence[str] | None = None) -> None:
    args = parse_args(args)
    split_family_hits(args.fasta, args.domtbl, args.length_threshold, args.hits)


if __name__ == "__main__":
    sys.exit(main())
