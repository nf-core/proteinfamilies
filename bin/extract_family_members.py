#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## Modified version for extracting family members.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import sys
import os
import gzip
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import tempfile
import shutil


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--fasta_folder",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Input folder with amino acid FASTA files.",
    )
    parser.add_argument(
        "-t",
        "--num_threads",
        type=int,
        default=4,
        help="Number of threads to use for parallel processing.",
    )
    parser.add_argument(
        "-o",
        "--out_tsv",
        required=True,
        metavar="FILE",
        type=str,
        help="Output TSV file with family_id and member_id.",
    )
    return parser.parse_args(args)


def extract_ids(filepath):
    """
    Extract all sequence IDs from a FASTA file.

    Args:
        filepath (str): Path to the FASTA (.fa or .fa.gz) file.

    Returns:
        list of str: List of sequence IDs.
    """
    open_func = gzip.open if filepath.endswith(".gz") else open
    ids = []
    with open_func(filepath, "rt") as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def process_fasta_file(filename, fasta_folder, tmp_dir):
    """
    Process a single FASTA file and write family_id, member_id pairs to a temp TSV.

    Args:
        filename (str): Filename of the FASTA file.
        fasta_folder (str): Directory containing the FASTA files.
        tmp_dir (str): Temporary directory for worker outputs.

    Returns:
        str: Path to the temporary TSV file for this family.
    """
    filepath = os.path.join(fasta_folder, filename)
    family_id = os.path.splitext(os.path.splitext(filename)[0])[0]
    member_ids = extract_ids(filepath)

    tmp_path = os.path.join(tmp_dir, f"{family_id}.tsv")
    with open(tmp_path, "w") as tmp_out:
        for mid in member_ids:
            tmp_out.write(f"{family_id}\t{mid}\n")

    return tmp_path


def extract_family_members(fasta_folder, num_threads, out_tsv):
    """
    Process all FASTA files in a folder in parallel and write family/member IDs to a TSV.
    """
    all_files = sorted(os.listdir(fasta_folder))
    tmp_dir = tempfile.mkdtemp(prefix="fam_members_")

    try:
        with ProcessPoolExecutor(max_workers=num_threads) as executor:
            tmp_files = list(
                executor.map(
                    partial(process_fasta_file, fasta_folder=fasta_folder, tmp_dir=tmp_dir),
                    all_files
                )
            )

        # Concatenate temp files into final TSV
        with open(out_tsv, "w") as out:
            out.write("family_id\tmember_id\n")
            for tmpf in tmp_files:
                with open(tmpf, "r") as tf:
                    shutil.copyfileobj(tf, out)

    finally:
        # Cleanup temp dir
        shutil.rmtree(tmp_dir)


def main(args=None):
    args = parse_args(args)
    extract_family_members(args.fasta_folder, args.num_threads, args.out_tsv)


if __name__ == "__main__":
    sys.exit(main())