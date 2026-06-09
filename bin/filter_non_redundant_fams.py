#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.
"""
Copies non-redundant family files from an input folder to the current working directory,
skipping any file whose family ID appears in the --redundant_ids list.
"""

import os
import sys
import argparse
import shutil


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_folder",
        required=True,
        metavar="FOLDER",
        type=str,
        help="All family files (hmm, seed_msa, full_msa or fasta folder).",
    )
    parser.add_argument(
        "-r",
        "--redundant_ids",
        required=True,
        metavar="FILE",
        type=str,
        help="Text file with one redundant family ID per line.",
    )
    return parser.parse_args(args)


def read_redundant_ids(filepath):
    """
    Args:
        filepath (str): Path to the file listing redundant family IDs.

    Returns:
        set[str]: Unique non-empty family IDs to exclude.
    """
    with open(filepath) as f:
        return set(line.strip() for line in f if line.strip())


def filter_files(input_dir, redundant_ids):
    """
    Copy non-redundant files to './' (Nextflow work dir). Family ID is the filename
    prefix before the first dot, so 'sample_1.hmm.gz' → 'sample_1'.

    Args:
        input_dir (str): Directory containing family files to copy.
        redundant_ids (set[str]): Family IDs that should be skipped.
    """
    for file in os.listdir(input_dir):
        fam_id = file.split(".")[0]
        if fam_id not in redundant_ids:
            shutil.copy(
                os.path.join(input_dir, file),
                os.path.join("./", file)
            )


def main(args=None):
    args = parse_args(args)
    redundant_ids = read_redundant_ids(args.redundant_ids)

    filter_files(args.input_folder, redundant_ids)


if __name__ == "__main__":
    sys.exit(main())
