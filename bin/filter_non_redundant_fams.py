#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

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
    with open(filepath) as f:
        return set(line.strip() for line in f if line.strip())


def filter_files(input_dir, redundant_ids):
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
