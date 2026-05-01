#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import sys
import os
import argparse


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="Merge selected seed MSA files based on a provided list of family IDs."
    )
    parser.add_argument(
        "-l",
        "--list",
        required=True,
        metavar="LIST",
        type=str,
        help="Comma-separated list of family IDs to merge, e.g. 'mgnifams_test_3,mgnifams_test_5'.",
    )
    parser.add_argument(
        "-f",
        "--folder",
        required=True,
        metavar="DIR",
        type=str,
        help="Folder containing seed MSA files (e.g. seed_msa/).",
    )
    parser.add_argument(
        "-o",
        "--out_file",
        required=True,
        metavar="FILE",
        type=str,
        help="Path to output merged alignment file.",
    )
    return parser.parse_args(args)


def merge_selected_alignments(family_ids, folder, out_file):
    """Merge all alignment files in 'folder' whose basenames (without extension)
    match one of the family IDs in 'family_ids'."""
    merged_contents = []

    for fam in family_ids:
        fname = next((os.path.join(folder, f) for f in os.listdir(folder) if os.path.splitext(f)[0] == fam), None)
        if os.path.exists(fname):
            with open(fname, "r") as f:
                content = f.read().strip()
                if content:
                    merged_contents.append(content)
        else:
            print(f"[WARNING] File not found: {fname}", file=sys.stderr)

    if merged_contents:
        with open(out_file, "w") as out:
            out.write("\n".join(merged_contents) + "\n")
        print(f"[INFO] Merged {len(merged_contents)} files into {out_file}")
    else:
        print("[WARNING] No files merged (none matched the provided list).", file=sys.stderr)


def main(args=None):
    args = parse_args(args)

    # parse list input
    if args.list.startswith("[") and args.list.endswith("]"):
        args.list = args.list[1:-1]
    family_ids = [x.strip() for x in args.list.split(",") if x.strip()]

    merge_selected_alignments(family_ids, args.folder, args.out_file)


if __name__ == "__main__":
    sys.exit(main())
