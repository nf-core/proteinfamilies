#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.
"""
Concatenates seed MSA files for a given list of family IDs into a single merged alignment
file, used as the starting point for rebuilding a merged family HMM. The same membership
is also written as a headerless representative<TAB>member TSV, the cluster format the
iterative family-generation algorithm rebuilds families from.
"""

import sys
import os
import gzip
import argparse
from typing import Sequence


def parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
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
    parser.add_argument(
        "-c",
        "--out_clusters",
        required=True,
        metavar="FILE",
        type=str,
        help="Path to output headerless representative<TAB>member TSV of the merged membership.",
    )
    return parser.parse_args(args)


def alignment_stem(filename: str) -> str:
    """
    Strip the alignment extension from a seed MSA filename, so it can be matched against a
    family ID. Handles the compressed alignments (e.g. '.fas.gz') that the iterative
    family-generation algorithm produces alongside the plain ones of the standard algorithm.

    Args:
        filename (str): Basename of a seed MSA file.

    Returns:
        str: The family ID the file belongs to.
    """
    if filename.endswith(".gz"):
        filename = filename[: -len(".gz")]
    return os.path.splitext(filename)[0]


def keep_non_empty_records(content: str) -> list[str]:
    """
    Split a FASTA alignment into records, dropping any that contain no residues.

    Gap-trimming (ClipKIT) removes columns but never rows, so a member whose residues all
    fell inside the trimmed columns survives as an all-gap record. Such records carry no
    sequence and abort the downstream realignment (FAMSA), so they are filtered out here.

    Args:
        content (str): Contents of a FASTA-format alignment file.

    Returns:
        list[str]: Records that still hold at least one residue, as FASTA blocks.
    """
    kept = []

    for record in content.split(">")[1:]:
        header, _, sequence = record.partition("\n")
        # '-' and '.' are both gap characters in the alignment formats handled here
        if sequence.replace("-", "").replace(".", "").strip():
            kept.append(f">{header}\n{sequence.strip()}")

    return kept


def merge_selected_alignments(family_ids: list[str], folder: str, out_file: str, out_clusters: str) -> None:
    """
    Merge seed alignments whose basename matches one of the requested family IDs.

    Args:
        family_ids (list[str]): Family IDs whose alignments should be concatenated.
        folder (str): Directory containing per-family seed alignments.
        out_file (str): Output path for the merged alignment file.
        out_clusters (str): Output path for the representative<TAB>member TSV.
    """
    merged_contents = []
    merged_files = 0
    dropped_records = 0

    for fam in family_ids:
        fname = next((os.path.join(folder, f) for f in os.listdir(folder) if alignment_stem(f) == fam), None)
        if fname and os.path.exists(fname):
            open_func = gzip.open if fname.endswith(".gz") else open
            with open_func(fname, "rt") as f:
                content = f.read().strip()
                if content:
                    kept = keep_non_empty_records(content)
                    dropped_records += content.count(">") - len(kept)
                    if kept:
                        merged_contents.extend(kept)
                        merged_files += 1
        else:
            print(f"[WARNING] No seed alignment found for family: {fam}", file=sys.stderr)

    if dropped_records:
        print(f"[WARNING] Dropped {dropped_records} all-gap record(s) while merging.", file=sys.stderr)

    if merged_contents:
        with open(out_file, "w") as out:
            # Literal concatenation of FASTA blocks — valid for any FASTA-format alignment.
            out.write("\n".join(merged_contents) + "\n")
        # The merged families become one cluster; its representative is the first member,
        # which is itself listed as a member, as in MMseqs2 cluster output.
        names = [record[1:].split()[0] for record in merged_contents]
        with open(out_clusters, "w") as out:
            for name in names:
                out.write(f"{names[0]}\t{name}\n")
        print(f"[INFO] Merged {merged_files} files into {out_file}")
    else:
        print("[WARNING] No files merged (none matched the provided list).", file=sys.stderr)


def main(args: Sequence[str] | None = None) -> None:
    args = parse_args(args)

    family_ids = [x.strip() for x in args.list.split(",") if x.strip()]

    merge_selected_alignments(family_ids, args.folder, args.out_file, args.out_clusters)


if __name__ == "__main__":
    sys.exit(main())
