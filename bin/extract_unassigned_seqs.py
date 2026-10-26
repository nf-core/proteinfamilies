#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.
"""
Writes the input sequences that no updated family holds, so they can go to family creation.

A name `seq/start-end` is the slice start..end of protein `seq`; a name without a range is the
whole protein. An input sequence is assigned if some family member of the same protein lies
inside it, so a member cut from it (a hit, in parent coordinates) assigns it, while a member
from another part of the protein does not.

Family files and the input may be plain or gzipped FASTA. Writes a gzipped FASTA.
"""

import sys
import gzip
import argparse
import math
import re
from pathlib import Path
from typing import Sequence

RANGE = re.compile(r"^(.+)/(\d+)-(\d+)$")


def parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Input sequences in FASTA format, plain or gzipped.",
    )
    parser.add_argument(
        "-m",
        "--families",
        default=[],
        metavar="FILE",
        nargs="*",
        type=str,
        help="Family member FASTA files, plain or gzipped (none if no family was updated).",
    )
    parser.add_argument(
        "-o",
        "--out_fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Output gzipped FASTA of the input sequences no family holds.",
    )
    return parser.parse_args(args)


def open_text(path: str):
    return (gzip.open if path.endswith(".gz") else open)(path, "rt")


def region(name: str) -> tuple[str, int, float]:
    """The protein a name belongs to and the residue range it covers (whole protein if no range)."""
    match = RANGE.match(name)
    if not match:
        return name, 1, math.inf
    start, end = sorted((int(match[2]), int(match[3])))
    return match[1], start, end


def extract_unassigned(fasta: str, families: Sequence[str], out_fasta: str) -> None:
    members: dict[str, list[tuple[float, float]]] = {}
    for family in families:
        with open_text(family) as f:
            for line in f:
                if line.startswith(">"):
                    protein, start, end = region(line[1:].split(maxsplit=1)[0])
                    members.setdefault(protein, []).append((start, end))

    # Streamed: the input can be far larger than memory
    written = 0
    keep = False
    with open_text(fasta) as f, gzip.open(out_fasta, "wt") as out:
        for line in f:
            if line.startswith(">"):
                protein, start, end = region(line[1:].split(maxsplit=1)[0])
                keep = not any(start <= s and e <= end for s, e in members.get(protein, []))
                written += keep
            if keep:
                out.write(line if line.endswith("\n") else line + "\n")
    print(f"Written {written} unassigned input sequences to {out_fasta}")


def main(args: Sequence[str] | None = None) -> None:
    args = parse_args(args)
    extract_unassigned(args.fasta, args.families, args.out_fasta)


if __name__ == "__main__":
    sys.exit(main())
