#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.
"""
Recalculates the `name/start-end` row coordinates (Pfam naming convention) of a trimmed
FASTA-format MSA.

Trimming removes alignment columns but keeps row names, so a row would claim residues it no
longer holds. The leading and trailing `trim` runs of the trimming log give the end columns
removed; each row's residues in them shift its range: `seq` -> `seq/(1+left)-(len-right)` and
`seq/s-e` -> `seq/(s+left)-(e-right)`. Interior removals (not ends-only trimming) are not
reflected in the range. Rows left with no residues are dropped.

Writes the renamed trimmed MSA and its degapped sequences, so both always match.
Plain string operations only (no Biopython): all per-row counting runs in C via str.count.
"""

import sys
import argparse
from typing import Sequence


def parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-u",
        "--untrimmed",
        required=True,
        metavar="FILE",
        type=str,
        help="Untrimmed MSA in FASTA format.",
    )
    parser.add_argument(
        "-t",
        "--trimmed",
        required=True,
        metavar="FILE",
        type=str,
        help="Trimmed MSA in FASTA format, rows in the same order as the untrimmed MSA.",
    )
    parser.add_argument(
        "-l",
        "--log",
        required=True,
        metavar="FILE",
        type=str,
        help="Trimming log with one '<position> <keep|trim> ...' line per untrimmed column (ClipKIT --log).",
    )
    parser.add_argument(
        "-o",
        "--out_msa",
        required=True,
        metavar="FILE",
        type=str,
        help="Output trimmed MSA with recalculated coordinates.",
    )
    parser.add_argument(
        "-f",
        "--out_fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Output degapped sequences of the output MSA.",
    )
    return parser.parse_args(args)


def read_fasta(path: str) -> list[tuple[str, str]]:
    """
    Read a (possibly line-wrapped) FASTA file in one pass.

    Returns:
        list of (header, sequence) tuples, header without '>'.
    """
    with open(path, "r") as f:
        content = f.read()
    records = []
    for record in content.split(">")[1:]:
        header, _, seq = record.partition("\n")
        records.append((header.rstrip(), seq.replace("\n", "").replace("\r", "")))
    return records


def end_trim_runs(log: str) -> tuple[int, int, int]:
    """
    Count the trimmed columns at each end of the untrimmed alignment.

    Returns:
        (left, right, width): leading and trailing `trim` run lengths, total columns.
    """
    with open(log, "r") as f:
        kept = [line.split(maxsplit=2)[1] == "keep" for line in f if line.strip()]
    width = len(kept)
    if True not in kept:  # everything trimmed
        return width, 0, width
    left = kept.index(True)
    right = kept[::-1].index(True)
    return left, right, width


def residues(seq: str, start: int = 0, end: int | None = None) -> int:
    """Number of non-gap characters in seq[start:end], counted in C."""
    end = len(seq) if end is None else end
    return (end - start) - seq.count("-", start, end) - seq.count(".", start, end)


def split_range(name: str) -> tuple[str, int | None, int | None]:
    """Split `seq/s-e` into (seq, s, e); names without a valid range return (name, None, None)."""
    base, sep, rng = name.rpartition("/")
    if sep:
        start, dash, end = rng.partition("-")
        if dash and start.isdigit() and end.isdigit():
            return base, int(start), int(end)
    return name, None, None


def recalculate(
    untrimmed: str, trimmed: str, log: str, out_msa: str, out_fasta: str
) -> None:
    left_cols, right_cols, width = end_trim_runs(log)
    before = read_fasta(untrimmed)
    after = read_fasta(trimmed)
    if len(before) != len(after):
        sys.exit(f"Row count differs: {len(before)} untrimmed vs {len(after)} trimmed.")

    msa_lines, fasta_lines = [], []
    for (header, seq), (trimmed_header, trimmed_seq) in zip(before, after):
        name, _, description = header.partition(" ")
        if trimmed_header.partition(" ")[0] != name:
            sys.exit(
                f"Row order differs: '{name}' untrimmed vs '{trimmed_header}' trimmed."
            )
        if len(seq) != width:
            sys.exit(f"Row '{name}' has {len(seq)} columns, the log has {width}.")

        degapped = trimmed_seq.replace("-", "").replace(".", "")
        if not degapped:  # no residues left
            continue

        base, start, end = split_range(name)
        if start is None:
            start, end = 1, residues(seq)
        start += residues(seq, 0, left_cols)
        end -= residues(seq, width - right_cols, width)

        new_header = f"{base}/{start}-{end}" + (
            f" {description}" if description else ""
        )
        msa_lines.append(f">{new_header}\n{trimmed_seq}")
        fasta_lines.append(f">{new_header}\n{degapped}")

    with open(out_msa, "w") as f:
        f.write("\n".join(msa_lines) + "\n" if msa_lines else "")
    with open(out_fasta, "w") as f:
        f.write("\n".join(fasta_lines) + "\n" if fasta_lines else "")


def main(args: Sequence[str] | None = None) -> None:
    args = parse_args(args)
    recalculate(args.untrimmed, args.trimmed, args.log, args.out_msa, args.out_fasta)


if __name__ == "__main__":
    sys.exit(main())
