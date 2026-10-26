#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.
"""
Pools the members of existing family MSAs with a sample's input sequences, so that updating
families searches both again.

MSAs may be Stockholm or aligned FASTA, plain or gzipped. Members are degapped and
upper-cased (Stockholm insert columns are lower case). A name `seq/start-end` is the slice
start..end of protein `seq`; a name without a range is the whole protein. A member is dropped
if its region is contained in an input sequence of the same protein (the input is the newer
copy) or in another kept member of it, so exact and nested duplicates collapse, while partial
overlaps and separate regions (e.g. two domains) are kept.

Writes the uncompressed pool (HMMER cannot search a gzip stream): the input sequences
first, then the kept members. Optionally also writes every member of each MSA, degapped, as
`<family>.faa.gz`: the FASTA of a family that passes the update as given.
"""

import sys
import gzip
import argparse
import math
import re
from pathlib import Path
from typing import Iterator, Sequence

GAPS = str.maketrans("", "", "-.~")
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
        "--msas",
        required=True,
        metavar="FILE",
        nargs="+",
        type=str,
        help="Existing family MSAs (Stockholm or aligned FASTA, plain or gzipped), or folders of them.",
    )
    parser.add_argument(
        "-o",
        "--out_fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Output pool of input sequences and kept MSA members, in FASTA format.",
    )
    parser.add_argument(
        "--out_members",
        metavar="DIR",
        type=Path,
        help="Optional output folder for each MSA's degapped members, as <family>.faa.gz.",
    )
    return parser.parse_args(args)


def open_text(path: Path):
    return (gzip.open if path.name.endswith(".gz") else open)(path, "rt")


def region(name: str) -> tuple[str, int, float]:
    """The protein a name belongs to and the residue range it covers (whole protein if no range)."""
    match = RANGE.match(name)
    if not match:
        return name, 1, math.inf
    start, end = sorted((int(match[2]), int(match[3])))
    return match[1], start, end


def contained(start: float, end: float, ranges: list[tuple[float, float]]) -> bool:
    return any(s <= start and end <= e for s, e in ranges)


def msa_rows(text: str) -> Iterator[tuple[str, str]]:
    """Yield (name, aligned sequence) rows of a Stockholm or aligned FASTA MSA."""
    if text.lstrip().startswith("# STOCKHOLM"):
        rows: dict[str, list[str]] = {}  # interleaved blocks: one line per row and block
        for line in text.splitlines():
            if line.strip() and not line.startswith(("#", "//")):
                name, seq = line.split(maxsplit=1)
                rows.setdefault(name, []).append(seq.strip())
        for name, chunks in rows.items():
            yield name, "".join(chunks)
    else:
        for record in text.split(">")[1:]:
            header, _, seq = record.partition("\n")
            yield header.split(maxsplit=1)[0], seq.replace("\n", "").replace("\r", "")


def pool_members(fasta: str, msas: Sequence[str], out_fasta: str, out_members: Path | None = None) -> None:
    files = []
    for msa in map(Path, msas):
        files.extend(sorted(p for p in msa.iterdir() if p.is_file()) if msa.is_dir() else [msa])

    with open(out_fasta, "w") as out:
        # Streamed: the input can be far larger than memory, only its names are kept
        input_ranges: dict[str, list[tuple[float, float]]] = {}
        line = "\n"
        with open_text(Path(fasta)) as f:
            for line in f:
                if line.startswith(">"):
                    protein, start, end = region(line[1:].split(maxsplit=1)[0])
                    input_ranges.setdefault(protein, []).append((start, end))
                out.write(line)
        if not line.endswith("\n"):
            out.write("\n")

        # protein -> [(start, end, file order, name, residues)]
        members: dict[str, list[tuple[int, float, int, str, str]]] = {}
        order = 0
        for path in files:
            with open_text(path) as f:
                rows = [(name, seq.translate(GAPS).upper()) for name, seq in msa_rows(f.read())]
            rows = [(name, residues) for name, residues in rows if residues]
            for name, residues in rows:
                protein, start, end = region(name)
                members.setdefault(protein, []).append((start, end, order, name, residues))
                order += 1
            if out_members:
                family = Path(path.name.removesuffix(".gz")).stem
                # mtime=0: identical members give identical files across runs
                with gzip.GzipFile(out_members / f"{family}.faa.gz", "wb", mtime=0) as out_family:
                    out_family.write("".join(f">{name}\n{residues}\n" for name, residues in rows).encode())

        kept = []
        for protein, group in members.items():
            kept_ranges: list[tuple[float, float]] = []
            # Leftmost and widest first: a region can then only be contained in one already kept
            for start, end, order, name, residues in sorted(group, key=lambda m: (m[0], -m[1], m[2])):
                if contained(start, end, input_ranges.get(protein, [])) or contained(start, end, kept_ranges):
                    continue
                kept_ranges.append((start, end))
                kept.append((order, name, residues))
        for _, name, residues in sorted(kept):  # file order
            out.write(f">{name}\n{residues}\n")
    print(f"Pooled {len(kept)} existing family members with the input sequences.")


def main(args: Sequence[str] | None = None) -> None:
    args = parse_args(args)
    if args.out_members:
        args.out_members.mkdir(parents=True, exist_ok=True)
    pool_members(args.fasta, args.msas, args.out_fasta, args.out_members)


if __name__ == "__main__":
    sys.exit(main())
