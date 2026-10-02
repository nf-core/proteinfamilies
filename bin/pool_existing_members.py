#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.
"""
Pools the members of existing family MSAs with a sample's input sequences, so that updating
families searches both again.

MSAs may be Stockholm or aligned FASTA, plain or gzipped. Members are degapped and
upper-cased (Stockholm insert columns are lower case). A member is dropped if the input holds
the same sequence (same name once a trailing `/start-end` range is removed), since the input
is the newer copy, or if an earlier member has the same name.

Writes the uncompressed pool (HMMER cannot search a gzip stream): the input sequences
first, then the kept members.
"""

import sys
import gzip
import argparse
import re
from pathlib import Path
from typing import Iterator, Sequence

GAPS = str.maketrans("", "", "-.~")
RANGE = re.compile(r"/\d+-\d+$")


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
    return parser.parse_args(args)


def open_text(path: Path):
    return (gzip.open if path.name.endswith(".gz") else open)(path, "rt")


def base_name(name: str) -> str:
    """The sequence a `name/start-end` slice is cut from."""
    return RANGE.sub("", name)


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


def pool_members(fasta: str, msas: Sequence[str], out_fasta: str) -> None:
    files = []
    for msa in map(Path, msas):
        files.extend(sorted(p for p in msa.iterdir() if p.is_file()) if msa.is_dir() else [msa])

    with open(out_fasta, "w") as out:
        # Streamed: the input can be far larger than memory, only its names are kept
        input_bases, line = set(), "\n"
        with open_text(Path(fasta)) as f:
            for line in f:
                if line.startswith(">"):
                    input_bases.add(base_name(line[1:].split(maxsplit=1)[0]))
                out.write(line)
        if not line.endswith("\n"):
            out.write("\n")

        seen = set()
        for path in files:
            with open_text(path) as f:
                rows = msa_rows(f.read())
            for name, seq in rows:
                residues = seq.translate(GAPS).upper()
                if residues and name not in seen and base_name(name) not in input_bases:
                    seen.add(name)
                    out.write(f">{name}\n{residues}\n")
    print(f"Pooled {len(seen)} existing family members with the input sequences.")


def main(args: Sequence[str] | None = None) -> None:
    args = parse_args(args)
    pool_members(args.fasta, args.msas, args.out_fasta)


if __name__ == "__main__":
    sys.exit(main())
