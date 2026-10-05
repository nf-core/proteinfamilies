#!/usr/bin/env python3
"""
Split existing family HMMs into one gzipped `<NAME>.hmm.gz` file per model.

The input is either a `.tar.gz` archive of HMM files or a single HMM library (e.g. the
pipeline's `<id>.lib.gz` or Pfam-A.hmm.gz). Files may be plain or gzipped and may hold several
models each. Every model becomes its own family, named by its NAME line, so families that
pass an update unchanged are gzipped like rebuilt and created ones.
"""

import argparse
import gzip
import io
import sys
import tarfile
from collections.abc import Iterator, Sequence
from pathlib import Path


def parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", required=True, type=Path, help="HMM .tar.gz archive or HMM library (plain or gzipped)")
    parser.add_argument("-o", "--outdir", required=True, type=Path, help="Output folder for the <NAME>.hmm.gz files")
    return parser.parse_args(args)


def lines(stream: io.BufferedReader) -> Iterator[str]:
    """Text lines of a stream, decompressed if it is gzipped (by magic bytes, not extension)."""
    data = gzip.GzipFile(fileobj=stream) if stream.peek(2)[:2] == b"\x1f\x8b" else stream
    for line in data:
        yield line.decode()


def input_lines(path: Path) -> Iterator[str]:
    """Lines of every HMM file in the input: each file member of a .tar.gz, or the library itself."""
    if path.name.endswith(".tar.gz"):
        with tarfile.open(path, "r|gz") as tar:
            for member in tar:
                if member.isfile():
                    yield from lines(tar.extractfile(member))
    else:
        with path.open("rb") as handle:
            yield from lines(handle)


def models(text: Iterator[str]) -> Iterator[tuple[str, list[str]]]:
    """(NAME, lines) of each model, ending at its `//` line."""
    model: list[str] = []
    name = None
    for line in text:
        model.append(line)
        if line.startswith("NAME"):
            name = line.split()[1]
        elif line.rstrip() == "//":
            if name is None:
                sys.exit("ERROR: an existing HMM has no NAME line.")
            yield name, model
            model, name = [], None
    if any(line.strip() for line in model):
        sys.exit(f"ERROR: existing HMM {name or '(no NAME)'} is truncated: no closing '//' line.")


def main(args: Sequence[str] | None = None) -> None:
    args = parse_args(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    seen = set()
    for name, model in models(input_lines(args.input)):
        if "/" in name:
            sys.exit(f"ERROR: existing HMM NAME '{name}' contains '/'.")
        if name in seen:
            sys.exit(f"ERROR: existing HMM NAME '{name}' is given more than once.")
        seen.add(name)
        # mtime=0: identical models give identical files across runs
        with gzip.GzipFile(args.outdir / f"{name}.hmm.gz", "wb", mtime=0) as out:
            out.write("".join(model).encode())
    if not seen:
        sys.exit(f"ERROR: no HMMs found in {args.input}.")


if __name__ == "__main__":
    main()
