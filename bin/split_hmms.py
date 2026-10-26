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
import re
import sys
import tarfile
from collections.abc import Iterator, Sequence
from pathlib import Path


def parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", required=True, type=Path, help="HMM .tar.gz archive or HMM library (plain or gzipped)")
    parser.add_argument("-o", "--outdir", required=True, type=Path, help="Output folder for the <NAME>.hmm.gz files")
    return parser.parse_args(args)


def lines(stream: io.BufferedReader, source: str) -> Iterator[str]:
    """Text lines of a stream, decompressed if it is gzipped (by magic bytes, not extension)."""
    data = gzip.GzipFile(fileobj=stream) if stream.peek(2)[:2] == b"\x1f\x8b" else stream
    try:
        for line in data:
            yield line.decode()
    except (UnicodeDecodeError, gzip.BadGzipFile, EOFError):
        sys.exit(f"ERROR: {source} is not a plain or gzipped text HMM file.")


def input_files(path: Path) -> Iterator[tuple[str, Iterator[str]]]:
    """(source, lines) of every HMM file in the input: each file member of a .tar.gz, or the library itself."""
    if path.name.endswith(".tar.gz"):
        with tarfile.open(path, "r|gz") as tar:
            for member in tar:
                # macOS tar adds AppleDouble `._<file>` metadata members, which are not HMMs
                if member.isfile() and not Path(member.name).name.startswith("._"):
                    source = f"{path.name}:{member.name}"
                    yield source, lines(tar.extractfile(member), source)
    else:
        with path.open("rb") as handle:
            yield path.name, lines(handle, path.name)


def models(text: Iterator[str], source: str) -> Iterator[tuple[str, list[str]]]:
    """(NAME, lines) of each model of one file, ending at its `//` line."""
    model: list[str] = []
    name = None
    for line in text:
        model.append(line)
        if line.startswith("NAME"):
            name = line.split()[1]
        elif line.rstrip() == "//":
            if name is None:
                sys.exit(f"ERROR: an existing HMM in {source} has no NAME line.")
            yield name, model
            model, name = [], None
    if any(line.strip() for line in model):
        sys.exit(f"ERROR: existing HMM {name or '(no NAME)'} in {source} is truncated or not an HMM: no closing '//' line.")


def main(args: Sequence[str] | None = None) -> None:
    args = parse_args(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    seen: dict[str, str] = {}  # NAME -> source
    for source, text in input_files(args.input):
        for name, model in models(text, source):
            # names become file names and shell arguments in later processes
            if not re.fullmatch(r"[A-Za-z0-9._-]+", name):
                sys.exit(f"ERROR: existing HMM NAME '{name}' in {source} may only contain letters, digits, '.', '_' and '-'.")
            if name in seen:
                sys.exit(f"ERROR: existing HMM NAME '{name}' is given more than once ({seen[name]}, {source}).")
            seen[name] = source
            # mtime=0: identical models give identical files across runs
            with gzip.GzipFile(args.outdir / f"{name}.hmm.gz", "wb", mtime=0) as out:
                out.write("".join(model).encode())
    if not seen:
        sys.exit(f"ERROR: no HMMs found in {args.input}.")


if __name__ == "__main__":
    main()
