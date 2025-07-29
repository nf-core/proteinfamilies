#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import sys
import os
import gzip
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from functools import partial


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--fasta_folder",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Input folder with amino acid FASTA files.",
    )
    parser.add_argument(
        "-t",
        "--num_threads",
        type=int,
        default=1,
        help="Number of threads to use for parallel processing."
    )
    parser.add_argument(
        "-m",
        "--metadata",
        required=True,
        metavar="FILE",
        type=str,
        help="Output CSV file with family ids, sizes and representative sequences.",
    )
    parser.add_argument(
        "-o",
        "--out_fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Output FASTA file with representative sequences.",
    )
    return parser.parse_args(args)


def parse_first_fasta(filepath):
    open_func = gzip.open if filepath.endswith(".gz") else open
    header = None
    seq_lines = []
    flag = False

    with open_func(filepath, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">") and not flag:
                flag = True
                if header is None:
                    header = line[1:].split()[0]  # first header only
                else:
                    break  # stop reading after second header
            elif header:
                seq_lines.append(line)

    # get count of total proteins on second pass
    with open_func(filepath, "rt") as f:
        size = sum(1 for line in f if line.startswith(">"))

    if header and seq_lines:
        sequence = "".join(seq_lines).replace("-", "").replace(".", "").upper()
        return header, sequence, size
    return None, None, size


def process_file(filename, fasta_folder):
    filepath = os.path.join(fasta_folder, filename)
    family_name = os.path.splitext(os.path.splitext(filename)[0])[0]

    header, sequence, size = parse_first_fasta(filepath)
    if header and sequence:
        return (
            family_name,           # Sample Name
            family_name,           # Family ID
            size,                  # Family size
            len(sequence),         # Representative Length
            header,                # Representative ID
            sequence               # Sequence
        )
    return None


def extract_first_sequences(fasta_folder, num_threads, metadata_file, out_fasta):
    all_files = sorted(os.listdir(fasta_folder))

    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(partial(process_file, fasta_folder=fasta_folder), all_files))

    with open(out_fasta, "w") as fasta_out, open(metadata_file, "w", newline="") as csv_out:
        # Write custom metadata lines to the metadata file
        csv_out.write(
            '# id: "family_metadata"\n'
            '# section_name: "Family Metadata"\n'
            '# description: "Family metadata table containing family ids and sizes along with representative sequences, ids and lengths."\n'
            '# format: "csv"\n'
            '# plot_type: "table"\n'
        )
        csv_writer = csv.writer(csv_out, quoting=csv.QUOTE_NONNUMERIC)
        csv_writer.writerow([
            "Sample Name",
            "Family Id",
            "Size",
            "Representative Length",
            "Representative Id",
            "Sequence",
        ])

        for res in results:
            if res:
                sample, fam_id, size, length, rep_id, seq = res
                csv_writer.writerow([sample, fam_id, size, length, rep_id, seq])
                fasta_out.write(f">{rep_id}\n{seq}\n")


def main(args=None):
    args = parse_args(args)
    extract_first_sequences(args.fasta_folder, args.num_threads, args.metadata, args.out_fasta)


if __name__ == "__main__":
    sys.exit(main())
