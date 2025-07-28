#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.

import sys
import os
import gzip
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ProcessPoolExecutor
from functools import partial


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--full_msa_folder",
        required=True,
        metavar="FOLDER",
        type=str,
        help="Input folder with fasta full alignments.",
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
        help="Name of the output csv file with family ids, sizes and representative sequences.",
    )
    parser.add_argument(
        "-o",
        "--out_fasta",
        required=True,
        metavar="FILE",
        type=str,
        help="Name of the output fasta file with family representative sequences.",
    )
    return parser.parse_args(args)


def process_file(filename, msa_folder):
    filepath = os.path.join(msa_folder, filename)
    format = "stockholm" if filepath.endswith(".sto.gz") else "fasta"

    open_func = gzip.open if filepath.endswith(".gz") else open
    with open_func(filepath, "rt") as fasta_file:
        parser = SeqIO.parse(fasta_file, format)
        try:
            # Only read first sequence
            first_record = next(parser)
            # Then count remaining records to get total family size without storing in memory
            family_size = 1 + sum(1 for _ in parser)
        except StopIteration:
            return None  # skip empty alignments

    if first_record:
        # Remove gaps from the sequence, and convert all to upper case
        cleaned_sequence = (
            str(first_record.seq).replace("-", "").replace(".", "").upper()
        )
        # Modify the ID to only include the part before the first space
        cleaned_id = first_record.id.split(" ")[0]
        # Create a new SeqRecord with the cleaned sequence and ID
        cleaned_record = SeqRecord(
            Seq(cleaned_sequence), id=cleaned_id, description=""
        )
        family_name = os.path.splitext(os.path.splitext(filename)[0])[0]
        return (family_name, family_size, len(cleaned_sequence), cleaned_id, cleaned_sequence, cleaned_record)

    return None


def extract_first_sequences(msa_folder, num_threads, metadata_file, out_fasta):
    all_files = sorted(os.listdir(msa_folder))

    # Process in parallel
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(partial(process_file, msa_folder=msa_folder), all_files))

    # Open the output FASTA file in write mode
    with open(out_fasta, "w") as fasta_out, open(
        metadata_file, "w", newline=""
    ) as csv_out:
        # Write custom metadata lines to the metadata file
        csv_out.write(
            '# id: "family_metadata"\n'
            '# section_name: "Family Metadata"\n'
            '# description: "Family metadata table containing family ids and sizes along with representative sequences, ids and lengths."\n'
            '# format: "csv"\n'
            '# plot_type: "table"\n'
        )
        csv_writer = csv.writer(csv_out, quoting=csv.QUOTE_NONNUMERIC)
        # Write the CSV header
        csv_writer.writerow(
            [
                "Sample Name",
                "Family Id",
                "Size",
                "Representative Length",
                "Representative Id",
                "Sequence",
            ]
        )

        for res in results:
            if res:
                family_name, size, rep_len, rep_id, sequence, record = res
                csv_writer.writerow([family_name, family_name, size, rep_len, rep_id, sequence])
                SeqIO.write(record, fasta_out, "fasta")


def main(args=None):
    args = parse_args(args)
    extract_first_sequences(args.full_msa_folder, args.num_threads, args.metadata, args.out_fasta)


if __name__ == "__main__":
    sys.exit(main())
