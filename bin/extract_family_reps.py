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


def extract_first_sequences(msa_folder, metadata_file, out_fasta):
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

        # Iterate over all files in the MSA folder
        for filename in os.listdir(msa_folder):
            filepath = os.path.join(msa_folder, filename)
            # Parse the MSA fasta file and extract the first sequence
            with (
                gzip.open(filepath, "rt")
                if filepath.endswith(".gz")
                else open(filepath, "r")
            ) as fasta_file:
                format = "stockholm" if filepath.endswith(".sto.gz") else "fasta"
                records = list(SeqIO.parse(fasta_file, format))
                family_size = len(records)
                if records:
                    first_record = records[0]
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
                    # Write the cleaned sequence to the FASTA file
                    SeqIO.write(cleaned_record, fasta_out, "fasta")
                    # Write the mapping to the CSV file
                    family_name = os.path.splitext(os.path.splitext(filename)[0])[0]
                    csv_writer.writerow(
                        [
                            family_name,
                            family_name,
                            family_size,
                            len(cleaned_sequence),
                            cleaned_id,
                            cleaned_sequence,
                        ]
                    )


def main(args=None):
    args = parse_args(args)
    extract_first_sequences(args.full_msa_folder, args.metadata, args.out_fasta)


if __name__ == "__main__":
    sys.exit(main())
