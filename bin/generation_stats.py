#!/usr/bin/env python

## Originally written by Evangelos Karatzas and released under the MIT license.
## See git repository (https://github.com/nf-core/proteinfamilies) for full license text.
"""
Summarises the families the standard algorithm built (or updated) for one report sample into
MultiQC custom content that mirrors the sections of the MultiQC mgnifam module: general stats
(Families In, % Successful), Family outcomes, Full MSA size, Model length and, for updates,
Model length change. Counts describe the generation step, before redundancy removal and merging.
"""

import sys
import os
import gzip
import json
import argparse
from collections import Counter
from typing import Sequence

PARENT = {
    "parent_id": "family_generation",
    "parent_name": "Family generation",
    "parent_description": (
        "Families as the standard algorithm built them (one sample per input, `<id>`) or updated "
        "them (`<id>_updated`), before redundancy removal and merging. The iterative algorithm's "
        "families are reported per chunk by the mgnifam section."
    ),
}


def parse_args(args: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample", required=True, help="Report sample name.")
    parser.add_argument(
        "-n",
        "--families_in",
        required=True,
        type=int,
        help="Clusters (creation) or existing models (update) that entered the algorithm.",
    )
    parser.add_argument(
        "-f",
        "--fasta_folder",
        required=True,
        help="Folder with one member FASTA per successful family (.faa or .faa.gz).",
    )
    parser.add_argument(
        "-m",
        "--hmm_folder",
        required=True,
        help="Folder with one HMM per successful family (.hmm or .hmm.gz).",
    )
    parser.add_argument(
        "-e",
        "--existing_hmm_folder",
        help="Update only: folder with the existing HMMs, matched to the new ones by family name.",
    )
    parser.add_argument(
        "-d",
        "--discard",
        action="append",
        default=[],
        metavar="REASON=COUNT",
        help="Families discarded for REASON; repeatable.",
    )
    parser.add_argument(
        "-r",
        "--remainder_reason",
        help="Reason for every input family that was neither successful nor given with --discard.",
    )
    parser.add_argument("-o", "--prefix", required=True, help="Output file prefix.")
    return parser.parse_args(args)


def family_files(folder: str | None) -> dict[str, str]:
    """
    Family name -> file path for the files of a folder; the name is the file name without
    '.gz' and its last extension. A missing folder holds no families.
    """
    if not folder or not os.path.isdir(folder):
        return {}
    return {
        os.path.splitext(name.removesuffix(".gz"))[0]: os.path.join(folder, name)
        for name in sorted(os.listdir(folder))
    }


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def count_sequences(path: str) -> int:
    with open_text(path) as f:
        return sum(1 for line in f if line.startswith(">"))


def model_length(path: str) -> int:
    """Match states of the (first) model in an HMMER3 file, from its LENG line."""
    with open_text(path) as f:
        for line in f:
            if line.startswith("LENG"):
                return int(line.split()[1])
    raise ValueError(f"No LENG line in {path}")


def capitalise(reason: str) -> str:
    return reason[:1].upper() + reason[1:]


def histogram(values: list[int]) -> dict[int, int]:
    return dict(sorted(Counter(values).items()))


def linegraph(section: str, name: str, description: str, helptext: str, xlab: str, sample: str, data: dict) -> dict:
    return {
        **PARENT,
        "id": f"family_generation_{section}",
        "section_name": name,
        "description": description,
        "helptext": helptext,
        "plot_type": "linegraph",
        "pconfig": {
            "id": f"family-generation-{section.replace('_', '-')}-plot",
            "title": f"Family generation: {name}",
            "xlab": xlab,
            "ylab": "Families",
            "y_decimals": False,
        },
        "data": {sample: data},
    }


def write(path: str, content: dict) -> None:
    with open(path, "w") as f:
        json.dump(content, f, indent=2)
        f.write("\n")


def main(args: Sequence[str] | None = None) -> None:
    args = parse_args(args)

    fastas = family_files(args.fasta_folder)
    hmms = family_files(args.hmm_folder)
    existing = family_files(args.existing_hmm_folder)
    successful = len(fastas)

    outcomes = Counter({"Successful": successful})
    for item in args.discard:
        reason, _, count = item.rpartition("=")
        outcomes[capitalise(reason)] += int(count)
    # Keep the bar summing to the input count, so a discard without a reason stays visible
    remainder = args.families_in - sum(outcomes.values())
    if remainder > 0:
        outcomes[capitalise(args.remainder_reason or "unaccounted")] += remainder

    general = {"families_input": args.families_in}
    # No input families: leave the percentage out rather than divide by zero
    if args.families_in:
        general["successful_pct"] = 100 * successful / args.families_in
    write(f"{args.prefix}_general_mqc.json", {
        "id": "family_generation_general",
        "plot_type": "generalstats",
        "headers": {
            "families_input": {
                "title": "Families In",
                "description": "Clusters (creation) or existing models (update) entering the standard algorithm",
                "scale": "Blues",
                "format": "{:,.0f}",
            },
            "successful_pct": {
                "title": "% Successful",
                "description": "Percentage of input families that produced a family",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
            },
        },
        "data": {args.sample: general},
    })

    write(f"{args.prefix}_outcomes_mqc.json", {
        **PARENT,
        "id": "family_generation_outcomes",
        "section_name": "Family outcomes",
        "description": "Families that were generated or updated successfully, and the reasons the rest were discarded.",
        "helptext": (
            "Each input cluster (creation) or existing model (update) ends in exactly one outcome: it "
            "becomes a successful family, or it is discarded with the reason shown here. A discarded "
            "existing family passes through unchanged."
        ),
        "plot_type": "bargraph",
        "pconfig": {"id": "family-generation-outcomes-plot", "title": "Family generation: Family outcomes", "ylab": "Families"},
        "data": {args.sample: dict(outcomes)},
    })

    # Histograms only for samples with families, as the mgnifam module plots no empty chunk
    if not fastas:
        return
    write(f"{args.prefix}_full_msa_size_mqc.json", linegraph(
        "full_msa_size", "Full MSA size", "Number of families per full alignment size.",
        "Each successful family's full alignment holds every sequence its final model recruited. "
        "Small alignments mark narrow families; very large ones can mark a family that absorbed "
        "unrelated sequences.",
        "Sequences in the full alignment", args.sample,
        histogram([count_sequences(path) for path in fastas.values()]),
    ))
    lengths = {family: model_length(hmms[family]) for family in fastas if family in hmms}
    write(f"{args.prefix}_model_length_mqc.json", linegraph(
        "model_length", "Model length", "Number of families per model length.",
        "The length of each successful family's final HMM, in match states.",
        "Match states in the family model", args.sample,
        histogram(list(lengths.values())),
    ))
    if args.existing_hmm_folder is not None:
        write(f"{args.prefix}_model_length_change_mqc.json", linegraph(
            "model_length_change", "Model length change",
            "Number of families per change in model length during the update.",
            "Final model length minus the existing model length, for every updated family. With "
            "`--skip_update_refinement` models are not rebuilt, so every family sits at 0. Otherwise a "
            "model can shrink as columns are trimmed, or grow when new members support extra columns.",
            "Match states gained or lost by the update", args.sample,
            histogram([length - model_length(existing[family]) for family, length in lengths.items() if family in existing]),
        ))


if __name__ == "__main__":
    sys.exit(main())
