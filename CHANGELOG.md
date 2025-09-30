# nf-core/proteinfamilies: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.4.0dev - [date]

### `Added`

- [#118](https://github.com/nf-core/proteinfamilies/pull/118)
  - Added preprint citation to the repo. (by @vagkaratzas)
  - Added separate metro map files for dark and light browser modes. (by @vagkaratzas)
  - Added new local module `EXTRACT_FAMILY_MEMBERS` that outputs a 2-col TSV file with final family identifiers and all their member sequence identifiers. (by @vagkaratzas)
- [#117](https://github.com/nf-core/proteinfamilies/pull/117)
  - Added `SEQKIT_SEQ` for optional sequence preprocessing in the quality check subworkflow. (by @vagkaratzas)
  - Added `SEQKIT_REPLACE` for optional sequence name parsing in the quality check subworkflow. (by @vagkaratzas)
  - Added `SEQKIT_RMDUP` for optional removal of duplicate names and sequences in the quality check subworkflow. (by @vagkaratzas)

### `Changed`

- [#118](https://github.com/nf-core/proteinfamilies/pull/118)
  - Swapped the local `CHECK_QUALITY` subworkflow with the new nf-core one `FAA_SEQFU_SEQKIT`. (by @vagkaratzas)
  - Based on protein family reproducibility benchmarks (i.e., computationally reproducing manually curated protein family resources), the `cluster_seq_identity` and `cluster_coverage` parameter default values have been updated to `0.3` and `0.5` (down from `0.5` and `0.9`) respectively. (by @vagkaratzas)
- [#117](https://github.com/nf-core/proteinfamilies/pull/117) - Swapped the local `SEQKIT_STATS` and the local `SEQKIT_STATS_TO_MQC` modules with the `SEQFU_STATS` one, which runs a bit faster and produces a MultiQC-ready output without the need for manual parsing. (by @vagkaratzas)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| seqfu   | -                | 1.20.3      |
| multiqc | 1.30             | 1.31        |

## v1.3.1 - [2025/09/22]

### `Fixed`

- [#112](https://github.com/nf-core/proteinfamilies/pull/112) - Fixed a bug in `EXTRACT_FAMILY_REPS`, where all sequences were pasted into the family representative one, and updated the relevant local nf-test. (by @vagkaratzas)

### `Changed`

- [#106](https://github.com/nf-core/proteinfamilies/pull/106) - Swapped the local `EXECUTE_CLUSTERING` subworkflow with the new nf-core `MMSEQS_FASTA_CLUSTER` one. (by @vagkaratzas)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| multiqc | 1.29             | 1.30        |

## v1.3.0 - [2025/08/06]

### `Changed`

- [#104](https://github.com/nf-core/proteinfamilies/pull/104) - Pulling `params` from local subworkflows into main workflow.
- [#103](https://github.com/nf-core/proteinfamilies/pull/103) - Parallelized execution for the `EXTRACT_FAMILY_REPS` local module and changed its input from `full_msa` to `fasta`.
- [#100](https://github.com/nf-core/proteinfamilies/pull/100) - `CAT_CAT` module replaced with `FIND_CONCATENATE` to avoid large scale `Argument list too long` errors.
- [#98](https://github.com/nf-core/proteinfamilies/pull/98) - nf-core tools template update to 3.3.2.

### `Added`

- [#105](https://github.com/nf-core/proteinfamilies/pull/105) - `CHECK_QUALITY` subworkflow added at the start of the pipeline.
  It utilizes the `seqkit/stats` nf-core module to generate a `MultiQC`-ready report with statistics for the input amino acid sequences.
  The metro-map has been updated to reflect this change.

## v1.2.0 - [2025/06/13]

### `Added`

- [#93](https://github.com/nf-core/proteinfamilies/pull/93)
  - Added nf-test and `meta.yml` file for local subworkflow `GENERATE_FAMILIES`.
  - Added nf-test and `meta.yml` file for local subworkflow `REMOVE_REDUNDANCY`.
  - Added nf-test and `meta.yml` file for local subworkflow `UPDATE_FAMILIES`.
- [#88](https://github.com/nf-core/proteinfamilies/pull/88)
  - Added nf-test and `meta.yml` file for local module `BRANCH_HITS_FASTA`.
  - Added nf-test and `meta.yml` file for local module `FILTER_NON_REDUNDANT_FAMS`.
  - Added nf-test and `meta.yml` file for local module `IDENTIFY_REDUNDANT_FAMS`.
  - Added nf-test and `meta.yml` file for local module `EXTRACT_FAMILY_REPS`.
  - Added the default pipeline end-to-end nf-test.

### `Changed`

- [#81](https://github.com/nf-core/proteinfamilies/pull/81) - nf-core tools template update to 3.3.1.

### `Fixed`

- [#80](https://github.com/nf-core/proteinfamilies/pull/80) - Fixed a bug where, due to a missing check for equal family sizes, non-redundant families were erroneously marked as redundant through transitive relationships and were removed

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| multiqc | 1.28             | 1.29        |

## v1.1.1 - [2025/05/17]

### `Changed`

- [#77](https://github.com/nf-core/proteinfamilies/pull/77) - Default branch changed from `master` to `main`.
- [#73](https://github.com/nf-core/proteinfamilies/pull/73) - Changed the fasta parsing library of the `CHUNK_CLUSTERS` local module, from `pyfastx` back to the latest version of `biopython`, and parallelized its writing mechanism, achieving decreased execution time.

### `Dependencies`

| Tool      | Previous version | New version |
| --------- | ---------------- | ----------- |
| biopython | 1.84             | 1.85        |
| pyfastx   | 2.2.0            |             |

### `Removed`

- [#73](https://github.com/nf-core/proteinfamilies/pull/73) - Deprecated `pyfastx` module version of `CHUNK_CLUSTERS`, since it was struggling performance-wise with larger datasets.

## v1.1.0 - [2025/05/06]

### `Added`

- [#69](https://github.com/nf-core/proteinfamilies/pull/69)
  - Added the `hhsuite/reformat` nf-core module to reformat `.sto` alignments to `.fas` when in-family sequence redundancy is not removed.
  - Added the option to save intermediate and final family fasta files throughout the workflow with various `save` parameters.
- [#58](https://github.com/nf-core/proteinfamilies/pull/58) - Added nf-test and `meta.yml` file for local module `REMOVE_REDUNDANCY_SEQS` (Hackathon 2025)
- [#56](https://github.com/nf-core/proteinfamilies/pull/56) - Added nf-test and `meta.yml` file for local module `FILTER_RECRUITED` (Hackathon 2025)
- [#55](https://github.com/nf-core/proteinfamilies/pull/55) - Added nf-test and `meta.yml` file for local module `CHUNK_CLUSTERS` (Hackathon 2025)
- [#54](https://github.com/nf-core/proteinfamilies/pull/54) - Added nf-test for local subworkflow `ALIGN_SEQUENCES` (Hackathon 2025)
- [#53](https://github.com/nf-core/proteinfamilies/pull/53) - Added nf-test for local subworkflow `EXECUTE_CLUSTERING` (Hackathon 2025)
- [#51](https://github.com/nf-core/proteinfamilies/pull/51) - Added nf-test and `meta.yml` file for local module `CALCULATE_CLUSTER_DISTRIBUTION` (Hackathon 2025)
- [#34](https://github.com/nf-core/proteinfamilies/pull/34) - Added the `EXTRACT_UNIQUE_CLUSTER_REPS` module, that calculates initial `MMseqs` clustering metadata, for each sample, to print with `MultiQC` (Id,Cluster Size,Number of Clusters)

### `Fixed`

- [#69](https://github.com/nf-core/proteinfamilies/pull/69) - Fixed a bug where redundant family alignments were not published properly, if intra-family redundancy removal mechanism was switched off [#68](https://github.com/nf-core/proteinfamilies/pull/68)
- [#65](https://github.com/nf-core/proteinfamilies/pull/65) - Fixed a bug in `CHUNK_CLUSTERS`, where pipeline would crash if the module filtered out all clusters, due to a high membership threshold [#64](https://github.com/nf-core/proteinfamilies/pull/64)
- [#35](https://github.com/nf-core/proteinfamilies/pull/35) - Fixed a bug in `remove_redundant_fams.py`, where comparison was between strings instead of integers to keep larger family
- [#33](https://github.com/nf-core/proteinfamilies/pull/33) - Fixed an always-true condition at the `filter_non_redundant_hmms.py` script, by adding missing parentheses
- [#29](https://github.com/nf-core/proteinfamilies/pull/29) - Fixed `hmmalign` empty input crash error, by preventing the `FILTER_RECRUITED` module from creating an empty output .fasta.gz file, when there are no remaining sequences after filtering the `hmmsearch` results [#28](https://github.com/nf-core/proteinfamilies/issues/28)

### `Changed`

- [#69](https://github.com/nf-core/proteinfamilies/pull/69)
  - Changed the publish directory architecture for HMMs, seed MSAs, full MSAs and family FASTA files, to make it more intuitive.
  - `REMOVE_REDUNDANT_FAMS` local module converted to `IDENTIFY_REDUNDANT_FAMS` to extract redundant family ids which will then be used downstream.
  - `FILTER_NON_REDUNDANT_HMMS` local module converted to `FILTER_NON_REDUNDANT_FAMS` and reused four times (HMM, seed MSA, full MSA, FASTA).
  - Changed the output format of the `EXTRACT_FAMILY_REPS` and `REMOVE_REDUNDANT_SEQS` local modules from `.fa` to `.faa`.
  - Metro map updated with new `hhsuite/reformat` module.
- [#57](https://github.com/nf-core/proteinfamilies/pull/57) - slight improvements of `nextflow_schema.json` (Hackathon 2025)
- [#57](https://github.com/nf-core/proteinfamilies/pull/57) - slight improtmenets of `assets/schema_input.json` (Hackathon 2025)
- [#34](https://github.com/nf-core/proteinfamilies/pull/34) - Swapped the `SeqIO` python library with `pyfastx` for the `CHUNK_CLUSTERS` module, quartering its duration
- [#32](https://github.com/nf-core/proteinfamilies/pull/32) - Updated `ClipKIT` 2.4.0 -> 2.4.1, that now also allows ends-only trimming, to completely replace the custom `CLIP_ENDS` module. Users can now also define its output format by setting the `--clipkit_out_format` parameter (default: `clipkit`)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| ClipKIT | 2.4.0            | 2.4.1       |
| pyfastx |                  | 2.2.0       |
| hhsuite |                  | 3.3.0       |
| multiqc | 1.27             | 1.28        |

### `Deprecated`

- [#32](https://github.com/nf-core/proteinfamilies/pull/32) - Deprecated `CLIP_ENDS` module and `--clipping_tool` parameter. The only option now is `ClipKIT`, covering both previous modes, via setting `--trim_ends_only`

## v1.0.0 - [2025/02/05]

Initial release of nf-core/proteinfamilies, created with the [nf-core](https://nf-co.re/) template

### `Added`

- Amino acid sequence clustering (mmseqs)
- Multiple sequence alignment (famsa, mafft, clipkit)
- Hidden Markov Model generation (hmmer)
- Between families redundancy removal (hmmer)
- In-family sequence redundancy removal (mmseqs)
- Family updating (hmmer, seqkit, mmseqs, famsa, mafft, clipkit)
- Family statistics presentation (multiqc)
