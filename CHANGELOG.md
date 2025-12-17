# nf-core/proteinfamilies: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.2.0dev - [date]

### `Added`

- [#142](https://github.com/nf-core/proteinfamilies/pull/142) - Added the `cmaple` module for optional phylogenetic tree inference for final family full MSAs. (by @vagkaratzas)
- [#140](https://github.com/nf-core/proteinfamilies/pull/140) - Using the new workflow output syntax to publish the downstream `nf-core/proteinannotator` samplesheet. (by @vagkaratzas)

### `Changed`

- [#142](https://github.com/nf-core/proteinfamilies/pull/142) - Updated metro-maps, citations, README.md and output.md to include `cmaple` phylogenetic trees. (by @vagkaratzas)
- [#141](https://github.com/nf-core/proteinfamilies/pull/141) - Aligning the default output directory of the new output syntax with the publishing directory of the pipeline. (by @vagkaratzas)
- [#140](https://github.com/nf-core/proteinfamilies/pull/140)
  - Updated metro-map to also depict the optional creation of downstream samplesheets for `nf-core/proteinfold` and `nf-core/proteinannotator`. Also swapped `seqkit/rmdup` and `seqkit/replace` modules to their proper execution sequence. (by @vagkaratzas)
  - Updated the `test_full` profile time and memory requirements to avoid AWS failures on release. (by @vagkaratzas)

### `Fixed`

- [#142](https://github.com/nf-core/proteinfamilies/pull/142) - Fixed a bug in `REMOVE_REDUNDANCY` subworkflow, where the combination of these skip flags `--skip_sequence_redundancy_removal true`, `--skip_additional_sequence_recruiting true` and `--skip_additional_sequence_recruiting false`, would execute `HHSUITE_REFORMAT_FILTERED` to reformat Stockholm alignments while they were already in fasta (or clipkit) format. (by @vagkaratzas)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| multiqc | 1.32             | 1.33        |
| cmaple  | -                | 1.1.0       |

## v2.1.0 - [2025/11/25]

### `Added`

- [#133](https://github.com/nf-core/proteinfamilies/pull/133) - Using the new workflow output syntax to publish the downstream `nf-core/proteinfold` samplesheet. (by @vagkaratzas)
- [#132](https://github.com/nf-core/proteinfamilies/pull/132) - Added optimized memory and time resources for `test` and `test_full` profiles. (by @vagkaratzas)

### `Changed`

- [#136](https://github.com/nf-core/proteinfamilies/pull/136) - Based on protein family reproducibility benchmarks, `cluster` is now the default MMseqs2 mode, due to its increased sensitivity compared to `linclust`. (by @vagkaratzas)
- [#135](https://github.com/nf-core/proteinfamilies/pull/135) - nf-core tools template update to 3.5.1. (by @vagkaratzas)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| mmseqs  | 17.b804f         | 18.8cc5c    |
| multiqc | 1.31             | 1.32        |

## v2.0.0 - [2025/10/14]

### `Added`

- [#124](https://github.com/nf-core/proteinfamilies/pull/124)
  - Added new subworkflow `MERGE_FAMILIES` that can optionally merge similar (but not redundant) generated protein families. (by @vagkaratzas)
  - Added new functionality to the local module `IDENTIFY_REDUNDANT_FAMS` which now also detects and outputs the identifiers of similar families that can optionally be merged downstream. These identifiers are written to _"/remove_redundancy/&lt;samplename&gt;/similar_fam_ids.txt"_, and the corresponding family pairwise similarity scores to _"/remove_redundancy/&lt;samplename&gt;/similarities.csv"_. (by @vagkaratzas)
  - Added new local module `POOL_SIMILAR_COMPONENTS` that generates family clusters, from a family-similarity edgelist. (by @vagkaratzas)
  - Added new local module `MERGE_SEEDS` that merges seed alignments of similar families, before restarting the family generation subworkflow. (by @vagkaratzas)
- [#118](https://github.com/nf-core/proteinfamilies/pull/118)
  - Added preprint citation to the repo. (by @vagkaratzas)
  - Added separate metro-map files for dark and light browser modes. (by @vagkaratzas)
  - Added new local module `EXTRACT_FAMILY_MEMBERS` which outputs a two-column TSV file containing the final family identifiers and their corresponding member sequence identifiers. The file is saved at _"/family_reps/&lt;samplename&gt;/&lt;samplename&gt;.tsv"_. (by @vagkaratzas)
- [#117](https://github.com/nf-core/proteinfamilies/pull/117)
  - Added `SEQKIT_SEQ` for optional sequence preprocessing in the quality check subworkflow. (by @vagkaratzas)
  - Added `SEQKIT_REPLACE` for optional sequence name parsing in the quality check subworkflow. (by @vagkaratzas)
  - Added `SEQKIT_RMDUP` for optional removal of duplicate names and sequences in the quality check subworkflow. (by @vagkaratzas)

### `Changed`

- [#128](https://github.com/nf-core/proteinfamilies/pull/128) - nf-core tools template update to 3.4.1.
- [#124](https://github.com/nf-core/proteinfamilies/pull/124)
  - Conditional workflow flags switched to their `skip` opposites; `--trim_msa` to `--skip_msa_trimming`, `--recruit_sequences_with_models` to `--skip_additional_sequence_recruiting`, `--remove_family_redundancy` to `--skip_family_redundancy_removal`, `--remove_sequence_redundancy` to `--skip_sequence_redundancy_removal`. (by @vagkaratzas)
- [#118](https://github.com/nf-core/proteinfamilies/pull/118)
  - Swapped the local `CHECK_QUALITY` subworkflow with the new nf-core one `FAA_SEQFU_SEQKIT`. (by @vagkaratzas)
  - Based on protein family reproducibility benchmarks (i.e., computationally reproducing manually curated protein family resources), the `cluster_seq_identity` and `cluster_coverage` parameter default values have been updated to `0.3` and `0.5` (down from `0.5` and `0.9`) respectively. (by @vagkaratzas)
- [#117](https://github.com/nf-core/proteinfamilies/pull/117) - Swapped the local `SEQKIT_STATS` and the local `SEQKIT_STATS_TO_MQC` modules with the `SEQFU_STATS` one, which runs a bit faster and produces a MultiQC-ready output without the need for manual parsing. (by @vagkaratzas)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| seqfu   | -                | 1.20.3      |
| multiqc | 1.30             | 1.31        |

### `Deprecated`

- [#124](https://github.com/nf-core/proteinfamilies/pull/124) - Deprecated `--trim_msa`, `--recruit_sequences_with_models`, `--remove_family_redundancy` and `--remove_sequence_redundancy`. (by @vagkaratzas)

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
  - Metro-map updated with new `hhsuite/reformat` module.
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
