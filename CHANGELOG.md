# nf-core/proteinfamilies: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - [2025/05/06]


### `Added`

- [#69](https://github.com/nf-core/proteinfamilies/pull/69) - Added the `hhsuite/reformat` nf-core module to reformat `.sto` alignments to `.fas` when in-family sequence redundancy is not removed.
  Also added the option to save intermediate and final family fasta files throughout the workflow with various `save` parameters.
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

- [#69](https://github.com/nf-core/proteinfamilies/pull/69) - Changed the publish directory architecture for HMMs, seed MSAs, full MSAs and family FASTA files, to make it more intuitive.
  `REMOVE_REDUNDANT_FAMS` local module converted to `IDENTIFY_REDUNDANT_FAMS` to extract redundant family ids which will then be used downstream.
  `FILTER_NON_REDUNDANT_HMMS` local module converted to `FILTER_NON_REDUNDANT_FAMS` and reused four times (HMM, seed MSA, full MSA, FASTA).
  Changed the output format of the `EXTRACT_FAMILY_REPS` and `REMOVE_REDUNDANT_SEQS` local modules from `.fa` to `.faa`.
  Metro map updated with new `hhsuite/reformat` module.
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
