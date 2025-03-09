# nf-core/proteinfamilies: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0 - [yyyy/mm/dd]

### `Added`

- [#34](https://github.com/nf-core/proteinfamilies/pull/34) - Added the `EXTRACT_UNIQUE_CLUSTER_REPS` module, that calculates initial `MMseqs` clustering metadata, for each sample, to print with `MultiQC` (Id,Cluster Size,Number of Clusters)

### `Fixed`

- [#33](https://github.com/nf-core/proteinfamilies/pull/33) - Fixed an always-true condition at the `filter_non_redundant_hmms.py` script, by adding missing parentheses
- [#29](https://github.com/nf-core/proteinfamilies/pull/29) - Fixed `hmmalign` empty input crash error, by preventing the `FILTER_RECRUITED` module from creating an empty output .fasta.gz file, when there are no remaining sequences after filtering the `hmmsearch` results [#28](https://github.com/nf-core/proteinfamilies/issues/28)

### `Changed`

- [#34](https://github.com/nf-core/proteinfamilies/pull/34) - Swapped the `SeqIO` python library with `pyfastx` for the `CHUNK_CLUSTERS` module, quartering its duration
- [#32](https://github.com/nf-core/proteinfamilies/pull/32) - Updated `ClipKIT` 2.4.0 -> 2.4.1, that now also allows ends-only trimming, to completely replace the custom `CLIP_ENDS` module. Users can now also define its output format by setting the `--clipkit_out_format` parameter (default: `clipkit`)

### `Dependencies`

| Tool    | Previous version | New version |
| ------- | ---------------- | ----------- |
| ClipKIT | 2.4.0            | 2.4.1       |
| pyfastx |                  | 2.2.0       |

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
