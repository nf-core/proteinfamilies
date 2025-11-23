# nf-core/proteinfamilies: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

Quality check:

- [SeqFu](#seqfu) for input amino acid sequences quality check (QC)
- [SeqKit](#seqkit) for preprocessing input amino acid sequences (i.e., gap removal, convert to upper case, validate, filter by length, replace special characters such as `/`, and remove duplicate sequences)

Initial clustering:

- [MMseqs2](#mmseqs2) initial clustering of input amino acid sequences and filtering with membership threshold

Multiple sequence alignment:

- [FAMSA](#famsa) aligner option. Best speed and sensitivity option to build seed multiple sequence alignment for the families
- [mafft](#mafft) aligner option. Fast but not as sensitive as FAMSA to build seed multiple sequence alignment for the families
- [ClipKIT](#clipkit) to optionally clip gapped portions of the multiple sequence alignment (MSA)

Generating family models:

- [hmmer](#hmmer) to build the family HMM (hmmbuild) and to optionally 'fish' additional sequences from the input fasta file (hmmsearch), with given thresholds, into the family and also build the family full MSA (hmmalign)

Removing redundancy:

- [hmmer](#hmmer-for-redundancy-removal) to match family representative sequences against other family models in order to keep non redundant and/or merge similar ones
- [MMseqs2](#mmseqs2-for-redundancy-removal) to strictly cluster the sequences within each of the remaining families, in order to still capture the evolutionary diversity within a family, but without keeping all the almost identical sequences
- [FAMSA](#famsa-for-redundancy-removal) aligner option. Re-align full MSA with final set of sequences
- [mafft](#mafft-for-redundancy-removal) aligner option. Re-align full MSA with final set of sequences
- [HH-suite3](#hh-suite3) to reformat raw full `.sto` MSAs to `.fas`, in case the user did not select to remove in-family sequence redundancy, that automatically re-aligns sequences with either FAMSA or mafft

Updating families:

- [untar](#untar) to decompress tarballs of existing hmms and msas
- [hmmer](#hmmer-for-updating-families) to match input sequences to existing families with hmmsearch as well as for rebuilding models with newly recruited sequences with hmmbuild
- [SeqKit](#seqkit-for-updating-families) to extract fasta formatted family sequences from their MSA files
- [MMseqs2](#mmseqs2-for-updating-families) to strictly cluster the sequences within each of the families to update
- [FAMSA](#famsa-for-updating-families) aligner option. Re-align full MSA with final set of sequences
- [mafft](#mafft-for-updating-families) aligner option. Re-align full MSA with final set of sequences
- [ClipKIT](#clipkit-for-updating-families) to optionally clip gapped portions of the multiple sequence alignment (MSA)

Reporting:

- [Extract family representatives](#extract-family-representatives) to produce the final metadata file along with a fasta of all family representative sequences (can be used downstream for structural prediction).
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

Downstream pipelines:

- [nf-core/proteinfold](#nf-coreproteinfold) downstream samplesheet generation from final family representative sequences

### SeqFu

<details markdown="1">
<summary>Output files</summary>

- `qc/`
  - `<samplename>/`
    - `<samplename>_before.tsv`: Statistics for the input amino acid sequences before preprocessing
    - `<samplename>_before_mqc.txt`: Statistics for the input amino acid sequences in MultiQC-ready format before preprocessing
    - `<samplename>_after.tsv`: (optional) Statistics for the input amino acid sequences after preprocessing
    - `<samplename>_after_mqc.txt`: (optional) Statistics for the input amino acid sequences in MultiQC-ready format after preprocessing
    - `<samplename>.log`: (optional) Output file with count of duplicate sequences that were found and removed

</details>

The `seqfu` module is used for statistics generation of input amino acid sequences, both before and after preprocessing.

[SeqFu](https://github.com/telatin/seqfu2) is a cross-platform compiled suite of tools to manipulate and inspect FASTA and FASTQ files.

### SeqKit

<details markdown="1">
<summary>Output files</summary>

- `qc/`
  - `<samplename>/`
    - `<samplename>.<suffix>`: Updated preprocessed input fasta file

</details>

The `seqkit` module is used for initial preprocessing of the input amino acid sequences
as well as during the `update_families` mode to extract sequences from family MSA, into intermediate fasta files (`seqkit` output folder).
The intermediate `update_families/fasta` folder contains the aggregation of existing family sequences along with their newly matching ones,
that will together produce the updated family MSA.

[SeqKit](https://github.com/shenwei356/seqkit) is a cross-platform and ultrafast toolkit for FASTA/Q file manipulation.

### MMseqs2

<details markdown="1">
<summary>Output files</summary>

- `mmseqs/`
  - `initial_clustering/`
    - `mmseqs_createtsv/`
      - `<samplename>.tsv`: tab-separated table containing 2 columns; the first one with the cluster representative sequences, and the second with the cluster members
    - `mmseqs_createdb/`
      - `<samplename>/`
        - `*`: (optional) mmseqs format db of fasta sequences. Can be turned on with --save_mmseqs_db
    - `mmseqs_linclust/`
      - `<samplename>/`
        - `*`: (optional) mmseqs format clustered db. Can be turned on with --save_mmseqs_clustering
    - `mmseqs_cluster/`
      - `<samplename>/`
        - `*`: (optional) mmseqs format clustered db. Can be turned on with --save_mmseqs_clustering
    - `<samplename>_clustering_distribution_mqc.csv`: CSV file with initial clustering metadata, from each sample, to print with MultiQC (column headers: Id,Cluster Size,Number of Clusters)
- `fasta/`
  - `mmseqs_initial_clustering_filtered/`
    - `<samplename>/`
      - `chunked_fasta/`
        - `*.faa`: (optional) fasta files with amino acid sequences of each cluster above the membership threshold

</details>

The `mmseqs_createtsv/<samplename>.tsv` contains the mmseqs clustering of sequences, which will then be filtered by size and split into chunks for further parallel processing.
The optionally saved `chunked_fasta` folder contains these fasta files of sequences for each cluster.
These per cluster fasta files act as input to produce downstream families in the next steps of the pipeline.
The original mmseqs db and the clustered mmseqs db can be optional saved to the output folder, but they won't be further utilised in this pipeline.

[MMseqs2](https://github.com/soedinglab/MMseqs2) clusters amino acid fasta files via either the 'cluster' or the 'linclust' algorithms.

### FAMSA aligner

<details markdown="1">
<summary>Output files</summary>

- `seed_msa/`
  - `raw/`
    - `famsa_align/`
      - `<samplename>/`
        - `<samplename>_*.aln`: fasta files with aligned amino acid sequences
  - `filtered/`
    - `<samplename>/`
      - `<samplename>_*.*`: filtered seed alignments after family redundancy removal
- `remove_redundancy/`
  - `merge_families/`
    - `seed_msa/`
      - `raw/`
        - `famsa_align/`
          - `<samplename>/`
            - `<samplename>_*.aln`: fasta files with aligned amino acid sequences from merged families

</details>

This folder contains the generated seed MSA family files, if `famsa` was chosen as the `--alignment_tool`.
These MSA files only contain the original sequences of each cluster as calculated by mmseqs.

[FAMSA](https://github.com/refresh-bio/FAMSA) is a progressive algorithm for large-scale multiple sequence alignments.

### mafft aligner

<details markdown="1">
<summary>Output files</summary>

- `seed_msa/`
  - `raw/`
    - `mafft_align/`
      - `<samplename>/`
        - `<samplename>_*.fas`: fasta files with aligned amino acid sequences
  - `filtered/`
    - `<samplename>/`
      - `<samplename>_*.*`: filtered seed alignments after family redundancy removal
- `remove_redundancy/`
  - `merge_families/`
    - `seed_msa/`
      - `raw/`
        - `mafft_align/`
          - `<samplename>/`
            - `<samplename>_*.*`: fasta files with aligned amino acid sequences from merged families

</details>

This folder contains the generated seed MSA family files, if `mafft` was chosen as the `--alignment_tool`.
These MSA files only contain the original sequences of each cluster as calculated by mmseqs.

[mafft](https://github.com/GSLBiotech/mafft) is a fast but not very sensitive multiple sequence alignment tool.

### ClipKIT

<details markdown="1">
<summary>Output files</summary>

- `seed_msa/`
  - `raw/`
    - `clipkit/`
      - `<samplename>/`
        - `<samplename>_*.clipkit`: gap-clipped fasta files of aligned amino acid sequences
  - `filtered/`
    - `<samplename>/`
      - `<samplename>_*.*`: filtered seed alignments after family redundancy removal
- `remove_redundancy/`
  - `merge_families/`
    - `seed_msa/`
      - `raw/`
        - `clipkit/`
          - `<samplename>/`
            - `<samplename>_*.clipkit`: gap-clipped fasta files of aligned amino acid sequences from merged families

</details>

If the `--skip_msa_trimming` parameter was set to `false`, then `clipkit` runs, and according to the `--gap_threshold` parameter,
gaps (above that threshold, across all aligned sequences) are either removed only at the ends of the MSA if `trim_ends_only` is set to `true`, or throughout the alignment otherwise.
Results are stored in the `seed_msa/raw` folder.

[ClipKIT](https://github.com/JLSteenwyk/ClipKIT) is a fast and flexible alignment trimming tool that keeps phylogenetically informative sites and removes others.

### hmmer

<details markdown="1">
<summary>Output files</summary>

- `hmmer/`
  - `hmmsearch/`
    - `<samplename>/`
      - `<samplename>_*.domtbl.gz`: (optional) hmmsearch results along parameters info. Can be turned on with `--save_hmmsearch_results`
      - `<samplename>_*.txt.gz`: (optional) hmmsearch execution log. Can be turned on with `--save_hmmsearch_results`
- `hmm/`
  - `filtered/`
    - `<samplename>/`
      - `<samplename>_*.hmm.gz`: filtered non-redundant compressed hmm model for the family
  - `raw/`
    - `<samplename>/`
      - `<samplename>_*.hmm.gz`: compressed hmm model for the family
      - `<samplename>_*.hmmbuild.txt`: (optional) hmmbuild execution log
- `full_msa/`
  - `raw/`
    - `hmmer_hmmalign/`
      - `<samplename>/`
        - `<samplename>_*.sto.gz`: compressed family full MSA produced by hmmalign (before checking for redundancy)
  - `filtered/`
    - `hmmsearch/`
      - `<samplename>/`
        - `<samplename>_*.*`: filtered full alignments after family redundancy removal
- `fasta/`
  - `hmmsearch_filtered_recruited/`
    - `<samplename>/`
      - `<samplename>_*.fasta.gz`: (optional) filtered fasta sequences after hmmsearch and applied thresholds
  - `non_redundant_family_filtered/`
    - `<samplename>/`
      - `<samplename>_*.fasta.gz`: (optional) filtered full alignment sequences after family redundancy removal in fasta format
- `remove_redundancy/`
  - `merge_families/`
    - `hmmer/`
      - `hmmsearch/`
        - `<samplename>/`
          - `<samplename>_*.domtbl.gz`: (optional) hmmsearch results along parameters info. Can be turned on with `--save_hmmsearch_results`
          - `<samplename>_*.txt.gz`: (optional) hmmsearch execution log. Can be turned on with `--save_hmmsearch_results`
    - `hmm/`
      - `raw/`
        - `<samplename>/`
          - `<samplename>_*.hmm.gz`: compressed hmm model for the merged family
          - `<samplename>_*.hmmbuild.txt`: (optional) hmmbuild execution log
    - `full_msa/`
      - `raw/`
        - `hmmer_hmmalign/`
          - `<samplename>/`
            - `<samplename>_*.sto.gz`: compressed merged family full MSA produced by hmmalign (after checking for redundancy)
    - `fasta/`
      - `hmmsearch_filtered_recruited/`
        - `<samplename>/`
          - `<samplename>_*.fasta.gz`: (optional) filtered fasta sequences of merged families after hmmsearch and applied thresholds

</details>

The `hmm/raw` folder contains all originally created family HMMs. These models will be used downstream to recruit additional sequences in families, to compute
full MSAs if `--skip_additional_sequence_recruiting` is set to `false`, and/or to remove among-family redundancies if `--skip_family_redundancy_removal` is set to `false`.
When `--skip_family_redundancy_removal` is set to `false`, the `hmm/filtered` folder will also be produced with the filtered subset of the original raw HMMs.
The HMMs (raw or filtered) can also be used in the `update_families` execution mode of the pipeline,
along with the families' respective full MSAs, to recruit sequences from a new input fasta file into the families, updating both family HMM and full MSA files.

[hmmer](https://github.com/EddyRivasLab/hmmer) is a fast and flexible alignment trimming tool that keeps phylogenetically informative sites and removes others.

### hmmer for redundancy removal

<details markdown="1">
<summary>Output files</summary>

- `remove_redundancy/`
  - `<samplename>/`
    - `redundant_fam_ids.txt`: redundant family identifiers that are being dropped
    - `similar_fam_ids.txt`: similar family identifiers that are being dropped (their merged versions will be used instead)
    - `similarities.csv`: CSV file containing pairwise family similarities above user-defined threshold
  - `hmmer/`
    - `concatenated/`
      - `<samplename>.hmm.gz`: (optional) concatenated compressed hmm model for all families in a given sample (pre redundancy removal)
    - `hmmsearch/`
      - `<samplename>/`
        - `<samplename>_*.domtbl.gz`: (optional) hmmsearch results of family reps against families' HMMs
  - `family_reps/`
    - `<samplename>/`
      - `<samplename>_meta_mqc.csv`: (optional) CSV file with metadata (column headers: Sample Name,Family Id,Size,Representative Length,Representative Id,Sequence)
      - `<samplename>_reps.faa`: (optional) fasta file of all family representative sequences (one sequence per family)
  - `merge_families/`
    - `<samplename>/`
      - `pooled_components.txt`: comma separated clusters of similar family ids
      - `<merged_id>.fas`: (optional) merged seed alignment of each pooled component
  - `skipped_ids/`
    - `<samplename>.txt`: (optional) concatenated redundant and similar (single) family ids that are filtered out

</details>

If one of `--skip_family_redundancy_removal` or `--skip_family_merging` is set to `false`, the `hmmer/hmmsearch` module is used
to identify family representative sequences that are identical or similar (respectively) to other family HMMs.
In case of redundancy, the smaller sized families are flagged for removal.
If `--skip_family_merging` is set to `false`, and if `hmmsearch_family_similarity_length_threshold` is correctly set
lower than `hmmsearch_family_redundancy_length_threshold` (or `skip_family_redundancy_removal` is set to `true`), then similar family seed alignments can be merged
and go through the `generate_families` subworkflow once more.
Most `remove_redundancy` outputs are optional folders that contain intermediate pipeline results, and therefore are not saved in the output results by default.

[hmmer](https://github.com/EddyRivasLab/hmmer) is a fast and flexible alignment trimming tool that keeps phylogenetically informative sites and removes others.

### MMseqs2 for redundancy removal

<details markdown="1">
<summary>Output files</summary>

- `mmseqs/`
  - `redundancy_clustering/`
    - `mmseqs_createtsv/`
      - `<samplename>/`
        - `<samplename>_*.tsv`: tab-separated table containing 2 columns; the first one with the cluster representative sequences, and the second with the cluster members
    - `mmseqs_createdb/`
      - `<samplename>/`
        - `*`: (optional) mmseqs format db of fasta sequences
    - `mmseqs_linclust/`
      - `<samplename>/`
        - `*`: (optional) mmseqs format clustered db
    - `mmseqs_cluster/`
      - `<samplename>/`
        - `*`: (optional) mmseqs format clustered db
- `fasta/`
  - `non_redundant_sequences_filtered/`
    - `<samplename>/`
      - `<samplename>_reps.faa`: (optional) fasta file of all family representative sequences (one sequence per family)

</details>

If `--skip_sequence_redundancy_removal` is set to `false`, the mmseqs clustering subworkflow will be executed
to very strictly cluster (`--cluster_seq_identity_for_redundancy` = 0.97, `cluster_coverage_for_redundancy` = 0.97,
`cluster_cov_mode_for_redundancy` = 0 -meaning both strands) in-family sequences, keeping only cluster representatives
before recalculating the family MSAs.

[MMseqs2](https://github.com/soedinglab/MMseqs2) clusters amino acid fasta files via either the 'cluster' or the 'linclust' algorithms.

### FAMSA for redundancy removal

<details markdown="1">
<summary>Output files</summary>

- `full_msa/`
  - `filtered/`
    - `famsa_align/`
      - `<samplename>/`
        - `<samplename>_*.aln`: family full MSA (after checking for sequence redundancy)

</details>

If `--skip_sequence_redundancy_removal` is set to `false`, then the full MSAs will be recalculated after in-family sequence redundancy is removed.
If the `--alignment_tool` is `famsa`, then this `famsa_align` folder will be created, containing the final full MSA files.

[FAMSA](https://github.com/refresh-bio/FAMSA) is a progressive algorithm for large-scale multiple sequence alignments.

### mafft for redundancy removal

<details markdown="1">
<summary>Output files</summary>

- `full_msa/`
  - `filtered/`
    - `mafft_align/`
      - `<samplename>/`
        - `<samplename>_*.fas`: family full MSA (after checking for sequence redundancy)

</details>

If `--skip_sequence_redundancy_removal` is set to `false`, then the full MSAs will be recalculated after in-family sequence redundancy is removed.
If the `--alignment_tool` is `mafft`, then this `mafft_align` folder will be created, containing the final full MSA files.

[mafft](https://github.com/GSLBiotech/mafft) is a fast but not very sensitive multiple sequence alignment tool.

### HH-suite3

<details markdown="1">
<summary>Output files</summary>

- `full_msa/`
  - `filtered/`
    - `hhsuite_reformat/`
      - `<samplename>/`
        - `<samplename>_*.fas.gz`: reformatted filtered full MSA files
  - `raw/`
    - `hhsuite_reformat/`
      - `<samplename>/`
        - `<samplename>_*.fas.gz`: reformatted raw full MSA files

</details>

If `--skip_sequence_redundancy_removal` is set to `true`, then either the raw (if `--skip_family_redundancy_removal` is set to `true`) or the filtered (if `--skip_family_redundancy_removal` is set to `false`) full `.sto` MSAs will be reformatted to `.fas`.

[HH-suite3](https://github.com/soedinglab/hh-suite) is an open-source software package for sensitive protein sequence searching based on the pairwise alignment of hidden Markov models (HMMs).

### untar

<details markdown="1">
<summary>Output files</summary>

- `untar/`
  - `hmm/`
    - `<samplename>/`
      - `<family_name>.{hmm.gz,hmm}`: (optional) decompressed input hmm tarball
  - `msa/`
    - `<samplename>/`
      - `<family_id>.{aln,fas}`: (optional) decompressed input msa tarball

</details>

### hmmer for updating families

<details markdown="1">
<summary>Output files</summary>

- `update_families/`
  - `hmmer/`
    - `concatenated/`
      - `<samplename>.hmm.gz`: (optional) concatenated compressed HMM models for all families in a given sample, to be used as input for hmmsearch, to determine which families will be updated with new sequences
    - `hmmsearch/`
      - `<samplename>/`
        - `<samplename>.domtbl.gz`: (optional) hmmsearch results of input fasta file against existing families' HMMs
  - `hmm/`
    - `<samplename>/`
      - `<family_id>.hmm.gz`: (optional) compressed family HMM after the update
      - `<family_id>.hmmbuild.txt`: (optional) hmmbuild execution log
  - `branch_fasta/`
    - `hits/`
      - `<family_id>.fasta`: (optional) subset of the input FASTA with hit sequences for each existing family
    - `<samplename>.fasta.gz`: (optional) FASTA file that contains all remaining non-hit input sequences, which will be passed to normal execution mode to create new families
  - `family_reps/`
    - `<samplename>/`
      - `<samplename>_meta_mqc.csv`: CSV file with metadata (column headers: Sample Name,Family Id,Size,Representative Length,Representative Id,Sequence)
      - `<samplename>_reps.faa`: fasta file of all family representative sequences (one sequence per family)
      - `<samplename>.tsv`: 2-column TSV file with family ids and all sequence member ids

</details>

The `update_families` execution mode is run if paths to `existing_hmms_to_update` and `existing_msas_to_update` are provided in the input samplesheet.csv.
The `hmmer/hmmsearch` module is used to match new incoming sequences in the existing family models.
In case of hits, the new sequences are reclustered along their matching family existing ones, and new models are build with `hmmer/hmmbuild`
in the `update_families/hmmer/hmmbuild` folder, from the respective new MSAs.

[hmmer](https://github.com/EddyRivasLab/hmmer) is a fast and flexible alignment trimming tool that keeps phylogenetically informative sites and removes others.

### SeqKit for updating families

<details markdown="1">
<summary>Output files</summary>

- `seqkit/`
  - `<samplename>_<family_id>.fastq`: (optional) fasta formatted family sequences from full MSA with gaps removed
- `update_families/`
  - `fasta/`
    - `<samplename>_<family_id>.fastq`: (optional) concatenated family fasta with newly recruited sequences

</details>

The `seqkit` module is mainly used during the `update_families` mode
to extract sequences from family MSA, into intermediate fasta files (`seqkit` output folder).
The intermediate `update_families/fasta` folder contains the aggregation of existing family sequences along with their newly matching ones,
that will together produce the updated family MSA.

[SeqKit](https://github.com/shenwei356/seqkit) is a cross-platform and ultrafast toolkit for FASTA/Q file manipulation.

### MMseqs2 for updating families

<details markdown="1">
<summary>Output files</summary>

- `mmseqs/`
  - `update_families/`
    - `mmseqs_createtsv/`
      - `<samplename>/`
        - `<family_id>.tsv`: tab-separated table containing 2 columns; the first one with the cluster representative sequences, and the second with the cluster members
    - `mmseqs_createdb/`
      - `<family_id>/`
        - `*`: (optional) mmseqs format db of fasta sequences
    - `mmseqs_linclust/`
      - `<family_id>/`
        - `*`: (optional) mmseqs format clustered db
    - `mmseqs_cluster/`
      - `<family_id>/`
        - `*`: (optional) mmseqs format clustered db
    - `non_redundant_sequences/`
      - `<samplename>/`
        - `<samplename>_reps.faa`: (optional) fasta file of all family representative sequences (one sequence per family)

</details>

Similarly to the in-family sequence redundancy removal mechanism, the mmseqs suite is used to strictly cluster
existing family sequences along newly recruited ones, keeping a non redundant set.
The new family representative sequences can now be found in the intermediate `mmseqs/non_redundant_sequences` folder.

[MMseqs2](https://github.com/soedinglab/MMseqs2) clusters amino acid fasta files via either the 'cluster' or the 'linclust' algorithms.

### FAMSA for updating families

<details markdown="1">
<summary>Output files</summary>

- `update_families/`
  - `full_msa/`
    - `famsa_align/`
    - `<samplename>/`
      - `<family_id>.aln`: family full MSA (after updating with new sequences)

</details>

In the `update_families` mode, if new sequences are added in an existing family,
and after (optionally) removing in-family sequence redundacny, if `--skip_sequence_redundancy_removal` is set to `false`,
then the family MSA is recalculated.
If the `--alignment_tool` is `famsa`, then this `famsa_align` folder will be created, containing the updated family MSA files.

[FAMSA](https://github.com/refresh-bio/FAMSA) is a progressive algorithm for large-scale multiple sequence alignments.

### mafft for updating families

<details markdown="1">
<summary>Output files</summary>

- `update_families/`
  - `full_msa/`
    - `mafft_align/`
    - `<samplename>/`
      - `<family_id>.fas`: family full MSA (after updating with new sequences)

</details>

In the `update_families` mode, if new sequences are added in an existing family,
and after (optionally) removing in-family sequence redundacny, if `--skip_sequence_redundancy_removal` is set to `false`,
then the family MSA is recalculated.
If the `--alignment_tool` is `mafft`, then this `mafft_align` folder will be created, containing the updated family MSA files.

[mafft](https://github.com/GSLBiotech/mafft) is a fast but not very sensitive multiple sequence alignment tool.

### ClipKIT for updating families

<details markdown="1">
<summary>Output files</summary>

- `update_families/`
  - `full_msa/`
    - `clipkit/`
      - `<samplename>/`
        - `<family_id>.clipkit`: gap-clipped fasta files of aligned amino acid sequences

</details>

If the `--skip_msa_trimming` parameter was set to `false`, then `clipkit` runs, and according to the `--gap_threshold` parameter,
gaps (above that threshold, across all aligned sequences) are either removed only at the ends of the MSA if `trim_ends_only` is set to `true`, or throughout the alignment otherwise.
Results are stored in the `update_families/full_msa` folder.

[ClipKIT](https://github.com/JLSteenwyk/ClipKIT) is a fast and flexible alignment trimming tool that keeps phylogenetically informative sites and removes others.

### Extract family representatives

<details markdown="1">
<summary>Output files</summary>

- `family_reps/`
  - `<samplename>/`
    - `<samplename>_meta_mqc.csv`: CSV file with metadata to print with MultiQC (column headers: Sample Name,Family Id,Size,Representative Length,Representative Id,Sequence)
    - `<samplename>_reps.faa`: fasta file of all family representative sequences (one sequence per family)
    - `<samplename>.tsv`: 2-column TSV file with family ids and all sequence member ids
- `update_families/`
  - `family_reps/`
    - `<samplename>/`
      - `<samplename>_meta_mqc.csv`: CSV file with metadata to print with MultiQC (column headers: Sample Name,Family Id,Size,Representative Length,Representative Id,Sequence)
      - `<samplename>_reps.faa`: fasta file of all family representative sequences (one sequence per family)

</details>

The final report of the nf-core/proteinfamilies pipeline.
The `*_meta_mqc.csv` file are used to report family metadata and statistics in the browser, via the MultiQC software.
The `*_reps.faa` protein fasta file contains all family representative sequence in one place.
This file can be further used as input in other pipelines such as nf-core/proteinfold for structural prediction
or in fasta annotation pipelines.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

Custom output MultiQC data includes a metadata file (`multiqc_data/multiqc_family_metadata.txt`) with family information such as: Sample,Family Id,Size,Representative Length,Representative Id,Sequence

This custom metadata is presented as a data table in the MultiQC report file.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

### nf-core/proteinfold

<details markdown="1">
<summary>Output files</summary>

- `proteinfold/`
  - `<samplename>/`
    - `<samplename>_reps.faa`: A copy of the amino acid fasta file with all family representative sequences.
  - `samplesheet.csv`: Downstream samplesheet to be used as the `nf-core/proteinfold` input.

</details>

[nf-core/proteinfold](https://nf-co.re/proteinfold) is a bioinformatics best-practice analysis pipeline for protein 3D structure prediction.
The samplesheet contains two columns; `id` and `fasta`, where `id` is the sequence identifier, and `fasta` the path to the sequence file.

Example samplesheet:

```csv title="samplesheet.csv"
id,fasta
T1024,https://raw.githubusercontent.com/nf-core/test-datasets/proteinfold/testdata/sequences/T1024.fasta
T1026,https://raw.githubusercontent.com/nf-core/test-datasets/proteinfold/testdata/sequences/T1026.fasta
```

An `nf-core/proteinfold` run command would look something like this:

```
nextflow run proteinfold -profile singularity,gpu --input /path/to/proteinfamilies/results/proteinfold/samplesheet.csv --outdir result --split_fasta --use_gpu true --mode alphafold2 --alphafold2_mode split_msa_prediction --alphafold2_db '/path/to/alphafold_db' --alphafold2_params_link '/path/to/alphafold_db/' --foldseek_search easysearch --foldseek_db pdb --foldseek_db_path '/path/to/foldseek/8-ef4e960/pdb/'
```

For more information, visit the [usage page](https://nf-co.re/proteinfold/dev/docs/usage) of the `nf-core/proteinfold` pipeline.
