# nf-core/proteinfamilies: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/proteinfamilies/usage](https://nf-co.re/proteinfamilies/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

**nf-core/proteinfamilies** is a bioinformatics pipeline that generates protein families from amino acid sequences and/or updates existing families with new sequences.
It takes a protein fasta file as input, clusters the sequences and then generates protein family Hidden Markov Models (HMMs) along with their multiple sequence alignments (MSAs).
Optionally, existing family HMMs (and their seed and/or full MSAs) can be given in order to update those families with new sequences in case of matching hits.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 2 mandatory and 3 optional columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

```csv
id,fasta,existing_hmms,existing_seed_msas,existing_full_msas
CONTROL_REP1,amino_acid_sequences_input.faa,,,
CONTROL_REP2,amino_acid_sequences_extra.faa.gz,existing_hmms.tar.gz,,
CONTROL_REP3,amino_acid_sequences_extra.faa.gz,existing_hmms.tar.gz,existing_seed_msas.tar.gz,existing_full_msas.tar.gz
```

| Column               | Description                                                                                                                                                                                                                                                                                                                         |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `id`                 | Custom sample name. Only letters, digits, dots (`.`), underscores (`_`) and dashes (`-`) are allowed, since the sample name is used to build output paths and to parse family names back out of them.                                                                                                                               |
| `fasta`              | Full path to amino acid fasta file. Allowed extensions are ".faa", ".fasta" and ".fa", with or without a following ".gz" for gzipped files.                                                                                                                                                                                         |
| `existing_hmms`      | (Optional) Full path to a ".tar.gz" archive of HMM files, or to one HMM library (".hmm", ".hmm.gz", ".lib" or ".lib.gz", e.g. this pipeline's `<id>.lib.gz`). A sample with existing HMMs is updated: its sequences are searched against these families, and only the sequences without hits go on to create new families.          |
| `existing_seed_msas` | (Optional, needs `existing_hmms`) Full path to a ".tar.gz" archive with seed MSAs (aligned FASTA, as this pipeline writes them; optionally gzipped), each named after its family's HMM `NAME` (as file name stem).                                                                                                                  |
| `existing_full_msas` | (Optional, needs `existing_hmms`) Full path to a ".tar.gz" archive with full MSAs (aligned FASTA, as this pipeline writes them, or Stockholm; optionally gzipped), each named after its family's HMM `NAME`. Their members are pooled with the input sequences and searched again, so families keep the old members that still hit. |

Input sequences named `<sequence>/<start>-<end>` (Pfam convention) are treated as slices of `<sequence>`: family members cut from them are named in the parent sequence's coordinates (a hit on residues 3-180 of `seqA/10-200` becomes `seqA/12-189`). Any other name is taken as a full protein.

### Updating existing families

Each model in `existing_hmms` is an existing family, identified by its `NAME` (case-sensitive; letters, digits, `.`, `_` and `-` only, as in Pfam); HMM files may hold one or more models, and their file names do not matter. Seed and full MSAs must be named after that `NAME`, without extensions (e.g. `fam_1.aln` for `NAME  fam_1`). Not every family needs an MSA, but every MSA file needs a family.

> [!WARNING]
> hmmsearch reports hits by HMM `NAME`, while MSAs are matched to their family by file name. The pipeline stops if two models share a `NAME`, if two MSA files give the same family, or if an MSA file is not named after an existing HMM `NAME`. It also stops if an existing family is named like the families this run creates for the sample (`<id>_<number>...`, e.g. after a previous run with the same `id`): use a new `id` for the update run, such as `<id>_r2`.
> HMM and MSA archives written by nf-core/proteinfamilies already follow these rules and can be used as they are.

The input sequences, together with the members of any `existing_full_msas` (gaps removed; a member `seq/<start>-<end>` is skipped if its region lies inside an input sequence of the same protein `seq`, where a name without a range is the whole protein, or inside another member; partial overlaps are kept), are searched against the existing HMMs.
Each family's hits are then rebuilt like a newly created family: optionally made non-redundant, aligned and trimmed into a new seed MSA, built into a new HMM, and used to recruit the new full MSA from the same pool (the new seed MSA serves as the full MSA with `--skip_recruiting`).
With `--skip_update_refinement`, the existing HMMs are kept instead: each one aligns its hits into the new full MSA (hmmalign), and its `existing_seed_msas` file, if given, passes through unchanged. Seed MSAs are never searched, so sequences only found in a seed MSA must also be in the `fasta` or in a full MSA to stay in their family.
Families without any hit, or whose rebuilt HMM recruits nothing, pass through as given (their existing HMM, seed and full MSA) and are listed with the reason in `update_families/passed_through_families/<id>_passed_through_existing_families.tsv`.
Input sequences that end up in no updated family go to family creation (a sample without any hit sends all of them); members of existing full MSAs that no family holds anymore are dropped and never create new families.
Updated families then go through redundancy removal together with the families created for the sample. An updated family is never removed: a created family redundant with it is, and two redundant updated families are both kept. The removed created family's sequences are dropped whatever its size, as in any redundant pair: an update keeps the existing families and adds to them. To keep every sequence instead, pool the old and new sequences and create the families anew. With `--family_redundancy_removal created_only`, updated families bypass this check, so created families may duplicate them. Similar families can be merged, but never two updated families: a merge holds at most one, so curated families keep their identity (to combine existing families, pool their sequences and create the families anew). When a pool links several updated families, directly or through created families, those updated families stay unmerged and only the pool's created families are merged. A merge holding an updated family recruits from the update pool (input sequences plus existing full MSA members) and keeps that family's name (`--merged_family_name new` names it like other merges instead). With `--family_merging created_only`, updated families are left out of merging (`none` turns redundancy removal or merging off for every family). With `--skip_update_refinement`, updated families are never merged (as with `--family_merging created_only`), since a merge rebuilds the family's HMM. Like created families, updated families then lose their redundant members (sequence redundancy removal, off with `--skip_sequence_redundancy_removal`), and their full MSAs are re-aligned from the remaining members (`--alignment_tool`). Passed-through families skip both: their members and full MSA stay as given.

Every run writes each sample's final families to `archives/<id>/<id>_{hmms,seed_msas,full_msas}.tar.gz` (see [output](output.md#archives-of-final-families)). To update them later, give the three archives in the existing columns, with the new sequences in `fasta` and a new `id` (the families created for `<id>` are named `<id>_<number>`, so reusing it would clash):

```csv title="samplesheet.csv"
id,fasta,existing_hmms,existing_seed_msas,existing_full_msas
s1_r2,new_sequences.faa.gz,results/archives/s1/s1_hmms.tar.gz,results/archives/s1/s1_seed_msas.tar.gz,results/archives/s1/s1_full_msas.tar.gz
```

### Migrating from v2 to v3

- **Samplesheet:** rename the columns `sample` → `id`, `existing_hmms_to_update` → `existing_hmms` and `existing_msas_to_update` → `existing_full_msas`, and add an `existing_seed_msas` column (may be left empty). MSA archives are now optional; a row with MSAs but no HMMs is rejected.
- **Existing HMMs** may be a `.tar.gz` archive or one HMM library; families are identified by HMM `NAME`, and MSA files must be named after it (see [Updating existing families](#updating-existing-families)).
- **Parameters:** `--skip_msa_trimming` is now `--skip_seed_msa_trimming`; `--clipkit_out_format`, `--save_update_families_pre_clipped_fasta` and `--save_update_families_clipped_fasta` are removed (updated families use the same `save_*` parameters as created ones); `--skip_update_refinement` is new.
- **Outputs:** `clipkit/` folders are now `trimmed/` (FASTA `.aln`). Updated families are published like created ones under `update_families/{seed_msa,hmm,full_msa}/raw/<tool>/<id>/` instead of `update_families/full_msa/<tool>/` and `update_families/fasta/`, then go through redundancy removal with the created families, so their final files sit next to them (e.g. `hmm/filtered/<id>/`, `full_msa/filtered/<tool>/<id>/`). Family representatives of all final families are in `family_reps/<id>/` (no more `update_families/family_reps/`). Existing families the update did not return pass through and are listed in `update_families/passed_through_families/`, and every sample's final families are archived under `archives/<id>/` for later updates.

Other parameters removed or renamed in v3 (a v2 command line with an old name stops with an unknown-parameter error):

| v2                                                     | v3                                                                                   |
| ------------------------------------------------------ | ------------------------------------------------------------------------------------ |
| `--hmmsearch_write_target`, `--hmmsearch_write_domain` | Removed: the per-domain table is always written, the per-target table was never used |
| `--cluster_seq_identity`                               | `--clustering_min_seq_identity`                                                      |
| `--cluster_coverage`                                   | `--clustering_min_coverage`                                                          |
| `--cluster_cov_mode`                                   | `--clustering_cov_mode`                                                              |
| `--cluster_size_threshold`                             | `--clustering_min_cluster_size`                                                      |
| `--clusters_per_chunk`                                 | `--iterative_clusters_per_chunk`                                                     |
| `--trim_ends_only`                                     | `--seed_msa_trimming_ends_only`                                                      |
| `--gap_threshold`                                      | `--seed_msa_trimming_max_gap_fraction`                                               |
| `--skip_additional_sequence_recruiting`                | `--skip_recruiting`                                                                  |
| `--hmmsearch_evalue_cutoff`                            | `--search_evalue_cutoff`                                                             |
| `--hmmsearch_query_length_threshold`                   | `--recruit_min_model_coverage`                                                       |
| `--hmmsearch_family_redundancy_length_threshold`       | `--family_redundancy_min_model_coverage`                                             |
| `--hmmsearch_family_similarity_length_threshold`       | `--family_similarity_min_model_coverage`                                             |
| `--cluster_seq_identity_for_redundancy`                | `--seq_redundancy_min_seq_identity`                                                  |
| `--cluster_coverage_for_redundancy`                    | `--seq_redundancy_min_coverage`                                                      |
| `--cluster_cov_mode_for_redundancy`                    | `--seq_redundancy_cov_mode`                                                          |

## Parameter specifications

Here we provide guidance regarding some parameter choices.

- `clustering_tool` ["cluster", "linclust"]: The mmseqs algorithm used for clustering.
  The `cluster` option is slower but more sensitive, and is recommended where there are sufficient compute resources available and a more sensitive search is called for.
  It tends to produce fewer and larger clusters than `linclust`.
  The `linclust` option is less sensitive, but extremely fast for clustering larger datasets.
- `clustering_cov_mode` [0, 1, 2]: The default bidirectional value for coverage mode (`clustering_cov_mode` = 0) automatically sets the MMseqs2 clustering mode to greedy cluster set.
  However, users can opt to override this parameter either indirectly, by changing the coverage mode, or directly, by setting the `--cluster-mode` argument in the modules configuration file.
- `alignment_tool` ["famsa", "mafft"]: Multiple Sequence Alignment (MSA) options.
  The `famsa` option is generally recommended as the best time-memory-accuracy combination.
  The `mafft` option offers various alignment strategies, but in general is slower and less sensitive than `famsa`.
- `seed_msa_trimming_ends_only`: Flag to either clip seed MSA gaps throughout the alignment, or only at the ends.
  Only used if `skip_seed_msa_trimming` is off. Full MSAs are never trimmed.
  The pipeline authors strongly recommend keeping `seed_msa_trimming_ends_only` on (default): gaps inside the sequences may still carry evolutionary significance, and only end trimming keeps row coordinates correct.

> [!WARNING]
> Trimmed MSA rows that lost residues are renamed `<sequence>/<start>-<end>` to the residues they still hold, recalculated from the residues removed at the alignment ends; rows that lost none keep their name.
> With `--seed_msa_trimming_ends_only false`, residues removed from interior columns are **not** reflected, so a row's range spans more residues than the row contains and no longer maps back to its exact source residues.
> Only turn it off if you need interior trimming and do not rely on row coordinates downstream.

## Family generation algorithms

`family_generation_algorithm` ["standard", "iterative"] selects how clusters become family models. Both produce the same kinds of outputs and feed the same redundancy removal, merging and reporting steps, so the choice is invisible downstream.

### `standard` (default)

Each cluster is chunked into its own FASTA file and processed by a chain of tools orchestrated by Nextflow: FAMSA or mafft aligns the cluster into a seed MSA, ClipKIT optionally trims it, `hmmbuild` builds the family HMM, and `hmmsearch` optionally recruits further members from the sample's sequence pool in a single pass, which `hmmalign` then aligns into the full MSA.

### `iterative`

Whole chunks of clusters, `iterative_clusters_per_chunk` at a time, are handed to [mgnifam](https://github.com/vagkaratzas/mgnifam), which builds a family from each cluster by looping: build an HMM, recruit members from the sequence pool, realign the expanded membership, and repeat up to three times or, until the family converges or the cluster is discarded. A single task therefore emits many families.

mgnifam performs each step in-process with its own libraries rather than by calling the pipeline's tools:

| Step               | Library                                          |
| ------------------ | ------------------------------------------------ |
| Alignment          | [pyfamsa](https://github.com/althonos/pyfamsa)   |
| Trimming           | [pytrimal](https://github.com/althonos/pytrimal) |
| HMM build & search | [pyhmmer](https://github.com/althonos/pyhmmer)   |

Because of that, the parameters below are honoured only by the `standard` algorithm. They are ignored when the `iterative` path creates or merges families, which always behaves as stated:

| Parameter                     | Behaviour of the `iterative` algorithm                                                             |
| ----------------------------- | -------------------------------------------------------------------------------------------------- |
| `alignment_tool`              | Always FAMSA, through pyfamsa                                                                      |
| `skip_seed_msa_trimming`      | Trimming is always applied, through pytrimal                                                       |
| `seed_msa_trimming_ends_only` | ClipKIT is not used; pytrimal trims by column gap occupancy (`seed_msa_trimming_max_gap_fraction`) |
| `skip_recruiting`             | Recruitment is always performed, and repeated until convergence                                    |
| `save_hmmsearch_results`      | Searching is in-process, so no hmmsearch report files exist                                        |

> [!NOTE]
> Updating existing families (samplesheet entries with existing HMMs) always runs the `standard` update path, whichever algorithm is selected, so the `standard` parameters above (e.g. `alignment_tool`, `skip_seed_msa_trimming`, `seed_msa_trimming_ends_only`, `seed_msa_trimming_max_gap_fraction`, `skip_recruiting`) apply to updated families.

The parameters both algorithms share are mapped onto their mgnifam equivalents:

| Parameter                            | mgnifam option                    |
| ------------------------------------ | --------------------------------- |
| `clustering_min_cluster_size`        | applied while chunking clusters   |
| `min_seq_length`                     | `--discard_min_rep_length`        |
| `max_seq_length`                     | `--discard_max_rep_length`        |
| `seq_redundancy_min_seq_identity`    | `--max_seq_identity`              |
| `seed_msa_trimming_max_gap_fraction` | `--max_gap_occupancy`             |
| `search_evalue_cutoff`               | `--recruit_evalue_cutoff`         |
| `recruit_min_model_coverage`         | `--recruit_hit_length_percentage` |

mgnifam's remaining options (`--discard_min_starting_membership`, `--max_seed_seqs`, `--batch_size`, `--prefetch_targets`) keep their tool defaults and can be set through `ext.args`, as described in [Custom Tool Arguments](#custom-tool-arguments).

`iterative_clusters_per_chunk` (default 1000) trades parallelism against scheduling overhead: smaller chunks give more tasks, better load balancing and finer-grained `-resume`, at the cost of more per-task startup. Set `save_iterative_family_metadata` to publish mgnifam's family roster, metadata, converged, successful and discarded records, its per-family HMM consensus match states (rf), representative sequence fasta, and its log.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/proteinfamilies --input ./samplesheet.csv --outdir ./results  -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/proteinfamilies -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/proteinfamilies
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/proteinfamilies releases page](https://github.com/nf-core/proteinfamilies/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
