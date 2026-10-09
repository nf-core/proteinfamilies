---
title: Contributing
markdownPlugin: checklist
---

# `nf-core/proteinfamilies`: Contributing guidelines

Hi there!
Thanks for taking an interest in improving nf-core/proteinfamilies.

This page describes the recommended nf-core way to contribute to both nf-core/proteinfamilies and nf-core pipelines in general, including:

- [General contribution guidelines](#general-contribution-guidelines): common procedures or guides across all nf-core pipelines.
- [Pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines): procedures or guides specific to the development conventions of nf-core/proteinfamilies.

> [!NOTE]
> If you need help using or modifying nf-core/proteinfamilies, ask on the nf-core Slack [#proteinfamilies](https://nfcore.slack.com/channels/proteinfamilies) channel ([join our Slack here](https://nf-co.re/join/slack)).

## General contribution guidelines

### Contribution quick start

To contribute code to any nf-core pipeline:

- [ ] Ensure you have Nextflow, nf-core tools, and nf-test installed. See the [nf-core/tools repository](https://github.com/nf-core/tools) for instructions.
- [ ] Check whether a GitHub [issue](https://github.com/nf-core/proteinfamilies/issues) about your idea already exists. If an issue does not exist, create one so that others are aware you are working on it.
- [ ] [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [nf-core/proteinfamilies repository](https://github.com/nf-core/proteinfamilies) to your GitHub account.
- [ ] Create a branch on your forked repository and make your changes following [pipeline conventions](#pipeline-contribution-conventions) (if applicable).
- [ ] To fix major bugs, name your branch `patch` and follow the [patch release](#patch-release) process.
- [ ] Update relevant documentation within the `docs/` folder, use nf-core/tools to update `nextflow_schema.json`, and update `CITATIONS.md`.
- [ ] Run and/or update tests. See [Testing](#testing) for more information.
- [ ] [Lint](#lint-tests) your code with nf-core/tools.
- [ ] Submit a pull request (PR) against the `dev` branch and request a review.

If you are not used to this workflow with Git, see the [GitHub documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or [Git resources](https://try.github.io/) for more information.

## Use of AI and LLMs

The nf-core stance on the use of AI and LLMs is that humans are still ultimately responsible for their submitted code, regardless of the tools they use.

If you’re using AI tools, try to stick by these guidelines:

- Keep PRs as small and focused as possible
- Avoid any unnecessary changes, such as moving or refactoring code (unless that is the explicit intention of the PR)
- Review all generated code yourself before opening a PR, and ensure that you understand it
- Engage with the community review process and expect to make revisions

For more detail, see the [blog post](https://nf-co.re/blog/2026/statement-on-ai) for a statement from the nf-core/core team.

### Getting help

For further information and help, see the [nf-core/proteinfamilies documentation](https://nf-co.re/proteinfamilies/usage) or ask on the nf-core [#proteinfamilies](https://nfcore.slack.com/channels/proteinfamilies) Slack channel ([join our Slack here](https://nf-co.re/join/slack)).

### GitHub Codespaces

You can contribute to nf-core/proteinfamilies without installing a local development environment on your machine by using [GitHub Codespaces](https://github.com/codespaces).

[GitHub Codespaces](https://github.com/codespaces) is an online developer environment that runs in your browser, complete with VS Code and a terminal.
Most nf-core repositories include a devcontainer configuration, which creates a GitHub Codespaces environment specifically for Nextflow development.
The environment includes pre-installed nf-core tools, Nextflow, and a few other helpful utilities via a Docker container.

To get started, open the repository in [Codespaces](https://github.com/nf-core/proteinfamilies/codespaces).

### Testing

Once you have made your changes, run the pipeline with nf-test to test them locally.
For additional information, use the `--verbose` flag to view the Nextflow console log output.

```bash
nf-test test --tag test --profile +docker --verbose
```

If you have added new functionality, ensure you update the test assertions in the `.nf.test` files in the `tests/` directory.
Update the snapshots with the following command:

```bash
nf-test test --tag test --profile +docker --verbose --update-snapshots
```

When you create a pull request with changes, GitHub Actions will run automatic tests.
Pull requests are typically reviewed when these tests are passing.

Two types of tests are typically run:

#### Lint tests

nf-core has a [set of guidelines](https://nf-co.re/docs/specifications/overview) which all pipelines must follow.
To enforce these, run linting with nf-core/tools:

```bash
nf-core pipelines lint <pipeline_directory>
```

If you encounter failures or warnings, follow the linked documentation printed to screen.
For more information about linting tests, see [nf-core/tools API documentation](https://nf-co.re/docs/nf-core-tools/api_reference/latest/pipeline_lint_tests/actions_awsfulltest).

#### Pipeline tests

Each nf-core pipeline should be set up with a minimal set of test data.
GitHub Actions runs the pipeline on this data to ensure it runs through and exits successfully.
If there are any failures then the automated tests fail.
These tests are run with the latest available version of Nextflow and the minimum required version specified in the pipeline code.

### Patch release

> [!WARNING]
> Only in the unlikely event of a release that contains a critical bug.

- [ ] Create a new branch `patch` on your fork based on `upstream/main` or `upstream/master`.
- [ ] Fix the bug and use nf-core/tools to bump the version to the next semantic version, for example, `1.2.3` → `1.2.4`.
- [ ] Open a Pull Request from `patch` directly to `main`/`master` with the changes.

### Pipeline contribution conventions

nf-core semi-standardises how you write code and other contributions to make the nf-core/proteinfamilies code and processing logic more understandable for new contributors and to ensure quality.

#### Add a new pipeline step

To contribute a new step to the pipeline, follow the general nf-core coding procedure.
Please also refer to the [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines):

- [ ] Define the corresponding [input channel](#channel-naming-schemes) into your new process from the expected previous process channel.
- [ ] Install a module with nf-core/tools, or write a local module (see [default processes resource requirements](#default-processes-resource-requirements)), and add it to the target `<workflow>.nf`.
- [ ] Define the output channel if needed. Mix the version output channel into `ch_versions` and relevant files into `ch_multiqc`.
- [ ] Add new or updated parameters to the `params` block of `main.nf` with a type and a [default value](#default-parameter-values).
- [ ] Add new or updated parameters and relevant help text to `nextflow_schema.json` with [nf-core/tools](#default-parameter-values).
- [ ] Add validation for relevant parameters to the pipeline utilisation section of `utils_nfcore_\_pipeline/main.nf` subworkflow.
- [ ] Perform local tests to validate that the new code works as expected.
  - [ ] If applicable, add a new test in the `tests` directory.
- [ ] Update `usage.md`, `output.md`, and `citation.md` as appropriate.
- [ ] [Lint](#lint-tests) the code with nf-core/tools.
- [ ] Update any diagrams or pipeline images as necessary.
- [ ] Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name cleanup, and module plots are in the appropriate order.
- [ ] If applicable, create a [MultiQC](https://seqera.io/multiqc/) module.
- [ ] Add a description of the output files and, if relevant, images from the MultiQC report to `docs/output.md`.

To update the minimum required Nextflow version, see the [Nextflow version bumping](#nextflow-version-bumping) section below. For more information about pipeline contributions, see [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines).

#### Channel naming schemes

Use the following naming schemes for channels to make the channel flow easier to understand:

- Initial process channel: `ch_output_from_<process>`
- Intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

#### Default parameter values

Parameters should be declared with a type and default value in the `params` block of `main.nf`, so that values given on the command line are cast (e.g. `min_seq_length: Integer = 30`).
Only parameters the config reads itself (e.g. `outdir`, `publish_dir_mode`, `save_intermediates`) get their default in the `params` scope of `nextflow.config`; booleans among them are also typed, without a default, in `main.nf`.
Every parameter should also be documented in the pipeline JSON schema.

To update `nextflow_schema.json`, run:

```bash
nf-core pipelines schema build
```

The schema builder interface that loads in your browser should automatically update the defaults in the parameter documentation.

#### Default processes resource requirements

If you write a local module, specify a default set of resource requirements for the process.

Sensible defaults for process resource requirements (CPUs, memory, time) should be defined in `conf/base.config`.
Specify these with generic `withLabel:` selectors, so they can be shared across multiple processes and steps of the pipeline.

nf-core provides a set of standard labels that you should follow where possible, as seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config).
These labels define resource defaults for single-core processes, modules that require a GPU, and different levels of multi-core configurations with increasing memory requirements.

Values assigned within these labels can be dynamically passed to a tool using the `${task.cpus}` and `${task.memory}` Nextflow variables in the `script:` block of a module (see an example in the [modules repository](https://github.com/nf-core/modules/blob/bd1b6a40f55933d94b8c9ca94ec8c1ea0eaf4b82/modules/nf-core/samtools/bam2fq/main.nf#L30)).

#### Nextflow version bumping

If you use a new feature from core Nextflow, bump the minimum required Nextflow version in the pipeline with:

```bash
nf-core pipelines bump-version --nextflow . <min_nf_version>
```

#### Images and figures guidelines

If you update images or graphics, follow the nf-core [style guidelines](https://nf-co.re/docs/community/brand/workflow-schematics).

## Pipeline specific contribution guidelines

When contributing to nf-core/proteinfamilies, please keep the following pipeline-specific checks in mind:

- Keep changes to protein family generation, updating, merging, and redundancy removal covered by the relevant `nf-test` tests under `tests/`, `modules/local/`, and `subworkflows/local/`.
- If you change input requirements, output files, or parameter behaviour, update `nextflow_schema.json`, `docs/usage.md`, `docs/output.md`, and the affected test profiles in `conf/test*.config`.
- Keep local helper scripts in `bin/` deterministic, streaming-friendly where possible, and compatible with compressed FASTA/HMMER/MMseqs2 outputs used by the existing modules.
- Use nf-core modules for shared tools where possible; keep pipeline-specific logic in `modules/local/` or `subworkflows/local/` with matching `meta.yml`, tests, and snapshots.
- For changes that affect the family creation or update workflows, run at least the smallest relevant profile, such as `test_minimal`, `test`, `test_update`, or `test_merge`, before opening a pull request.
- Check that generated family identifiers, member tables, representative sequences, HMM libraries, family archives, and the downstream samplesheet remain stable or document any intentional changes clearly in the pull request.

### Parameter naming

New parameters follow the naming used across the pipeline, so users can guess a name from the step it controls:

- Value parameters are named `<stage>_<attribute>`, after the pipeline stage they control, not after the tool that runs it: `clustering_*`, `iterative_*`, `seed_msa_trimming_*`, `search_*`, `recruit_*`, `family_redundancy_*`, `family_similarity_*`, `seq_redundancy_*` (e.g. `--recruit_min_model_coverage`, not `--hmmsearch_query_length_threshold`). Tools may change; the names should not have to.
- Thresholds say which bound they set, `min_*` or `max_*` (e.g. `--clustering_min_seq_identity`, `--seed_msa_trimming_max_gap_fraction`); E-value thresholds end in `_evalue_cutoff`.
- Steps that run by default are turned off with `skip_<step>` (default `false`); optional steps are turned on with `run_<step>` (default `false`). Avoid double negatives such as a `skip_*` parameter that defaults to `true`.
- A policy with more than two states is one enum parameter (e.g. `--family_redundancy_removal all|created_only|none`, `--deduplicate_by name|sequence`), not several booleans that can contradict each other. List its options in a comment after its default in `main.nf` (e.g. `family_merging: String = 'all' // ['all', 'created_only', 'none']`).
- `--save_intermediates` is the only `save_*` parameter (see below); do not add per-step `save_*` flags.
- Keep `nextflow_schema.json` groups in pipeline order (preprocessing, clustering, family generation, seed MSA trimming, recruiting and search, family redundancy, sequence redundancy, phylogeny), and the `params` block of `main.nf` in the same order.

### Outputs

- Final results are published at paths that do not depend on parameters: `qc/<id>/`, `clustering/<id>/`, `families/<id>/` (HMM library, seed MSA, full MSA and FASTA archives, members, representatives, reports), `families/samplesheet.csv` and `phylogeny/<id>/`. A final file is absent only when it would have no content, and `docs/output.md` lists each case.
- Every other file is an intermediate: publish it under `intermediates/` with `enabled: params.save_intermediates` in `conf/modules.config` (the default `publishDir` already does this for processes without their own entry).
- Do not publish the same result twice (e.g. a loose copy of files already in an archive).
