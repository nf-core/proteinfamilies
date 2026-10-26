#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/proteinfamilies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/proteinfamilies
    Website: https://nf-co.re/proteinfamilies
    Slack  : https://nfcore.slack.com/channels/proteinfamilies
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROTEINFAMILIES         } from './workflows/proteinfamilies'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_proteinfamilies_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Typed, so that values given on the command line are cast; parameters the config reads
// itself (output, publishing, logging, institutional and test data options) stay in nextflow.config
params {
    // Input options
    input: Path?

    // Preprocessing
    skip_preprocessing: Boolean
    min_seq_length: Integer = 30
    max_seq_length: Integer = 5000
    deduplicate_by: String = 'name' // ['name', 'sequence']
    // Clustering
    clustering_tool: String = 'cluster' // ['linclust', 'cluster']
    clustering_min_seq_identity: Float = 0.3
    clustering_min_coverage: Float = 0.5
    clustering_cov_mode: Integer = 0
    clustering_min_cluster_size: Integer = 25
    // Family generation
    family_generation_algorithm: String = 'standard' // ['standard', 'iterative']
    iterative_clusters_per_chunk: Integer = 1000
    alignment_tool: String = 'famsa' // ['famsa', 'mafft']
    skip_update_refinement: Boolean
    // Seed MSA trimming
    skip_seed_msa_trimming: Boolean
    seed_msa_trimming_ends_only: Boolean = true
    seed_msa_trimming_max_gap_fraction: Float = 0.5
    // Recruiting and search
    skip_recruiting: Boolean
    search_evalue_cutoff: Float = 0.001
    recruit_min_model_coverage: Float = 0.9
    // Family redundancy
    family_redundancy_removal: String = 'all' // ['all', 'created_only', 'none']
    family_merging: String = 'all' // ['all', 'created_only', 'none']
    merged_family_name: String = 'existing' // ['existing', 'new']
    family_redundancy_min_model_coverage: Float = 1
    family_similarity_min_model_coverage: Float = 0.9
    // Sequence redundancy
    skip_sequence_redundancy_removal: Boolean
    seq_redundancy_min_seq_identity: Float = 0.9
    seq_redundancy_min_coverage: Float = 0.9
    seq_redundancy_cov_mode: Integer = 0
    // Phylogeny
    run_phylogenetic_inference: Boolean

    // MultiQC options
    multiqc_config: Path?
    multiqc_title: String?
    multiqc_logo: Path?
    max_multiqc_email_size: String = '25.MB'
    multiqc_methods_description: Path?

    // Boilerplate options
    email: String?
    email_on_fail: String?
    plaintext_email: Boolean
    help_full: Boolean
    show_hidden: Boolean
    version: Boolean
    validate_params: Boolean = true

    // Defaults in nextflow.config, which reads them; typed here so command line values are cast
    save_intermediates: Boolean
    monochrome_logs: Boolean
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_PROTEINFAMILIES {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    PROTEINFAMILIES (
        samplesheet,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    family_reps             = PROTEINFAMILIES.out.family_reps
    passed_through_families = PROTEINFAMILIES.out.passed_through_families
    merged_families         = PROTEINFAMILIES.out.merged_families
    multiqc_report          = PROTEINFAMILIES.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PROTEINFAMILIES (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_PROTEINFAMILIES.out.multiqc_report
    )

    protein_reps_samplesheet = NFCORE_PROTEINFAMILIES.out.family_reps
        .map { meta, file ->
            [
                id: meta.id,
                fasta: file
            ]
        }

    publish:
    family_reps                  = protein_reps_samplesheet
    passed_through_families      = NFCORE_PROTEINFAMILIES.out.passed_through_families
    merged_families              = NFCORE_PROTEINFAMILIES.out.merged_families
}

output {
    // Per-sample reports are named <id>_<report>.tsv, so the sample folder comes from the file name
    passed_through_families {
        path { report -> "families/${report.name - '_passed_through_existing_families.tsv'}/" }
        mode params.publish_dir_mode
    }

    merged_families {
        path { report -> "families/${report.name - '_merged_families.tsv'}/" }
        mode params.publish_dir_mode
    }

    // Input samplesheet for nf-core/proteinfold and nf-core/proteinannotator
    family_reps {
        path { sample -> "families/${sample.id}/" }
        mode params.publish_dir_mode
        index {
            path 'families/samplesheet.csv'
            header true
            sep ','
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
