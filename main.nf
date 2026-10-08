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
