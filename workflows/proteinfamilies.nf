/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { UPDATE_FAMILIES    } from '../subworkflows/local/update_families'
include { EXECUTE_CLUSTERING } from '../subworkflows/local/execute_clustering'
include { GENERATE_FAMILIES  } from '../subworkflows/local/generate_families'
include { REMOVE_REDUNDANCY  } from '../subworkflows/local/remove_redundancy'

//
// MODULE: Local to the pipeline
//
include { CALCULATE_CLUSTER_DISTRIBUTION } from "../modules/local/calculate_cluster_distribution/main"
include { CHUNK_CLUSTERS                 } from '../modules/local/chunk_clusters/main'
include { EXTRACT_FAMILY_REPS            } from '../modules/local/extract_family_reps/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROTEINFAMILIES {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_family_reps   = Channel.empty()

    ch_samplesheet_for_create = Channel.empty()
    ch_samplesheet_for_update = Channel.empty()

    ch_branch_result = ch_samplesheet
        .branch { _meta, _fasta, existing_hmms_to_update, existing_msas_to_update ->
            to_create: !existing_hmms_to_update?.size() && !existing_msas_to_update?.size()
            to_update: existing_hmms_to_update?.size() && existing_msas_to_update?.size()
        }

    /************************************/
    /* Splitting the samplesheet into 2 */
    /* - Entries to create new families */
    /*   (they only have sequences)     */
    /* - Entries to update existing     */
    /*   families (existing HMM models  */
    /*   and MSAs)                      */
    /************************************/
    ch_samplesheet_for_create = ch_branch_result.to_create
        .map { meta, fasta, _existing_hmms, _existing_msas ->
            [meta, fasta]
        }
    ch_samplesheet_for_update = ch_branch_result.to_update

    // Updating existing families
    if (ch_branch_result.to_update) {
        UPDATE_FAMILIES( ch_samplesheet_for_update )
        ch_versions = ch_versions.mix( UPDATE_FAMILIES.out.versions )

        ch_family_reps = ch_family_reps.mix( UPDATE_FAMILIES.out.updated_family_reps )
        ch_samplesheet_for_create = ch_samplesheet_for_create.mix( UPDATE_FAMILIES.out.no_hit_seqs )
    }

    // Creating new families
    // Clustering
    EXECUTE_CLUSTERING( ch_samplesheet_for_create, params.clustering_tool )
    ch_versions = ch_versions.mix( EXECUTE_CLUSTERING.out.versions )

    CALCULATE_CLUSTER_DISTRIBUTION( EXECUTE_CLUSTERING.out.clusters )
    ch_versions = ch_versions.mix( CALCULATE_CLUSTER_DISTRIBUTION.out.versions )

    CHUNK_CLUSTERS( EXECUTE_CLUSTERING.out.clusters, EXECUTE_CLUSTERING.out.seqs, params.cluster_size_threshold )
    ch_versions = ch_versions.mix( CHUNK_CLUSTERS.out.versions )

    // Multiple sequence alignment
    GENERATE_FAMILIES( ch_samplesheet_for_create, CHUNK_CLUSTERS.out.fasta_chunks )
    ch_versions = ch_versions.mix( GENERATE_FAMILIES.out.versions )

    // Remove redundant sequences and families
    REMOVE_REDUNDANCY( GENERATE_FAMILIES.out.seed_msa, GENERATE_FAMILIES.out.full_msa, GENERATE_FAMILIES.out.fasta, GENERATE_FAMILIES.out.hmm )
    ch_versions = ch_versions.mix( REMOVE_REDUNDANCY.out.versions )

    // Post-processing
    ch_full_msa = REMOVE_REDUNDANCY.out.full_msa
        .map { meta, aln -> [ [id: meta.id], aln ] }
        .groupTuple(by: 0)

    EXTRACT_FAMILY_REPS( ch_full_msa )
    ch_versions = ch_versions.mix( EXTRACT_FAMILY_REPS.out.versions )
    ch_family_reps = ch_family_reps.mix( EXTRACT_FAMILY_REPS.out.map )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    ch_multiqc_files = ch_multiqc_files.mix(CALCULATE_CLUSTER_DISTRIBUTION.out.mqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_family_reps.collect { it[1] }.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
