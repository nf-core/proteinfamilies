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
include { FAA_SEQFU_SEQKIT     } from '../subworkflows/nf-core/faa_seqfu_seqkit/main'
include { UPDATE_FAMILIES      } from '../subworkflows/local/update_families'
include { MMSEQS_FASTA_CLUSTER } from '../subworkflows/nf-core/mmseqs_fasta_cluster'
include { GENERATE_FAMILIES    } from '../subworkflows/local/generate_families'
include { REMOVE_REDUNDANCY    } from '../subworkflows/local/remove_redundancy'

//
// MODULE: Local to the pipeline
//
include { CALCULATE_CLUSTER_DISTRIBUTION } from '../modules/local/calculate_cluster_distribution/main'
include { CHUNK_CLUSTERS                 } from '../modules/local/chunk_clusters/main'
include { EXTRACT_FAMILY_MEMBERS         } from '../modules/local/extract_family_members/main'
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

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_family_reps   = channel.empty()

    ch_samplesheet_for_create = channel.empty()
    ch_samplesheet_for_update = channel.empty()

    ch_input_for_qc = ch_samplesheet
        .map { meta, fasta, _existing_hmms_to_update, _existing_msas_to_update ->
            [ meta, fasta ]
        }

    FAA_SEQFU_SEQKIT( ch_input_for_qc, params.skip_preprocessing )
    ch_versions = ch_versions.mix( FAA_SEQFU_SEQKIT.out.versions )

    // Replace input fasta and join back in samplesheet to ensure in sync in case of multiple sequence files
    ch_samplesheet_updated = ch_samplesheet
        .combine(FAA_SEQFU_SEQKIT.out.fasta, by: 0)
        .map {
            meta, _fasta, existing_hmms, existing_msas, updated_fasta ->
            [ meta, updated_fasta, existing_hmms, existing_msas ]
        }

    ch_branch_result = ch_samplesheet_updated
        .branch { _meta, _updated_fasta, existing_hmms_to_update, existing_msas_to_update ->
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
        .map { meta, updated_fasta, _existing_hmms, _existing_msas ->
            [meta, updated_fasta]
        }
    ch_samplesheet_for_update = ch_branch_result.to_update

    // Updating existing families
    if (ch_branch_result.to_update) {
        UPDATE_FAMILIES (
            ch_samplesheet_for_update,
            params.hmmsearch_query_length_threshold,
            params.skip_sequence_redundancy_removal,
            params.clustering_tool,
            params.alignment_tool,
            params.skip_msa_trimming,
            params.clipkit_out_format
        )
        ch_versions = ch_versions.mix( UPDATE_FAMILIES.out.versions )

        ch_family_reps = ch_family_reps.mix( UPDATE_FAMILIES.out.updated_family_reps )
        ch_samplesheet_for_create = ch_samplesheet_for_create.mix( UPDATE_FAMILIES.out.no_hit_seqs )
    }

    // Creating new families
    // Clustering
    MMSEQS_FASTA_CLUSTER (
        ch_samplesheet_for_create,
        params.clustering_tool
    )
    ch_versions = ch_versions.mix( MMSEQS_FASTA_CLUSTER.out.versions )

    CALCULATE_CLUSTER_DISTRIBUTION( MMSEQS_FASTA_CLUSTER.out.clusters )
    ch_versions = ch_versions.mix( CALCULATE_CLUSTER_DISTRIBUTION.out.versions )

    CHUNK_CLUSTERS( MMSEQS_FASTA_CLUSTER.out.clusters, MMSEQS_FASTA_CLUSTER.out.seqs, params.cluster_size_threshold )
    ch_versions = ch_versions.mix( CHUNK_CLUSTERS.out.versions )

    ch_fasta_chunks = CHUNK_CLUSTERS.out.fasta_chunks
        .transpose()
        .map { meta, file_path ->
            [ [id: meta.id, chunk: file(file_path, checkIfExists: true).baseName.tokenize('_').last()], file_path ]
        }

    // Multiple sequence alignments, model building and sequence recruiting
    GENERATE_FAMILIES (
        ch_samplesheet_for_create,
        ch_fasta_chunks,
        params.alignment_tool,
        params.skip_msa_trimming,
        params.clipkit_out_format,
        params.hmmsearch_write_target,
        params.hmmsearch_write_domain,
        params.skip_additional_sequence_recruiting,
        params.hmmsearch_query_length_threshold
    )
    ch_versions = ch_versions.mix( GENERATE_FAMILIES.out.versions )

    // Remove redundant sequences and families
    REMOVE_REDUNDANCY (
        ch_samplesheet_for_create,
        GENERATE_FAMILIES.out.seed_msa,
        GENERATE_FAMILIES.out.full_msa,
        GENERATE_FAMILIES.out.fasta,
        GENERATE_FAMILIES.out.hmm,
        params.skip_family_redundancy_removal,
        params.skip_family_merging,
        params.hmmsearch_family_redundancy_length_threshold,
        params.hmmsearch_family_similarity_length_threshold,
        params.skip_sequence_redundancy_removal,
        params.clustering_tool,
        params.alignment_tool,
        params.skip_msa_trimming,
        params.clipkit_out_format,
        params.hmmsearch_write_target,
        params.hmmsearch_write_domain,
        params.skip_additional_sequence_recruiting,
        params.hmmsearch_query_length_threshold
    )
    ch_versions = ch_versions.mix( REMOVE_REDUNDANCY.out.versions )

    // Post-processing
    ch_fasta = REMOVE_REDUNDANCY.out.fasta
        .map { meta, aln -> [ [id: meta.id], aln ] }
        .groupTuple(by: 0)

    EXTRACT_FAMILY_MEMBERS( ch_fasta )
    ch_versions = ch_versions.mix( EXTRACT_FAMILY_MEMBERS.out.versions )

    EXTRACT_FAMILY_REPS( ch_fasta )
    ch_versions = ch_versions.mix( EXTRACT_FAMILY_REPS.out.versions )
    ch_family_reps = ch_family_reps.mix( EXTRACT_FAMILY_REPS.out.map )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'proteinfamilies_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    ch_multiqc_files = ch_multiqc_files.mix(FAA_SEQFU_SEQKIT.out.multiqc_files.collect{it[1]}.ifEmpty([]))
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
    family_reps    = EXTRACT_FAMILY_REPS.out.fasta
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
