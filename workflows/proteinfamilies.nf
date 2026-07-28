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

include { FAA_SEQFU_SEQKIT                                 } from '../subworkflows/nf-core/faa_seqfu_seqkit/main'
include { UPDATE_FAMILIES                                  } from '../subworkflows/local/update_families'
include { MMSEQS_FASTA_CLUSTER                             } from '../subworkflows/nf-core/mmseqs_fasta_cluster'
include { CALCULATE_CLUSTER_DISTRIBUTION                   } from '../modules/local/calculate_cluster_distribution/main'
include { CHUNK_CLUSTERS                                   } from '../modules/local/chunk_clusters/main'
include { GENERATE_FAMILIES                                } from '../subworkflows/local/generate_families'
include { REMOVE_REDUNDANCY                                } from '../subworkflows/local/remove_redundancy'
include { FIND_CONCATENATE as FIND_CONCATENATE_HMM_LIBRARY } from '../modules/nf-core/find/concatenate'
include { CMAPLE                                           } from '../modules/nf-core/cmaple/main'
include { EXTRACT_FAMILY_MEMBERS                           } from '../modules/local/extract_family_members/main'
include { EXTRACT_FAMILY_REPS                              } from '../modules/local/extract_family_reps/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Two-path pipeline: samples providing existing HMMs+MSAs go through UPDATE_FAMILIES to recruit
// new sequences into those families; all other samples take the de-novo path
// (cluster → align → HMM build). Sequences not assigned to any existing family during the
// update path are forwarded to the de-novo path so nothing is discarded.
workflow PROTEINFAMILIES {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_multiqc_files = channel.empty()
    ch_samplesheet_for_create = channel.empty()
    ch_samplesheet_for_update = channel.empty()
    ch_family_reps            = channel.empty()

    ch_input_for_qc = ch_samplesheet
        .map { meta, fasta, _existing_hmms_to_update, _existing_msas_to_update ->
            [ meta, fasta ]
        }

    FAA_SEQFU_SEQKIT( ch_input_for_qc, params.skip_preprocessing )

    // Replace input fasta and join back in samplesheet to ensure in sync in case of multiple sequence files
    ch_samplesheet_updated = ch_samplesheet
        .combine(FAA_SEQFU_SEQKIT.out.fasta, by: 0)
        .map {
            meta, _fasta, existing_hmms, existing_msas, updated_fasta ->
            [ meta, updated_fasta, existing_hmms, existing_msas ]
        }

    // ?.size() is Groovy's null-safe operator: absent samplesheet columns yield null (falsy).
    // Both HMMs and MSAs must be present (non-null, non-zero file size) to route to to_update.
    ch_branch_result = ch_samplesheet_updated
        .branch { _meta, _updated_fasta, existing_hmms_to_update, existing_msas_to_update ->
            to_create: !existing_hmms_to_update?.size() && !existing_msas_to_update?.size()
            to_update: existing_hmms_to_update?.size() && existing_msas_to_update?.size()
        }

    // Entries with existing models go to UPDATE_FAMILIES; entries with sequences only go to de-novo creation.
    ch_samplesheet_for_create = ch_branch_result.to_create
        .map { meta, updated_fasta, _existing_hmms, _existing_msas ->
            [meta, updated_fasta]
        }
    ch_samplesheet_for_update = ch_branch_result.to_update

    // Updating existing families
    UPDATE_FAMILIES (
        ch_samplesheet_for_update,
        params.hmmsearch_query_length_threshold,
        params.skip_sequence_redundancy_removal,
        params.clustering_tool,
        params.alignment_tool,
        params.skip_msa_trimming,
        params.clipkit_out_format,
        params.save_update_families_pre_clipped_fasta,
        params.save_update_families_clipped_fasta
    )

    ch_family_reps = ch_family_reps.mix( UPDATE_FAMILIES.out.updated_family_reps )
    // Sequences not assigned to any existing family during update feed the de-novo creation path.
    ch_samplesheet_for_create = ch_samplesheet_for_create.mix( UPDATE_FAMILIES.out.no_hit_seqs )

    // Creating new families
    // Clustering
    MMSEQS_FASTA_CLUSTER (
        ch_samplesheet_for_create,
        params.clustering_tool
    )

    CALCULATE_CLUSTER_DISTRIBUTION( MMSEQS_FASTA_CLUSTER.out.clusters )

    CHUNK_CLUSTERS( MMSEQS_FASTA_CLUSTER.out.clusters, MMSEQS_FASTA_CLUSTER.out.seqs, params.cluster_size_threshold )

    // tokenize('_').last() extracts the numeric suffix from filenames like 'sample_1.faa.gz' as the chunk ID.
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

    // Collect all final HMMs per sample and concatenate into a .lib.gz library.
    // Strip chunk from meta (keep only id) so all family HMMs within a sample are grouped together.
    ch_hmm_for_library = UPDATE_FAMILIES.out.hmm
        .map { meta, model -> [ [id: meta.id], model ] }
        .mix(
            REMOVE_REDUNDANCY.out.hmm
                .map { meta, model -> [ [id: meta.id], model ] }
        )
        .groupTuple()

    FIND_CONCATENATE_HMM_LIBRARY( ch_hmm_for_library )

    // Infer Phylogenetic relations of full MSAs
    if (!params.skip_phylogenetic_inference) {
        CMAPLE (
            REMOVE_REDUNDANCY.out.full_msa
                .map { meta, file -> [ meta, file, [] ] }
        )
    }

    // Post-processing
    ch_fasta = REMOVE_REDUNDANCY.out.fasta
        .map { meta, aln -> [ [id: meta.id], aln ] }
        .groupTuple(by: 0)

    EXTRACT_FAMILY_MEMBERS( ch_fasta )

    EXTRACT_FAMILY_REPS( ch_fasta )
    ch_family_reps = ch_family_reps.mix( EXTRACT_FAMILY_REPS.out.map )

    //
    // Collate and save software versions
    //
    // channel.topic("versions") replaced legacy ch_versions by catching version emissions from any
    // process that publishes to the 'versions' topic. distinct() prevents duplicates when many
    // parallel tasks emit the same tool version. Branched into Path instances (versions.yml
    // files from legacy modules) and [process, tool, version] tuples (from topic-aware modules).
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

    def ch_collated_versions = softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'proteinfamilies_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    ch_multiqc_files = ch_multiqc_files.mix(FAA_SEQFU_SEQKIT.out.multiqc_files.collect { file -> file[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CALCULATE_CLUSTER_DISTRIBUTION.out.mqc.collect { file -> file[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_family_reps.collect { file -> file[1] }.ifEmpty([]))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'proteinfamilies'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    family_reps    = EXTRACT_FAMILY_REPS.out.fasta
    multiqc_report = MULTIQC.out.report.map { _meta, report -> report } // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
