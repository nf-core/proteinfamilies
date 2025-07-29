include { UNTAR as UNTAR_HMM            } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_MSA            } from '../../../modules/nf-core/untar/main'
include { validateMatchingFolders       } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { FIND_CONCATENATE as CAT_HMM   } from '../../../modules/nf-core/find/concatenate/main'
include { HMMER_HMMSEARCH               } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { BRANCH_HITS_FASTA             } from '../../../modules/local/branch_hits_fasta'
include { SEQKIT_SEQ                    } from '../../../modules/nf-core/seqkit/seq/main'
include { FIND_CONCATENATE as CAT_FASTA } from '../../../modules/nf-core/find/concatenate/main'
include { EXECUTE_CLUSTERING            } from '../../../subworkflows/local/execute_clustering'
include { REMOVE_REDUNDANT_SEQS         } from '../../../modules/local/remove_redundant_seqs/main'
include { ALIGN_SEQUENCES               } from '../../../subworkflows/local/align_sequences'
include { CLIPKIT                       } from '../../../modules/nf-core/clipkit/main'
include { HMMER_HMMBUILD                } from '../../../modules/nf-core/hmmer/hmmbuild/main'
include { EXTRACT_FAMILY_REPS           } from '../../../modules/local/extract_family_reps/main'

workflow UPDATE_FAMILIES {
    take:
    ch_samplesheet_for_update // channel: [meta, sequences, existing_hmms_to_update, existing_msas_to_update]

    main:
    ch_versions            = Channel.empty()
    ch_updated_family_reps = Channel.empty()
    ch_no_hit_seqs         = Channel.empty()

    ch_input_for_untar = ch_samplesheet_for_update
        .multiMap { meta, _fasta, existing_hmms_to_update, existing_msas_to_update ->
            hmm: [ meta, existing_hmms_to_update ]
            msa: [ meta, existing_msas_to_update ]
        }

    UNTAR_HMM( ch_input_for_untar.hmm )
    ch_versions = ch_versions.mix( UNTAR_HMM.out.versions )

    UNTAR_MSA( ch_input_for_untar.msa )
    ch_versions = ch_versions.mix( UNTAR_MSA.out.versions )

    // check that the HMMs and the MSAs match
    // join to ensure in sync
    ch_folders_to_validate = UNTAR_HMM.out.untar
        .join(UNTAR_MSA.out.untar)
        .multiMap { meta, folder1, folder2 ->
            ch_hmm_folder: [meta, folder1]
            ch_msa_folder: [meta, folder2]
        }
    validateMatchingFolders(ch_folders_to_validate.ch_hmm_folder, ch_folders_to_validate.ch_msa_folder)

    // Squeeze the HMMs into a single file
    CAT_HMM( UNTAR_HMM.out.untar.map { meta, folder -> [meta, file("${folder.toUriString()}/*", checkIfExists: true)] } )
    ch_versions = ch_versions.mix( CAT_HMM.out.versions )

    // Prep the sequences to search against the HMM concatenated model of families
    ch_input_for_hmmsearch = CAT_HMM.out.file_out
        .combine(ch_samplesheet_for_update, by: 0)
        .map { meta, concatenated_hmm, fasta, _existing_hmms_to_update, _existing_msas_to_update -> [meta, concatenated_hmm, fasta, false, false, true] }

    HMMER_HMMSEARCH( ch_input_for_hmmsearch )
    ch_versions = ch_versions.mix( HMMER_HMMSEARCH.out.versions )

    ch_input_for_branch_hits = HMMER_HMMSEARCH.out.domain_summary
        .join(ch_samplesheet_for_update)
        .multiMap { meta, domtbl, fasta, _existing_hmms_to_update, _existing_msas_to_update ->
            domtbl: [ meta, domtbl ]
            fasta: [ meta, fasta ]
        }

    // Branch hit families/fasta proteins from non hit fasta proteins
    BRANCH_HITS_FASTA ( ch_input_for_branch_hits.fasta, ch_input_for_branch_hits.domtbl, params.hmmsearch_query_length_threshold )
    ch_versions = ch_versions.mix( BRANCH_HITS_FASTA.out.versions )
    ch_no_hit_seqs = BRANCH_HITS_FASTA.out.non_hit_fasta

    ch_hits_fasta = BRANCH_HITS_FASTA.out.hits
        .transpose()
        .map { meta, file ->
            [[id: meta.id, family: file.getSimpleName()], file]
        }

    ch_family_msas = UNTAR_MSA.out.untar
        .map { meta, folder ->
            [meta, file("${folder.toUriString()}/*", checkIfExists: true)]
        }
        .transpose()
        .map { meta, file ->
            [[id: meta.id, family: file.getSimpleName()], file]
        }

    // Keep fasta with family sequences by removing gaps
    SEQKIT_SEQ( ch_family_msas )
    ch_versions = ch_versions.mix(SEQKIT_SEQ.out.versions)

    // Match newly recruited sequences with existing ones for each family
    ch_input_for_cat = SEQKIT_SEQ.out.fastx
        .combine(ch_hits_fasta, by: 0)
        .map { meta, family_fasta, new_fasta  ->
            [meta, [family_fasta, new_fasta]]
        }

    // Aggregate each family's MSA sequences with the newly recruited ones
    CAT_FASTA( ch_input_for_cat )
    ch_versions = ch_versions.mix( CAT_FASTA.out.versions )
    ch_fasta = CAT_FASTA.out.file_out

    if (params.remove_sequence_redundancy) {
        // Strict clustering to remove redundancy
        EXECUTE_CLUSTERING( ch_fasta, params.clustering_tool )
        ch_versions = ch_versions.mix( EXECUTE_CLUSTERING.out.versions )

        REMOVE_REDUNDANT_SEQS( EXECUTE_CLUSTERING.out.clusters, EXECUTE_CLUSTERING.out.seqs )
        ch_versions = ch_versions.mix( REMOVE_REDUNDANT_SEQS.out.versions )
        ch_fasta = REMOVE_REDUNDANT_SEQS.out.fasta
    }

    ALIGN_SEQUENCES( ch_fasta, params.alignment_tool )
    ch_versions = ch_versions.mix( ALIGN_SEQUENCES.out.versions )
    ch_msa = ALIGN_SEQUENCES.out.alignments

    if (params.trim_msa) {
        CLIPKIT( ch_msa, params.clipkit_out_format )
        ch_versions = ch_versions.mix( CLIPKIT.out.versions )
        ch_msa = CLIPKIT.out.clipkit
    }

    HMMER_HMMBUILD( ch_msa, [] )
    ch_versions = ch_versions.mix( HMMER_HMMBUILD.out.versions )

    ch_fasta = ch_fasta
        .map { meta, faa -> [ [id: meta.id], faa ] }
        .groupTuple(by: 0)
    EXTRACT_FAMILY_REPS( ch_fasta )
    ch_versions = ch_versions.mix( EXTRACT_FAMILY_REPS.out.versions )
    ch_updated_family_reps = ch_updated_family_reps.mix( EXTRACT_FAMILY_REPS.out.map )

    emit:
    versions            = ch_versions
    no_hit_seqs         = ch_no_hit_seqs
    updated_family_reps = ch_updated_family_reps
}
