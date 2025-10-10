/*
    REMOVAL OF REDUNDANT SEQUENCES AND FAMILIES
*/

include { EXTRACT_FAMILY_REPS                                        } from '../../../modules/local/extract_family_reps/main'
include { FIND_CONCATENATE as FIND_CONCATENATE_HMMS                  } from '../../../modules/nf-core/find/concatenate'
include { HMMER_HMMSEARCH                                            } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { IDENTIFY_REDUNDANT_FAMS                                    } from '../../../modules/local/identify_redundant_fams/main'
include { MERGE_FAMILIES                                             } from '../../../subworkflows/local/merge_families/main'
include { FIND_CONCATENATE as FIND_CONCATENATE_SKIP_IDS              } from '../../../modules/nf-core/find/concatenate'
include { FILTER_NON_REDUNDANT_FAMS as FILTER_NON_REDUNDANT_HMM      } from '../../../modules/local/filter_non_redundant_fams/main'
include { FILTER_NON_REDUNDANT_FAMS as FILTER_NON_REDUNDANT_SEED_MSA } from '../../../modules/local/filter_non_redundant_fams/main'
include { FILTER_NON_REDUNDANT_FAMS as FILTER_NON_REDUNDANT_FULL_MSA } from '../../../modules/local/filter_non_redundant_fams/main'
include { FILTER_NON_REDUNDANT_FAMS as FILTER_NON_REDUNDANT_FASTA    } from '../../../modules/local/filter_non_redundant_fams/main'
include { MMSEQS_FASTA_CLUSTER                                       } from '../../../subworkflows/nf-core/mmseqs_fasta_cluster'
include { REMOVE_REDUNDANT_SEQS                                      } from '../../../modules/local/remove_redundant_seqs/main'
include { ALIGN_SEQUENCES                                            } from '../../../subworkflows/local/align_sequences'
include { HHSUITE_REFORMAT as HHSUITE_REFORMAT_FILTERED              } from '../../../modules/nf-core/hhsuite/reformat/main'
include { HHSUITE_REFORMAT as HHSUITE_REFORMAT_RAW                   } from '../../../modules/nf-core/hhsuite/reformat/main'

workflow REMOVE_REDUNDANCY {
    take:
    sequences                                    // tuple val(meta), path(fasta)
    seed_msa                                     // tuple val(meta), path(clipkit)
    full_msa                                     // tuple val(meta), path(sto.gz)
    fasta                                        // tuple val(meta), path(fasta.gz)
    hmm                                          // tuple val(meta), path(hmm.gz)
    skip_family_redundancy_removal               // boolean
    skip_family_merging                          // boolean
    hmmsearch_family_redundancy_length_threshold // number [0.0, 1.0]
    hmmsearch_family_similarity_length_threshold // number [0.0, 1.0]
    skip_sequence_redundancy_removal             // boolean
    clustering_tool                              // string ["linclust", "cluster"]
    alignment_tool                               // string ["famsa", "mafft"]
    skip_msa_trimming                            // boolean
    clipkit_out_format                           // string (default: clipkit)
    hmmsearch_write_target                       // boolean
    hmmsearch_write_domain                       // boolean
    skip_additional_sequence_recruiting          // boolean
    hmmsearch_query_length_threshold             // number [0.0, 1.0]

    main:
    ch_versions        = Channel.empty()
    ch_merged_seed_msa = Channel.empty()
    ch_merged_full_msa = Channel.empty()
    ch_merged_fasta    = Channel.empty()
    ch_merged_hmm      = Channel.empty()

    if (!skip_family_redundancy_removal || !skip_family_merging) {
        ch_fasta = fasta
            .map { meta, faa -> [[id: meta.id], faa] }
            .groupTuple(by: 0)
        EXTRACT_FAMILY_REPS( ch_fasta )
        ch_versions = ch_versions.mix( EXTRACT_FAMILY_REPS.out.versions )

        ch_hmm = hmm
            .map { meta, model -> [[id: meta.id], model] }
            .groupTuple(by: 0)
        FIND_CONCATENATE_HMMS( ch_hmm )
        ch_versions = ch_versions.mix( FIND_CONCATENATE_HMMS.out.versions )

        ch_input_for_hmmsearch = FIND_CONCATENATE_HMMS.out.file_out
            .combine(EXTRACT_FAMILY_REPS.out.fasta, by: 0)
            .map { meta, model, seqs -> [meta, model, seqs, false, false, true] }

        HMMER_HMMSEARCH( ch_input_for_hmmsearch )
        ch_versions = ch_versions.mix( HMMER_HMMSEARCH.out.versions )

        // Join to ensure in sync
        ch_input_for_redundant_fam_identification = EXTRACT_FAMILY_REPS.out.map
            .join(HMMER_HMMSEARCH.out.domain_summary)
            .multiMap { meta, map, domtbl ->
                map: [meta, map]
                domtbl: [meta, domtbl]
            }

        IDENTIFY_REDUNDANT_FAMS( ch_input_for_redundant_fam_identification.map, \
            ch_input_for_redundant_fam_identification.domtbl, \
            hmmsearch_family_redundancy_length_threshold, hmmsearch_family_similarity_length_threshold )
        ch_versions = ch_versions.mix( IDENTIFY_REDUNDANT_FAMS.out.versions )

        ch_seed_msa = seed_msa
            .map { meta, fas -> [[id: meta.id], fas] }
            .groupTuple(by: 0)

        if (!skip_family_merging) {
            MERGE_FAMILIES( IDENTIFY_REDUNDANT_FAMS.out.similarities, ch_seed_msa, \
            sequences, alignment_tool, skip_msa_trimming, clipkit_out_format, \
            hmmsearch_write_target, hmmsearch_write_domain, \
            skip_additional_sequence_recruiting, hmmsearch_query_length_threshold )
            ch_versions = ch_versions.mix( MERGE_FAMILIES.out.versions )

            ch_merged_seed_msa = MERGE_FAMILIES.out.seed_msa
            ch_merged_full_msa = MERGE_FAMILIES.out.full_msa
            ch_merged_fasta    = MERGE_FAMILIES.out.fasta
            ch_merged_hmm      = MERGE_FAMILIES.out.hmm

            ch_seed_msa = seed_msa
                .mix(ch_merged_seed_msa)
                .map { meta, fas -> [[id: meta.id], fas] }
                .groupTuple(by: 0)

            ch_hmm = hmm
                .mix(ch_merged_hmm)
                .map { meta, model -> [[id: meta.id], model] }
                .groupTuple(by: 0)
        }

        fasta = fasta
            .mix(ch_merged_fasta)
            .map { meta, fas -> [[id: meta.id], fas] }
            .groupTuple(by: 0)

        full_msa = full_msa
            .mix(ch_merged_full_msa)
            .map { meta, fas -> [[id: meta.id], fas] }
            .groupTuple(by: 0)

        ch_skip_ids = IDENTIFY_REDUNDANT_FAMS.out.redundant_ids
            .concat(IDENTIFY_REDUNDANT_FAMS.out.similar_ids)
            .groupTuple(by: 0)
        FIND_CONCATENATE_SKIP_IDS( ch_skip_ids )
        ch_versions = ch_versions.mix( FIND_CONCATENATE_SKIP_IDS.out.versions )

        // Join to ensure in sync
        ch_input_for_fam_removal = FIND_CONCATENATE_SKIP_IDS.out.file_out
            .join(fasta)
            .join(ch_hmm)
            .join(ch_seed_msa)
            .join(full_msa)
            .multiMap { meta, ids, seq, model, seed, full ->
                ids: [meta, ids]
                seq: [meta, seq]
                model: [meta, model]
                seed: [meta, seed]
                full: [meta, full]
            }

        FILTER_NON_REDUNDANT_HMM( ch_input_for_fam_removal.model, ch_input_for_fam_removal.ids )
        ch_versions = ch_versions.mix( FILTER_NON_REDUNDANT_HMM.out.versions )

        FILTER_NON_REDUNDANT_SEED_MSA( ch_input_for_fam_removal.seed, ch_input_for_fam_removal.ids )
        ch_versions = ch_versions.mix( FILTER_NON_REDUNDANT_SEED_MSA.out.versions )

        FILTER_NON_REDUNDANT_FULL_MSA( ch_input_for_fam_removal.full, ch_input_for_fam_removal.ids )
        ch_versions = ch_versions.mix( FILTER_NON_REDUNDANT_FULL_MSA.out.versions )

        full_msa = FILTER_NON_REDUNDANT_FULL_MSA.out.filtered
            .transpose()
            .map { meta, file ->
                def filename = file.getSimpleName()
                def chunk = filename.split("${meta.id}_", 2)[1]  // Split by meta.id_ and take remainder, to also match merged ids
                [[id: meta.id, chunk: chunk], file]
            }

        FILTER_NON_REDUNDANT_FASTA( ch_input_for_fam_removal.seq, ch_input_for_fam_removal.ids  )
        ch_versions = ch_versions.mix( FILTER_NON_REDUNDANT_FASTA.out.versions )

        fasta = FILTER_NON_REDUNDANT_FASTA.out.filtered
            .transpose()
            .map { meta, file ->
                def filename = file.getSimpleName()
                def chunk = filename.split("${meta.id}_", 2)[1]  // Split by meta.id_ and take remainder, to also match merged ids
                [[id: meta.id, chunk: chunk], file]
            }
    }

    if (!skip_sequence_redundancy_removal) {
        MMSEQS_FASTA_CLUSTER( fasta, clustering_tool )
        ch_versions = ch_versions.mix( MMSEQS_FASTA_CLUSTER.out.versions )

        REMOVE_REDUNDANT_SEQS( MMSEQS_FASTA_CLUSTER.out.clusters, MMSEQS_FASTA_CLUSTER.out.seqs )
        ch_versions = ch_versions.mix( REMOVE_REDUNDANT_SEQS.out.versions )
        fasta = REMOVE_REDUNDANT_SEQS.out.fasta

        ALIGN_SEQUENCES( REMOVE_REDUNDANT_SEQS.out.fasta, alignment_tool )
        ch_versions = ch_versions.mix( ALIGN_SEQUENCES.out.versions )
    } else if (!skip_family_redundancy_removal) {
        HHSUITE_REFORMAT_FILTERED( full_msa, "sto", "fas" )
        ch_versions = ch_versions.mix( HHSUITE_REFORMAT_FILTERED.out.versions )
    } else if (!skip_additional_sequence_recruiting) {
        // if both skip_family_redundancy_removal and skip_sequence_redundancy_removal true, 
        // and skip_additional_sequence_recruiting also true (sequences not recruited, so not hmmalign .sto)
        HHSUITE_REFORMAT_RAW( full_msa, "sto", "fas" )
        ch_versions = ch_versions.mix( HHSUITE_REFORMAT_RAW.out.versions )
    }

    emit:
    versions = ch_versions
    fasta = fasta
}
