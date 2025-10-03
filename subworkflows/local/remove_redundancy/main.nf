/*
    REMOVAL OF REDUNDANT SEQUENCES AND FAMILIES
*/

include { EXTRACT_FAMILY_REPS                                        } from '../../../modules/local/extract_family_reps/main'
include { FIND_CONCATENATE                                           } from '../../../modules/nf-core/find/concatenate'
include { HMMER_HMMSEARCH                                            } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { IDENTIFY_REDUNDANT_FAMS                                    } from '../../../modules/local/identify_redundant_fams/main'
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
    seed_msa                                     // tuple val(meta), path(clipkit)
    full_msa                                     // tuple val(meta), path(sto.gz)
    fasta                                        // tuple val(meta), path(fasta.gz)
    hmm                                          // tuple val(meta), path(hmm.gz)
    remove_family_redundancy                     // boolean
    hmmsearch_family_redundancy_length_threshold // number [0.0, 1.0]
    hmmsearch_family_similarity_length_threshold // number [0.0, 1.0]
    remove_sequence_redundancy                   // boolean
    clustering_tool                              // string ["linclust", "cluster"]
    alignment_tool                               // string ["famsa", "mafft"]

    main:
    ch_versions = Channel.empty()

    if (remove_family_redundancy) {
        ch_fasta = fasta
            .map { meta, faa -> [[id: meta.id], faa] }
            .groupTuple(by: 0)
        EXTRACT_FAMILY_REPS( ch_fasta )
        ch_versions = ch_versions.mix( EXTRACT_FAMILY_REPS.out.versions )

        ch_hmm = hmm
            .map { meta, model -> [[id: meta.id], model] }
            .groupTuple(by: 0)
        FIND_CONCATENATE( ch_hmm )
        ch_versions = ch_versions.mix( FIND_CONCATENATE.out.versions )

        ch_input_for_hmmsearch = FIND_CONCATENATE.out.file_out
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

        fasta = fasta
            .map { meta, fas -> [[id: meta.id], fas] }
            .groupTuple(by: 0)

        ch_seed_msa = seed_msa
            .map { meta, fas -> [[id: meta.id], fas] }
            .groupTuple(by: 0)

        full_msa = full_msa
            .map { meta, fas -> [[id: meta.id], fas] }
            .groupTuple(by: 0)

        // Join to ensure in sync
        ch_input_for_fam_removal = IDENTIFY_REDUNDANT_FAMS.out.redundant_ids
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
                [[id: meta.id, chunk: file.getSimpleName().split('_')[-1]], file]
            }

        FILTER_NON_REDUNDANT_FASTA( ch_input_for_fam_removal.seq, ch_input_for_fam_removal.ids  )
        ch_versions = ch_versions.mix( FILTER_NON_REDUNDANT_FASTA.out.versions )

        fasta = FILTER_NON_REDUNDANT_FASTA.out.filtered
            .transpose()
            .map { meta, file ->
                [[id: meta.id, chunk: file.getSimpleName().split('_')[-1]], file]
            }
    }

    if (remove_sequence_redundancy) {
        MMSEQS_FASTA_CLUSTER( fasta, clustering_tool )
        ch_versions = ch_versions.mix( MMSEQS_FASTA_CLUSTER.out.versions )

        REMOVE_REDUNDANT_SEQS( MMSEQS_FASTA_CLUSTER.out.clusters, MMSEQS_FASTA_CLUSTER.out.seqs )
        ch_versions = ch_versions.mix( REMOVE_REDUNDANT_SEQS.out.versions )
        fasta = REMOVE_REDUNDANT_SEQS.out.fasta

        ALIGN_SEQUENCES( REMOVE_REDUNDANT_SEQS.out.fasta, alignment_tool )
        ch_versions = ch_versions.mix( ALIGN_SEQUENCES.out.versions )
    } else if (remove_family_redundancy) {
        HHSUITE_REFORMAT_FILTERED( full_msa, "sto", "fas" )
        ch_versions = ch_versions.mix( HHSUITE_REFORMAT_FILTERED.out.versions )
    } else { // both remove_family_redundancy and remove_sequence_redundancy false, different publish dir
        HHSUITE_REFORMAT_RAW( full_msa, "sto", "fas" )
        ch_versions = ch_versions.mix( HHSUITE_REFORMAT_RAW.out.versions )
    }

    emit:
    versions = ch_versions
    fasta = fasta
}
