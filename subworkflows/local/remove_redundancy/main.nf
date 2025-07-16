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
include { EXECUTE_CLUSTERING                                         } from '../../../subworkflows/local/execute_clustering'
include { REMOVE_REDUNDANT_SEQS                                      } from '../../../modules/local/remove_redundant_seqs/main'
include { ALIGN_SEQUENCES                                            } from '../../../subworkflows/local/align_sequences'
include { HHSUITE_REFORMAT as HHSUITE_REFORMAT_FILTERED              } from '../../../modules/nf-core/hhsuite/reformat/main'
include { HHSUITE_REFORMAT as HHSUITE_REFORMAT_RAW                   } from '../../../modules/nf-core/hhsuite/reformat/main'

workflow REMOVE_REDUNDANCY {
    take:
    seed_msa // tuple val(meta), path(fas)
    full_msa // tuple val(meta), path(fas)
    fasta    // tuple val(meta), path(fasta)
    hmm      // tuple val(meta), path(hmm)

    main:
    ch_versions = Channel.empty()

    if (params.remove_family_redundancy) {
        ch_msa = full_msa
            .map { meta, aln -> [[id: meta.id], aln] }
            .groupTuple(by: 0)
        EXTRACT_FAMILY_REPS( ch_msa )
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
            ch_input_for_redundant_fam_identification.domtbl, params.hmmsearch_family_length_threshold )
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

    if (params.remove_sequence_redundancy) {
        EXECUTE_CLUSTERING( fasta, params.clustering_tool )
        ch_versions = ch_versions.mix( EXECUTE_CLUSTERING.out.versions )

        REMOVE_REDUNDANT_SEQS( EXECUTE_CLUSTERING.out.clusters, EXECUTE_CLUSTERING.out.seqs )
        ch_versions = ch_versions.mix( REMOVE_REDUNDANT_SEQS.out.versions )

        ALIGN_SEQUENCES( REMOVE_REDUNDANT_SEQS.out.fasta, params.alignment_tool )
        ch_versions = ch_versions.mix( ALIGN_SEQUENCES.out.versions )
        full_msa = ALIGN_SEQUENCES.out.alignments
    } else if (params.remove_family_redundancy) {
        HHSUITE_REFORMAT_FILTERED( full_msa, "sto", "fas" )
        ch_versions = ch_versions.mix( HHSUITE_REFORMAT_FILTERED.out.versions )
        full_msa = HHSUITE_REFORMAT_FILTERED.out.msa
    } else { // both remove_family_redundancy and remove_sequence_redundancy false, different publish dir
        HHSUITE_REFORMAT_RAW( full_msa, "sto", "fas" )
        ch_versions = ch_versions.mix( HHSUITE_REFORMAT_RAW.out.versions )
        full_msa = HHSUITE_REFORMAT_RAW.out.msa
    }

    emit:
    versions = ch_versions
    full_msa = full_msa
}
