/*
    REMOVAL OF REDUNDANT SEQUENCES AND FAMILIES

    Two independent redundancy stages applied in order:
    1. Family-level: concatenates all HMMs, searches family representatives against them,
       merges similar families (if enabled), and removes fully redundant ones.
    2. Sequence-level: clusters all sequences within each family and removes duplicates,
       then re-aligns the remaining members.
    Either or both stages can be skipped via parameters.

    Created and updated families enter together, keyed [id, family] by file stem. Updated
    families (existing NAMEs, never `<id>_<digit>...`) are never dropped: a created family
    redundant with an updated one is. Families without a seed MSA cannot merge.
*/

include { EXTRACT_FAMILY_REPS                                        } from '../../../modules/local/extract_family_reps/main'
include { FIND_CONCATENATE as FIND_CONCATENATE_HMMS                  } from '../../../modules/nf-core/find/concatenate'
include { HMMER_HMMSEARCH                                            } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { IDENTIFY_REDUNDANT_FAMS                                    } from '../../../modules/local/identify_redundant_fams/main'
include { MERGE_FAMILIES                                             } from '../../../subworkflows/local/merge_families/main'
include { FIND_CONCATENATE as FIND_CONCATENATE_SKIP_IDS              } from '../../../modules/nf-core/find/concatenate'
include { fileStem                                                   } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { isCreatedFamily                                            } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
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
    sequences                                    // tuple val(meta), path(faa): the sequences families were created from
    update_pool                                  // tuple val(meta), path(fasta): update samples' searched pool
    seed_msa                                     // tuple val(meta), path({aln,fas}), meta [id, family]
    full_msa                                     // tuple val(meta), path({sto.gz,aln,fas}), meta [id, family]
    fasta                                        // tuple val(meta), path(faa.gz), meta [id, family]
    hmm                                          // tuple val(meta), path(hmm.gz), meta [id, family]
    skip_family_redundancy_removal               // boolean
    skip_family_merging                          // boolean
    hmmsearch_family_redundancy_length_threshold // number [0.0, 1.0]
    hmmsearch_family_similarity_length_threshold // number [0.0, 1.0]
    skip_sequence_redundancy_removal             // boolean
    clustering_tool                              // string ["linclust", "cluster"]
    family_generation_algorithm                  // string ["standard", "iterative"]
    alignment_tool                               // string ["famsa", "mafft"]
    skip_seed_msa_trimming                       // boolean
    hmmsearch_write_target                       // boolean
    hmmsearch_write_domain                       // boolean
    skip_additional_sequence_recruiting          // boolean
    hmmsearch_query_length_threshold             // number [0.0, 1.0]
    merged_family_name                           // string ["existing", "new"]

    main:
    ch_merged_seed_msa = channel.empty()
    ch_merged_full_msa = channel.empty()
    ch_merged_fasta    = channel.empty()
    ch_merged_hmm      = channel.empty()
    ch_merged_families = channel.empty()
    ch_output_hmm      = channel.empty()

    // FAMILY REDUNDANCY REMOVAL MECHANISM
    // Block runs if either feature is enabled — both share the same HMM-search infrastructure.
    if (!skip_family_redundancy_removal || !skip_family_merging) {
        ch_fasta    = perSample(fasta)
        ch_hmm      = perSample(hmm)
        ch_seed_msa = perSample(seed_msa)
        ch_full_msa = perSample(full_msa)

        EXTRACT_FAMILY_REPS( ch_fasta )

        FIND_CONCATENATE_HMMS( ch_hmm )

        ch_input_for_hmmsearch = FIND_CONCATENATE_HMMS.out.file_out
            .combine(EXTRACT_FAMILY_REPS.out.fasta, by: 0)
            .map { meta, model, seqs -> [meta, model, seqs, false, false, true] }

        HMMER_HMMSEARCH( ch_input_for_hmmsearch )

        // Per sample: the updated families, and the families without a seed MSA to merge from
        ch_family_roles = hmm
            .map { meta, model -> [[id: meta.id], fileStem(model)] }
            .groupTuple()
            .join(seed_msa.map { meta, seed -> [[id: meta.id], fileStem(seed)] }.groupTuple(), remainder: true)
            .map { meta, families, seeded ->
                [meta, families.findAll { family -> !isCreatedFamily(meta.id, family) }.sort(), (families - (seeded ?: [])).sort()]
            }

        // Join to ensure in sync
        ch_input_for_redundant_fam_identification = EXTRACT_FAMILY_REPS.out.map
            .join(HMMER_HMMSEARCH.out.domain_summary)
            .join(ch_family_roles)
            .multiMap { meta, map, domtbl, updated, seedless ->
                map: [meta, map]
                domtbl: [meta, domtbl]
                roles: [meta, updated, seedless]
            }

        IDENTIFY_REDUNDANT_FAMS (
            ch_input_for_redundant_fam_identification.map,
            ch_input_for_redundant_fam_identification.domtbl,
            ch_input_for_redundant_fam_identification.roles,
            hmmsearch_family_redundancy_length_threshold,
            hmmsearch_family_similarity_length_threshold
        )

        if (!skip_family_merging) {
            // A merge recruits from the sequences its families were built from: created families
            // from `sequences`, updated ones from their sample's update pool
            ch_merge_sequences = sequences
                .map { meta, faa -> [[id: meta.id, pool: 'create'], faa] }
                .mix(update_pool.map { meta, faa -> [[id: meta.id, pool: 'update'], faa] })

            MERGE_FAMILIES (
                IDENTIFY_REDUNDANT_FAMS.out.similarities,
                ch_seed_msa,
                ch_merge_sequences,
                family_generation_algorithm,
                alignment_tool,
                skip_seed_msa_trimming,
                hmmsearch_write_target,
                hmmsearch_write_domain,
                skip_additional_sequence_recruiting,
                hmmsearch_query_length_threshold,
                merged_family_name
            )

            ch_merged_seed_msa = MERGE_FAMILIES.out.seed_msa
            ch_merged_full_msa = MERGE_FAMILIES.out.full_msa
            ch_merged_fasta    = MERGE_FAMILIES.out.fasta
            ch_merged_hmm      = MERGE_FAMILIES.out.hmm
            ch_merged_families = MERGE_FAMILIES.out.merged_families
        }

        // if --skip_family_redundancy_removal true, redundant_ids is returned empty by the script
        ch_skip_ids = IDENTIFY_REDUNDANT_FAMS.out.redundant_ids
        // will only remove similar families (e.g., _1 and _7) if merging them (i.e., will keep _1_7)
        if (!skip_family_merging) {
            ch_skip_ids = ch_skip_ids.concat( IDENTIFY_REDUNDANT_FAMS.out.similar_ids )
        }
        ch_skip_ids = ch_skip_ids.groupTuple(by: 0)

        FIND_CONCATENATE_SKIP_IDS( ch_skip_ids )

        // Join to ensure in sync. Merged families are kept as they are, next to the filtered
        // originals: a merge may take the name of an updated family it replaces. A sample whose
        // families have no seed MSA (updated without refinement) has no seed MSAs to filter.
        ch_input_for_fam_removal = FIND_CONCATENATE_SKIP_IDS.out.file_out
            .join(ch_fasta)
            .join(ch_hmm)
            .join(ch_full_msa)
            .join(ch_seed_msa, remainder: true)
            .join(perSample(ch_merged_fasta), remainder: true)
            .join(perSample(ch_merged_hmm), remainder: true)
            .join(perSample(ch_merged_seed_msa), remainder: true)
            .join(perSample(ch_merged_full_msa), remainder: true)
            .multiMap { meta, ids, seq, model, full, seed, merged_seq, merged_model, merged_seed, merged_full ->
                ids: [meta, ids]
                seq: [meta, seq, merged_seq ?: []]
                model: [meta, model, merged_model ?: []]
                seed: [meta, seed ?: [], merged_seed ?: []]
                full: [meta, full, merged_full ?: []]
            }

        FILTER_NON_REDUNDANT_HMM( ch_input_for_fam_removal.model, ch_input_for_fam_removal.ids )
        ch_output_hmm = FILTER_NON_REDUNDANT_HMM.out.filtered
            .transpose()   // unpack [meta, [f1,f2,...]] → individual [meta, file] tuples

        ch_seeds_for_removal = ch_input_for_fam_removal.seed
            .join(ch_input_for_fam_removal.ids)
            .filter { _meta, seed, merged_seed, _ids -> seed || merged_seed }
            .multiMap { meta, seed, merged_seed, ids ->
                seed: [meta, seed, merged_seed]
                ids: [meta, ids]
            }
        FILTER_NON_REDUNDANT_SEED_MSA( ch_seeds_for_removal.seed, ch_seeds_for_removal.ids )
        seed_msa = FILTER_NON_REDUNDANT_SEED_MSA.out.filtered
            .transpose()

        FILTER_NON_REDUNDANT_FULL_MSA( ch_input_for_fam_removal.full, ch_input_for_fam_removal.ids )

        full_msa = FILTER_NON_REDUNDANT_FULL_MSA.out.filtered
            .transpose()
            .map { meta, file -> [[id: meta.id, family: fileStem(file)], file] }

        FILTER_NON_REDUNDANT_FASTA( ch_input_for_fam_removal.seq, ch_input_for_fam_removal.ids  )

        fasta = FILTER_NON_REDUNDANT_FASTA.out.filtered
            .transpose()
            .map { meta, file -> [[id: meta.id, family: fileStem(file)], file] }
    } else {
        ch_output_hmm = hmm  // raw individual [meta(id,family), file] tuples
    }
    // END FAMILY REDUNDANCY REMOVAL MECHANISM

    if (!skip_sequence_redundancy_removal) {
        // SEQUENCE REDUNDANCY REMOVAL MECHANISM
        MMSEQS_FASTA_CLUSTER( fasta, clustering_tool ) // fasta channel contains all sequences of full MSA

        REMOVE_REDUNDANT_SEQS( MMSEQS_FASTA_CLUSTER.out.clusters, MMSEQS_FASTA_CLUSTER.out.seqs )
        fasta = REMOVE_REDUNDANT_SEQS.out.fasta

        // Full MSAs are never trimmed, so the fasta keeps matching them
        full_msa = ALIGN_SEQUENCES( REMOVE_REDUNDANT_SEQS.out.fasta, alignment_tool, true ).alignments
        // END SEQUENCE REDUNDANCY REMOVAL MECHANISM
    } else {
        // REFORMATTING FULL MSA
        // Recruited full MSAs (hmmalign, Stockholm) become aligned FASTA; seed MSAs reused as full
        // MSAs (no recruiting) already are. Two module aliases are required because Nextflow
        // prevents calling the same import more than once in a workflow: HHSUITE_REFORMAT_FILTERED
        // for full MSAs after filtering/merging, HHSUITE_REFORMAT_RAW when neither ran.
        ch_full_msa_format = full_msa.branch { _meta, msa ->
            stockholm: msa.name.endsWith('.sto.gz')
            fasta: true
        }
        if (!skip_family_redundancy_removal || !skip_family_merging) {
            ch_reformatted = HHSUITE_REFORMAT_FILTERED( ch_full_msa_format.stockholm, "sto", "fas" ).msa
        } else { // did not go through filtering processes
            ch_reformatted = HHSUITE_REFORMAT_RAW( ch_full_msa_format.stockholm, "sto", "fas" ).msa
        }
        full_msa = ch_full_msa_format.fasta.mix(ch_reformatted)
        // END REFORMATTING FULL MSA
    }

    emit:
    seed_msa        = seed_msa
    fasta           = fasta
    full_msa        = full_msa
    hmm             = ch_output_hmm
    merged_families = ch_merged_families // [meta, [[merged_id, 'member,...'], ...]], samples with merges only
}

// One [[id], [files]] per sample
def perSample(ch_files) {
    ch_files
        .map { meta, file -> [[id: meta.id], file] }
        .groupTuple(by: 0)
}
