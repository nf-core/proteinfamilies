/*
    UPDATE EXISTING FAMILIES HMM AND MSA

    Assigns new sequences to existing families by searching against a concatenated HMM
    library. Hit sequences are merged with each family's existing members, re-aligned,
    and used to rebuild the family HMM. Sequences matching no family are emitted as
    no_hit_seqs for downstream de-novo family creation.
*/

include { UNTAR as UNTAR_HMM            } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_MSA            } from '../../../modules/nf-core/untar/main'
include { validateMatchingFolders       } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { FIND_CONCATENATE as CAT_HMM   } from '../../../modules/nf-core/find/concatenate/main'
include { GUNZIP                        } from '../../../modules/nf-core/gunzip/main'
include { HMMER_HMMSEARCH               } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { BRANCH_HITS_FASTA             } from '../../../modules/local/branch_hits_fasta'
include { fileStem                      } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { SEQKIT_SEQ                    } from '../../../modules/nf-core/seqkit/seq/main'
include { FIND_CONCATENATE as CAT_FASTA } from '../../../modules/nf-core/find/concatenate/main'
include { MMSEQS_FASTA_CLUSTER          } from '../../../subworkflows/nf-core/mmseqs_fasta_cluster'
include { REMOVE_REDUNDANT_SEQS         } from '../../../modules/local/remove_redundant_seqs/main'
include { ALIGN_SEQUENCES               } from '../../../subworkflows/local/align_sequences'
include { HMMER_HMMBUILD                } from '../../../modules/nf-core/hmmer/hmmbuild/main'
include { EXTRACT_FAMILY_MEMBERS        } from '../../../modules/local/extract_family_members/main'
include { EXTRACT_FAMILY_REPS           } from '../../../modules/local/extract_family_reps/main'

workflow UPDATE_FAMILIES {
    take:
    ch_samplesheet_for_update                   // channel: [meta, sequences, existing_hmms_to_update, existing_msas_to_update]
    hmmsearch_query_length_threshold            // number [0.0, 1.0]
    skip_sequence_redundancy_removal            // boolean
    clustering_tool                             // string ["linclust", "cluster"]
    alignment_tool                              // string ["famsa", "mafft"]
    skip_seed_msa_trimming                      // boolean

    main:
    ch_updated_family_reps = channel.empty()
    ch_no_hit_seqs         = channel.empty()

    ch_input_for_untar = ch_samplesheet_for_update
        .multiMap { meta, _fasta, existing_hmms_to_update, existing_msas_to_update ->
            hmm: [ meta, existing_hmms_to_update ]
            msa: [ meta, existing_msas_to_update ]
        }

    UNTAR_HMM( ch_input_for_untar.hmm )

    UNTAR_MSA( ch_input_for_untar.msa )

    // Validate that each HMM archive and MSA archive contain matching files. A mismatch
    // would cause silent per-family key-join failures in the combine steps below.
    // join ensures HMM/MSA tarballs are processed in sync per sample.
    ch_folders_to_validate = UNTAR_HMM.out.untar
        .join(UNTAR_MSA.out.untar)
        .multiMap { meta, folder1, folder2 ->
            ch_hmm_folder: [meta, folder1]
            ch_msa_folder: [meta, folder2]
        }
    validateMatchingFolders(ch_folders_to_validate.ch_hmm_folder, ch_folders_to_validate.ch_msa_folder)

    // Squeeze the HMMs into a single file
    CAT_HMM( UNTAR_HMM.out.untar.map { meta, folder -> [meta, file("${folder.toUriString()}/*", checkIfExists: true)] } )

    // HMMER rewinds the target database for every query HMM; a gzip stream cannot rewind, so
    // only the first HMM of the library would be searched against a gzipped FASTA.
    ch_branched_sequences = ch_samplesheet_for_update
        .map { meta, fasta, _existing_hmms_to_update, _existing_msas_to_update -> [meta, fasta] }
        .branch { _meta, fasta ->
            compressed  : fasta.name.endsWith('.gz')
            uncompressed: true
        }

    GUNZIP( ch_branched_sequences.compressed )

    ch_sequences = ch_branched_sequences.uncompressed.mix( GUNZIP.out.gunzip )

    // Prep the sequences to search against the HMM concatenated model of families
    ch_input_for_hmmsearch = CAT_HMM.out.file_out
        .join(ch_sequences)
        .map { meta, concatenated_hmm, fasta -> [meta, concatenated_hmm, fasta, false, false, true] }

    HMMER_HMMSEARCH( ch_input_for_hmmsearch )

    ch_input_for_branch_hits = HMMER_HMMSEARCH.out.domain_summary
        .join(ch_samplesheet_for_update)
        .multiMap { meta, domtbl, fasta, _existing_hmms_to_update, _existing_msas_to_update ->
            domtbl: [ meta, domtbl ]
            fasta: [ meta, fasta ]
        }

    // Branch hit families/fasta proteins from non hit fasta proteins
    BRANCH_HITS_FASTA ( ch_input_for_branch_hits.fasta, ch_input_for_branch_hits.domtbl, hmmsearch_query_length_threshold )
    ch_no_hit_seqs = BRANCH_HITS_FASTA.out.non_hit_fasta

    // Both channels use [id, family] meta so they can be combined by key to pair each
    // newly recruited sequence file with its corresponding family's MSA.
    ch_hits_fasta = BRANCH_HITS_FASTA.out.hits
        .transpose()
        .map { meta, file ->
            [[id: meta.id, family: fileStem(file)], file]
        }

    ch_family_msas = UNTAR_MSA.out.untar
        .map { meta, folder ->
            [meta, file("${folder.toUriString()}/*", checkIfExists: true)]
        }
        .transpose()
        .map { meta, file ->
            [[id: meta.id, family: fileStem(file)], file]
        }

    // Keep fasta with family sequences by removing gaps
    SEQKIT_SEQ( ch_family_msas )

    // Match newly recruited sequences with existing ones for each family
    ch_input_for_cat = SEQKIT_SEQ.out.fastx
        .combine(ch_hits_fasta, by: 0)
        .map { meta, family_fasta, new_fasta  ->
            [meta, [family_fasta, new_fasta]]
        }

    // Aggregate each family's MSA sequences with the newly recruited ones
    CAT_FASTA( ch_input_for_cat )
    ch_fasta = CAT_FASTA.out.file_out

    if (!skip_sequence_redundancy_removal) {
        // Strict clustering to remove redundancy
        MMSEQS_FASTA_CLUSTER( ch_fasta, clustering_tool )

        REMOVE_REDUNDANT_SEQS( MMSEQS_FASTA_CLUSTER.out.clusters, MMSEQS_FASTA_CLUSTER.out.seqs )
        ch_fasta = REMOVE_REDUNDANT_SEQS.out.fasta
    }

    // The (trimmed) alignment builds the HMM and is the family's MSA; its rows are the family fasta
    ALIGN_SEQUENCES( ch_fasta, alignment_tool, skip_seed_msa_trimming )
    ch_fasta = ALIGN_SEQUENCES.out.sequences

    HMMER_HMMBUILD( ALIGN_SEQUENCES.out.alignments, [] )

    // Strip family from meta and group by sample ID so EXTRACT_FAMILY_MEMBERS/REPS
    // receive all families for a sample together.
    ch_fasta = ch_fasta
        .map { meta, faa -> [ [id: meta.id], faa ] }
        .groupTuple(by: 0)

    EXTRACT_FAMILY_MEMBERS( ch_fasta )

    EXTRACT_FAMILY_REPS( ch_fasta )
    ch_updated_family_reps = ch_updated_family_reps.mix( EXTRACT_FAMILY_REPS.out.map )

    emit:
    no_hit_seqs         = ch_no_hit_seqs
    updated_family_reps = ch_updated_family_reps
    hmm                 = HMMER_HMMBUILD.out.hmm
}
