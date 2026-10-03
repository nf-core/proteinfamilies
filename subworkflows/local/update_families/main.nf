/*
    UPDATE EXISTING FAMILIES

    Assigns new sequences to existing families by searching them, together with the members of
    any existing full MSAs, against a concatenated library of the existing HMMs. Each family's
    hits are then rebuilt like a newly created family (GENERATE_FAMILIES): aligned and trimmed
    into a new seed MSA, built into a new HMM and, unless skipped, used to recruit the full MSA
    from the same sequence pool. With skip_update_refinement, the existing HMMs are kept and
    only align their hits into new full MSAs. Input sequences matching no family are emitted as no_hit_seqs
    for downstream de-novo family creation.
*/

include { UNTAR as UNTAR_HMM            } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_SEED_MSA       } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_FULL_MSA       } from '../../../modules/nf-core/untar/main'
include { validateHmmNames              } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { validateMsaStems              } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { FIND_CONCATENATE as CAT_HMM   } from '../../../modules/nf-core/find/concatenate/main'
include { GUNZIP                        } from '../../../modules/nf-core/gunzip/main'
include { POOL_EXISTING_MEMBERS         } from '../../../modules/local/pool_existing_members/main'
include { HMMER_HMMSEARCH               } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { BRANCH_HITS_FASTA             } from '../../../modules/local/branch_hits_fasta'
include { fileStem                      } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { MMSEQS_FASTA_CLUSTER          } from '../../../subworkflows/nf-core/mmseqs_fasta_cluster'
include { REMOVE_REDUNDANT_SEQS         } from '../../../modules/local/remove_redundant_seqs/main'
include { GENERATE_FAMILIES             } from '../../../subworkflows/local/generate_families'
include { HMMER_HMMALIGN                } from '../../../modules/nf-core/hmmer/hmmalign/main'
include { EXTRACT_FAMILY_MEMBERS        } from '../../../modules/local/extract_family_members/main'
include { EXTRACT_FAMILY_REPS           } from '../../../modules/local/extract_family_reps/main'

workflow UPDATE_FAMILIES {
    take:
    ch_samplesheet_for_update           // channel: [meta, sequences, existing_hmms, existing_seed_msas, existing_full_msas]; MSAs may be []
    hmmsearch_query_length_threshold    // number [0.0, 1.0]
    skip_sequence_redundancy_removal    // boolean
    clustering_tool                     // string ["linclust", "cluster"]
    alignment_tool                      // string ["famsa", "mafft"]
    skip_seed_msa_trimming              // boolean
    hmmsearch_write_target              // boolean
    hmmsearch_write_domain              // boolean
    skip_additional_sequence_recruiting // boolean
    skip_update_refinement              // boolean: keep the existing HMMs, only rebuild the full MSAs

    main:
    ch_updated_family_reps = channel.empty()

    ch_input_for_untar = ch_samplesheet_for_update
        .multiMap { meta, _fasta, existing_hmms, existing_seed_msas, existing_full_msas ->
            hmm: [ meta, existing_hmms ]
            seed_msa: [ meta, existing_seed_msas ]
            full_msa: [ meta, existing_full_msas ]
        }

    UNTAR_HMM( ch_input_for_untar.hmm )
    UNTAR_SEED_MSA( ch_input_for_untar.seed_msa.filter { _meta, archive -> archive } )
    UNTAR_FULL_MSA( ch_input_for_untar.full_msa.filter { _meta, archive -> archive } )

    // Families are matched by HMM NAME (hmmsearch) and by file stem (MSAs): both must agree
    validateHmmNames( UNTAR_HMM.out.untar )
    validateMsaStems( UNTAR_HMM.out.untar, UNTAR_SEED_MSA.out.untar )
    validateMsaStems( UNTAR_HMM.out.untar, UNTAR_FULL_MSA.out.untar )

    // Squeeze the HMMs into a single file
    CAT_HMM( UNTAR_HMM.out.untar.map { meta, folder -> [meta, file("${folder.toUriString()}/*", checkIfExists: true)] } )

    // The searched pool: the input sequences, plus the members of existing full MSAs, so that
    // the families keep the old members that still hit. HMMER rewinds the target database for
    // every query HMM and a gzip stream cannot rewind, so the pool is always uncompressed.
    ch_input_fasta = ch_samplesheet_for_update
        .map { meta, fasta, _existing_hmms, _existing_seed_msas, _existing_full_msas -> [meta, fasta] }

    POOL_EXISTING_MEMBERS( ch_input_fasta.join(UNTAR_FULL_MSA.out.untar) )

    ch_branched_sequences = ch_input_fasta
        .join(UNTAR_FULL_MSA.out.untar, remainder: true)
        .filter { _meta, _fasta, full_msas -> !full_msas }
        .map { meta, fasta, _full_msas -> [meta, fasta] }
        .branch { _meta, fasta ->
            compressed  : fasta.name.endsWith('.gz')
            uncompressed: true
        }

    GUNZIP( ch_branched_sequences.compressed )

    ch_pool = ch_branched_sequences.uncompressed
        .mix( GUNZIP.out.gunzip )
        .mix( POOL_EXISTING_MEMBERS.out.fasta )

    ch_input_for_hmmsearch = CAT_HMM.out.file_out
        .join(ch_pool)
        .map { meta, concatenated_hmm, pool -> [meta, concatenated_hmm, pool, false, false, true] }

    HMMER_HMMSEARCH( ch_input_for_hmmsearch )

    // Hits are cut from the pool, but only input sequences can be non-hits
    ch_input_for_branch_hits = HMMER_HMMSEARCH.out.domain_summary
        .join(ch_input_fasta)
        .join(POOL_EXISTING_MEMBERS.out.fasta, remainder: true)
        .multiMap { meta, domtbl, fasta, pool ->
            domtbl: [ meta, domtbl ]
            fasta: [ meta, fasta, pool ?: [] ]
        }

    // Branch hit families from input sequences without hits
    BRANCH_HITS_FASTA ( ch_input_for_branch_hits.fasta, ch_input_for_branch_hits.domtbl, hmmsearch_query_length_threshold )

    // Families without any hit are kept unchanged: their existing HMM passes through, and they
    // are listed per sample. A sample without any hit emits no hits at all.
    ch_zero_hit_hmms = UNTAR_HMM.out.untar
        .join(BRANCH_HITS_FASTA.out.hits, remainder: true)
        .map { meta, folder, hits ->
            def hit_families = [hits].flatten().findAll().collect { hit -> fileStem(hit) }
            [ meta, folder.listFiles().toList().findAll { hmm -> !(fileStem(hmm) in hit_families) } ]
        }

    // [id, family] meta, as for created families' chunks
    ch_hits = BRANCH_HITS_FASTA.out.hits
        .transpose()
        .map { meta, file ->
            [[id: meta.id, family: fileStem(file)], file]
        }

    if (skip_update_refinement) {
        // The existing HMMs stay as they are and align their own hits into the new full MSAs;
        // a provided seed MSA passes through unchanged
        ch_existing_hmms = UNTAR_HMM.out.untar
            .flatMap { meta, folder -> folder.listFiles().collect { hmm -> [[id: meta.id, family: fileStem(hmm)], hmm] } }

        ch_input_for_hmmalign = ch_hits
            .join(ch_existing_hmms)
            .multiMap { meta, seqs, hmm ->
                seq: [ meta, seqs ]
                hmm: [ hmm ]
            }

        HMMER_HMMALIGN( ch_input_for_hmmalign.seq, ch_input_for_hmmalign.hmm )

        ch_seed_msa = UNTAR_SEED_MSA.out.untar
            .flatMap { meta, folder -> folder.listFiles().collect { msa -> [[id: meta.id, family: fileStem(msa)], msa] } }
            .join(ch_hits)
            .map { meta, msa, _hits -> [ meta, msa ] }
        ch_full_msa = HMMER_HMMALIGN.out.sto
        ch_fasta    = ch_hits
        ch_hmm      = ch_existing_hmms // zero-hit families included
    } else {
        ch_fasta = ch_hits
        if (!skip_sequence_redundancy_removal) {
            // Strict clustering to remove redundancy
            MMSEQS_FASTA_CLUSTER( ch_fasta, clustering_tool )

            REMOVE_REDUNDANT_SEQS( MMSEQS_FASTA_CLUSTER.out.clusters, MMSEQS_FASTA_CLUSTER.out.seqs )
            ch_fasta = REMOVE_REDUNDANT_SEQS.out.fasta
        }

        // Rebuild each family like a created one: new seed MSA and HMM, full MSA recruited from the pool
        GENERATE_FAMILIES(
            ch_pool,
            ch_fasta,
            alignment_tool,
            skip_seed_msa_trimming,
            hmmsearch_write_target,
            hmmsearch_write_domain,
            skip_additional_sequence_recruiting,
            hmmsearch_query_length_threshold
        )

        ch_seed_msa = GENERATE_FAMILIES.out.seed_msa
        ch_full_msa = GENERATE_FAMILIES.out.full_msa
        ch_fasta    = GENERATE_FAMILIES.out.fasta
        ch_hmm      = GENERATE_FAMILIES.out.hmm
            .mix( ch_zero_hit_hmms.transpose().map { meta, hmm -> [[id: meta.id, family: fileStem(hmm)], hmm] } )
    }

    // Strip family from meta and group by sample ID so EXTRACT_FAMILY_MEMBERS/REPS
    // receive all families for a sample together.
    ch_fasta_per_sample = ch_fasta
        .map { meta, faa -> [ [id: meta.id], faa ] }
        .groupTuple(by: 0)

    EXTRACT_FAMILY_MEMBERS( ch_fasta_per_sample )

    EXTRACT_FAMILY_REPS( ch_fasta_per_sample )
    ch_updated_family_reps = ch_updated_family_reps.mix( EXTRACT_FAMILY_REPS.out.map )

    emit:
    seed_msa            = ch_seed_msa
    full_msa            = ch_full_msa
    fasta               = ch_fasta
    hmm                 = ch_hmm
    zero_hit_families   = ch_zero_hit_hmms.map { meta, hmms -> [meta, hmms.collect { hmm -> fileStem(hmm) }.sort()] }
    no_hit_seqs         = BRANCH_HITS_FASTA.out.non_hit_fasta
    updated_family_reps = ch_updated_family_reps
}
