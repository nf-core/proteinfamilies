/*
    UPDATE EXISTING FAMILIES

    PREPARE (engine-agnostic): untar and validate the existing families, and pool the input
    sequences with the members of any existing full MSAs.

    ENGINE (standard, below; mgnifam update_families from v3.1.0): the existing HMMs search the
    pool. Each family's hits are rebuilt like a newly created family (GENERATE_FAMILIES): aligned
    and trimmed into a new seed MSA, built into a new HMM and, unless skipped, used to recruit the
    full MSA from the pool. With skip_update_refinement, the existing HMMs only align their hits
    into new full MSAs. An engine returns a family by emitting its full MSA, plus a new HMM, seed
    MSA and family FASTA where it built them. Pooled full-MSA members without hits are dropped,
    they never go to de-novo family creation.

    FINALISE (engine-agnostic): returned families are emitted with the engine's files, plus the
    provided seed and HMM where the engine kept them, for redundancy removal next to created
    families. Every other family passes through as given (its full MSA's degapped members as the
    family FASTA) and is listed with a reason. Only input FASTA sequences (never pooled full-MSA
    members) that no returned family holds are emitted as no_hit_seqs for de-novo creation.
*/

include { SPLIT_HMMS                     } from '../../../modules/local/split_hmms/main'
include { UNTAR as UNTAR_SEED_MSA       } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_FULL_MSA       } from '../../../modules/nf-core/untar/main'
include { validateHmmNames              } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { validateMsaStems              } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { FIND_CONCATENATE as CAT_HMM   } from '../../../modules/nf-core/find/concatenate/main'
include { GUNZIP                        } from '../../../modules/nf-core/gunzip/main'
include { POOL_EXISTING_MEMBERS         } from '../../../modules/local/pool_existing_members/main'
include { HMMER_HMMSEARCH               } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { SPLIT_FAMILY_HITS             } from '../../../modules/local/split_family_hits'
include { fileStem                      } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { MMSEQS_FASTA_CLUSTER          } from '../../../subworkflows/nf-core/mmseqs_fasta_cluster'
include { REMOVE_REDUNDANT_SEQS         } from '../../../modules/local/remove_redundant_seqs/main'
include { GENERATE_FAMILIES             } from '../../../subworkflows/local/generate_families'
include { HMMER_HMMALIGN                } from '../../../modules/nf-core/hmmer/hmmalign/main'
include { EXTRACT_UNASSIGNED_SEQS       } from '../../../modules/local/extract_unassigned_seqs/main'

workflow UPDATE_FAMILIES {
    take:
    ch_samplesheet_for_update           // channel: [meta, sequences, existing_hmms, existing_seed_msas, existing_full_msas]; MSAs may be []
    recruit_min_model_coverage          // number [0.0, 1.0]
    skip_sequence_redundancy_removal    // boolean
    clustering_tool                     // string ["linclust", "cluster"]
    alignment_tool                      // string ["famsa", "mafft"]
    skip_seed_msa_trimming              // boolean
    skip_recruiting                     // boolean
    skip_update_refinement              // boolean: keep the existing HMMs, only rebuild the full MSAs

    main:
    ch_input_for_untar = ch_samplesheet_for_update
        .multiMap { meta, _fasta, existing_hmms, existing_seed_msas, existing_full_msas ->
            hmm: [ meta, existing_hmms ]
            seed_msa: [ meta, existing_seed_msas ]
            full_msa: [ meta, existing_full_msas ]
        }

    // One <NAME>.hmm.gz per model, from an archive or a library
    SPLIT_HMMS( ch_input_for_untar.hmm )
    UNTAR_SEED_MSA( ch_input_for_untar.seed_msa.filter { _meta, archive -> archive } )
    UNTAR_FULL_MSA( ch_input_for_untar.full_msa.filter { _meta, archive -> archive } )

    // Families are named by HMM NAME; MSAs are matched to them by file stem
    validateHmmNames( SPLIT_HMMS.out.hmms )
    validateMsaStems( SPLIT_HMMS.out.hmms, UNTAR_SEED_MSA.out.untar )
    validateMsaStems( SPLIT_HMMS.out.hmms, UNTAR_FULL_MSA.out.untar )

    // Squeeze the HMMs into a single file
    CAT_HMM( SPLIT_HMMS.out.hmms.map { meta, folder -> [meta, file("${folder.toUriString()}/*", checkIfExists: true)] } )

    // Provided family files, [id, family] meta
    ch_existing_hmm      = familyFiles( SPLIT_HMMS.out.hmms )
    ch_existing_seed_msa = familyFiles( UNTAR_SEED_MSA.out.untar )
    ch_existing_full_msa = familyFiles( UNTAR_FULL_MSA.out.untar )

    // The searched pool: the input sequences, plus the members of existing full MSAs, so that
    // the families keep the old members that still hit. HMMER rewinds the target database for
    // every query HMM and a gzip stream cannot rewind, so the pool is always uncompressed.
    ch_input_fasta = ch_samplesheet_for_update
        .map { meta, fasta, _existing_hmms, _existing_seed_msas, _existing_full_msas -> [meta, fasta] }

    POOL_EXISTING_MEMBERS( ch_input_fasta.join(UNTAR_FULL_MSA.out.untar) )
    // The degapped members of each existing full MSA, the FASTA of a family that passes through
    ch_existing_fasta = familyFiles( POOL_EXISTING_MEMBERS.out.members )

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

    // ENGINE (standard)
    HMMER_HMMSEARCH( ch_input_for_hmmsearch )

    // Hits are cut from the pool
    ch_input_for_split_hits = HMMER_HMMSEARCH.out.domain_summary
        .join(ch_pool)
        .multiMap { meta, domtbl, pool ->
            domtbl: [ meta, domtbl ]
            fasta: [ meta, pool ]
        }

    SPLIT_FAMILY_HITS ( ch_input_for_split_hits.fasta, ch_input_for_split_hits.domtbl, recruit_min_model_coverage )

    // [id, family] meta, as for created families' chunks
    ch_hits = SPLIT_FAMILY_HITS.out.hits
        .transpose()
        .map { meta, file -> [[id: meta.id, family: fileStem(file)], file] }

    if (skip_update_refinement) {
        // The existing HMMs stay as they are and align their own hits into the new full MSAs
        ch_input_for_hmmalign = ch_hits
            .join(ch_existing_hmm)
            .multiMap { meta, seqs, hmm ->
                seq: [ meta, seqs ]
                hmm: [ hmm ]
            }

        HMMER_HMMALIGN( ch_input_for_hmmalign.seq, ch_input_for_hmmalign.hmm )

        ch_engine_seed_msa = channel.empty()
        ch_engine_full_msa = HMMER_HMMALIGN.out.sto
        ch_engine_fasta    = ch_hits
        ch_engine_hmm      = channel.empty()
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
            skip_recruiting,
            recruit_min_model_coverage
        )

        // A family whose new HMM recruits nothing has no full MSA, so it is not returned
        ch_engine_full_msa = GENERATE_FAMILIES.out.full_msa
        ch_engine_seed_msa = GENERATE_FAMILIES.out.seed_msa.join(ch_engine_full_msa).map { meta, seed, _full -> [meta, seed] }
        ch_engine_hmm      = GENERATE_FAMILIES.out.hmm.join(ch_engine_full_msa).map { meta, hmm, _full -> [meta, hmm] }
        ch_engine_fasta    = GENERATE_FAMILIES.out.fasta
    }
    // Families with hits that the engine did not return
    ch_engine_reasons = ch_hits
        .join(ch_engine_full_msa, remainder: true)
        .filter { _meta, hits, full_msa -> hits && !full_msa }
        .map { meta, _hits, _full_msa -> [meta, 'no recruits'] }
    // END ENGINE

    // FINALISE
    ch_seed_msa = byReturned(passThrough(ch_existing_seed_msa, ch_engine_seed_msa), ch_engine_full_msa)
    ch_full_msa = byReturned(passThrough(ch_existing_full_msa, ch_engine_full_msa), ch_engine_full_msa)
    ch_hmm      = byReturned(passThrough(ch_existing_hmm, ch_engine_hmm), ch_engine_full_msa)
    ch_fasta    = byReturned(passThrough(ch_existing_fasta, ch_engine_fasta), ch_engine_full_msa)

    // Existing families the engine did not return pass through, listed per sample
    ch_passed_through_families = ch_existing_hmm
        .join(ch_engine_full_msa, remainder: true)
        .filter { _meta, hmm, full_msa -> hmm && !full_msa }
        .map { meta, _hmm, _full_msa -> [meta, 'no hits'] }
        .join(ch_engine_reasons, remainder: true)
        .filter { _meta, default_reason, _reason -> default_reason }
        .map { meta, default_reason, reason -> [[id: meta.id], [meta.family, reason ?: default_reason]] }
        .groupTuple()
    ch_passed_through_families = SPLIT_HMMS.out.hmms
        .join(ch_passed_through_families, remainder: true)
        .map { meta, _folder, passed_through -> [meta, (passed_through ?: []).sort { family_reason -> family_reason[0] }] }

    // Input sequences not in any returned family go to family creation
    EXTRACT_UNASSIGNED_SEQS(
        ch_input_fasta
            .join(ch_fasta.returned.map { meta, faa -> [ [id: meta.id], faa ] }.groupTuple(by: 0), remainder: true)
            .map { meta, fasta, family_fastas -> [meta, fasta, family_fastas ?: []] }
    )

    emit:
    seed_msa                = ch_seed_msa.returned
    full_msa                = ch_full_msa.returned
    fasta                   = ch_fasta.returned
    hmm                     = ch_hmm.returned
    passed_through_seed_msa = ch_seed_msa.passed_through
    passed_through_full_msa = ch_full_msa.passed_through
    passed_through_fasta    = ch_fasta.passed_through
    passed_through_hmm      = ch_hmm.passed_through
    passed_through_families = ch_passed_through_families // [meta, [[family, reason], ...]], [] if every family was returned
    no_hit_seqs             = EXTRACT_UNASSIGNED_SEQS.out.fasta
    pool                    = ch_pool
}

// One [[id, family], file] per file of each sample's folder
def familyFiles(ch_folders) {
    ch_folders.flatMap { meta, folder -> folder.listFiles().collect { file -> [[id: meta.id, family: fileStem(file)], file] } }
}

// A provided family file ([id, family] meta) unless the engine built a new one for that family
def passThrough(ch_existing_files, ch_engine_files) {
    ch_engine_files.mix(
        ch_existing_files
            .join(ch_engine_files, remainder: true)
            .filter { _meta, existing, engine -> existing && !engine }
            .map { meta, existing, _engine -> [meta, existing] }
    )
}

// A family's files, branched by whether the engine returned the family (emitted its full MSA)
def byReturned(ch_files, ch_engine_full_msa) {
    ch_files
        .join(ch_engine_full_msa.map { meta, _full_msa -> [meta, true] }, remainder: true)
        .filter { _meta, file, _returned -> file }
        .branch { meta, file, returned ->
            returned: returned
                return [meta, file]
            passed_through: true
                return [meta, file]
        }
}
