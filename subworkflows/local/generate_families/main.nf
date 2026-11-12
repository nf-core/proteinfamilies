/*
    FAMILY MODEL GENERATION

    Builds per-cluster HMMs from seed alignments. Two FASTA channels serve distinct roles:
      ch_fasta   — individual cluster chunks, each aligned into a seed MSA
      sequences  — the full per-sample sequence pool, searched with each cluster HMM to
                   recruit additional members beyond the initial cluster (unless
                   skip_recruiting is true, in which case the seed
                   MSA doubles as the final full MSA and its rows become the family fasta).
*/

include { ALIGN_SEQUENCES  } from '../../../subworkflows/local/align_sequences'
include { HMMER_HMMBUILD   } from '../../../modules/nf-core/hmmer/hmmbuild/main'
include { HMMER_HMMSEARCH  } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { FILTER_RECRUITED } from '../../../modules/local/filter_recruited/main'
include { HMMER_HMMALIGN   } from '../../../modules/nf-core/hmmer/hmmalign/main'

workflow GENERATE_FAMILIES {
    take:
    sequences                           // tuple val(meta), path(fasta)
    ch_fasta                            // tuple val(meta), path(fasta)
    alignment_tool                      // string ["famsa", "mafft"]
    skip_seed_msa_trimming              // boolean
    skip_recruiting                     // boolean
    recruit_min_model_coverage          // number [0.0, 1.0]

    main:
    ch_seed_msa = channel.empty()
    ch_full_msa = channel.empty()
    ch_hmm      = channel.empty()

    ALIGN_SEQUENCES( ch_fasta, alignment_tool, skip_seed_msa_trimming )
    ch_seed_msa = ALIGN_SEQUENCES.out.alignments

    HMMER_HMMBUILD( ch_seed_msa, [] )
    ch_hmm = HMMER_HMMBUILD.out.hmm

    // Combine on a chunk-free [id] key (plus the pool, for merged families) so each cluster's HMM
    // matches the sample sequence pool; the original meta rides along as an extra element and
    // the key is dropped after.
    ch_input_for_hmmsearch = ch_hmm
        .map { meta, hmm -> [ meta.subMap('id', 'pool'), meta, hmm ] }
        .combine(sequences, by: 0)
        .map { _id, meta, hmm, seqs -> [ meta, hmm, seqs, false, false, true ] }

    if (!skip_recruiting) {
        HMMER_HMMSEARCH( ch_input_for_hmmsearch )

        // Combine with same id to ensure in sync
        ch_input_for_filter_recruited = HMMER_HMMSEARCH.out.domain_summary
            .map { meta, domtbl -> [ meta.subMap('id', 'pool'), meta, domtbl ] }
            .combine(sequences, by: 0)
            .map { _id, meta, domtbl, seqs -> [ meta, domtbl, seqs ] }

        FILTER_RECRUITED( ch_input_for_filter_recruited, recruit_min_model_coverage )
        ch_fasta = FILTER_RECRUITED.out.fasta

        // Join to ensure in sync
        ch_input_for_hmmalign = ch_fasta
            .join(ch_hmm)
            .multiMap { meta, seqs, hmms ->
                seq: [ meta, seqs ]
                hmm: [ hmms ]
            }

        HMMER_HMMALIGN( ch_input_for_hmmalign.seq, ch_input_for_hmmalign.hmm )
        ch_full_msa = HMMER_HMMALIGN.out.sto
    } else {
        // Seed MSA serves as the final full MSA when additional sequence recruiting is skipped,
        // and the fasta must hold exactly its (possibly trimmed and renamed) rows.
        ch_full_msa = ch_seed_msa
        ch_fasta    = ALIGN_SEQUENCES.out.sequences
    }

    emit:
    seed_msa = ch_seed_msa
    full_msa = ch_full_msa
    fasta    = ch_fasta
    hmm      = ch_hmm
}
