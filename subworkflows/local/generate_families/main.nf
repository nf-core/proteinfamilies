/*
    FAMILY MODEL GENERATION
*/

include { ALIGN_SEQUENCES  } from '../../../subworkflows/local/align_sequences'
include { CLIPKIT          } from '../../../modules/nf-core/clipkit/main'
include { HMMER_HMMBUILD   } from '../../../modules/nf-core/hmmer/hmmbuild/main'
include { HMMER_HMMSEARCH  } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { FILTER_RECRUITED } from '../../../modules/local/filter_recruited/main'
include { HMMER_HMMALIGN   } from '../../../modules/nf-core/hmmer/hmmalign/main'

workflow GENERATE_FAMILIES {
    take:
    sequences                        // tuple val(meta), path(fasta)
    ch_fasta                         // tuple val(meta), path(fasta)
    alignment_tool                   // string ["famsa", "mafft"]
    skip_msa_trimming                // boolean
    clipkit_out_format               // string (default: clipkit)
    hmmsearch_write_target           // boolean
    hmmsearch_write_domain           // boolean
    recruit_sequences_with_models    // boolean
    hmmsearch_query_length_threshold // number [0.0, 1.0]

    main:
    ch_versions = Channel.empty()
    ch_seed_msa = Channel.empty()
    ch_full_msa = Channel.empty()
    ch_hmm      = Channel.empty()

    ALIGN_SEQUENCES( ch_fasta, alignment_tool )
    ch_versions = ch_versions.mix( ALIGN_SEQUENCES.out.versions )
    ch_seed_msa = ALIGN_SEQUENCES.out.alignments

    if (!skip_msa_trimming) {
        CLIPKIT( ch_seed_msa, clipkit_out_format )
        ch_versions = ch_versions.mix( CLIPKIT.out.versions )
        ch_seed_msa = CLIPKIT.out.clipkit
    }

    HMMER_HMMBUILD( ch_seed_msa, [] )
    ch_versions = ch_versions.mix( HMMER_HMMBUILD.out.versions )
    ch_hmm = HMMER_HMMBUILD.out.hmm

    // Combine with same id to ensure in sync
    ch_input_for_hmmsearch = ch_hmm
        .map { meta, hmm -> [ [id: meta.id], meta, hmm ] }
        .combine(sequences, by: 0)
        .map { id, meta, hmm, seqs -> [ meta, hmm, seqs, false, hmmsearch_write_target, hmmsearch_write_domain ] }

    if (recruit_sequences_with_models) {
        HMMER_HMMSEARCH( ch_input_for_hmmsearch )
        ch_versions = ch_versions.mix( HMMER_HMMSEARCH.out.versions )

        // Combine with same id to ensure in sync
        ch_input_for_filter_recruited = HMMER_HMMSEARCH.out.domain_summary
            .map { meta, domtbl -> [ [id: meta.id], meta, domtbl ] }
            .combine(sequences, by: 0)
            .map { id, meta, domtbl, seqs -> [ meta, domtbl, seqs ] }

        FILTER_RECRUITED( ch_input_for_filter_recruited, hmmsearch_query_length_threshold )
        ch_versions = ch_versions.mix( FILTER_RECRUITED.out.versions )
        ch_fasta = FILTER_RECRUITED.out.fasta

        // Join to ensure in sync
        ch_input_for_hmmalign = ch_fasta
            .join(ch_hmm)
            .multiMap { meta, seqs, hmms ->
                seq: [ meta, seqs ]
                hmm: [ hmms ]
            }

        HMMER_HMMALIGN( ch_input_for_hmmalign.seq, ch_input_for_hmmalign.hmm )
        ch_versions = ch_versions.mix( HMMER_HMMALIGN.out.versions )
        ch_full_msa = HMMER_HMMALIGN.out.sto
    } else {
        ch_full_msa = ch_seed_msa
    }

    emit:
    versions = ch_versions
    seed_msa = ch_seed_msa
    full_msa = ch_full_msa
    fasta    = ch_fasta
    hmm      = ch_hmm
}
