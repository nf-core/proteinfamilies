/*
    MULTIPLE SEQUENCE ALIGNMENT
*/

include { FAMSA_ALIGN } from '../../../modules/nf-core/famsa/align/main'
include { MAFFT_ALIGN } from '../../../modules/nf-core/mafft/align/main'

workflow ALIGN_SEQUENCES {
    take:
    sequences      // tuple val(meta), path(fasta)
    alignment_tool // string: MSA tool

    main:
    ch_versions   = Channel.empty()
    ch_alignments = Channel.empty()

    if (alignment_tool == 'famsa') {
        alignment_res = FAMSA_ALIGN( sequences, [[:],[]], false )
        ch_versions = ch_versions.mix( FAMSA_ALIGN.out.versions.first() )
        ch_alignments = alignment_res.alignment
    } else { // fallback: mafft
        alignment_res = MAFFT_ALIGN( sequences, [[:], []], [[:], []], [[:], []], [[:], []], [[:], []], false )
        ch_versions = ch_versions.mix( MAFFT_ALIGN.out.versions.first() )
        ch_alignments = alignment_res.fas
    }

    emit:
    versions   = ch_versions
    alignments = ch_alignments
}
