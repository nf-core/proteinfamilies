/*
    MULTIPLE SEQUENCE ALIGNMENT

    Dispatches to FAMSA or MAFFT based on alignment_tool. Any value other than 'famsa'
    falls back to MAFFT.
*/

include { FAMSA_ALIGN } from '../../../modules/nf-core/famsa/align/main'
include { MAFFT_ALIGN } from '../../../modules/nf-core/mafft/align/main'

workflow ALIGN_SEQUENCES {
    take:
    sequences      // tuple val(meta), path(fasta)
    alignment_tool // string: MSA tool

    main:
    ch_alignments = channel.empty()

    if (alignment_tool == 'famsa') {
        alignment_res = FAMSA_ALIGN( sequences, [[:],[]], false )
        ch_alignments = alignment_res.alignment
    } else { // fallback: mafft
        // [[:],[]] placeholders for optional MAFFT inputs (addfragments, seed alignment, query,
        // gap-open penalties, gap-extend penalties) that are not used in this pipeline.
        alignment_res = MAFFT_ALIGN( sequences, [[:], []], [[:], []], [[:], []], [[:], []], [[:], []], false )
        ch_alignments = alignment_res.fas
    }

    emit:
    alignments = ch_alignments
}
