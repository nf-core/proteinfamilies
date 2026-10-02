/*
    MULTIPLE SEQUENCE ALIGNMENT

    Dispatches to FAMSA or MAFFT based on alignment_tool. Any value other than 'famsa'
    falls back to MAFFT. Unless skip_trimming, the alignment is then trimmed with ClipKIT
    and its rows' name/start-end coordinates recalculated to the residues they still hold.
    The emitted sequences always match the emitted alignments.
*/

include { FAMSA_ALIGN             } from '../../../modules/nf-core/famsa/align/main'
include { MAFFT_ALIGN             } from '../../../modules/nf-core/mafft/align/main'
include { CLIPKIT                 } from '../../../modules/nf-core/clipkit/main'
include { RECALCULATE_COORDINATES } from '../../../modules/local/recalculate_coordinates/main'

workflow ALIGN_SEQUENCES {
    take:
    sequences      // tuple val(meta), path(fasta)
    alignment_tool // string: MSA tool
    skip_trimming  // boolean

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

    ch_sequences = sequences
    if (!skip_trimming) {
        // ClipKIT writes FASTA (the input format); its 'clipkit' extension never clashes with the aligners' .aln/.fas
        CLIPKIT( ch_alignments, 'clipkit', [] )

        // The trimmed rows are rebuilt from the untrimmed MSA and the log's keep columns, so CLIPKIT.out.clipkit is unused
        RECALCULATE_COORDINATES( ch_alignments.join(CLIPKIT.out.log) )
        ch_alignments = RECALCULATE_COORDINATES.out.alignment
        ch_sequences  = RECALCULATE_COORDINATES.out.fasta
    }

    emit:
    alignments = ch_alignments // tuple val(meta), path(msa)
    sequences  = ch_sequences  // tuple val(meta), path(fasta): the degapped rows of alignments
}
