/*
    INFER SEQUENCE PHYLOGENY OF FULL ALIGNMENTS
*/

include { GUNZIP } from '../../../modules/nf-core/gunzip/main'
include { CMAPLE } from '../../../modules/nf-core/cmaple/main'
workflow INFER_PHYLOGENY {
    take:
    full_msa                             // tuple val(meta), path({aln,fas,clipkit,fas.gz})
    skip_sequence_redundancy_removal     // boolean
    skip_additional_sequence_recruiting  // boolean
    skip_family_redundancy_removal       // boolean

    main:

    ch_versions = channel.empty()

    if (skip_sequence_redundancy_removal
        && !(
            skip_additional_sequence_recruiting
            && skip_family_redundancy_removal
        ) // if not already in uncompressed state
    ) {
        full_msa    = GUNZIP( full_msa ).gunzip // cmaple does not work with .gz alignments
        ch_versions = ch_versions.mix( GUNZIP.out.versions.first() )
    }

    CMAPLE (
        full_msa
            .map { meta, file -> [ meta, file, [] ] }
    )
    ch_versions = ch_versions.mix( CMAPLE.out.versions.first() )

    emit:
    treefile   = CMAPLE.out.treefile
    cmaple_log = CMAPLE.out.log
    versions   = ch_versions
}
