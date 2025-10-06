/*
    FAMILY MERGING
*/

include { POOL_SIMILAR_COMPONENTS } from '../../../modules/local/pool_similar_components/main'
// include { MERGE_SEEDS             } from '../../../modules/local/merge_seeds/main'
// include { GENERATE_FAMILIES  } from '../../../subworkflows/local/generate_families'

workflow MERGE_FAMILIES {
    take:
    similarities // tuple val(meta), path(txt)
    seed_msa     // tuple val(meta), path(aln)

    main:
    ch_versions = Channel.empty()

    similarities.view()
    seed_msa.view()
    POOL_SIMILAR_COMPONENTS( similarities )
    ch_versions = ch_versions.mix( POOL_SIMILAR_COMPONENTS.out.versions )
    // MERGE_SEEDS( )
    // ch_versions = ch_versions.mix( MERGE_SEEDS.out.versions )

    // GENERATE_FAMILIES( )
    // ch_versions = ch_versions.mix( GENERATE_FAMILIES.out.versions )

    emit:
    versions = ch_versions
}
