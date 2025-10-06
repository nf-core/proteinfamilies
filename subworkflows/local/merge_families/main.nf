/*
    FAMILY MERGING
*/

include { POOL_SIMILAR_COMPONENTS } from '../../../modules/local/pool_similar_components/main'
include { MERGE_SEEDS             } from '../../../modules/local/merge_seeds/main'
// include { GENERATE_FAMILIES  } from '../../../subworkflows/local/generate_families'

workflow MERGE_FAMILIES {
    take:
    similarities // tuple val(meta), path(txt)
    seed_msa     // tuple val(meta), path(aln)

    main:
    ch_versions = Channel.empty()

    POOL_SIMILAR_COMPONENTS( similarities )
    ch_versions = ch_versions.mix( POOL_SIMILAR_COMPONENTS.out.versions )

    ch_pooled_components = POOL_SIMILAR_COMPONENTS.out.pooled_components
        .splitCsv( by:1 )
        .map { meta, components ->
            // Extract the suffix from each component (part after last underscore)
            def suffixes = components.collect { component ->
                component.split('_').last()
            }
            // Join the suffixes with underscores
            def combinedSuffix = suffixes.join('_')
            // Create new meta with combined ID
            def newMeta = meta + [id: "${meta.id}_${combinedSuffix}"]
            return [newMeta, components]
        }

    MERGE_SEEDS( ch_pooled_components, seed_msa.first() )
    ch_versions = ch_versions.mix( MERGE_SEEDS.out.versions )

    // GENERATE_FAMILIES( )
    // ch_versions = ch_versions.mix( GENERATE_FAMILIES.out.versions )

    emit:
    versions = ch_versions
}
