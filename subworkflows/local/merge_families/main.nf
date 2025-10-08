/*
    FAMILY MERGING
*/

include { POOL_SIMILAR_COMPONENTS } from '../../../modules/local/pool_similar_components/main'
include { MERGE_SEEDS             } from '../../../modules/local/merge_seeds/main'
include { GENERATE_FAMILIES       } from '../../../subworkflows/local/generate_families'

workflow MERGE_FAMILIES {
    take:
    similarities                     // tuple val(meta), path(txt)
    seed_msa                         // tuple val(meta), path(aln)
    sequences                        // tuple val(meta), path(fasta)
    alignment_tool                   // string ["famsa", "mafft"]
    skip_msa_trimming                // boolean
    clipkit_out_format               // string (default: clipkit)
    hmmsearch_write_target           // boolean
    hmmsearch_write_domain           // boolean
    recruit_sequences_with_models    // boolean
    hmmsearch_query_length_threshold // number [0.0, 1.0]

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
            // Keep original id, add new field merged_id
            def newMeta = meta + [merged_id: "${meta.id}_${combinedSuffix}"]
            return [newMeta, components]
        }

    MERGE_SEEDS( ch_pooled_components, seed_msa.first() )
    ch_versions = ch_versions.mix( MERGE_SEEDS.out.versions )

    GENERATE_FAMILIES( sequences, MERGE_SEEDS.out.merged_seed_msa, \
            alignment_tool, skip_msa_trimming, clipkit_out_format, \
            hmmsearch_write_target, hmmsearch_write_domain, \
            recruit_sequences_with_models, hmmsearch_query_length_threshold )
    ch_versions = ch_versions.mix( GENERATE_FAMILIES.out.versions )

    emit:
    versions = ch_versions
    seed_msa = GENERATE_FAMILIES.out.seed_msa
    full_msa = GENERATE_FAMILIES.out.full_msa
    fasta    = GENERATE_FAMILIES.out.fasta
    hmm      = GENERATE_FAMILIES.out.hmm
}
