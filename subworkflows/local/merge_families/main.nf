/*
    FAMILY MERGING

    Groups similar families into pools, merges their seed MSAs into a single combined
    seed, then rebuilds final family models via GENERATE_FAMILIES. The merged_id
    encodes which original families were combined (e.g., 'sample_1_7' from 'sample_1'
    and 'sample_7'); for very large pools it collapses to a stable hash so the
    resulting output filename stays within the filesystem's name-length limit.
*/

include { POOL_SIMILAR_COMPONENTS } from '../../../modules/local/pool_similar_components/main'
include { MERGE_SEEDS             } from '../../../modules/local/merge_seeds/main'
include { GENERATE_FAMILIES       } from '../../../subworkflows/local/generate_families'

workflow MERGE_FAMILIES {
    take:
    similarities                        // tuple val(meta), path(txt)
    seed_msa                            // tuple val(meta), path(aln)
    sequences                           // tuple val(meta), path(fasta)
    alignment_tool                      // string ["famsa", "mafft"]
    skip_msa_trimming                   // boolean
    clipkit_out_format                  // string (default: clipkit)
    hmmsearch_write_target              // boolean
    hmmsearch_write_domain              // boolean
    skip_additional_sequence_recruiting // boolean
    hmmsearch_query_length_threshold    // number [0.0, 1.0]

    main:

    POOL_SIMILAR_COMPONENTS( similarities )

    ch_pooled_components = POOL_SIMILAR_COMPONENTS.out.pooled_components
        .splitCsv( by:1 )
        .map { meta, components ->
            // Extract the suffix from each component (part after last underscore)
            def suffixes = components.collect { component ->
                component.split('_').last()
            }
            // Readable id encoding every combined family, e.g. 'sample_1_7'
            def readableId = "${meta.id}_${suffixes.join('_')}"
            // merged_id becomes the output-file prefix for every merged-family process, so it
            // must fit the filesystem's 255-byte name limit. Large pools (dozens of families)
            // would overflow it; in that case fall back to a short, stable hash of the members.
            def merged_id = readableId.length() <= 200
                ? readableId
                : "${meta.id}_${suffixes.size()}fams_${suffixes.join('_').md5().take(10)}"
            // Keep original id, add new field merged_id
            def newMeta = meta + [merged_id: merged_id]
            return [newMeta, components.join(',')]
        }

    // .first() converts seed_msa to a value channel so each pooled-component group can combine
    // with the full seed MSA collection; without it, a queue channel would be consumed after
    // the first pairing.
    MERGE_SEEDS( ch_pooled_components, seed_msa.first() )

    GENERATE_FAMILIES (
        sequences,
        MERGE_SEEDS.out.merged_seed_msa,
        alignment_tool,
        skip_msa_trimming,
        clipkit_out_format,
        hmmsearch_write_target,
        hmmsearch_write_domain,
        skip_additional_sequence_recruiting,
        hmmsearch_query_length_threshold
    )

    emit:
    seed_msa = GENERATE_FAMILIES.out.seed_msa
    full_msa = GENERATE_FAMILIES.out.full_msa
    fasta    = GENERATE_FAMILIES.out.fasta
    hmm      = GENERATE_FAMILIES.out.hmm
}
