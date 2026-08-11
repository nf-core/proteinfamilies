/*
    FAMILY MERGING

    Groups similar families into pools, merges their seed MSAs into a single combined
    seed, then rebuilds final family models via GENERATE_FAMILIES. The merged_id
    encodes which original families were combined (e.g., 'sample_1_7' from 'sample_1'
    and 'sample_7'); for very large pools it collapses to a stable hash so the
    resulting output filename stays within the filesystem's name-length limit.
*/

include { POOL_SIMILAR_COMPONENTS       } from '../../../modules/local/pool_similar_components/main'
include { MERGE_SEEDS                   } from '../../../modules/local/merge_seeds/main'
include { GENERATE_FAMILIES             } from '../../../subworkflows/local/generate_families'
include { GENERATE_FAMILIES_ITERATIVELY } from '../../../subworkflows/local/generate_families_iteratively'

workflow MERGE_FAMILIES {
    take:
    similarities                        // tuple val(meta), path(txt)
    seed_msa                            // tuple val(meta), path(aln)
    sequences                           // tuple val(meta), path(fasta)
    family_generation_algorithm         // string ["standard", "iterative"]
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
            // Extract each component's family suffix, the part after the sample id. Splitting
            // on the last underscore instead would collapse the compound suffixes the
            // iterative algorithm produces ('2_1' and '3_1' would both become '1').
            // split() takes a regex, which is safe here because assets/schema_input.json
            // restricts sample names to letters, digits, dots, underscores and dashes.
            def suffixes = components.collect { component ->
                component.split("${meta.id}_", 2)[1]
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

    // Each pooled group is paired with its own sample's seed collection. combine() by sample id
    // repeats that collection for every group of the sample, without letting one sample's groups
    // consume another sample's seeds.
    ch_input_for_merge_seeds = ch_pooled_components
        .map { meta, components -> [ [id: meta.id], meta, components ] }
        .combine(seed_msa, by: 0)
        .multiMap { id, meta, components, seeds ->
            components: [ meta, components ]
            seed_msa  : [ id, seeds ]
        }

    MERGE_SEEDS( ch_input_for_merge_seeds.components, ch_input_for_merge_seeds.seed_msa )

    // The iterative algorithm rebuilds a family from a cluster rather than from a seed
    // alignment, so it takes the merged membership MERGE_SEEDS writes alongside the seed.
    if (family_generation_algorithm == 'iterative') {
        GENERATE_FAMILIES_ITERATIVELY( sequences, MERGE_SEEDS.out.clusters )
        ch_families = GENERATE_FAMILIES_ITERATIVELY.out
    } else {
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
        ch_families = GENERATE_FAMILIES.out
    }

    emit:
    seed_msa = ch_families.seed_msa
    full_msa = ch_families.full_msa
    fasta    = ch_families.fasta
    hmm      = ch_families.hmm
}
