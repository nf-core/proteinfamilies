/*
    FAMILY MERGING

    Groups similar families into pools, merges their seed MSAs into a single combined
    seed, then rebuilds final family models via GENERATE_FAMILIES. The merged_id
    encodes which original families were combined (e.g., 'sample_1_7' from 'sample_1'
    and 'sample_7'); for very large pools it collapses to a stable hash so the
    resulting output filename stays within the filesystem's name-length limit. With
    merged_family_name 'existing', a merge holding an updated family (at most one, see
    POOL_SIMILAR_COMPONENTS) takes its name instead.
*/

include { POOL_SIMILAR_COMPONENTS       } from '../../../modules/local/pool_similar_components/main'
include { fileStem                      } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { isCreatedFamily               } from '../../../subworkflows/local/utils_nfcore_proteinfamilies_pipeline'
include { MERGE_SEEDS                   } from '../../../modules/local/merge_seeds/main'
include { GENERATE_FAMILIES             } from '../../../subworkflows/local/generate_families'
include { GENERATE_FAMILIES_ITERATIVELY } from '../../../subworkflows/local/generate_families_iteratively'

workflow MERGE_FAMILIES {
    take:
    similarities                        // tuple val(meta), path(csv), val(updated_families)
    seed_msa                            // tuple val(meta), path(aln)
    sequences                           // tuple val(meta), path(fasta), meta [id, pool: 'create' or 'update']
    family_generation_algorithm         // string ["standard", "iterative"]
    alignment_tool                      // string ["famsa", "mafft"]
    skip_seed_msa_trimming              // boolean
    skip_additional_sequence_recruiting // boolean
    hmmsearch_query_length_threshold    // number [0.0, 1.0]
    merged_family_name                  // string ["existing", "new"]

    main:

    POOL_SIMILAR_COMPONENTS( similarities )

    ch_pooled_components = POOL_SIMILAR_COMPONENTS.out.pooled_components
        .splitCsv( by:1 )
        .map { meta, components ->
            // Extract each created component's family suffix, the part after the sample id.
            // Splitting on the last underscore instead would collapse the compound suffixes the
            // iterative algorithm produces ('2_1' and '3_1' would both become '1'). Updated
            // families keep their whole name and come last, so the merged name still starts
            // like a created one ('<id>_<digit>...').
            def created = components.findAll { component -> isCreatedFamily(meta.id, component) }
            def suffixes = created.collect { component -> component.substring(meta.id.length() + 1) } + (components - created)
            // Readable id encoding every combined family, e.g. 'sample_1_7'
            def readableId = "${meta.id}_${suffixes.join('_')}"
            // merged_id becomes the output-file prefix for every merged-family process, so it
            // must fit the filesystem's 255-byte name limit. Large pools (dozens of families)
            // would overflow it; in that case fall back to a short, stable hash of the members.
            def newId = readableId.length() <= 200
                ? readableId
                : "${meta.id}_${suffixes.size()}fams_${suffixes.join('_').md5().take(10)}"
            // A merge holding an updated family (at most one) keeps its name (standard algorithm
            // only: mgnifam names every family it builds `<prefix>_<n>`)
            def updated = (components - created).sort()
            def merged_id = merged_family_name == 'existing' && family_generation_algorithm == 'standard' && updated
                ? updated[0]
                : newId
            // Keep original id, add new field merged_id, and the pool to recruit from: a merge
            // holding an updated family recruits from its sample's update pool
            def newMeta = meta + [merged_id: merged_id, pool: updated ? 'update' : 'create']
            return [newMeta, components.join(',')]
        }

    // Pair each pool with its own sample's seeds, keeping only the pool's members: staging every
    // seed of the sample into every task exhausts the head job's heap on large samples.
    // Seeds are keyed by family ID once per sample, so each pool is a lookup, not a scan.
    ch_input_for_merge_seeds = ch_pooled_components
        .map { meta, components -> [ [id: meta.id], meta, components ] }
        .combine(seed_msa.map { id, seeds -> [ id, seeds.collectEntries { seed -> [ fileStem(seed), seed ] } ] }, by: 0)
        .multiMap { id, meta, components, seedsByFamily ->
            components: [ meta, components ]
            seed_msa  : [ id, seedsByFamily.subMap(components.split(',')).values() as List ]
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
            skip_seed_msa_trimming,
            skip_additional_sequence_recruiting,
            hmmsearch_query_length_threshold
        )
        ch_families = GENERATE_FAMILIES.out
    }

    emit:
    seed_msa        = ch_families.seed_msa
    full_msa        = ch_families.full_msa
    fasta           = ch_families.fasta
    hmm             = ch_families.hmm
    pooled_ids      = POOL_SIMILAR_COMPONENTS.out.pooled_ids
    merged_families = ch_pooled_components.map { meta, components -> [[id: meta.id], [meta.merged_id, components]] }.groupTuple() // [meta, [[merged_id, 'member,...'], ...]]
}
