/*
    CLUSTER CHUNKING AND FAMILY MODEL GENERATION

    Splits a sample's clusters into chunks and turns each chunk into family models with the
    selected algorithm. Chunking lives here because its form follows the algorithm: the
    standard one aligns one cluster per task and so takes a FASTA per cluster, while the
    iterative one builds many families per task and so takes a TSV of clusters plus the
    sequence pool to recruit from.

    Both algorithms emit the same per-family channels, so everything downstream is unaware
    of which one ran.
*/

include { CHUNK_CLUSTERS               } from '../../../modules/local/chunk_clusters/main'
include { GENERATE_FAMILIES            } from '../../../subworkflows/local/generate_families'
include { GENERATE_FAMILIES_ITERATIVELY } from '../../../subworkflows/local/generate_families_iteratively'

workflow CHUNK_AND_GENERATE_FAMILIES {
    take:
    sequences                           // tuple val(meta), path(fasta)
    clusters                            // tuple val(meta), path(tsv)
    family_generation_algorithm         // string ["standard", "iterative"]
    cluster_size_threshold              // integer
    clusters_per_chunk                  // integer
    alignment_tool                      // string ["famsa", "mafft"]
    skip_msa_trimming                   // boolean
    clipkit_out_format                  // string (default: clipkit)
    hmmsearch_write_target              // boolean
    hmmsearch_write_domain              // boolean
    skip_additional_sequence_recruiting // boolean
    hmmsearch_query_length_threshold    // number [0.0, 1.0]

    main:
    def iterative = family_generation_algorithm == 'iterative'

    CHUNK_CLUSTERS( clusters, sequences, cluster_size_threshold, iterative ? 'tsv' : 'fasta', clusters_per_chunk )

    // tokenize('_').last() extracts the numeric suffix from filenames like 'sample_1.faa.gz' as the chunk ID.
    ch_chunks = ( iterative ? CHUNK_CLUSTERS.out.cluster_chunks : CHUNK_CLUSTERS.out.fasta_chunks )
        .transpose()
        .map { meta, file_path ->
            [ [id: meta.id, chunk: file(file_path, checkIfExists: true).baseName.tokenize('_').last()], file_path ]
        }

    if (iterative) {
        GENERATE_FAMILIES_ITERATIVELY( sequences, ch_chunks )
        ch_families = GENERATE_FAMILIES_ITERATIVELY.out
    } else {
        GENERATE_FAMILIES (
            sequences,
            ch_chunks,
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
