process CHUNK_CLUSTERS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8372f6241b480332d91bc00a88ec8c72c8f7fcc9994177a5dd67a07007cd6e32/data' :
        'community.wave.seqera.io/library/biopython:1.85--6f761292fa9881b4' }"

    input:
    tuple val(meta) , path(clustering)
    tuple val(meta2), path(sequences)
    val(size_threshold)
    val(out_format)
    val(clusters_per_chunk)

    output:
    tuple val(meta), path("chunked_fasta/*"), emit: fasta_chunks  , optional: true
    tuple val(meta), path("chunked_tsv/*")  , emit: cluster_chunks, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('biopython'), eval("python -c \"import importlib.metadata; print(importlib.metadata.version('biopython'))\""), emit: versions_biopython, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The tsv format only ever reads sequence names out of the clustering file,
    // so the sequence pool is left compressed and untouched.
    def is_compressed = sequences.getName().endsWith(".gz") && out_format == "fasta" ? true : false
    def fasta_name    = sequences.name.replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $sequences > $fasta_name
    fi

    chunk_clusters.py \\
        --clustering ${clustering} \\
        --sequences ${fasta_name} \\
        --threshold ${size_threshold} \\
        --threads $task.cpus \\
        --out_format ${out_format} \\
        --clusters_per_chunk ${clusters_per_chunk} \\
        --out_folder chunked_${out_format}
    """

    stub:
    def extension = out_format == "tsv" ? "tsv" : "fasta"
    """
    mkdir -p chunked_${out_format}
    touch chunked_${out_format}/1.${extension}
    """
}
