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

    output:
    tuple val(meta), path("chunked_fasta/*"), emit: fasta_chunks, optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def is_compressed = sequences.getName().endsWith(".gz") ? true : false
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
        --out_folder chunked_fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p chunked_fasta
    touch chunked_fasta/1.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """
}
