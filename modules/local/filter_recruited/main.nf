process FILTER_RECRUITED {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8372f6241b480332d91bc00a88ec8c72c8f7fcc9994177a5dd67a07007cd6e32/data' :
        'community.wave.seqera.io/library/biopython:1.85--6f761292fa9881b4' }"

    input:
    tuple val(meta), path(domtbl), path(fa)
    val(length_threshold)

    output:
    tuple val(meta), path("${prefix}.fasta.gz"), emit: fasta, optional: true
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_recruited.py \\
        --domtbl ${domtbl} \\
        --fasta ${fa} \\
        --length_threshold ${length_threshold} \\
        --out_fasta ${prefix}.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """
}
