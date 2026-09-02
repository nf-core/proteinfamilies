process BRANCH_HITS_FASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8372f6241b480332d91bc00a88ec8c72c8f7fcc9994177a5dd67a07007cd6e32/data' :
        'community.wave.seqera.io/library/biopython:1.85--6f761292fa9881b4' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(domtbl)
    val(length_threshold)

    output:
    tuple val(meta), path("hits/*")    , emit: hits
    tuple val(meta), path("*.fasta.gz"), emit: non_hit_fasta
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions
    tuple val("${task.process}"), val('biopython'), eval("python -c \"import importlib.metadata; print(importlib.metadata.version('biopython'))\""), emit: versions_biopython, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    branch_hits_fasta.py \\
        --fasta ${fasta} \\
        --domtbl ${domtbl} \\
        --length_threshold ${length_threshold} \\
        --hits hits \\
        --non_hits ${prefix}.fasta.gz
    """

    stub:
    """
    mkdir -p hits
    touch hits/test.fasta
    echo "" | gzip > test.fasta.gz
    """
}
