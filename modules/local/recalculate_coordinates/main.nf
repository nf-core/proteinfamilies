process RECALCULATE_COORDINATES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d0fee0217685dc2501570e1dd12076f9466175d0a37335ca424390cffec5fa1/data' :
        'community.wave.seqera.io/library/python:3.13.1--d00663700fcc8bcf' }"

    input:
    // Untrimmed MSA staged in a subfolder, so it never clashes with the outputs (e.g. FAMSA's ${prefix}.aln)
    tuple val(meta), path(untrimmed, stageAs: "untrimmed/*"), path(trim_log)

    output:
    tuple val(meta), path("${prefix}.aln"), emit: alignment
    tuple val(meta), path("${prefix}.faa"), emit: fasta
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    recalculate_coordinates.py \\
        --untrimmed ${untrimmed} \\
        --log ${trim_log} \\
        --out_msa ${prefix}.aln \\
        --out_fasta ${prefix}.faa
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aln
    touch ${prefix}.faa
    """
}
