process EXTRACT_UNASSIGNED_SEQS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d0fee0217685dc2501570e1dd12076f9466175d0a37335ca424390cffec5fa1/data' :
        'community.wave.seqera.io/library/python:3.13.1--d00663700fcc8bcf' }"

    input:
    tuple val(meta), path(fasta), path(families, stageAs: "families/*") // families: [] if no family was updated

    output:
    tuple val(meta), path("${prefix}.faa.gz"), emit: fasta
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_unassigned"
    def families_arg = families ? "--families ${families}" : ''
    """
    extract_unassigned_seqs.py \\
        --fasta ${fasta} \\
        ${families_arg} \\
        --out_fasta ${prefix}.faa.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_unassigned"
    """
    echo "" | gzip > ${prefix}.faa.gz
    """
}
