process SPLIT_HMMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d0fee0217685dc2501570e1dd12076f9466175d0a37335ca424390cffec5fa1/data' :
        'community.wave.seqera.io/library/python:3.13.1--d00663700fcc8bcf' }"

    input:
    tuple val(meta), path(hmms) // .tar.gz archive of HMM files, or one HMM library

    output:
    tuple val(meta), path("${prefix}"), emit: hmms
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_hmms"
    """
    split_hmms.py \\
        --input ${hmms} \\
        --outdir ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_hmms"
    """
    mkdir ${prefix}
    echo "" | gzip > ${prefix}/${meta.id}_stub.hmm.gz
    """
}
