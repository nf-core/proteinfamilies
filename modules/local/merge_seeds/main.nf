process MERGE_SEEDS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d0fee0217685dc2501570e1dd12076f9466175d0a37335ca424390cffec5fa1/data' :
        'community.wave.seqera.io/library/python:3.13.1--d00663700fcc8bcf' }"

    input:
    tuple val(meta) , val(similarities)
    tuple val(meta2), path(seed_msa, stageAs: "seed_msa/*")

    output:
    tuple val(meta), path("${prefix}.fas")         , emit: merged_seed_msa
    tuple val(meta), path("${prefix}_clusters.tsv"), emit: clusters
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_seeds.py \\
        --list "${similarities}" \\
        --folder seed_msa \\
        --out_file ${prefix}.fas \\
        --out_clusters ${prefix}_clusters.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fas
    touch ${prefix}_clusters.tsv
    """
}
