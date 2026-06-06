process POOL_SIMILAR_COMPONENTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/64/64288a2d9e2c21c4f3cd9f8cb8476c04d323de98dbd8a310265077566f9bb75a/data' :
        'community.wave.seqera.io/library/networkx_pandas_python:ab2ee9a9e2c80a69' }"

    input:
    tuple val(meta), path(similarities)

    output:
    tuple val(meta), path("pooled_components.txt"), emit: pooled_components
    tuple val("${task.process}"), val('python'),   eval("python --version 2>&1 | sed 's/Python //'"),                                                               emit: versions_python,  topic: versions
    tuple val("${task.process}"), val('pandas'),   eval("python -c \"import importlib.metadata; print(importlib.metadata.version('pandas'))\""),                    emit: versions_pandas,  topic: versions
    tuple val("${task.process}"), val('networkx'), eval("python -c \"import importlib.metadata; print(importlib.metadata.version('networkx'))\""),                  emit: versions_networkx, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pool_similar_components.py \\
        --input_csv ${similarities} \\
        --out_file pooled_components.txt
    """

    stub:
    """
    touch pooled_components.txt
    """
}
