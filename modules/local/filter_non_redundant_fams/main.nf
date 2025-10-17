process FILTER_NON_REDUNDANT_FAMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d0fee0217685dc2501570e1dd12076f9466175d0a37335ca424390cffec5fa1/data' :
        'community.wave.seqera.io/library/python:3.13.1--d00663700fcc8bcf' }"

    input:
    tuple val(meta) , path(files, stageAs: "input_folder/*")
    tuple val(meta2), path(redundant_ids)

    output:
    tuple val(meta), path("*.${extension}"), emit: filtered
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    extension = files instanceof List ? files[0].extension : files.extension
    """
    filter_non_redundant_fams.py \\
        --input_folder input_folder  \\
        --redundant_ids ${redundant_ids}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    extension = files[0].extension
    """
    touch test.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
