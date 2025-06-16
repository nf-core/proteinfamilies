process IDENTIFY_REDUNDANT_FAMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4a/4af7c28266e09eb3af86e75b1c61e9876d2827acddfe793cef5a704a1e98fb9b/data' :
        'community.wave.seqera.io/library/pandas_python:f5ae1ede7f6e9d73' }"

    input:
    tuple val(meta) , path(mapping)
    tuple val(meta2), path(domtbl)
    val(length_threshold)

    output:
    tuple val(meta), path("redundant_fam_ids.txt"), emit: redundant_ids
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    identify_redundant_fams.py \\
        --mapping ${mapping} \\
        --domtbl ${domtbl} \\
        --length_threshold ${length_threshold} \\
        --out_file redundant_fam_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
    END_VERSIONS
    """

    stub:
    """
    touch redundant_fam_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import importlib.metadata; print(importlib.metadata.version('pandas'))")
    END_VERSIONS
    """
}
