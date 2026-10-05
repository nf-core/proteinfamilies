process FILTER_NON_REDUNDANT_FAMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d0fee0217685dc2501570e1dd12076f9466175d0a37335ca424390cffec5fa1/data' :
        'community.wave.seqera.io/library/python:3.13.1--d00663700fcc8bcf' }"

    input:
    tuple val(meta) , path(files, stageAs: "input_folder/*"), path(kept, stageAs: "kept/*")
    tuple val(meta2), path(redundant_ids)

    output:
    tuple val(meta), path(pattern), emit: filtered
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Created and updated families' MSAs may differ in format (e.g. aligned FASTA and Stockholm)
    def extensions = ([files] + [kept]).flatten().collect { file -> file.extension }.unique().sort()
    pattern = extensions.size() == 1 ? "*.${extensions[0]}" : "*.{${extensions.join(',')}}"
    """
    filter_non_redundant_fams.py \\
        --input_folder input_folder  \\
        --kept_folder kept \\
        --redundant_ids ${redundant_ids}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = [files].flatten()[0].extension
    pattern = "*.${extension}"
    """
    touch ${prefix}_1.${extension}
    """
}
