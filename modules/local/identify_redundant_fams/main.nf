process IDENTIFY_REDUNDANT_FAMS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

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
