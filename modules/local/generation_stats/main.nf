process GENERATION_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7d0fee0217685dc2501570e1dd12076f9466175d0a37335ca424390cffec5fa1/data' :
        'community.wave.seqera.io/library/python:3.13.1--d00663700fcc8bcf' }"

    input:
    tuple val(meta), val(families_in), val(discards), path(fasta, stageAs: "fasta/*"), path(hmm, stageAs: "hmm/*"), path(existing_hmm, stageAs: "existing_hmm/*")
    val(remainder_reason)

    output:
    tuple val(meta), path("${prefix}_*_mqc.json"), emit: multiqc
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def discard_args = discards.collect { reason, count -> "--discard '${reason}=${count}'" }.join(' ')
    """
    generation_stats.py \\
        --sample ${prefix} \\
        --families_in ${families_in} \\
        --fasta_folder fasta \\
        --hmm_folder hmm \\
        ${existing_hmm ? '--existing_hmm_folder existing_hmm' : ''} \\
        ${discard_args} \\
        ${remainder_reason ? "--remainder_reason '${remainder_reason}'" : ''} \\
        --prefix ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_general_mqc.json
    touch ${prefix}_outcomes_mqc.json
    """
}
