process GENERATE_FASTA_INDEX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hb6cb901_4' :
        'quay.io/biocontainers/hmmer:3.4--hb6cb901_4' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta}.ssi"), emit: ssi
    tuple val("${task.process}"), val('easel'), eval("esl-sfetch -h | sed '2!d;s/^# Easel *//;s/ .*//'"), emit: versions_easel, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // No ext.prefix: Easel writes the index beside its input and records the sequence
    // file name inside it, so the index only resolves for a FASTA of this exact name.
    """
    esl-sfetch --index ${fasta}
    """

    stub:
    """
    touch ${fasta}.ssi
    """
}
