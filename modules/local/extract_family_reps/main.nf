process EXTRACT_FAMILY_REPS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8372f6241b480332d91bc00a88ec8c72c8f7fcc9994177a5dd67a07007cd6e32/data' :
        'community.wave.seqera.io/library/biopython:1.85--6f761292fa9881b4' }"

    input:
    tuple val(meta), path(aln, stageAs: "aln/*")

    output:
    tuple val(meta), path("${prefix}_reps.faa")    , emit: fasta
    tuple val(meta), path("${prefix}_meta_mqc.csv"), emit: map
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_family_reps.py \\
        --full_msa_folder aln \\
        --metadata ${prefix}_meta_mqc.csv \\
        --out_fasta ${prefix}_reps.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_reps.faa
    touch ${prefix}_meta_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """
}
