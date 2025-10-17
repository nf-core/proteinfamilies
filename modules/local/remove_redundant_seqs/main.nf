process REMOVE_REDUNDANT_SEQS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8372f6241b480332d91bc00a88ec8c72c8f7fcc9994177a5dd67a07007cd6e32/data' :
        'community.wave.seqera.io/library/biopython:1.85--6f761292fa9881b4' }"

    input:
    tuple val(meta) , path(clustering)
    tuple val(meta2), path(sequences)

    output:
    tuple val(meta), path("${prefix}.faa"), emit: fasta
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def is_compressed = sequences.getName().endsWith(".gz") ? true : false
    def fasta_name    = sequences.name.replace(".gz", "")
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $sequences > $fasta_name
    fi

    remove_redundant_seqs.py \\
        --clustering ${clustering} \\
        --sequences ${fasta_name} \\
        --out_fasta ${prefix}.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """
}
