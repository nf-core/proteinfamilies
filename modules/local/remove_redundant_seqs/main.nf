process REMOVE_REDUNDANT_SEQS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/eb3700531c7ec639f59f084ab64c05e881d654dcf829db163539f2f0b095e09d/data' :
        'community.wave.seqera.io/library/biopython:1.84--3318633dad0031e7' }"

    input:
    tuple val(meta) , path(clustering)
    tuple val(meta2), path(sequences)

    output:
    tuple val(meta), path("${prefix}_reps.faa"), emit: fasta
    path "versions.yml"                        , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """
}
