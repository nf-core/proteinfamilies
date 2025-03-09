process CHUNK_CLUSTERS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/43/438649357e4dcaab676b1bff95e3aace1decb36a658d6257869a641155867e0c/data' :
        'community.wave.seqera.io/library/pip_pyfastx:c1d255a74c4291f8' }"

    input:
    tuple val(meta) , path(clustering)
    tuple val(meta2), path(sequences)
    val(size_threshold)

    output:
    tuple val(meta), path("chunked_fasta/*"), emit: fasta_chunks
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def is_compressed = sequences.getName().endsWith(".gz") ? true : false
    def fasta_name    = sequences.name.replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $sequences > $fasta_name
    fi

    chunk_clusters.py \\
        --clustering ${clustering} \\
        --sequences ${fasta_name} \\
        --threshold ${size_threshold} \\
        --out_folder chunked_fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import importlib.metadata; print(importlib.metadata.version('biopython'))")
    END_VERSIONS
    """
}
