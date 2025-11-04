process CREATE_PROTEINFOLD_SAMPLESHEET {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(split_files)

    output:
    tuple val(meta), path("samplesheet_${prefix}.csv"), emit: samplesheet
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "id,fasta" > samplesheet_${prefix}.csv

    for fasta_file in *.fa; do
        if [ -f "\$fasta_file" ]; then
            seq_id=\$(head -n1 "\$fasta_file" | sed 's/^>//' )
            # Get the full path to the original file by following symlinks
            full_path=\$(readlink -f "\$fasta_file")

            echo "\$seq_id,\$full_path" >> samplesheet_${prefix}.csv
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "id,fasta" > samplesheet_${prefix}.csv
    echo "test_seq1,/path/to/test_seq1.fa" >> samplesheet_${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
