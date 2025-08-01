process SEQKIT_STATS_TO_MQC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("${prefix}_stats_mqc.csv"), emit: mqc
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<EOF > ${prefix}_stats_mqc.csv
    # id: "seqkit_stats_summary"
    # section_name: "SeqKit Statistics Summary"
    # description: "Statistics about the input protein FASTA file from SeqKit."
    # format: "csv"
    # plot_type: "table"
    File,File Format,Sequence Type,Number of Sequences,Total Length (aa),Minimum Length,Average Length,Maximum Length,Length Q1 (25th percentile),Length Q2 (Median),Length Q3 (75th percentile),Total Gaps,N50 Length,N50 Sequence Count
    EOF

    # Skip header line and extract relevant columns (1–14) from tab-separated input
    tail -n +2 ${stats} | cut -f1-14 --output-delimiter=, >> ${prefix}_stats_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_stats_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
