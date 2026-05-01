process CALCULATE_CLUSTER_DISTRIBUTION {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(clusters)

    output:
    tuple val(meta), path("${prefix}_clustering_distribution_mqc.csv"), emit: mqc
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    num_sequences=\$(wc -l < "${clusters}")

    cat <<EOF > ${prefix}_clustering_distribution_mqc.csv
    # id: "${prefix}_cluster_distribution"
    # section_name: "Initial Clustering Sizes Distribution for Sample ${prefix}"
    # description: "A total of <b> \$num_sequences </b> input sequences were processed, with seed clusters identified using the MMseqs suite. The table below provides a breakdown of cluster sizes, showing the number of times each unique cluster size was observed."
    # format: "csv"
    # plot_type: "table"
    Id,Cluster Size,Number of Clusters
    EOF

    awk '{count[\$1]++} END {for (c in count) size[count[c]]++} END {for (s in size) print s "," s "," size[s]}' ${clusters} | sort -t, -k3,3nr --parallel=${task.cpus} >> ${prefix}_clustering_distribution_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_clustering_distribution_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
