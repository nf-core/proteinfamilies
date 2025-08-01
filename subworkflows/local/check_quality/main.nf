/*
    AMINO ACID SEQUENCE QUALITY CHECK
*/

include { SEQKIT_STATS        } from '../../../modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS_TO_MQC } from '../../../modules/local/seqkit_stats_to_mqc/main'

workflow CHECK_QUALITY {
    take:
    fasta // tuple val(meta), path(fasta)
    
    main:
    ch_versions = Channel.empty()

    SEQKIT_STATS( fasta )
    ch_versions = ch_versions.mix( SEQKIT_STATS.out.versions )

    SEQKIT_STATS_TO_MQC( SEQKIT_STATS.out.stats )
    ch_versions = ch_versions.mix( SEQKIT_STATS_TO_MQC.out.versions )

    emit:
    versions         = ch_versions
    seqkit_stats_mqc = SEQKIT_STATS_TO_MQC.out.mqc
}
