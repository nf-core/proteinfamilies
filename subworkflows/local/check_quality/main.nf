/*
    AMINO ACID SEQUENCE QUALITY CHECK
*/

include { SEQFU_STATS } from '../../../modules/nf-core/seqfu/stats/main'

workflow CHECK_QUALITY {
    take:
    fasta // tuple val(meta), path(fasta)

    main:
    ch_multiqc_files = Channel.empty()
    ch_versions      = Channel.empty()

    SEQFU_STATS( fasta )
    ch_multiqc_files = ch_multiqc_files.mix( SEQFU_STATS.out.multiqc )
    ch_versions      = ch_versions.mix( SEQFU_STATS.out.versions )

    emit:
    versions         = ch_versions
    multiqc_files    = ch_multiqc_files
}
