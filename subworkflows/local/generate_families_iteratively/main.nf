/*
    ITERATIVE FAMILY MODEL GENERATION

    Alternative to GENERATE_FAMILIES that hands a whole chunk of clusters to mgnifam, which
    loops build HMM -> recruit members -> realign per cluster until each family converges or
    is discarded, and so emits many families per task instead of one.

    Its outputs are reshaped here into the same per-family channels GENERATE_FAMILIES emits:
      seed_msa — Stockholm seeds reformatted to aligned FASTA, the form MERGE_SEEDS reads
      full_msa — left as gzipped Stockholm, as HMMER_HMMALIGN produces in the standard path
      fasta    — full MSA members with their gaps removed, mgnifam emits no such file itself
      hmm      — as produced
*/

include { GUNZIP                                          } from '../../../modules/nf-core/gunzip/main'
include { HMMER_ESLSFETCHINDEX                            } from '../../../modules/nf-core/hmmer/eslsfetchindex/main'
include { MGNIFAM_GENERATEFAMILIES                        } from '../../../modules/nf-core/mgnifam/generatefamilies/main'
include { HHSUITE_REFORMAT as HHSUITE_REFORMAT_SEED       } from '../../../modules/nf-core/hhsuite/reformat/main'
include { HHSUITE_REFORMAT as HHSUITE_REFORMAT_FULL       } from '../../../modules/nf-core/hhsuite/reformat/main'
include { SEQKIT_SEQ as SEQKIT_SEQ_ITERATIVE_MSA_TO_FASTA } from '../../../modules/nf-core/seqkit/seq/main'

workflow GENERATE_FAMILIES_ITERATIVELY {
    take:
    sequences       // tuple val(meta), path(fasta)
    clusters_chunks // tuple val(meta), path(tsv)

    main:
    // Easel cannot seek within a gzip stream, so mgnifam needs the sequence pool uncompressed.
    ch_branched_sequences = sequences
        .branch { _meta, fasta ->
            compressed  : fasta.name.endsWith('.gz')
            uncompressed: true
        }

    GUNZIP( ch_branched_sequences.compressed )

    ch_sequences = ch_branched_sequences.uncompressed.mix( GUNZIP.out.gunzip )

    // Indexed once per sample: every chunk of a sample searches the same sequence pool, and
    // mgnifam would otherwise rebuild an identical index per chunk.
    HMMER_ESLSFETCHINDEX( ch_sequences )

    ch_indexed_sequences = ch_sequences.join( HMMER_ESLSFETCHINDEX.out.ssi )

    // Combine on a chunk-free [id] key so each chunk matches the full sample sequence pool;
    // the original meta rides along as an extra element and the key is dropped after.
    ch_input_for_mgnifam = clusters_chunks
        .map { meta, tsv -> [ [id: meta.id], meta, tsv ] }
        .combine(ch_indexed_sequences, by: 0)
        .map { _id, meta, tsv, seqs, ssi -> [ meta, tsv, seqs, ssi ] }

    MGNIFAM_GENERATEFAMILIES( ch_input_for_mgnifam )

    // One task emits every family of its chunk. Family files are named <prefix>_<n>, with
    // prefix set per chunk in modules.config, so the sample id splits each name into the
    // per-family chunk key the downstream subworkflows expect.
    ch_seed_msa = MGNIFAM_GENERATEFAMILIES.out.seed_msa.transpose().map { meta, file -> [ [id: meta.id, chunk: file.simpleName.split("${meta.id}_", 2)[1]], file ] }
    ch_full_msa = MGNIFAM_GENERATEFAMILIES.out.full_msa.transpose().map { meta, file -> [ [id: meta.id, chunk: file.simpleName.split("${meta.id}_", 2)[1]], file ] }
    ch_hmm      = MGNIFAM_GENERATEFAMILIES.out.hmm.transpose().map { meta, file -> [ [id: meta.id, chunk: file.simpleName.split("${meta.id}_", 2)[1]], file ] }

    HHSUITE_REFORMAT_SEED( ch_seed_msa, "sto", "fas" )

    // mgnifam has no per-family member FASTA to emit, so it is recovered from the full MSA
    HHSUITE_REFORMAT_FULL( ch_full_msa, "sto", "fas" )

    SEQKIT_SEQ_ITERATIVE_MSA_TO_FASTA( HHSUITE_REFORMAT_FULL.out.msa )

    emit:
    seed_msa = HHSUITE_REFORMAT_SEED.out.msa
    full_msa = ch_full_msa
    fasta    = SEQKIT_SEQ_ITERATIVE_MSA_TO_FASTA.out.fastx
    hmm      = ch_hmm
}
