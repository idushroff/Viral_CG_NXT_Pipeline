#!/usr/bin/env nextflow

/* ------------------------------------------------------------
 UPDATE
*/ 

nextflow.enable.dsl=2

params.ref_fasta = "${workflow.projectDir}/data/nCoV-2019.reference.fasta"
params.bed_file  = "${workflow.projectDir}/data/nCoV-2019.bed"
params.out_dir   = "${workflow.projectDir}/results"

params.ivar_min_qual    = 20
params.ivar_min_depth   = 30
params.ivar_freq_thresh = 0.6


workflow {
    Channel
        .fromFilePairs('data/fastq/SRR*_{1,2}.fastq.gz')
        .set { read_pairs }

    // read_pairs.view { "✅ Found pair: ${it}" }

    ALIGN_AND_SORT(read_pairs)
        .set { sorted_bam }

    TRIM_IVAR(sorted_bam)
        .set { trimmed_bam }

    CONSENSUS_IVAR(trimmed_bam)
    .set{consensus_fasta}
    
    NEXTCLADE_QC(consensus_fasta)
}





process NEXTCLADE_QC {
    tag "$sample_id"
    container = 'ghcr.io/nextstrain/nextclade:3.0.0'
    publishDir "${params.out_dir}/nextclade", mode: 'copy'

    input:
    tuple val(sample_id), file(fasta)

    output:
    tuple val(sample_id), file("${sample_id}_nextclade.tsv"), file("${sample_id}_aligned.fasta")

    script:
    """
    nextclade run \\
        --input-fasta ${fasta} \\
        --input-dataset "https://github.com/nextstrain/nextclade_data/releases/latest/download/sars-cov-2" \\
        --output-tsv ${sample_id}_nextclade.tsv \\
        --output-fasta ${sample_id}_aligned.fasta
    """
}