#!/usr/bin/env nextflow

/*
 * FIXES:
 * - Correct handling of the BED file by using a split channel to pair with all BAMs
 * - Ensure multiple consensus FASTAs are produced
 * - Minor clarity/structure improvements
 */

nextflow.enable.dsl = 2

// Parameters
params.ref_fasta          = "${workflow.projectDir}/data/nCoV-2019.reference.fasta"
params.bed_file           = "${workflow.projectDir}/data/nCoV-2019.bed"
params.out_dir            = "${workflow.projectDir}/results"

params.ivar_min_qual      = 20
params.ivar_min_depth     = 30
params.ivar_freq_thresh   = 0.6

workflow {

    // Pair FASTQ files
    Channel
        .fromFilePairs('data/fastq/SRR*_{1,2}.fastq.gz')
        .set { read_pairs }

    // Workflow execution
    ALIGN_AND_SORT(read_pairs)
        .set { sorted_bams }

    bed_channel = Channel.fromPath(params.bed_file)

    TRIM_IVAR(sorted_bams.combine(bed_channel))
        .set { trimmed_bams }

    CONSENSUS_IVAR(trimmed_bams)
        .set{consensus_fasta}
    
    NEXTCLADE_QC(consensus_fasta)
}

// Process 1: Align and sort reads
process ALIGN_AND_SORT {
    tag "$sample_id"
    container = 'ghcr.io/tgen/containers/bwa_mem2_samtools:2.2.1-23080315'

    input:
    tuple val(sample_id), file(reads)

    output:
    tuple val(sample_id), file("${sample_id}.sorted.bam")

    script:
    """
    zcat ${reads[0]} > ${sample_id}_1.fastq
    zcat ${reads[1]} > ${sample_id}_2.fastq

    bwa-mem2 index ${params.ref_fasta}
    bwa-mem2 mem -t 4 ${params.ref_fasta} ${sample_id}_1.fastq ${sample_id}_2.fastq | \
        samtools view -bS - | samtools sort -o ${sample_id}.sorted.bam
    """
}

// Process 2: Trim primers using iVar
process TRIM_IVAR {
    tag "$sample_id"
    container = 'community.wave.seqera.io/library/ivar:1.4.4--89c11d667b4e27d1'

    input:
    tuple val(sample_id), file(bam), file(bed_file)

    output:
    tuple val(sample_id), file("${sample_id}.trimmed.bam")

    script:
    """
    samtools index ${bam}
    ivar trim -i ${bam} -b ${bed_file} -p ${sample_id}.trimmed -m ${params.ivar_min_depth}
    """
}

// Process 3: Generate consensus using iVar
process CONSENSUS_IVAR {
    tag "$sample_id"
    publishDir "${params.out_dir}/consensus", mode: 'copy'
    container = 'community.wave.seqera.io/library/ivar:1.4.4--89c11d667b4e27d1'

    input:
    tuple val(sample_id), file(trimmed_bam)

    output:
    file("${sample_id}.consensus.fasta")

    script:
    """
    samtools mpileup -A -d 0 --reference ${params.ref_fasta} ${trimmed_bam} | \
        ivar consensus -p ${sample_id}.consensus -q ${params.ivar_min_qual} -t ${params.ivar_freq_thresh} -m ${params.ivar_min_depth}
    
    mv ${sample_id}.consensus.fa ${sample_id}.consensus.fasta
    """
}

// Process 4: Nextclade QC
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