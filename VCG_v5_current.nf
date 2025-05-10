#!/usr/bin/env nextflow

/* ------------------------------------------------------------
 FIXED - hardcoded paths to relative paths
 FIXED - switched from of seqera container to github container with both 
 bwa and samtools to improve code efficiency
 NOTE - seqera containers is free therefore I've left it as is for IVAR
 FIXED - hard coded values in my processes (e.g., the quality score threshold
 in ivar consensus). Changed these to params (similarly to what I've done for
 ref_fast, bed_file, etc). 
 -- Added params.ivar_min_qual, params.ivar_min_depth, and params.ivar_freq_thresh at the top.
 -- Updated TRIM_IVAR and CONSENSUS_IVAR to use these new parameters.
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

    read_pairs.view { "✅ Found pair: ${it}" }

    ALIGN_AND_SORT(read_pairs)
        .set { sorted_bam }

    TRIM_IVAR(sorted_bam)
        .set { trimmed_bam }

    CONSENSUS_IVAR(trimmed_bam)
}


process ALIGN_AND_SORT {
    tag "$sample_id"

    // publishDir "${params.out_dir}/aligned", mode: 'copy'

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


/*
 * Step 3: Trim primers from BAM file using iVar and BED file
 */
process TRIM_IVAR {
    tag "$sample_id"
    
    container = 'community.wave.seqera.io/library/ivar:1.4.4--89c11d667b4e27d1'
    
    input:
    tuple val(sample_id), file(bam)  // Input sorted BAM

    output:
    tuple val(sample_id), file("${sample_id}.trimmed.bam")  // Output trimmed BAM

    script:
    """
    # Index the BAM file
    samtools index ${bam}

    # Trim primer sequences using BED coordinates
    ivar trim -i ${bam} -b ${params.bed_file} -p ${sample_id}.trimmed -m ${params.ivar_min_depth}
    """
}

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
