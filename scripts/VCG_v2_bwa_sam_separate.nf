#!/usr/bin/env nextflow

/* ------------------------------------------------------------
 NOTE - switched from quay.io instead of seqera containers
 NOTE - bwa-mem2 and samtools had to be separated as not in same seqera container 
 FIXED - decompressed reads now in the same process as alignment
*/

// Enable Nextflow DSL2 syntax
nextflow.enable.dsl=2

// Define input parameters: reference genome, BED file, and output directory
params.ref_fasta   = 'data/nCoV-2019.reference.fasta'
params.bed_file    = 'data/nCoV-2019.bed'
params.out_dir     = 'results'

/*
 * Main workflow definition
 */
workflow {
    // Create a channel of paired FASTQ files
    Channel
        .fromFilePairs('data/fastq/SRR*_{1,2}.fastq.gz')
        .set { read_pairs }

    // Print the paired FASTQ file paths to help debug file matching
    read_pairs.view { "✅ Found pair: ${it}" }

    // Step 1: Align reads using BWA-MEM2
    ALIGN_BWA(read_pairs)
        .set { sam_files }

    // Step 2: Convert SAM to sorted BAM using SAMtools
    SAM_TO_SORTED_BAM(sam_files)
        .set { sorted_bam }

    // Step 3: Trim primers from BAM files using iVar and BED file
    TRIM_IVAR(sorted_bam)
        .set { trimmed_bam }

    // Step 4: Generate consensus genome from trimmed BAM using iVar
    CONSENSUS_IVAR(trimmed_bam)
}

/*
 * Step 1: Align paired-end reads to the reference genome using BWA-MEM2
 */
process ALIGN_BWA {
    tag "$sample_id"  // Tags each job with sample ID for easier tracking

    //uncomment the line below if you would like to see these files in the results directory - WARNING: requires a lot of disk space
    //*publishDir "${params.out_dir}/trimmed", mode: 'copy' 

    // Use BWA-MEM2 container from Seqera community
    container = 'community.wave.seqera.io/library/bwa-mem2:2.2.1--1842774b9b0b4729'
    

    input:
    tuple val(sample_id), file(reads)  // Paired FASTQ files

    output:
    tuple val(sample_id), file("${sample_id}.sam")  // Output SAM file

    script:
    """
    // Decompress the gzipped FASTQ files
    zcat ${reads[0]} > ${sample_id}_1.fastq
    zcat ${reads[1]} > ${sample_id}_2.fastq

    // Index the reference genome (needed for BWA)
    bwa-mem2 index ${params.ref_fasta}

    // Perform alignment
    bwa-mem2 mem -t 4 ${params.ref_fasta} ${sample_id}_1.fastq ${sample_id}_2.fastq > ${sample_id}.sam
    """
}

/*
 * Step 2: Convert SAM to BAM and sort using SAMtools
 */
process SAM_TO_SORTED_BAM {
    tag "$sample_id"

    //uncomment the line below if you would like to see these files in the results directory - WARNING: requires a lot of disk space
    //*publishDir "${params.out_dir}/trimmed", mode: 'copy' 

    container = 'community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c'

    input:
    tuple val(sample_id), file(sam)  // Input SAM file

    output:
    tuple val(sample_id), file("${sample_id}.sorted.bam")  // Output sorted BAM file

    script:
    """
    // Convert SAM to BAM and sort
    samtools view -bS ${sam} | samtools sort -o ${sample_id}.sorted.bam
    """
}

/*
 * Step 3: Trim primers from BAM file using iVar and BED file
 */
process TRIM_IVAR {
    tag "$sample_id"
    
    //uncomment the line below if you would like to see these files in the results directory - WARNING: requires a lot of disk space    
    //*publishDir "${params.out_dir}/trimmed", mode: 'copy' 

    container = 'community.wave.seqera.io/library/ivar:1.4.4--89c11d667b4e27d1'

    input:
    tuple val(sample_id), file(bam)  // Input sorted BAM

    output:
    tuple val(sample_id), file("${sample_id}.trimmed.bam")  // Output trimmed BAM

    script:
    """
    // Index the BAM file
    samtools index ${bam}

    // Trim primer sequences using BED coordinates
    ivar trim -i ${bam} -b ${params.bed_file} -p ${sample_id}.trimmed -m 30
    """
}

/*
 * Step 4: Generate consensus genome using iVar
 */
process CONSENSUS_IVAR {
    tag "$sample_id"

    // Save output consensus genomes to the output directory
    publishDir "${params.out_dir}/consensus", mode: 'copy'

    container = 'community.wave.seqera.io/library/ivar:1.4.4--89c11d667b4e27d1'

    input:
    tuple val(sample_id), file(trimmed_bam)  // Input BAM after trimming

    output:
    file("${sample_id}.consensus.fasta")  // Final consensus FASTA

    script:
    """
    // Generate pileup and create consensus sequence
    samtools mpileup -A -d 0 --reference ${params.ref_fasta} ${trimmed_bam} |
        ivar consensus -p ${sample_id}.consensus -q 20 -t 0.6 -m 30

    // Rename output file for clarity
    mv ${sample_id}.consensus.fa ${sample_id}.consensus.fasta
    """
}
