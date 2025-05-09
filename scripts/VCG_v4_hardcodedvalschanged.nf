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

// Enable Nextflow DSL2 syntax
nextflow.enable.dsl=2

// Define input parameters: reference genome, BED file, and output directory

params.ref_fasta = "${workflow.projectDir}/data/nCoV-2019.reference.fasta"
params.bed_file  = "${workflow.projectDir}/data/nCoV-2019.bed"
params.out_dir   = "${workflow.projectDir}/results"

// New params for iVar settings
params.ivar_min_qual    = 20
params.ivar_min_depth   = 30
params.ivar_freq_thresh = 0.6

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

    // Step 1&2: Align reads using BWA-MEM2 & Convert SAM to sorted BAM using SAMtools
    ALIGN_AND_SORT(read_pairs)
        .set { sorted_bam }

    // Step 3: Trim primers from BAM files using iVar and BED file
    TRIM_IVAR(sorted_bam)
        .set { trimmed_bam }

    // Step 4: Generate consensus genome from trimmed BAM using iVar
    CONSENSUS_IVAR(trimmed_bam)
}

/*
 * Step 1 & 2 Combined: Align paired-end reads to the reference genome using BWA-MEM2
 * and convert SAM to sorted BAM using SAMtools
 */
process ALIGN_AND_SORT {
    tag "$sample_id"  // Tags each job with sample ID for easier tracking

    // Uncomment the line below if you would like to see these files in the results directory - WARNING: requires a lot of disk space
    // publishDir "${params.out_dir}/trimmed", mode: 'copy' 

    // Use container with both BWA-MEM2 and SAMtools
    container = 'ghcr.io/tgen/containers/bwa_mem2_samtools:2.2.1-23080315'

    input:
    tuple val(sample_id), file(reads)  // Paired FASTQ files

    output:
    tuple val(sample_id), file("${sample_id}.sorted.bam")  // Output sorted BAM file

    script:
    """

    // Decompress the gzipped FASTQ files
    zcat ${reads[0]} > ${sample_id}_1.fastq
    zcat ${reads[1]} > ${sample_id}_2.fastq

    bash -c '

    // Index the reference genome (needed for BWA)
    bwa-mem2 index ${params.ref_fasta}
    
    // Perform alignment and pipe directly to SAMtools for BAM conversion and sorting
    bwa-mem2 mem -t 4 ${params.ref_fasta} ${reads[0]}_1.fastq ${reads[1]}_2.fastq | samtools view -bS - | samtools sort -o ${sample_id}.sorted.bam
    """
}

/*
 * Step 3: Trim primers from BAM file using iVar and BED file
 */
process TRIM_IVAR {
    tag "$sample_id"
    
    //uncomment the line below if you would like to see these files in the results directory - WARNING: requires a lot of disk space    
    // publishDir "${params.out_dir}/trimmed", mode: 'copy' 

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
    ivar trim -i ${bam} -b ${params.bed_file} -p ${sample_id}.trimmed -m ${params.ivar_min_depth}
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
    samtools mpileup -A -d 0 --reference ${params.ref_fasta} ${trimmed_bam} | \
        ivar consensus -p ${sample_id}.consensus -q ${params.ivar_min_qual} -t ${params.ivar_freq_thresh} -m ${params.ivar_min_depth}

    // Rename output file for clarity
    mv ${sample_id}.consensus.fa ${sample_id}.consensus.fasta
    """
}