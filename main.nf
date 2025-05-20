#!/usr/bin/env nextflow

/*
 * NOTES: need to figure out how to add samtools flagstat in process 1. 
 latest update: updated with Check Coverage process 
 */

nextflow.enable.dsl = 2 // Enable Nextflow DSL2 syntax

// Parameters
params.ref_fasta          = "${workflow.projectDir}/data/nCoV-2019.reference.fasta"
params.out_dir            = "${workflow.projectDir}/results"

params.ivar_min_qual      = 20  //-q
params.ivar_min_depth     = 30  //-m
params.ivar_freq_thresh   = 0.6 //-t



/*
 * Main workflow definition
 */
workflow {

    // Pair FASTQ files
    Channel
        .fromFilePairs('data/fastq/SRR*_{1,2}.fastq.gz')
        .set { read_pairs }

    // Workflow execution
    ALIGN_AND_SORT(read_pairs) // Process 1: Align reads using BWA-MEM2 & Convert SAM to sorted BAM using SAMtools
        .set { sorted_bams }
    

    bed_channel = Channel.fromPath('data/nCoV-2019.bed')

    TRIM_IVAR(sorted_bams.combine(bed_channel)) // Process 2: Trim primers from BAM files using iVar and BED file
        .set { trimmed_bams }

    CHECK_COVERAGE(trimmed_bams) // Process 3:
        .set{coverage_reports}

    ref_file = Channel.fromPath('data/nCoV-2019.reference.fasta')

    CONSENSUS_IVAR(trimmed_bams.combine(ref_file)) // Process 4: Generate consensus genome from trimmed BAM using iVar
        .set{consensus_fasta}
    
}

// Process 1: Align paired-end reads to the reference genome using BWA-MEM2 and convert SAM to sorted BAM using SAMtools
process ALIGN_AND_SORT {
    
    tag "$sample_id" // Tags each job with sample ID for easier tracking
    
    container = 'ghcr.io/tgen/containers/bwa_mem2_samtools:2.2.1-23080315' // Use container with both BWA-MEM2 and SAMtools

    input:
    tuple val(sample_id), file(reads) // Paired FASTQ files

    output:
    tuple val(sample_id), file("${sample_id}.sorted.bam") // Output sorted BAM file

    // Decompress the gzipped FASTQ files
    // Index the reference genome (needed for BWA)
    // Perform alignment and pipe directly to SAMtools for BAM conversion and sorting
    // LAST line along with the modified output for this process gives me a file for each sorted.bam file tell me how well they are aligned to the reference. You can find these in the work folder. 
    
    script:
    """
    zcat ${reads[0]} > ${sample_id}_1.fastq
    zcat ${reads[1]} > ${sample_id}_2.fastq

    bwa-mem2 index ${params.ref_fasta}

    bwa-mem2 mem -t 4 ${params.ref_fasta} ${sample_id}_1.fastq ${sample_id}_2.fastq | \
        samtools view -bS - | samtools sort -o ${sample_id}.sorted.bam
    """
}

// Process 2: Trim primers from BAM file using iVar and BED file
process TRIM_IVAR {
    tag "$sample_id"
    container = 'community.wave.seqera.io/library/ivar:1.4.4--89c11d667b4e27d1'

    input:
    tuple val(sample_id), file(bam), file(bed_file) // Input sorted BAM combined to the Bed_file

    output:
    tuple val(sample_id), file("${sample_id}.trimmed.sorted.bam")  // Output trimmed BAM

    // Index the BAM file
    // Trim primer sequences using BED coordinates
    
    script:
    """
    samtools index ${bam}
    ivar trim -i ${bam} -b ${bed_file} -p ${sample_id}.trimmed -m ${params.ivar_min_depth}
    
    samtools sort -o ${sample_id}.trimmed.sorted.bam ${sample_id}.trimmed.bam
    """
}

// Process 3: 
process CHECK_COVERAGE {
    tag "$sample_id"
    publishDir "${params.out_dir}/coverage", mode: 'copy'
    
    input:
    tuple val(sample_id), file(trimmed_bam) // Input trimmed BAM file

    output:
    tuple val(sample_id), file("${sample_id}.coverage.txt") // Coverage report output
    
    // Run samtools depth on the sorted BAM file
    script:
    """
    samtools depth ${sample_id}.trimmed.sorted.bam > ${sample_id}.coverage.txt
    """
}


// Process 4: Generate consensus using iVar
process CONSENSUS_IVAR {
    tag "$sample_id"

    // Save output consensus genomes to the output directory
    publishDir "${params.out_dir}/consensus", mode: 'copy'
    container = 'community.wave.seqera.io/library/ivar:1.4.4--89c11d667b4e27d1'

    input:
    tuple val(sample_id), file(trimmed_bam), file(ref_file) // Input trimmed BAM after trimming and coverage report

    output:
    file("${sample_id}.consensus.fasta") // Final consensus FASTA

    script:
    """    
    samtools mpileup -A -d 0 --reference ${ref_file} ${trimmed_bam} | \
        ivar consensus -p ${sample_id}.consensus -q ${params.ivar_min_qual} -t ${params.ivar_freq_thresh} -m ${params.ivar_min_depth} > ${sample_id}.consensus.fa
    
    mv ${sample_id}.consensus.fa ${sample_id}.consensus.fasta
    """
}

