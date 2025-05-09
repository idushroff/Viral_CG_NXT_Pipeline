#!/usr/bin/env nextflow

/*
 FIXED - hardcoded paths to relative paths
 NOTE - used quay.io instead of seqera containers
 NOTE - bwa-mem2 and samtools used in the same process
 FIXED - originally the bwa-mem2 was bwa (this was my mistake)
 CHANGE - decompress reads can be done in the same process as alignment
*/

nextflow.enable.dsl=2

// Define input parameters: reference genome, BED file, and output directory
params.fastq_dir   = 'data/fastq'
params.ref_fasta   = 'data/nCoV-2019.reference.fasta'
params.bed_file    = 'data/nCoV-2019.bed'
params.out_dir     = 'results'

// ------------------------------------------------------------
// Decompress paired-end FASTQ files
// ------------------------------------------------------------
process decompress_reads {
    container 'quay.io/biocontainers/pigz:2.6--h27826a3_0'

    input:
    tuple path(read1_gz), path(read2_gz)

    output:
    tuple path("${read1_gz.baseName}.fastq"), path("${read2_gz.baseName}.fastq")

    script:
    """
    gzip -d -c ${read1_gz} > ${read1_gz.baseName}.fastq
    gzip -d -c ${read2_gz} > ${read2_gz.baseName}.fastq
    """
}

// ------------------------------------------------------------
// Align paired-end reads with BWA-MEM2 and sort BAM
// ------------------------------------------------------------
process align_reads {
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_7'

    input:
    tuple path(read1), path(read2)
    path reference_genome

    output:
    path "aligned.bam"

    script:
    """
    bwa-mem2 mem -t 4 ${reference_genome} ${read1} ${read2} | \
    samtools view -bS - | samtools sort -o aligned.bam
    """
}

// ------------------------------------------------------------
// Trim primers using iVar
// ------------------------------------------------------------
process trim_reads {
    container 'quay.io/biocontainers/ivar:1.4.2--h9ee0642_0'

    input:
    path bam_file
    path bed_file

    output:
    path "trimmed.bam"

    script:
    """
    ivar trim -i ${bam_file} -b ${bed_file} -p trimmed
    """
}

// ------------------------------------------------------------
// Generate consensus genome using iVar
// ------------------------------------------------------------
process generate_consensus {
    container 'quay.io/biocontainers/ivar:1.4.2--h9ee0642_0'

    input:
    path trimmed_bam
    path bed_file

    output:
    path "consensus.fasta"

    script:
    """
    ivar consensus -i ${trimmed_bam} -b ${bed_file} -p consensus
    """
}

// ------------------------------------------------------------
// Main workflow logic
// ------------------------------------------------------------
workflow {

    // Match R1 and R2 files based on naming pattern
    Channel
        .fromFilePairs("${params.fastq_dir}/*_{R1,R2}.fastq.gz", flat: true)
        .set { paired_reads }

    // Run processes in sequence
    paired_reads
        | decompress_reads
        | align_reads(file(params.reference_genome))
        | trim_reads(file(params.bed_file))
        | generate_consensus(file(params.bed_file))
}
