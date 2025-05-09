#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.bam_dir     = "${workflow.projectDir}/results/aligned"
params.bed_file    = "${workflow.projectDir}/data/nCoV-2019.bed"
params.out_dir   = "${workflow.projectDir}/results"

process DEBUG_BAM {
    tag "$sample_id"

    publishDir "debug_outputs", mode: 'copy'  // <-- Add this line

    container = 'community.wave.seqera.io/library/ivar:1.4.4--89c11d667b4e27d1'

    input:
    tuple val(sample_id), file(bam)

    output:
    file("${sample_id}.bam_debug.txt")

    script:
    """
    echo "=== BAM Header for ${bam} ===" > ${sample_id}.bam_debug.txt
    samtools view -H ${bam} >> ${sample_id}.bam_debug.txt

    echo "\\n=== First 10 Alignments ===" >> ${sample_id}.bam_debug.txt
    samtools view ${bam} | head -n 10 >> ${sample_id}.bam_debug.txt

    echo "\\n=== Contig Names from BAM Header ===" >> ${sample_id}.bam_debug.txt
    samtools view -H ${bam} | grep '^@SQ' | cut -f2 | sed 's/SN://' >> ${sample_id}.bam_debug.txt

    echo "\\n=== Contig Names from BED File ===" >> ${sample_id}.bam_debug.txt
    cut -f1 ${params.bed_file} | sort | uniq >> ${sample_id}.bam_debug.txt

    echo "\\n=== If these sets don't match exactly, ivar will fail to find primers ===" >> ${sample_id}.bam_debug.txt
    """
}


workflow {
    Channel
        .fromPath("${params.bam_dir}/*.sorted.bam")
        .map { bam_file -> 
            def sample_id = bam_file.getBaseName().replace('.sorted', '')
            tuple(sample_id, bam_file)
        }
        .set { bam_files }

    DEBUG_BAM(bam_files)
}
