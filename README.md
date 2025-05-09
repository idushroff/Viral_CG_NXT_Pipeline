Nextflow Viral Consensus Genome Pipeline

This pipeline is designed to process viral sequencing data, from raw paired-end FASTQ files to a consensus genome using tools such as BWA, iVar, and Samtools. It demonstrates an automated workflow with Nextflow, starting from sequence alignment to primer trimming and final consensus generation. This process can be particularly useful for analyzing viral tiled amplicon sequencing data.

Overview
This Nextflow pipeline consists of the following steps:

Download paired FASTQ files from NCBI (use #5 dataset from the benchmarking repository).

Align the reads to the reference genome using BWA-MEM2 and sort the resulting BAM file using Samtools.

Trim primers from the aligned BAM file using iVar and a BED file for primer coordinates.

Generate a consensus genome using iVar and output a FASTA file with the consensus sequence.

Prerequisites
Before running the pipeline, make sure the following software and tools are available:

Nextflow: Nextflow Installation Guide

Docker (for containerized execution of tools like BWA-MEM2, Samtools, and iVar): Docker Installation Guide

Reference Files:

Download the reference FASTA file from the NCBI or your preferred source.

Download the BED file that contains primer coordinates.

Obtain raw paired-end FASTQ files (can be downloaded from NCBI or another source).

Troubleshooting and Notes
Memory/CPU Requirements: Depending on your dataset size, you might need to adjust memory and CPU allocation for the processes. This can be done in your config file.

File Paths: Ensure that your file paths (for FASTQ, reference, and BED files) are correct and accessible by the pipeline.

License
This pipeline is licensed under the MIT License.