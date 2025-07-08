# 🧬 Nextflow Viral Consensus Genome Pipeline 

This Nextflow pipeline processes raw viral tiled amplicon sequencing data (paired-end FASTQ) to generate a consensus genome using tools such as BWA-MEM2, Samtools, and iVar.

It automates the entire workflow—from read alignment to primer trimming and consensus generation—and is suitable for analyzing tiled amplicon sequencing data such as SARS-CoV-2.

## 📋 Overview

The pipeline consists of the following steps:

1. Download raw paired-end FASTQ files

- For benchmarking, you can use dataset #5 from the NCBI benchmarking repository.

2. Align reads to reference genome

- Using BWA-MEM2 for alignment

- BAM file is sorted using Samtools

3. Primer trimming

- Performed with iVar using a BED file containing primer coordinates

4. Consensus genome generation

- With iVar

- Output: consensus genome in FASTA format

## 🛠️ Prerequisites

Before running the pipeline, ensure the following are installed:

- Nextflow
- Docker (used to containerize BWA-MEM2, Samtools, and iVar)

Required Input Files (can be found in the data folder) or you can aquire them from your preferred source:

- Reference genome in FASTA format
- BED file with primer coordinates
- Raw paired-end FASTQ files (_R1.fastq.gz and _R2.fastq.gz)

## ⚠️ Notes & Troubleshooting
Resource tuning: For large datasets, increase memory/CPU allocation in nextflow.config

File paths: Make sure all input files and directories are correctly specified and accessible

Docker: Ensure Docker is running before executing the pipeline

## 📁 Directory Structure Example

    project/
    ├── data/
    │   └── fastq/
    │       ├── sample1_R1.fastq.gz
    │       └── sample1_R2.fastq.gz
    │   └── reference/
    │       ├── nCoV-2019.reference.fasta
    │       ├── sars-cov-2-nonvoivoc.tsv
    │       └── nCoV-2019.bed
    ├── results/
    ├── README.md
    ├── nextflow.config
    ├── .gitignore
    └── main.nf    

## 📄 License
This project is licensed under the MIT License.