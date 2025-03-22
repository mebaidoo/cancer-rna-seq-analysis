# RNA-Seq Analysis of Yeast (GSE11209) using Python
This repository contains the complete RNA-seq analysis pipeline for Saccharomyces cerevisiae (yeast) using the dataset GSE11209. The analysis involves downloading raw RNA-seq FASTQ files, building a genome index, performing quality control (QC), trimming, and aligning the reads to the genome. Below is an overview of each step in the process.
## Overview
1. **Data Acquisition**:
   - Raw RNA-seq FASTQ files from the **GSE11209** dataset are downloaded using the **Entrez** library (via GSE, GSM, SRA IDs) and the **SRA Toolkit**.
   - Saccharomyces cerevisiae genome and annotation files (R64-1-1 version) are downloaded from Ensembl Fungi for the construction of a genome index.
2. **Building the Genome Index**:
   - The genome FASTA file and GTF annotation file are used to build a genome index using **STAR** aligner.
   - The index allows for efficient and accurate alignment of the RNA-seq reads to the yeast genome.
3. **Quality Control**:
   - Quality control of raw FASTQ files is performed using **FastQC**.
   - The output QC reports are summarized using **MultiQC**.
4. **Read Trimming**:
   - Reads are trimmed using **Trim Galore** to remove low-quality bases and adapter sequences from the RNA-seq reads.
   - Trimmed files are stored in a separate directory.
5. **Read Alignment**:
   - Trimmed reads are aligned to the Saccharomyces cerevisiae genome using **STAR** aligner.
   - Aligned reads are output in **BAM** format, sorted by coordinate.
