#!/bin/bash
# Preprocesing of metatrasncriptomics data

# KneadData is a tool designed to perform quality control on metagenomic sequencing data, especially data from microbiome experiments including decontamination of human reads
# https://huttenhower.sph.harvard.edu/kneaddata/ 
# How to run KneadData on Paired-End Inputs:
# kneaddata --input1 $INPUT1 --input2 $INPUT2 --reference-db $DATABASE --output $OUTPUT_DIR
kneaddata --input File_R1.fastq.gz --input File_R2.fastq.gz --reference-db /home/soporte/kneaddata --output /home/soporte/Escritorio/RNAseqAno2024/microbiome/microbiome --trimmomatic /home/soporte/Trimmomatic-0.39 --bypass-trf --remove-intermediate-output --threads 35

# MetaPhlAn is a computational tool for profiling the composition of microbial communities (Bacteria, Archaea and Eukaryotes) from metagenomic shotgun sequencing data (i.e. not 16S) with species-level. 
# https://huttenhower.sph.harvard.edu/metaphlan/
# How to run MetaPhlAn
metaphlan  File_kneaddata_paired_R1.fastq,File_kneaddata_paired_R2.fastq --input_type fastq --bowtie2out FileName.bz2 --add_viruses --nproc 35 > FileName_profile.txt

# HUMAnN 3.0 is the next iteration of HUMAnN, the HMP Unified Metabolic Analysis Network. HUMAnN is a method for efficiently and accurately profiling the abundance of microbial metabolic pathways and other molecular functions from metagenomic or metatranscriptomic sequencing data.
# https://huttenhower.sph.harvard.edu/humann/
# How to run HUMAnN 3.0
humann -i File_kneaddata_paired_R1.fastq -i File_kneaddata_paired_R2.fastq -o /home/soporte/Escritorio/RNAseqAno2023/microbiome/humann_results --nucleotide-database ~/humann/databases/chocophlan --protein-database ~/humann/databases/uniref --threads 35
