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
# Output: set of files .bz2 and _profile.txt with the MetaPhlAn profiles (identified taxa and their abundances)
# To merge MetaPhlAn profiles from multiple samples 
merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt


# HUMAnN 3.0 is the next iteration of HUMAnN, the HMP Unified Metabolic Analysis Network. HUMAnN is a method for efficiently and accurately profiling the abundance of microbial metabolic pathways and other molecular functions from metagenomic or metatranscriptomic sequencing data.
# https://huttenhower.sph.harvard.edu/humann/
# How to run HUMAnN 3.0
humann -i File_kneaddata_paired_R1.fastq -i File_kneaddata_paired_R2.fastq -o /home/soporte/Escritorio/RNAseqAno2023/microbiome/humann_results --nucleotide-database ~/humann/databases/chocophlan --protein-database ~/humann/databases/uniref --threads 35
#The results will be three main output files for each input file named:
$SAMPLE_genefamilies.tsv
$SAMPLE_pathabundance.tsv
$SAMPLE_pathcoverage.tsv

# Normalize the $SAMPLE_genefamilies.tsv abundance output files
#   Basic usage: 
#	$ humann_renorm_table --input $TABLE --units $CHOICE --output $TABLE2
#   $TABLE = gene/pathway table (tsv format)
#   $CHOICE = "relab" (relative abundance) or "cpm" (copies per million)
#   $TABLE2 = normalized gene/pathway table
#   Run with -h to see additional command line options
#!/bin/bash
for f in *_genefamilies.tsv
	do echo "Processing genefamilies $f file.."
humann_renorm_table --input "$f" --output "$f"_genefamilies_relab.tsv --units relab
done
for f in *_pathabundance.tsv
	do echo "Processing pathabundaces $f file.."
humann_renorm_table --input "$f" --output "$f"_pathabundance_relab.tsv --units relab
done

#El output es el nombre del archivo joined, --file_name es el tipo de file a joinear
$ humann_join_tables --input ~/Escritorio/RNAseqAno2023/microbiome/humann_results --output all_normalizedgenefamilies_relab.tsv --file_name _genefamilies_relab.tsv
$ humann_join_tables --input ~/Escritorio/RNAseqAno2023/microbiome/humann_results --output all_pathcoverage.tsv --file_name _pathcoverage.tsv
$ humann_join_tables --input ~/Escritorio/RNAseqAno2023/microbiome/humann_results --output all_normalizedpathabundance_relab.tsv --file_name _pathabundance_relab.tsv
