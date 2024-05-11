# ASCC-transcriptomics 
Bioinformatics code used for Transcriptomic Portrait of ASCC
![Screenshot from 2024-05-10 23-17-21](https://github.com/mabba777/ASCC-transcriptomics/assets/5058918/e5b2dc43-da64-40a3-84e2-00a9930be3e7)


# 0- Data access from Gene Expression Omnibus
The raw data have been submitted to NCBI GEO database with accession number GSE253560.
The preprocessed host transcriptome data and the metatranscriptome are available at the data_matrix folder

# 1- Metatranscriptomics analysis
For metatranscriptomic analysis, the obtained RNA sequencing data were processed using the Biobakery suite of tools: KneadData was used to separate the human and the non-human reads; taxonomic profiling was performed using MetaPhlAn to identify and quantify microbial taxa at species level present in the anal samples.
For determining the relative differential abundance and the multivariable association between subjects’ metadata and microbial features, we used the MaAsLin2 package from the bioBakery suite in R/Bioconductor.
Metatranscriptomic pathway analysis was conducted using the HMP Unified Metabolic Analysis Network 3 (HUMAnN3) pipeline to investigate potential variations in metabolic pathways.

# 2- Trascriptomics analysis of host gene expression

# 3- Immune profiling of preivasive and invise samples

# 4- Mutational analysis of ASCC based on their transcriptome
The preprocessed reads previously used for the transcriptomic analysis were aligned and mapped to the human genome reference GRCh38 using the Subjunc aligner algorithm provided by Rsubread R/Bioconductor package. Subjunc aligner was developed for aligning RNA-seq reads and for the detection of exon-exon junctions at the same time. The Subjunc mapping results (BAM files) were used for genomic variants detection using the exactSNP variant caller algorithm provided by Rsubread package. The VariantAnnotation R/Bioconductor package was subsequently used for SNPs and InDels filtering of the obtained VCF files based on quality (QUAL > 20) and coverage (DP>10) metrics. Identified variants were annotated, filtered and interpreted using OpenCRAVAT and their aggregated variant databases and resources (GnomAD, Cancer Genome Interpreter, Cancer Hotspots, CIVIC, Cosmic, SIFT, PolyPhen2) for the prediction of somatic mutations in cancer driver genes. 

