# RNAseq lncRNAs analisis con NGS data (Illumina Genome Analyzer) All-in-one 06/04/2021
# Análisis de trabscriptomas con Illumina RNAseq data: All-in-one
#Instalacion de los paquetes necesarios
BiocManager::install('Rfastp') #QC y procesamiento de fastq
BiocManager::install('Rsubread') #Alineamiento y mapeo
BiocManager::install('BSgenome.Hsapiens.UCSC.hg19') #genoma de referencia
BiocManager::install('limma') #Expression diferencial
BiocManager::install('edgeR') #Expresión diferencial
BiocManager::install('org.Hs.eg.db') #Anotación de genes
BiocManager::install('EnhancedVolcano') #Volcano plots
BiocManager::install('clusterProfiler') #Análisis de enriquecimiento funcional
BiocManager::install('EnhancedVolcano')

#Descargamos Genoma de referencia e indexado para Rsubread 
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome.Hsapiens.UCSC.hg19
genome <- BSgenome.Hsapiens.UCSC.hg19
seqlengths(genome)
genome$chr1  # same as genome[["chr1"]]
export(genome, "hg19.fasta")# Export as FASTA file.
library(Rsubread)
buildindex(basename="my_index_hg19",reference="hg19.fasta")

##Control de calidad de los archivos FASTQ y procesamiento de reads------------------------------------
library(Rfastp)
pe_read1 <- "Sample1_R1.fastq.gz"
pe_read2 <- "Sample2_R2.fastq.gz"
json_report <- rfastp(read1 = pe_read1, read2 = pe_read2,
                      outputFastq = "Sample1_filtered")

#Generate report tables/plots and quality curves
qcSummary(json_report)
curvePlot(json_report, curves = "quality_curves")

#Alineamiento de reads sobre genoma de referencia------------------------------------------------------
library(Rsubread)
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome.Hsapiens.UCSC.hg19
genome <- BSgenome.Hsapiens.UCSC.hg19
seqlengths(genome)
export(genome, "hg19.fasta")# Export as FASTA file.
buildindex(basename="my_index",reference="hg19.fasta")
#Alineamiento de los Fastq files para generar los Bam files 
align(index="my_index_hg19",readfile1="FH1005_R1.fastq.gz",readfile2="FH1005_R2.fastq.gz",output_file="FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=30)
align(index="my_index_hg19",readfile1="Mesri-11732-001-RNAseq_S1_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-001-RNAseq_S1_ME_L001_R2_001.fastq.gz",output_file="FH1005_L.bam",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-002-RNAseq_S2_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-002-RNAseq_S2_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-003-RNAseq_S3_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-003-RNAseq_S3_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-004-RNAseq_S4_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-004-RNAseq_S4_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-005-RNAseq_S5_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-005-RNAseq_S5_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-006-RNAseq_S6_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-006-RNAseq_S6_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-007-RNAseq_S7_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-007-RNAseq_S7_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-008-RNAseq_S8_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-008-RNAseq_S8_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-009-RNAseq_S9_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-009-RNAseq_S9_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-010-RNAseq_S10_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-010-RNAseq_S10_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-011-RNAseq_S11_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-011-RNAseq_S11_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-012-RNAseq_S12_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-012-RNAseq_S12_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-013-RNAseq_S13_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-013-RNAseq_S13_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)
align(index="my_index_hg19",readfile1="Mesri-11732-014-RNAseq_S14_ME_L001_R1_001.fastq.gz",readfile2="Mesri-11732-014-RNAseq_S14_ME_L001_R2_001.fastq.gz",output_file="nonQC_FH1005_L",type="rna",input_format="gzFASTQ",output_format="BAM",nthreads=35)


rm()
gc()

#Generación de la counts matrix
myfiles <- c("108731_5.bam","119203_4.bam","154979_2.bam","159770_4.bam","192373_2.bam","196782_2.bam","203403_5.bam","207178_9.bam","R307T.bam","R376T.bam","R710T.bam","R811T.bam","401044T_BT401044T.bam","401651T_BT401651T.bam",
             "402513T_BT402513T.bam","411611T.bam","533424_6.bam","641358_7.bam","642256_6.bam","648979_3.bam","706783_11.bam","723711_11.bam","750937_3.bam","90724_5.bam","92468_8.bam","A1_BNA1.bam","A2_BNA2.bam",
             "A3_BNA3.bam","A4_BNA4.bam","A5_BNA5.bam","A8_BNA8.bam","A9.bam","A10.bam","A11.bam","BNA8_CD10.bam","BNA8_EpCAM.bam","TLCCH_1.bam","TLDCIS_COM.bam","MCF10A_pLOCLinc885_2.bam","MCF10A_pLOCLinc885_1.bam",
             "MCF10A_pLOCLinc1024_2.bam","MCF10A_pLOCLinc1024_1.bam","MCF10A_pLOCempty_2.bam","MCF10A_pLOCempty_1.bam","MCF10A_pCDHHotair_2.bam","MCF10A_pCDHHotair_1.bam","MCF10A_pCDHempty_2.bam",
             "MCF10A_pCDHempty_1.bam","DCIScom_pLOCLinc885_2.bam","DCIScom_pLOCLinc885_1.bam","DCIScom_pLOCLinc1024_2.bam","DCIScom_pLOCLinc1024_1.bam","DCIScom_pLOCempty_2.bam","DCIScom_pLOCempty_1.bam",
             "DCIScom_pCDHHotair_2.bam","DCIScom_pCDHHotair_1.bam","DCIScom_pCDHempty_2.bam","DCIScom_pCDHempty_1.bam","4_DCISCOM_pLoc1011_2.bam","3_DCISCOM_pLoc1011_1.bam","2_DCISCOM_1024_2.bam","1_DCISCOM_1024_1.bam") #All samples

fc <- featureCounts(files=myfiles,annot.inbuilt="hg19",isPairedEnd=TRUE,nthreads=8)
write.table(fc$counts,file="Hotair_counts_matrix.csv",sep="\t")#Para salir a otros test estadsíticos

mygroups <- c(rep("MCF10A_empty",2),rep("MCF10A_Hotair",2),rep("DCIS_empty",2),rep("DCIS_Hotair",2))
d <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")], group = mygroups)
dim(d)
d$counts
d$samples


# Filtrado de datos: setear el N al numero de replicas del grupo de menor tamaño 
keep <- rowSums(cpm(d)>1) >= 2
d <- d[keep,] 
dim(d) 
d$samples$lib.size <- colSums(d$counts)

# Normalización
d <- calcNormFactors(d)
d$samples

# Exploración de datos
plotMDS(d)

# Estimación de dispersión y analisis estadistico
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
plotBCV(d)
et <- exactTest(d)
#Pair comparisons
# Comparación en MCF10A cells
et <- exactTest(d, pair=c("MCF10A_empty","MCF10A_Hotair"))
top <- topTags(et,n=40000)
top
write.table(top$table,"MCF10A_Hotair_stat.csv",sep="\t")
# Comparación en DCIS-COM cells
et <- exactTest(d, pair=c("DCIS_empty","DCIS_Hotair"))
top <- topTags(et,n=40000)
top
write.table(top$table,"DCIS_Hotair_stat.csv",sep="\t")

# LOG2 counts per million for top genes RECOMENDED FOR HCL
cpm <- cpm(d)
cpm <- cpm(d, prior.count=2, log=TRUE)
write.table(cpm,"Log2cpm.csv",sep="\t")

#Anotation----------------------------------------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("org.Mm.eg.db")  

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
fbids <- rownames(cpm)
cols = c("SYMBOL","GENENAME")
annots <- select(org.Hs.eg.db,keys=fbids,columns=cols,keytype="ENTREZID")
countTable <- cbind(ENTREZID=rownames(as.data.frame(cpm)),as.data.frame(cpm)) 
newTable = merge(countTable,annots,by="ENTREZID") 
write.table(newTable,file="normalized_CPM_annotated_data.csv",sep="\t")

#Unsupervised analysis: Ward Hierarchical Clustering----------------------------
expression <- t(cpm) #transposicion de exprs
distancia <- dist(expression, method = "euclidean") # distance matrix
fit <- hclust(distancia, method="ward.D2")
plot(fit)

#Plot the most basic volcano plot
mydata <- read.csv("MCF10A_Hotair_stat_edited.csv", header=TRUE, sep="\t")
attach(mydata)
plot(logFC,-log2(FDR))#Volcano plot elemental

library(EnhancedVolcano)
Volcano <- EnhancedVolcano(mydata,
                           lab = mydata$SYMBOL,
                           x = "logFC",
                           y = "PValue",
                           xlim = c(-10, 10),
                           title = 'MCF10A empty vs. Hotair',
                           pCutoff = 0.01,
                           FCcutoff = 1.0,
                           transcriptPointSize = 1.5,
                           transcriptLabSize = 3,
                           legendLabSize = 10,
                           legendIconSize = 3)
Volcano

#Selección de los DEGs (Differentially Expressed Genes) basado en críterios de significancia
FirstCriteria_Pvalue <- subset(mydata, PValue < 0.01, select = c("EntrezID","SYMBOL","GENENAME","PValue","logFC","FDR"))
SecondCriteria_FC <- subset(FirstCriteria_FDR, abs(logFC) >= 1, select = c("EntrezID","SYMBOL","GENENAME","PValue","logFC","FDR"))#Up- and Down-modulated in tumors
write.table(SecondCriteria_FC,file="selectedgenes.csv",sep="\t")


#Functional Enrichment
library(clusterProfiler)
de <- as.vector(FirstCriteria_Pvalue$EntrezID)
de <- as.vector(mydata$GeneID)
ego <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
head(ego)

library(enrichplot)
goplot(ego)#Gene Ontology (GO) is organized as a directed acyclic graph.

barplot(ego, showCategory=20) #Barplot del GO functional enrichment

dotplot(ego2, showCategory=30, orderBy = "x") #Dotplot del GO functional enrichment

#Gene-Concept Network
## remove redundent GO terms
ego2 <- simplify(ego)
cnetplot(ego2, foldChange=mydata$Overall.log2.Fold.Change)

cnetplot(ego2, foldChange=mydata$Overall.log2.Fold.Change,circular = TRUE, colorEdge = TRUE)

