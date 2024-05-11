#Genomics Variants Annotation
#Instalacion de paquetes
#BiocManager::install("VariantAnnotation")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38") #Last R version
#BiocManager::install("SNPlocs.Hsapiens.dbSNP150.GRCh38") #Previous R version
#BiocManager::install("PolyPhen.Hsapiens.dbSNP131")

library("VariantAnnotation")
library("BSgenome.Hsapiens.UCSC.hg38")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("SNPlocs.Hsapiens.dbSNP150.GRCh38")
library("SNPlocs.Hsapiens.dbSNP155.GRCh38")

#Lectura del VCF
setwd("~/Escritorio/RNAseqAno2023/merged")
library(VariantAnnotation)
vcf <- readVcf("PGA871_SNPs_filtered.vcf", "hg38")
header(vcf)
meta(header(vcf))$fileformat
#Extraemos las posiciones genomicas de las variantes
head(rowRanges(vcf))
rd <- rowRanges(vcf)
#Obtenemos los datos de dbSNP
all_snps <- SNPlocs.Hsapiens.dbSNP150.GRCh38
#Chequeo de problemas: 
#seqlevels in the two objects dosn't match
seqlevels(rd)
seqlevels(all_snps)
#seqlevelStyle in the two objects dosn't match
seqlevelsStyle(rd)
seqlevelsStyle(all_snps)
#genome in the tow objects dosen't match
genome(rd)
genome(all_snps)

#Unify seqnames and Process
tar_chr <- as.vector(seqnames(rd)@values)
tar_chr <- gsub("chr", "", tar_chr)
tar_chr[grepl(tar_chr, pattern = "M")] <- "MT"

#dbSNP annotation
vcf <- names(rowRanges(vcf))
my_snps <- snpsBySeqname(all_snps, c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"))
anno <- mcols(my_snps)$RefSNP_id
in_dbSNP <- vcf_rsids %in% chr22_rsids
my_snps[1:2]

all <- locateVariants(vcf, all_snps, AllVariants())

library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
vcf_rsids <- names(rowRanges(vcf))
chr22snps <- snpsBySeqname(SNPlocs.Hsapiens.dbSNP150.GRCh38, "22")
chr22_rsids <- mcols(chr22snps)$RefSNP_id
in_dbSNP <- vcf_rsids %in% chr22_rsids
table(in_dbSNP)




#Predict amino acid changes---------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
head(seqlevels(txdb))

#Sanity check to confirm we have matching seqlevels and genomes 
intersect(seqlevels(txdb), seqlevels(vcf))
unique(genome(txdb))
unique(genome(vcf))

coding <- predictCoding(vcf, txdb, seqSource = Hsapiens)
coding[1]

#Transform into data.frame
matA <- data.frame(Variant=names(coding),chromosome=seqnames(coding),start=start(coding),end=end(coding),QUAL=coding$QUAL,
                   ref_allele=as.character(coding$REF),alt_allele=unlist(lapply(lapply(coding$ALT,`[[`,1),as.character)),
                   GeneID=coding$GENEID,TxID=coding$TXID,Protein_posi=unlist(lapply(lapply(coding$PROTEINLOC,`[[`,1),as.integer)),
                   ref_AA=as.character(coding$REFAA),alt_AA=as.character(coding$VARAA),Type=coding$CONSEQUENCE)

matA$aaChange <- paste0("p.",matA$ref_AA,matA$Protein_posi,matA$alt_AA)
matA <- dplyr::select(matA,-Protein_posi,-ref_AA,-alt_AA)
matA[1:3, ]

## Gene name annotation get symbols ids from gene entrez
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
geneentrez <- matA$GeneID
genesymbol <- select(org.Hs.eg.db, keys=geneentrez, keytype="ENTREZID", columns="SYMBOL")
matA$Symbol <- genesymbol$SYMBOL

#Subseting by mutation types
seleccion <- c("nonsynonymous","nonsense")
relevants <- subset(matA, matA$Type %in% seleccion)


##Eliminar SYMBOLS = NA
relevantsNonNA <- subset(relevants, !is.na(Symbol))

relevantsNonNA$TxID <- NULL

relevantsNonNAUnicos <- unique(relevantsNonNA)

write.table(relevantsNonNAUnicos,"PGA871_SNV.csv", sep="\t")

#How many variations in coding region
var_in_coding <- data.frame(varName = names(vcf), in_coding = names(vcf) %in% 
                              matA$Variant, stringsAsFactors = FALSE)
table(var_in_coding$in_coding)
#How many types of mutations in coding 
taC <- table(matA$Type)
taC_dat <- as.data.frame(taC)
taC
#Mutation types in coding region
library("ggplot2")
ggplot(taC_dat,aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat='Identity')+
  labs(x="",y="Counts",fill="")+
  theme(legend.position = "none")



all <- locateVariants(vcf, txdb, AllVariants())
#How many types of mutations in coding 


