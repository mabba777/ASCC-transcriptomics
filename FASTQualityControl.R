#An Ultra-Fast All-in-One FASTQ preprocessor
#https://www.bioconductor.org/packages/devel/bioc/vignettes/Rfastp/inst/doc/Rfastp.html#installation
library(Rfastp)

#a normal QC run for single-end fastq 
se_json_report <- rfastp(read1 = se_read1, 
                         outputFastq = paste0(outputPrefix, "_se"), thread = 4)


#a normal QC run for paired-end fastq files:Multiple files (https://rpubs.com/frm4001/895346)
fastq.files <- list.files(pattern = ".fastq.gz$", full.names = FALSE)
for (i in seq(1, length(fastq.files), 2)){
  pair1 <- fastq.files[i]
  pair2 <- fastq.files[i+1]
  output_name <- paste0("filtered_", fastq.files[i])
  json_report <- rfastp(read1= pair1 , read2 = pair2, outputFastq = output_name, thread=20)
}


#a normal QC run for paired-end fastq files.
pe_read1 <- system.file("/data/Ahuja_RNAseq/extdata","hMSC-219-KS1_S10_L001_R1_001.fastq.gz")
pe_read2 <- system.file("extdata","hMSC-219-KS1_S10_L001_R2_001.fastq.gz")
outputPrefix <- tempfile(tmpdir = tempdir())
pe_json_report <- rfastp(read1 = pe_read1, read2 = pe_read2,
                         outputFastq = paste0(outputPrefix, "_pe"))

pe_read1 <- system.file("extdata","reads1.fastq.gz",package="Rfastp")
pe_read2 <- system.file("extdata","reads2.fastq.gz",package="Rfastp")
outputPrefix <- tempfile(tmpdir = tempdir())

#Generate report tables/plots
#A data frame for the summary.
dfsummary <- qcSummary(pe_json_report)
#a ggplot2 object of base quality plot.
p1 <- curvePlot(pe_json_report)
p1
#a ggplot2 object of GC Content plot.
p2 <- curvePlot(pe_json_report, curve="content_curves")
p2



#Este script es de utilidad para agregar diferetes reporte generfknvkxv
#install.packages("fastqcr")
library(fastqcr)
library(dplyr)
# Reporte de un unico FASTQC samples 
qc <- qc_read("Mesri-20220505-KS-219-GEX-3v3-1-20220531_S33_L001_R1_001.fastq.gz")
qc_plot(qc, "summary")
qc_plot(qc, "Basic statistics")
qc_plot(qc, "Per base sequence quality")
qc_plot(qc, "Per sequence quality scores")
qc_plot(qc, "Per base sequence content")
qc_plot(qc, "Per sequence GC content")
qc_plot(qc, "Per base N content")
qc_plot(qc, "Sequence length distribution")
qc_plot(qc, "Sequence duplication levels")
qc_plot(qc, "Overrepresented sequences")
qc_plot(qc, "Adapter content")
qc_plot(qc, "Kmer content")

# Para generar reportes de multiples muestras
qc <- qc_aggregate("~/Desktop/LanariExomeSeq/QC", progressbar = FALSE)#multiple samples path a carpatea fastqc.zip files
summary <-summary(qc)
stats <- qc_stats(qc)
fails <- qc_fails(qc)
warned <- qc_warned(qc)
problemas<- qc_problems(qc, "module", compact = FALSE)

write.table(summary, "summary.csv", sep="\t")
write.table(stats, "stats.csv", sep="\t")
write.table(problemas, "qc_problems.csv", sep="\t")

--------------------------------------------
  