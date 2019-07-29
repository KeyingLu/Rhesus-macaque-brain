

library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(dendextend)
library(gtools)
library(BiocParallel)
library(data.table)
library(Rtsne)


setwd("/mnt/data2/Rhesus_brain")
dir.create("revise0530/transcripts_count_res_with_TACO_gtf")


#############################################################
#                      loading Rhesus                       #
#############################################################
RhesusTransCount <- read.csv("revise0530/Rhesus_transcripts_counts_by_prepDE_with_TACO_minExpr5.5_assembly_gtf.csv", header=T, row.names = 1)
dim(RhesusTransCount) # 53661  408
RhesusTransCount <- RhesusTransCount[apply(RhesusTransCount, 1, function(x){sum(x>1)>3}),]
dim(RhesusTransCount) # 51406   408
save(RhesusTransCount, file="revise0530/transcripts_count_res_with_TACO_gtf/RhesusTransCount.RData")



#############################################################
#                           colData                         #
#############################################################
load("revise0530/gene_count_res/Rhesus_colData.RData")
save(Rhesus_colData, file="revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_colData.RData")



#############################################################
#                           DESeq2                          #
#############################################################
## run DESeq2 on R console
library(DESeq2)
library(BiocParallel)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_colData.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/RhesusTransCount.RData")


all(rownames(Rhesus_colData) %in% colnames(RhesusTransCount))
str(Rhesus_colData)
Rhesus_TransObject <- DESeqDataSetFromMatrix(RhesusTransCount, Rhesus_colData, ~SR+ID)
Rhesus_TransObject <- DESeq(Rhesus_TransObject, parallel = T, BPPARAM=MulticoreParam(10))
save(Rhesus_TransObject, file="revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_DESeq2_object.RData")













