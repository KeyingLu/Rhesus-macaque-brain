
# library(devtools)
# devtools::install_github("ctlab/linseed")



#######################################################
#            Single cell data (cibersort)            #
#######################################################
# /mnt/data1/Tools/R-3.6.0/bin/R
# library(devtools)
# install_github("rosedu1/deconvSeq")
# .libPaths("/home/forest/deconvseq/packrat/lib/x86_64-pc-linux-gnu/3.4.2")
library(deconvSeq)
library(reshape2)
setwd("/mnt/data2/Rhesus_brain")
load("SingleCell/Sestan.adultMonkeyNuclei.Psychencode.Rdata") # umi2 meta2
load("SingleCell/rhesus_fetal_scRNAseq.celltype.RData") # meta
load("revise0530/gene_count_res/RhesusGeneCount.RData")
load("revise0530/gene_count_res/Rhesus_norlCount.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


SigGenes <- read.table("SingleCell/rhesus_adult_snRNA-seq_celltype_signature.txt", 
                       stringsAsFactors = F, sep= "\t",
                       skip = 2)
colnames(SigGenes) <- SigGenes[1,]
SigGenes <- SigGenes[-1,-which(colnames(SigGenes)=="NA")]

top24_SigGenes <- list()
for(i in seq(2, dim(SigGenes)[2], 2)){
  CT <- colnames(SigGenes)[i]
  top24_SigGenes[[CT]] <- SigGenes[1:24,i-1]
}


top24_SigGenes <- melt(top24_SigGenes)
SigExpr <- umi2[as.character(top24_SigGenes$value), rownames(meta2)]
SigExpr <- SigExpr[unique(rownames(SigExpr)),]
dim(SigExpr) # 504 26933

SigExpr.mean <- sapply(levels(meta2$subtype), function(x){rowMeans(SigExpr[, meta2$subtype==x])})
dim(SigExpr.mean)
rownames(SigExpr.mean) <- sapply(rownames(SigExpr.mean), function(x){strsplit(x, "[|]")[[1]][2]})
SigExpr.mean <- data.frame(Gene_Name=rownames(SigExpr.mean), SigExpr.mean)
write.table(SigExpr.mean, "SingleCell/Signature_genes_Expr.mean.txt",
            row.names = F, col.names = T, quote = F, sep="\t")

norlCount <- Rhesus_norlCount
norlCount <- data.frame(Gene_Name=rownames(norlCount), norlCount)
# norlCount <- data.frame(Gene_Name=gtf_ensembl_gene[rownames(norlCount), "gene_name"], norlCount)

write.table(norlCount, "SingleCell/Rhesus_norlCount.txt",
            row.names = F, col.names = T, quote = F, sep="\t")



##
library(pheatmap)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dendsort)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")


Est.prop <- t(read.table("SingleCell/CIBERSORT.Output_Job1.txt",
                       header=T, row.names = 1, sep="\t"))
Est.prop <- Est.prop[1:21,]


colInfo <- Rhesus_colData[colnames(Est.prop ),c("SR_merge", "ID")]
head(colInfo)
str(colInfo)

ann_colors <- list(SR_merge = c(anno_colors$SR,anno_colors$merge)[unique(colInfo$SR_merge)],
                   ID = anno_colors$ID)


pheatmap(Est.prop,
         filename="revise0530/gene_count_res/CIBERSORT_Rhesus_deconv_prop_heatmap.pdf",
         col=colorRampPalette(rev(brewer.pal(11,"RdBu")[c(2:5,10)]))(100),
         show_rownames=T,
         show_colnames = F,
         annotation_col = colInfo,
         annotation_colors = ann_colors,
         cluster_rows = as.hclust(hclust(dist(Est.prop, method = "maximum"))))



Est.prop.melt <- melt(Est.prop)
colnames(Est.prop.melt) <- c("CellType", "Samples", "Proportions")

Rhesus_colData$SR

Est.prop.melt <- data.frame(Est.prop.melt, 
                            Rhesus_colData[as.character(Est.prop.melt$Samples), c("ID", "SR", "SR_merge", "Sex", "Age_Stage")])
pp <- Est.prop.melt[Est.prop.melt$CellType %in%
                      names(rowSums(Est.prop))[rowSums(Est.prop)>0],]

pp$CellType <- factor(as.character(pp$CellType))
pp$Region <- paste(pp$SR_merge, pp$SR, sep=":")



### 
ggplot(pp, aes(x=SR_merge, y=Proportions, fill=SR_merge, color=SR_merge)) +
  geom_boxplot(alpha=0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
  facet_wrap(~CellType, ncol=5, scales = "free_y") +
  scale_fill_manual(values = c(anno_colors$SR, anno_colors$merge)[unique(pp$SR_merge)]) +
  scale_color_manual(values = c(anno_colors$SR, anno_colors$merge)[unique(pp$SR_merge)])


ggsave("revise0530/gene_count_res/Rhesus_deconv_prop_boxplot.pdf",
       height = 12, width = 18)

p <- pp[, c("CellType", "Samples", "Proportions", "SR_merge")]
colnames(p)[4] <- "Region"
write.table(p, "SourceData/Fig.S2c.txt",
            col.names = T, row.names = F, quote=F, sep="\t")



Rhesus_norlCount <- counts(Rhesus_GeneObject, normalized=T)
save(Rhesus_norlCount, file="revise0530/gene_count_res/Rhesus_norlCount.RData")


#########
# correlation 
########
library(Rtsne)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")

df <- Est.prop
corr <- cor(df)
pp <- melt(corr)
pp <- pp[1:(dim(pp)[1]/2), ]
pp$Var1 <- as.character(pp$Var1)
pp$Var2 <- as.character(pp$Var2)
pp$ID1 <- Rhesus_colData[pp$Var1, "ID"]
pp$ID2 <- Rhesus_colData[pp$Var2, "ID"]
pp$SR1 <- Rhesus_colData[pp$Var1, "SR"]
pp$SR2 <- Rhesus_colData[pp$Var2, "SR"]
pp <- pp[pp$SR1==pp$SR2,]
pp <- pp[!(pp$ID1==pp$ID2),]
pp$SR_merge <-  Rhesus_colData[pp$Var1, "SR_merge"]
pp$SR <- paste(pp$SR_merge, pp$SR1, sep=":")
pp <- pp[order(pp$SR_merge, pp$SR1), ]
pp$SR1 <- factor(pp$SR1, levels = unique(pp$SR1))
pp$SR_merge2 <- pp$SR_merge
pp$SR_merge2[pp$SR_merge2!="cortex"] <- "Subcortical"


ggplot(pp, aes(x=SR1, y=value, color=SR1)) +
  geom_boxplot() +
  geom_jitter(size=0.5, color="grey") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
  scale_fill_manual(values = anno_colors$SR) +
  scale_color_manual(values = anno_colors$SR) +
  facet_wrap(~SR_merge2, ncol=1, scales = "free_x") +
  ylab("Correlation") +
  guides(fill=F, colour=F)

ggsave("revise0530/gene_count_res/deconv_boxplot_cor.pdf")


p <- pp[, c("value", "SR_merge2", "SR1")]
colnames(p) <- c("Expression", "Region1", "Region2")
write.table(p, "SourceData/Fig.S2c.txt",
            col.names = T, row.names = F, quote=F, sep="\t")


