library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(dendextend)
library(gtools)
library(BiocParallel)
library(data.table)
library(Rtsne)


setwd("/mnt/data2/Rhesus_brain")
dir.create("revise0530/gene_count_res")


##########################################################
#                     loading Rhesus                     #
##########################################################
Rhesus_Info <- read.table("revise0530/Rhesus_sample_info.txt", 
                          header = T, sep = "\t", stringsAsFactors = F)
rownames(Rhesus_Info) <- Rhesus_Info[,1]
Rhesus_Info[,3] <- gsub(" ", "_", Rhesus_Info[,3])
head(Rhesus_Info)

Rhesus_Info_2 <- read.table("revise0530/region_info.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
head(Rhesus_Info_2)

RhesusGeneCount <- read.table("revise0530/Rhesus_count.matrix", header=T, row.names = 1)
dim(RhesusGeneCount) # [1] 30807   408
order <- apply(RhesusGeneCount, 1, function(x){sum(x>1)>3})
RhesusGeneCount <- RhesusGeneCount[order,]
dim(RhesusGeneCount) # 23651   408

save(RhesusGeneCount, file="revise0530/gene_count_res/RhesusGeneCount.RData")





##########################################################
#                      loading gtf                       #
##########################################################
library(devtools)
library(rtracklayer)
gtf_ensembl <- rtracklayer::import("/mnt/data1/Ref/Rhesus_macaque_8.0.1/Macaca_mulatta.Mmul_8.0.1.91.chr.gtf")
gtf_ensembl <- as.data.frame(gtf_ensembl)
View(gtf_ensembl[1:10,])
gtf_ensembl_gene <- gtf_ensembl[which(gtf_ensembl$type=="gene"),
                                c(1,2,3,5,10,23,13)]
gtf_ensembl_trans <- gtf_ensembl[which(gtf_ensembl$type=="transcript"), ]
gtf_ensembl_exon <- gtf_ensembl[which(gtf_ensembl$type=="exon"), 
                                c(1, 2, 3, 24, 14, 5, 7,10,23)]
colnames(gtf_ensembl_exon)[c(1,3)] <- c("chrom", "stop")
rownames(gtf_ensembl_gene) <- gtf_ensembl_gene$gene_id
rownames(gtf_ensembl_trans) <- gtf_ensembl_trans$transcript_id
save(gtf_ensembl_gene, file="revise0530/gene_count_res/gtf_ensembl_gene.RData")
save(gtf_ensembl_trans, file="revise0530/gene_count_res/gtf_ensembl_trans.RData")
save(gtf_ensembl_exon, file="revise0530/gene_count_res/gtf_ensembl_exon.RData")



ref_coding_trans <- gtf_ensembl_trans[which(gtf_ensembl_trans$gene_biotype == "protein_coding"),]
ref_coding_trans <- ref_coding_trans[-which(is.na(ref_coding_trans$gene_name)),]
View(ref_coding_trans)
write.table(rownames(ref_coding_trans), "revise0530/gene_count_res/ref_coding_trans_id.txt",
            quote=F, sep="\t", row.names = F, col.names = F)
write.table(ref_coding_trans[, c("transcript_id", "gene_name")], "revise0530/gene_count_res/ref_coding_trans_id_gene_name.txt",
            quote=F, sep="\t", row.names = F, col.names = F)




##########################################################
#                    Rhesus colData                      #
##########################################################
Rhesus_Info <- read.table("revise0530/Rhesus_sample_info.txt", 
                          header = T, sep = "\t", stringsAsFactors = F)
rownames(Rhesus_Info) <- Rhesus_Info[,1]
Rhesus_Info[,3] <- gsub(" ", "_", Rhesus_Info[,3])
head(Rhesus_Info)

Rhesus_Info_2 <- read.table("revise0530/region_info.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)
head(Rhesus_Info_2)


Rhesus_RIN <- read.table("revise0530/RNA_RIN.txt", header = F, 
                         row.names = 1, stringsAsFactors = F)
rownames(Rhesus_RIN) <- paste("X", rownames(Rhesus_RIN), sep="")

map <- setNames(c("Rhesus_1", "Rhesus_2", "Rhesus_3", "Rhesus_4", "Rhesus_5", "Rhesus_6", "Rhesus_7", "Rhesus_8"),
                c("X2", "X3", "X4", "X6", "X7", "X8", "X10", "X11"))
a <- sapply(colnames(RhesusGeneCount), function(x) strsplit(x, "_")[[1]][1])
b <- sapply(colnames(RhesusGeneCount), function(x) strsplit(x, "_")[[1]][2])

Rhesus_colData <- Rhesus_Info[map[a],1:4]
Rhesus_colData <- data.frame(Rhesus_colData, Rhesus_Info_2[b,2:3])
rownames(Rhesus_colData) <- colnames(RhesusGeneCount)
head(Rhesus_colData)

colnames(Rhesus_colData)[c(1, 5:6)] <- c("ID", "SR", "BR")
head(Rhesus_colData)
str(Rhesus_colData)
Rhesus_colData$BR <- gsub(" ", "_", Rhesus_colData$BR)
SR <- as.character(Rhesus_colData$SR)
BR <- as.character(Rhesus_colData$BR)

testInfo <- SR
testInfo[BR %in% c("frontal_lobe","parietal_lobe","occipital_lobe","temporal_lobe","Limbic_cortex")] <- "cortex"
testInfo[SR == "OB"] <- "OB"
Rhesus_colData$SR_cortex <- testInfo


Rhesus_colData$SR_merge <- Rhesus_colData$SR_cortex
Rhesus_colData$SR_merge[which(Rhesus_colData$SR_merge=="DG"|Rhesus_colData$SR_merge=="CA1"|Rhesus_colData$SR_merge=="CA3")] <- "HIP"
Rhesus_colData$SR_merge[which(Rhesus_colData$SR_merge=="ACb"|Rhesus_colData$SR_merge=="CN"|Rhesus_colData$SR_merge=="PUT")] <- "STR"
Rhesus_colData$SR_merge[which(Rhesus_colData$SR_merge=="CBC"|Rhesus_colData$SR_merge=="CBV")] <- "CB"


Rhesus_colData$RIN <- Rhesus_RIN[rownames(Rhesus_colData), 1]

save(Rhesus_colData, file="revise0530/gene_count_res/Rhesus_colData.RData")



############################################################
#                     ann_colors                           #
############################################################
library(RColorBrewer)
library(wesanderson)
setwd("/mnt/data2/Rhesus_brain")

anno_colors <- list(SR=setNames(colorRampPalette(brewer.pal(12,"Paired"))(51), levels(factor(Rhesus_colData$SR))),
                    ID=setNames(c(wes_palette('Royal1'),wes_palette('GrandBudapest2')), levels(factor(Rhesus_colData$ID))),
                    BR=setNames(brewer.pal(12,"Paired")[1:10], levels(factor(Rhesus_colData$BR))),
                    merge = c(cortex=brewer.pal(7, "Reds")[6],
                              CB="#29AAE1",
                              HIP="#B6DCAE",
                              STR="#794FA3"),
                    Age_Stage = c(Mid="#3F60AC", Young="#8B2052"),
                    Sex = c(c(female="#BF645D", male="#F8CE7A")))

save(anno_colors, file="revise0530/gene_count_res/Rhesus_anno_colors.RData")




###########################################################
#                       DESeq2                            #
###########################################################
## run DESeq2 on R console
library(DESeq2)
library(BiocParallel)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/RhesusGeneCount.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
str(Rhesus_colData)


Rhesus_GeneObject <- DESeqDataSetFromMatrix(RhesusGeneCount, 
                                               Rhesus_colData, 
                                               ~SR+Age_Stage+Sex)
Rhesus_GeneObject <- DESeq(Rhesus_GeneObject, 
                              parallel = T, BPPARAM=MulticoreParam(40))
save(Rhesus_GeneObject, 
     file="revise0530/gene_count_res/Rhesus_DESeq2_object.RData")




###############################################################
#                         t-SNE                               #
###############################################################
library(Rtsne)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")

vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)


# PCA-based t-SNE
seed = 1:50
per = c(5, 10, 15, 20, 25, 30, 35)
for(i in seed){
  # for(perplexity in per){
    i=46
    set.seed(i)
    perplexity = 40
    tsne <- Rtsne::Rtsne(t(vsd_data), initial_dims=50,
                         max_iter = 5000,perplexity=perplexity)
    rownames(tsne$Y) <- colnames(vsd_data)
    colnames(tsne$Y) <- c("tsne1", "tsne2")
    Rhesus_colData <- data.frame(colData(Rhesus_GeneObject))
    p <- data.frame(tsne$Y, Rhesus_colData[rownames(tsne$Y),])
    head(p)
    title = paste("t-SNE : seed=", i, " perplexity=", perplexity, sep="")
    
    # ggplot(p, aes(x=tsne1, y=tsne2, colour=SR)) +
    ggplot(p, aes(x=tsne1, y=tsne2, colour=ID)) +
      geom_point(size=0.3) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face="bold"),
            axis.title = element_text(face="bold"),
            axis.text = element_text(face="bold"),
            legend.title = element_text(face="bold")) +
      scale_colour_manual(values=anno_colors$ID) +
      # scale_colour_manual(values=anno_colors$SR) +
      # scale_colour_manual(values=brewer.pal(12,"Paired")) +
      # geom_text(aes(label=rownames(p)), vjust=-1.2, size=0.5, colour="black") +
      # geom_text(aes(label=SR), vjust=-1.2, size=1, colour="black") +
      ggtitle(title) +
      guides(colour = guide_legend(override.aes = list(size=2))) +
      xlab("t-SNE1") +
      ylab("t-SNE2")
    
    # 
    # filename = paste("revise0530/gene_count_res/t_SNE_res/Rhesus_Gene_t-SNE_SR_seed", i, "_per", perplexity, ".pdf", sep="")
    # ggsave(filename,width = 7.5, height = 5)
    filename = paste("revise0530/gene_count_res/t_SNE_res/Rhesus_Gene_t-SNE_ID_seed", i,  "_per", perplexity, ".pdf", sep="")
    ggsave(filename,width = 6, height = 5)
  # }
}



###############################################################
#                         MDS plot                            #
###############################################################
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")

vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)


# d <- dist(t(vsd_data)) # euclidean distances between the rows
# fit <- cmdscale(d, k=2, eig = T) # k is the number of dim
# fit # view results
# colnames(fit) <- c("PC1", "PC2")

s <- svd(vsd_data-rowMeans(vsd_data))
PC1 <- s$d[1]*s$v[,1]
PC2 <- s$d[2]*s$v[,2]
var <- s$d^2/sum(s$d^2)


p <- data.frame(PC1, PC2, Rhesus_colData[colnames(vsd_data),])
head(p)

# ggplot(p, aes(x=PC1, y=PC2, colour=SR)) +
# ggplot(p, aes(x=PC1, y=PC2, colour=Age_Stage)) +
# ggplot(p, aes(x=PC1, y=PC2, colour=Sex)) +
ggplot(p, aes(x=PC1, y=PC2, colour=ID)) +
  geom_point(size=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.title = element_text(face="bold"),
        axis.text = element_text(face="bold"),
        legend.title = element_text(face="bold")) +
  scale_colour_manual(values=anno_colors$ID) +
  # scale_colour_manual(values=anno_colors$SR) +
  # scale_colour_manual(values=brewer.pal(12,"Paired")) +
  # geom_text(aes(label=rownames(p)), vjust=-1.2, size=0.5, colour="black") +
  geom_text(aes(label=SR), vjust=-1.2, size=1, colour="black") +
  # ggtitle(title) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  xlab("PC1") +
  ylab("PC2")


filename = paste("revise0530/gene_count_res/Rhesus_Gene_MDS_SR.pdf", sep="")
ggsave(filename,width = 7.5, height = 5)
filename = paste("revise0530/gene_count_res/Rhesus_Gene_MDS_Sex.pdf", sep="")
ggsave(filename,width = 5.5, height = 5)
filename = paste("revise0530/gene_count_res/Rhesus_Gene_MDS_Age_Stage.pdf", sep="")
ggsave(filename,width = 5.5, height = 5)
filename = paste("revise0530/gene_count_res/Rhesus_Gene_MDS_ID.pdf", sep="")
ggsave(filename,width = 5.5, height = 5)




###########################################################
#                   DESeq2 (exclude 6_36)                 #
###########################################################
## run DESeq2 on R console
# library(DESeq2)
# library(BiocParallel)
# options(stringAsFactors=FALSE)
# setwd("/mnt/data2/Rhesus_brain")
# load("revise0530/gene_count_res/RhesusGeneCount.RData")
# load("revise0530/gene_count_res/Rhesus_colData.RData")
# str(Rhesus_colData)
# 
# 
# RhesusGeneCount <- RhesusGeneCount[, -which(colnames(RhesusGeneCount)=="X6_36")]
# Rhesus_colData <- Rhesus_colData[colnames(RhesusGeneCount),]
# Rhesus_GeneObject_ex1 <- DESeqDataSetFromMatrix(RhesusGeneCount, 
#                                             Rhesus_colData, 
#                                             ~SR+Age_Stage+Sex)
# Rhesus_GeneObject_ex1 <- DESeq(Rhesus_GeneObject_ex1, 
#                            parallel = T, BPPARAM=MulticoreParam(40))
# save(Rhesus_GeneObject_ex1, 
#      file="revise0530/gene_count_res/Rhesus_GeneObject_ex1.RData")




###########################################################
#          DESeq2 (for cortex exclude chrX,Y)             #
###########################################################
## run DESeq2 on R console
library(DESeq2)
library(BiocParallel)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/RhesusGeneCount.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
str(Rhesus_colData)

####
genes <- gtf_ensembl_gene[-which(gtf_ensembl_gene$seqnames %in% c("X", "Y")),]$gene_id
RhesusGeneCount <- RhesusGeneCount[intersect(genes, rownames(RhesusGeneCount)), 
                                   rownames(Rhesus_colData)[Rhesus_colData$SR_cortex %in% "cortex"]]
order <- apply(RhesusGeneCount, 1, function(x){sum(x>1)>3})
RhesusGeneCount <- RhesusGeneCount[order, ]
dim(RhesusGeneCount) # 21306   248
###
Rhesus_colData <- Rhesus_colData[colnames(RhesusGeneCount),]
Rhesus_GeneObject_cortex <- DESeqDataSetFromMatrix(RhesusGeneCount, 
                                                Rhesus_colData, 
                                                ~ID+SR)
Rhesus_GeneObject_cortex <- DESeq(Rhesus_GeneObject_cortex,
                                  test = "LRT",
                                  reduced = ~SR,
                               parallel = T, BPPARAM=MulticoreParam(40))
save(Rhesus_GeneObject_cortex, 
     file="revise0530/gene_count_res/Rhesus_GeneObject_cortex.RData")




############################################################
#                    t-SNE for cortex                      #
############################################################
library(Rtsne)
library(DESeq2)
library(ggplot2)
library(wesanderson)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_GeneObject_cortex.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")

vsd <- varianceStabilizingTransformation(Rhesus_GeneObject_cortex, blind=FALSE)
vsd_data <- assay(vsd)
dim(vsd_data) # 21306   248


# Alternative: Distance-based t-SNE:
for(seed in 1:10){
  for(perplexity in seq(5, 40, 5)){
    # seed = 10
    set.seed(seed)
    # perplexity = 6
    tsne <- Rtsne::Rtsne(t(vsd_data), initial_dims=50, 
                         max_iter = 5000,perplexity=perplexity)
    rownames(tsne$Y) <- colnames(vsd_data)
    colnames(tsne$Y) <- c("tsne1", "tsne2")
    p <- data.frame(tsne$Y, Rhesus_colData[rownames(tsne$Y),])
    head(p)
    title = paste("t-SNE for cortex ", "seed:", seed, " ", "perplexity:", perplexity, sep="")
    
    ggplot(p, aes(x=tsne1, y=tsne2, colour=ID)) +
      geom_point(size=0.6) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      scale_colour_manual(values=anno_colors$ID) +
      # scale_colour_manual(values=anno_colors$BR) +
      # geom_text(aes(label=SR), vjust=-1.2, size=1, colour="black") +
      ggtitle(title) +
      guides(colour = guide_legend(override.aes = list(size=2)))
    
    filename = paste("revise0530/gene_count_res/t_SNE_res/Rhesus_t-SNE_ID_for_cortex_ex_sex_seed", seed, "_per", perplexity, ".pdf", sep="")
    ggsave(filename, width = 6, height = 5)
    # filename = paste("matrix0420/gene_count_res/Rhesus_t-SNE_BR_for_cortex_ex_sex_per", perplexity, ".pdf", sep="")
    # ggsave(filename, width = 6.5, height = 5)
    
  }
}





################################################################
#              Likelihood ratio test for SR                    #
################################################################
# library(DESeq2)
# library(BiocParallel)
# setwd("/mnt/data2/Rhesus_brain")
# load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
# load("revise0530/gene_count_res/Rhesus_colData.RData")
# load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
# 
# 
# # the p value is for the likelihood ratio test of all the variables and all the levels, 
# # while the log fold change is a single comparison from among those variables and levels.
# 
# dds <- DESeq(Rhesus_GeneObject, test="LRT", 
#              reduced=~Age_Stage+Sex,
#              parallel = T, BPPARAM=MulticoreParam(40))
# res <- results(dds)
# 
# 
# res <- data.frame(res)
# res <- res[order(res$padj),]
# res <- res[res$padj<0.05,] # 18437     6
# genes <- rownames(res)[1:5000]
# 
# 
# ########
# library(RColorBrewer)
# library(pheatmap)
# Rhesus_norCounts <- counts(Rhesus_GeneObject_ex1, normalized=TRUE)
# dim(Rhesus_norCounts) # 23651   407
# df <- Rhesus_norCounts
# 
# colInfo <- Rhesus_colData[colnames(df),c("SR_cortex", "ID")]
# head(colInfo)
# str(colInfo)
# 
# ann_colors <- list(SR_cortex = c(anno_colors$SR,anno_colors$merge)[unique(colInfo$SR_cortex)],
#                    ID = anno_colors$ID)
# 
# dfheat <- Rhesus_norCounts[genes, ]
# 
# 
# method1 = c("euclidean", "maximum", "manhattan", "canberra", "minkowski")
# method2 = c("ward.D", "ward.D2", "single",  "complete", "average")
# method3 = c("average","min")
# for(i in method1){
#   for(j in method2){
#     for(z in method3){
#       filename= paste("revise0530/gene_count_res/test_SR_profile_", i, "_", j, "_",z,".pdf", sep="")
#       title = paste("distance:", i, " hclust:",j, " type:",z, sep = "")
#       pheatmap(
#         dfheat,
#         col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
#         breaks = c(-4,seq(-3,3,length=20),4),
#         legend_breaks = c(-3,-1,1,3),
#         filename = filename,
#         scale="row",
#         main = title,
#         width = 10,
#         height = 10,
#         # clustering_method = "complete",
#         border_color=NA,
#         # fontsize_col = 0.5,
#         # annotation_row = geneCluster,
#         annotation_col = colInfo,
#         annotation_colors = ann_colors,
#         annotation_legend = T,
#         show_rownames = F,
#         show_colnames = F,
#         # cluster_rows = F,
#         cluster_cols = as.hclust(dendsort(hclust(dist(t(dfheat), method=i), method = j), type=z))
#       )
#     }
#   }
# }
# ##



################################################################
#         Likelihood ratio test for ID (in cortex)             #
################################################################
library(DESeq2)
library(BiocParallel)
library(ConsensusClusterPlus)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_GeneObject_cortex.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


####
resultsNames(Rhesus_GeneObject_cortex)
res <- results(Rhesus_GeneObject_cortex, contrast=c("ID", "Rhesus_3", "Rhesus_1"))
res <- results(Rhesus_GeneObject_cortex)
save(res, file="revise0530/gene_count_res/cortex_LRT_res.RData")


res <- data.frame(res)
res <- res[order(res$padj),]
res <- res[res$padj<0.05,]
dim(res) # 17411     6
genes <- rownames(res)[1:1000]



##### gene cluster 
Rhesus_norCounts <- counts(Rhesus_GeneObject_cortex, normalized=TRUE)
dim(Rhesus_norCounts) # 21306   248
df <- Rhesus_norCounts
top_df <- df[genes,]
ttop_df <- t(top_df)
ttop_df = scale(ttop_df)
title="revise0530/gene_count_res/ConsensusCluster_LRT_DifGenes_ID_in_cortex_top1000"
DIfGenes_Results_ID_in_cortex = ConsensusClusterPlus(ttop_df,maxK=15,reps=100,pItem=0.8,pFeature=1,
                                                     title=title,clusterAlg="km",
                                                     innerLinkage = "ward.D2",
                                                     finalLinkage = "ward.D2",
                                                     seed = 1,
                                                     distance="euclidean",
                                                     plot="png")
cal <- calcICL(DIfGenes_Results_ID_in_cortex)
cal$clusterConsensus
save(DIfGenes_Results_ID_in_cortex, 
     file="revise0530/gene_count_res/ConsensusCluster_LRT_DifGenes_ID_in_cortex_top1000.RData")



#################
library(pheatmap)
library(dendsort)
library(RColorBrewer)
load("revise0530/gene_count_res/ConsensusCluster_LRT_DifGenes_ID_in_cortex_top1000.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")

k=8
res <- DIfGenes_Results_ID_in_cortex[[k]]
geneCluster <- data.frame(geneCluster=factor(sort(res$consensusClass)))


colInfo <- Rhesus_colData[colnames(df),c("Age_Stage", "ID", "Sex")]
head(colInfo)
str(colInfo)

ann_colors <- list(Age_Stage = anno_colors$Age_Stage,
                   Sex = anno_colors$Sex,
                   ID = anno_colors$ID)


# for(i in seq(500, 5000, 500)){
for(i in c(1000)){
  dfheat <- top_df[rownames(geneCluster), ]
  filename= paste("revise0530/gene_count_res/ID_LRT_profile_top", i, ".pdf", sep="")
  pheatmap(
    dfheat,
    col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
    breaks = c(-4,seq(-3,3,length=20),4),
    legend_breaks = c(-3,-1,1,3),
    filename = filename,
    scale="row",
    # main = title,
    width = 10,
    height = 10,
    # clustering_method = "complete",
    border_color=NA,
    # fontsize_col = 0.5,
    # annotation_row = geneCluster,
    gaps_row = cumsum(table(geneCluster$geneCluster)),
    annotation_col = colInfo,
    annotation_colors = ann_colors,
    annotation_legend = T,
    show_rownames = F,
    show_colnames = F,
    cluster_rows = F
  )
}


method1 = c("euclidean", "maximum", "manhattan", "canberra", "minkowski")
method2 = c("ward.D", "ward.D2", "single",  "complete", "average")
method3 = c("average","min")
for(i in method1){
  for(j in method2){
    for(a in method3){
      # i = "euclidean"
      # j = "single"
      # a = "min"
      main = paste("dist:", i, " hclust:", j, " dendsort:", a, sep="")
      filename = paste("revise0530/gene_count_res/ID_LRT_profile_top1000_", i, "_",j, "_",a, ".pdf", sep="")
      dist_data <- dist(t(dfheat), method=i) 
      hc <- hclust(dist_data, method = j) # 
      dd <- dendsort(hc, isReverse = F, type=a) # 
      hc_sorted  <- as.hclust(dd)
      pheatmap(
      dfheat,
      col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
      breaks = c(-3,seq(-2,2,length=20),3),
      legend_breaks = c(-2,-1,1,2),
      filename = filename,
      scale="row",
      main = main,
      width = 10,
      height = 10,
      # clustering_method = "complete",
      border_color=NA,
      # fontsize_col = 0.5,
      # annotation_row = geneCluster,
      gaps_row = cumsum(table(geneCluster$geneCluster)),
      annotation_col = colInfo,
      annotation_colors = ann_colors,
      annotation_legend = T,
      show_rownames = F,
      show_colnames = F,
      cluster_rows = F,
      cluster_cols = hc_sorted 
      )
    }
  }
}

##



##### function 
####################
library(gProfileR)
setwd("/mnt/data2/Rhesus_brain")
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/ConsensusCluster_LRT_DifGenes_ID_in_cortex_top1000.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")

k=8
res <- DIfGenes_Results_ID_in_cortex[[k]]
geneCluster <- data.frame(geneCluster=factor(sort(res$consensusClass)))


geneCluster_function_list <- list()
for(cluster in 1:k){
  gprofiler_res <- gprofiler(
    query=rownames(geneCluster)[geneCluster$geneCluster==cluster], 
    organism = "mmulatta")
  gprofiler_res <- gprofiler_res[order(gprofiler_res[,3]),]
  gprofiler_res$gene_name <- gprofiler_res$intersection
  if(dim(gprofiler_res)[1]==0) next
  else{for(i in 1:dim(gprofiler_res)[1]){
    x <- gprofiler_res$gene_name[i]
    gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],5], collapse=",")
  }
    re <- paste("cluster", cluster, sep="")
    geneCluster_function_list[[re]] <- gprofiler_res
    filename=paste("revise0530/gene_count_res/LRT_ID_geneCluster", cluster,"_function.txt", sep="")
    write.table(gprofiler_res, filename,
                quote = F, sep="\t", row.names = F, col.names = T)
  }
}


######
gprofiler_res <- gprofiler(
  query=rownames(geneCluster), organism = "mmulatta")
gprofiler_res <- gprofiler_res[order(gprofiler_res[,3]),]
gprofiler_res$gene_name <- gprofiler_res$intersection
for(i in 1:dim(gprofiler_res)[1]){
  x <- gprofiler_res$gene_name[i]
  gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],"gene_name"], collapse=",")
}
filename=paste("revise0530/gene_count_res/LRT_ID_top1000_function.txt", sep="")
write.table(gprofiler_res, filename,
            quote = F, sep="\t", row.names = F, col.names = T)


gene_name <- sapply(gprofiler_res$gene_name, function(x){strsplit(x, ",")[[1]]})

### term plot
pp <- gprofiler_res[gprofiler_res$domain %in% c("BP", "keg", "MF"), c("p.value", "domain", "term.name")]
order = pp$term.name[order(pp$p.value, decreasing=T)]
pp$term.name <- factor(pp$term.name, levels = order)

ggplot(pp, aes(x=term.name, y=-log10(p.value))) +
  geom_bar(stat="identity", fill="#4D9DE0", width = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=1, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        # legend.position = c(0.8,0.8),
        strip.background = element_rect(fill="white", colour = "black", size=1) ,
        legend.title = element_text(face="bold")) +
  # geom_text(aes(label=symbol), vjust=0.5, hjust=0,size=4, colour="black") +
  # annotate("text", x=pp$NAME[2], y=-1.5,size=3,
  #          label="nominal P values: \n*P < 0.05 \n**P < 0.005 \n***P < 0.0005 \nNS: not significant") +
  # scale_x_discrete(limits=c(pp$term.name))+
  coord_flip() +
  xlab("") 
# facet_grid(.~module, scales = "free_x", space = "free_x") 


ggsave("revise0530/gene_count_res/LRT_for_ID_top1000_function_bar.pdf",
       height = 7, width = 6)



##### write down
####################
library(gProfileR)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/cortex_LRT_res.RData") # res
cortex_LRT_res <- res
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/ConsensusCluster_LRT_DifGenes_ID_in_cortex_top1000.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")


k=8
res <- DIfGenes_Results_ID_in_cortex[[k]]
geneCluster <- data.frame(geneCluster=factor(sort(res$consensusClass)))


cluster_info_list <- list()
for(cluster in 1:k){
  genes <- rownames(geneCluster)[geneCluster$geneCluster==cluster]
  gene_name <- gtf_ensembl_gene[genes, "gene_name"]
  sub_res <- cortex_LRT_res[genes, ]
  sub_res <- sub_res[order(sub_res$padj),]
  sub_res <- data.frame(gene_id=genes, gene_name, 
                        sub_res[, c("baseMean", "pvalue", "padj")],
                        Cluster=cluster)
  cluster_info_list[[cluster]] <- sub_res
  # filename=paste("revise0530/gene_count_res/LRT_ID_geneCluster", cluster,"_order.txt", sep="")
  # write.table(sub_res, filename,
  #               quote = F, sep="\t", row.names = F, col.names = T)
}

pp <- Reduce(rbind, cluster_info_list)
write.table(pp, "revise0530/gene_count_res/LRT_ID_geneCluster_info.txt",
              quote = F, sep="\t", row.names = F, col.names = T)


#################
## selected term genes
################
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_GeneObject_cortex.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/cortex_LRT_res.RData")
cortex_LRT_res <- res
filename=paste("revise0530/gene_count_res/LRT_ID_top1000_function.txt", sep="")
gf <- read.table(filename, sep="\t", header = T, stringsAsFactors = F)

##
pp <- gf[-which(gf$domain=="hp"), c("p.value", "domain", "term.name", "intersection", "gene_name")]

term_list <- list(term1="catalytic activity",
                  term2="oxidoreductase activity",
                  term3="MHC protein complex")
term_res_list <- list()
for(term in names(term_list)){
  term_name <- term_list[[term]]
  intersection <- pp[pp$term.name %in% term_name, "intersection"]
  intersection <- strsplit(intersection, ",")[[1]]
  term_res <- cortex_LRT_res[intersection,]
  term_res <- term_res[order(term_res$padj),]
  term_res$gene_name <- gtf_ensembl_gene[rownames(term_res), "gene_name"]
  term_res <- term_res[complete.cases(term_res),]
  if(term == "term2"){
    if(any(rownames(term_res) %in% rownames(term_res_list[[1]])[1:4])){
      term_res <- term_res[-which(rownames(term_res) %in% rownames(term_res_list[[1]])[1:4]),]}
  }
  if(term == "term3"){
    if(any(rownames(term_res) %in% rownames(term_res_list[[1]])[1:4] | rownames(term_res) %in% rownames(term_res_list[[2]])[1:4])){
      term_res <- term_res[-which(rownames(term_res) %in% rownames(term_res_list[[1]])[1:4] | rownames(term_res) %in% rownames(term_res_list[[2]])[1:4]),]}
  }
  term_res_list[[term_name]] <- term_res 
}

melt(lapply(term_res_list, function(x){rownames(x)[1:4]}))


genes <- melt(lapply(term_res_list, function(x){rownames(x)[1:4]}))
colnames(genes) <- c("gene", "term_name")
rownames(genes) <- as.character(genes$gene)
head(genes)
# df <- counts(Rhesus_GeneObject_cortex, normalized=T)
vsd <- varianceStabilizingTransformation(Rhesus_GeneObject_cortex, blind=FALSE)
vsd_data <- assay(vsd)
colData <- Rhesus_colData[colnames(vsd_data),]
df <- vsd_data[rownames(genes), ]
dim(df)

pp <- melt(df)
colnames(pp) <- c("gene", "sample_id", "expression")
pp$ID <- colData[as.character(pp$sample_id), "ID"]
pp$gene_name <- gtf_ensembl_gene[as.character(pp$gene), "gene_name"]
pp$term_name <- genes[as.character(pp$gene), "term_name"]
pp$term_name <- factor(pp$term_name, levels=unique(pp$term_name))
pp$gene_name <- factor(pp$gene_name, levels=unique(pp$gene_name))

head(pp)
str(pp)


ggplot(pp, aes(x=gene_name, y=expression, fill=ID, colour=ID)) +
  geom_violin(alpha=0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.7, color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(color = "black"),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white", size=1)) +
  scale_colour_manual(values = c(anno_colors$ID)) +
  scale_fill_manual(values = c(anno_colors$ID)) +
  facet_wrap(~term_name, nrow=1, scales="free") +
  # annotate("text", x = 4, y = 8, label = paste("Mid: ", Mid_term_name, "\nYoung: ",Young_term_name, sep="")) +
  xlab("") +
  guides(fill=F, colour=F) +
  ylab("expression(VST)") +
  ggtitle("")


ggsave(paste("revise0530/gene_count_res/Cortex_ID_term_top4_genes_violin.pdf", sep=""),
       width=12, height = 4)



#################################################################
#                  coefficients (individual)                    #
#################################################################
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_GeneObject_cortex.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")


###
res <- results(Rhesus_GeneObject_cortex)
res <- res[order(res$padj),]
res[1:10,]
mcols(Rhesus_GeneObject_cortex)
coef <- coef(Rhesus_GeneObject_cortex)[, ]



##
Age_table <- list()
for(r in names(Age_permutation_genes)){
  object <- Age_res[[r]]
  Mid <- Age_permutation_genes[[r]][["Mid"]][, c(2, 5)]
  Young <- Age_permutation_genes[[r]][["Young"]][, c(2, 5)]
  coef <- coef(object)[, c("Age_Stage_Young_vs_Mid", "Sex_male_vs_female")]
  if(dim(Mid)[1]>0){Mid <- data.frame(region=r, hyper="Mid", Mid, coef[rownames(Mid),, drop=F])}
  if(dim(Young)[1]>0){Young <- data.frame(region=r, hyper="Young", Young, coef[rownames(Young),, drop=F])}
  Age_table[[r]] <- rbind(Mid, Young)
}


Age_table <- Reduce(rbind, Age_table)
Age_table <- data.frame(gene=rownames(Age_table), gene_name=gtf_ensembl_gene[rownames(Age_table), "gene_name"], Age_table)
write.table(Age_table, "revise0530/gene_count_res/Age_table.txt",
            row.names = F, col.names = T, sep="\t", quote = F)




###################################################################
#      compare to expression in different tissue of Macaque       #
###################################################################
library(DESeq2)
library(BiocParallel)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/cortex_LRT_res.RData")
cortex_LRT_res <- res
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")

##
RhesusGeneCount <- read.table("revise0530/Rhesus_count.matrix", header=T, row.names = 1)
dim(RhesusGeneCount) # 30807   408
Rhesus_tissue <- read.table("Macaque_tissue_data/matrix/Rhesus_tissue_count.matrix",
                             row.names = 1, header=T)
dim(Rhesus_tissue) # 30807    10
##
GeneCount <- data.frame(RhesusGeneCount, Rhesus_tissue)
GeneCount <- GeneCount[apply(GeneCount, 1, function(x){sum(x>3)>1}), ]
dim(GeneCount) # 22274   418

colData <- sapply(colnames(Rhesus_tissue), function(x){strsplit(x, "[.]")[[1]][5]})
colData <- data.frame(tissue=colData)
colData <- rbind(data.frame(tissue=setNames(Rhesus_colData[,"SR_merge"], rownames(Rhesus_colData))), 
                 colData)

##
count <- GeneCount
Rhesus_tissue_object <- DESeqDataSetFromMatrix(count, 
                                       colData, 
                                       ~tissue)
Rhesus_tissue_object <- DESeq(Rhesus_tissue_object, 
                              parallel = T, 
                              BPPARAM=MulticoreParam(40))
save(Rhesus_tissue_object, file="revise0530/gene_count_res/Rhesus_tissue_objectata.RData")


tissue_norCounts <- counts(Rhesus_tissue_object, normalized=TRUE)
vsd <- varianceStabilizingTransformation(tissue_norCounts, blind=FALSE)
vsd_data <- assay(vsd)



### CV : coefficient of variation
sd <- rowSds(tissue_norCounts)
mean <- rowMeans(tissue_norCounts)
cv <- sd/mean

pdf("revise0530/gene_count_res/tissue_CV_hist.pdf")
hist(cv)
dev.off()


head(cv)
cv <- sort(cv, decreasing = T)
cv_top1000 <- cv[1:1000]


##
cortex_LRT_res <- data.frame(cortex_LRT_res)
cortex_LRT_res <- cortex_LRT_res[order(cortex_LRT_res$padj),]
cortex_LRT_res <- cortex_LRT_res[cortex_LRT_res$padj<0.05,]
dim(cortex_LRT_res) # 17411     6
genes <- rownames(cortex_LRT_res)[1:1000]

##
inter <- intersect(names(cv_top1000), genes)
gtf_ensembl_gene[inter, "gene_name"]



###################################################################
#        compare to expression in different tissue of human       #
###################################################################
library(CePa)
setwd("/mnt/data2/Rhesus_brain")

###
human_count <- read.gct("human/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct")
colnames(human_count) <- gsub("[.]", "_", colnames(human_count))
###
human_colData <- read.table("human/GTEx_v7_Annotations_SampleAttributesDS.txt", 
                            sep="\t", row.names = 1, header = T, quote = '')
rownames(human_colData) <- gsub("-", "_", rownames(human_colData))
human_colData <- human_colData[colnames(human_count),]
###
length(unique(sapply(colnames(human_count), function(x){strsplit(x, "_")[[1]][2]})))  # 714 persons
table(sapply(colnames(human_count), function(x){strsplit(x, "_")[[1]][2]}))
human_colData$sample_id <- sapply(rownames(human_colData), function(x){strsplit(x, "_")[[1]][2]})

unique(human_colData$SMTS)
test <- table(human_colData$SMTS, human_colData$sample_id)
test <- as.array(test)
# a <- sort(rowSums(test), decreasing=T)
# b <- sort(colSums(test), decreasing=T)
# test <- test[names(a), names(b)]



# tissues <- c("Breast", "Colon", "Esophagus", "Heart", "Kidney", "Liver", 
#   "Lung", "Muscle", "Ovary", "Pancreas", "Prostate", "Spleen",
#   "Stomach", "Testis", "Thyroid", "Uterus", "Vagina")

nTissue <- 15
colData <- human_colData
colData <- colData[-which(colData$SMTS=="Brain"),]
tissues <- sort(table(colData$SMTS), decreasing = T)
tissues <- names(tissues)[1:nTissue]
colData <- colData[colData$SMTS %in% tissues, ] 
test <- table(colData$SMTS, colData$sample_id)
test <- as.array(test)
samples <- colnames(test)[apply(test, 2, function(x){sum(x>0)>=nTissue})]
colData <- colData[colData$sample_id %in% samples, ] 
length(unique(colData$sample_id))# 10个 individual,
dim(colData) # 227  63
count <- human_count[, rownames(colData)]
rownames(count) <- as.vector(sapply(rownames(count), function(x){strsplit(x, "[.]")[[1]][1]}))
dim(count) # 56202   227


### DESeq
library(DESeq2)
library(BiocParallel)

count <- count[apply(count, 1, function(x){sum(x>1)>3}),]
dim(count) # 44822   227
human_object <- DESeqDataSetFromMatrix(count, 
                                       colData, 
                                       ~SMTS+sample_id)
human_object <- DESeq(human_object, 
                      test = "LRT",
                      reduced = ~SMTS,
                      parallel = T, 
                      BPPARAM=MulticoreParam(40))
save(human_object, 
     file="revise0530/gene_count_res/human_object_tissue.RData")

###########
res <- results(human_object)  
res <- res[complete.cases(res),]
res <- res[order(res$padj),]
res <- res[res$padj<0.05,]
dim(res) # 18343     6
human_tissue_LRT_res <- res
save(human_tissue_LRT_res, file="revise0530/gene_count_res/human_tissue_LRT_res.RData")

human_norCounts <- counts(human_object, normalized=TRUE)
vsd <- varianceStabilizingTransformation(human_object, blind=FALSE)
vsd_data <- assay(vsd)


#################
library(reshape2)
library(pheatmap)
library(RColorBrewer)
load("revise0530/gene_count_res/human_tissue_LRT_res.RData")
load("revise0530/gene_count_res/cortex_LRT_res.RData")
Rhesus_cortex_LRT_res <- res
load("revise0530/gene_count_res/human_object_tissue.RData")
load("revise0530/gene_count_res/ConsensusCluster_LRT_DifGenes_ID_in_cortex_top1000.RData")

##
orth_genes <- read.table("ANNOTATION/humanGRCH38_Rhesus_orth_genes_one2one.txt", stringsAsFactors = F)
colnames(orth_genes) <- c("Human_gene_id", "Rhesus_gene_id", "Homologous_type")
rownames(orth_genes) <- orth_genes$Rhesus_gene_id
head(orth_genes)
dim(orth_genes) # 19810     3

##
Rhesus_cortex_LRT_res <- Rhesus_cortex_LRT_res[complete.cases(Rhesus_cortex_LRT_res),]
Rhesus_cortex_LRT_res <- Rhesus_cortex_LRT_res[order(Rhesus_cortex_LRT_res$padj),]
Rhesus_cortex_LRT_res <- Rhesus_cortex_LRT_res[Rhesus_cortex_LRT_res$padj<0.05,]
dim(Rhesus_cortex_LRT_res) # 17411
sub_Rhesus_cortex_LRT_res <- Rhesus_cortex_LRT_res[rownames(Rhesus_cortex_LRT_res) %in% orth_genes$Rhesus_gene_id, ]
dim(sub_Rhesus_cortex_LRT_res) # 13897     6
rownames(sub_Rhesus_cortex_LRT_res) <- orth_genes[rownames(sub_Rhesus_cortex_LRT_res), "Human_gene_id"]
##
sub_human_tissue_LRT_res <- human_tissue_LRT_res[rownames(human_tissue_LRT_res) %in% orth_genes$Human_gene_id,]
dim(sub_human_tissue_LRT_res) # 9545     6
##
inter <- intersect(rownames(sub_Rhesus_cortex_LRT_res)[1:1000], 
                   rownames(sub_human_tissue_LRT_res)[1:1000])
length(inter) # 93


sub_norCounts <- human_norCounts[rownames(human_norCounts) %in% sub_orth_genes$Human_gene_id, ]
sub_vsd_data <- vsd_data[rownames(vsd_data) %in% sub_orth_genes$Human_gene_id, ]


##
## scale function
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

# scale_df <- scale_mat(df, "row")
scale_df <- scale_mat(sub_vsd_data, "row")
scale_df[scale_df>4] <- 3

ann_colors <- list(sample_id=colorRampPalette(brewer.pal(12,"Set3"))(15),
                   SMTS=colorRampPalette(brewer.pal(8,"Dark2"))(16))
names(ann_colors$sample_id) <- unique(colData$sample_id)
names(ann_colors$SMTS) <- unique(colData$SMTS)

pheatmap(scale_df,
        filename="revise0530/gene_count_res/ID_dif_genes_in_human_tissue_expression.pdf",
        col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
        breaks = c(-3,seq(-2,2,length=20),3),
        legend_breaks = c(-2,-1,1,2),
        show_rownames = F,
        show_colnames = F,
        fontsize = 8,
        annotation_colors = ann_colors,
        annotation_col=colData[, c("sample_id", "SMTS")])




###################################################################
#                  LRT for ID (exclude cortex)                    #
###################################################################
## run DESeq2 on R console
library(DESeq2)
library(BiocParallel)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/RhesusGeneCount.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
str(Rhesus_colData)

####
genes <- gtf_ensembl_gene[-which(gtf_ensembl_gene$seqnames %in% c("X", "Y")),]$gene_id
RhesusGeneCount <- RhesusGeneCount[intersect(genes, rownames(RhesusGeneCount)), 
                                   rownames(Rhesus_colData)[-which(Rhesus_colData$SR_cortex %in% "cortex")]]
order <- apply(RhesusGeneCount, 1, function(x){sum(x>1)>3})
RhesusGeneCount <- RhesusGeneCount[order, ]
dim(RhesusGeneCount) # 21507   160
###
Rhesus_colData <- Rhesus_colData[colnames(RhesusGeneCount),]
sub_Rhesus_GeneObject <- DESeqDataSetFromMatrix(RhesusGeneCount, 
                                                   Rhesus_colData, 
                                                   ~ID+SR)
sub_Rhesus_GeneObject <- DESeq(sub_Rhesus_GeneObject, 
                                  test="LRT",
                                  reduced=~SR,
                                  parallel = T, 
                                  BPPARAM=MulticoreParam(40))
save(sub_Rhesus_GeneObject, 
     file="revise0530/gene_count_res/sub_Rhesus_GeneObject.RData")


##
res0 <- results(sub_Rhesus_GeneObject)
res0 <- res0[complete.cases(res0),]
res0 <- res0[order(res0$padj),]
res0 <- res0[res0$padj<0.05,]
dim(res0) # 15461     6
genes0 <- rownames(res0)[1:1000]


##
load("revise0530/gene_count_res/cortex_LRT_res.RData")
resID <- res
resID <- data.frame(resID)
resID <- resID[order(resID$padj),]
resID <- resID[resID$padj<0.05,]
dim(resID) # 17411     6
genes <- rownames(resID)[1:1000]


##
inter <- intersect(rownames(res0), rownames(resID))
inter <- intersect(genes0, genes)
length(inter)



#################################################################
#                        HIP LRT test                           #
#################################################################
library(DESeq2)
library(pheatmap)
library(BiocParallel)
library(reshape2)
library(ConsensusClusterPlus)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/gene_count_res/RhesusGeneCount.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")

sub_RhesusGeneCount <- RhesusGeneCount[, rownames(Rhesus_colData)[Rhesus_colData$SR_merge=="HIP"]]
sub_RhesusGeneCount <- sub_RhesusGeneCount[apply(sub_RhesusGeneCount, 1, function(x){sum(x>1)>3}), ]
sub_colData <- Rhesus_colData[colnames(sub_RhesusGeneCount),]

Rhesus_GeneObject_HIP <- DESeqDataSetFromMatrix(sub_RhesusGeneCount, 
                                                sub_colData, 
                                                   ~ID+SR)
Rhesus_GeneObject_HIP <- DESeq(Rhesus_GeneObject_HIP, 
                               test = "LRT",
                               reduced = ~ID,
                               parallel = T, 
                               BPPARAM=MulticoreParam(40))
save(Rhesus_GeneObject_HIP, 
     file="revise0530/gene_count_res/Rhesus_GeneObject_HIP.RData")


### dif genes
library(dendsort)
library(RColorBrewer)
load("revise0530/gene_count_res/Rhesus_GeneObject_HIP.RData")

object <- Rhesus_GeneObject_HIP
region_hyper_genes <- list()
for(r1 in levels(object$SR)){
  hyper_genes <- list()
  for(r2 in levels(object$SR)[-which(levels(object$SR)==r1)]){
    res <- results(object, contrast = c("SR", r1, r2))
    res <- res[complete.cases(res),]
    res <- res[order(res$padj),]
    res <- res[res$padj<0.05 & res$log2FoldChange > 1,]
    hyper_genes[[r2]] <- rownames(res)
  }
  hyper_genes <- Reduce(intersect, hyper_genes)
  region_hyper_genes[[r1]] <- hyper_genes
}


genes <- melt(region_hyper_genes)
colnames(genes) <- c("gene_id", "region")
genes$gene_id <- as.character(genes$gene_id)
save(genes, file="revise0530/gene_count_res/HIP_genes.RData")


##### hyper genes expression
load("revise0530/gene_count_res/HIP_genes.RData")
Rhesus_norCounts <- counts(Rhesus_GeneObject_HIP, normalized=TRUE)
dim(Rhesus_norCounts) # 18858    24
df <- Rhesus_norCounts
dfheat <- df[genes$gene_id,]
dim(dfheat) #  177  24

###
colInfo <- Rhesus_colData[colnames(df),c("SR"), drop=F]
colInfo <- colInfo[order(colInfo$SR),,drop=F]
head(colInfo)
str(colInfo)

ann_colors <- list(SR = setNames(brewer.pal(9, "Pastel1")[1:3], unique(colInfo$SR)))
gene_order <- rownames(geneCluster)


## scale function
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

# scale_df <- scale_mat(df, "row")
scale_df <- scale_mat(dfheat, "row")
scale_df[scale_df>4] <- 3


filename= paste("revise0530/gene_count_res/LRT_log2FC_HIP_profile.pdf", sep="")
pheatmap(
  # dfheat,
  scale_df[, rownames(colInfo)],
  col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
  breaks = c(-3,seq(-2,2,length=20),3),
  legend_breaks = c(-2,-1,1,2),
  filename = filename,
  # scale="row",
  # main = title,
  # width = 10,
  # height = 10,
  # clustering_method = "complete",
  border_color=NA,
  # fontsize_col = 0.5,
  # annotation_row = geneCluster,
  annotation_col = colInfo,
  annotation_colors = ann_colors,
  gaps_row = cumsum(table(genes$region)),
  annotation_legend = T,
  show_rownames = F,
  show_colnames = F,
  cluster_rows = F,
  cluster_cols = F)




######################### function 
library(gProfileR)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/HIP_genes.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


gprofiler_res_list <- list()
for(r in unique(genes$region)){
  gprofiler_res <- gprofiler(
    query=genes$gene_id[genes$region==r], 
    organism = "mmulatta")
  gprofiler_res <- gprofiler_res[order(gprofiler_res[,3]),]
  gprofiler_res$gene_name <- gprofiler_res$intersection
  if(dim(gprofiler_res)[1]==0) next
  else{for(i in 1:dim(gprofiler_res)[1]){
    x <- gprofiler_res$gene_name[i]
    gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],"gene_name"], collapse=",")
  }
    gprofiler_res_list[[r]] <- gprofiler_res
    filename=paste("revise0530/gene_count_res/HIP_LRT_", r,"_function.txt", sep="")
    write.table(gprofiler_res, filename,
                quote = F, sep="\t", row.names = F, col.names = T)
  }
}


############
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/HIP_genes.RData")
load("revise0530/gene_count_res/Rhesus_GeneObject_CB.RData")

###
object <- Rhesus_GeneObject_HIP
res <- results(object)

dfheat <- dfheat[, rownames(colInfo)]
mean <- sapply(unique(colInfo$SR), function(x){rowMeans(dfheat[,colInfo$SR==x])})
head(mean)
mean <- data.frame(gene_name=gtf_ensembl_gene[rownames(mean), "gene_name"],
                   genes, mean)
mean <- data.frame(mean, res[rownames(mean), 
                             c("log2FoldChange", "pvalue", "padj")])

write.table(mean, "revise0530/gene_count_res/HIP_LRT_gene_expression.txt", 
            row.names = F, col.names = T, quote = F, sep="\t")



#################################################################
#                        STR LRT test                           #
#################################################################
library(DESeq2)
library(pheatmap)
library(BiocParallel)
library(reshape2)
library(ConsensusClusterPlus)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/gene_count_res/RhesusGeneCount.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")

sub_RhesusGeneCount <- RhesusGeneCount[, rownames(Rhesus_colData)[Rhesus_colData$SR_merge=="STR"]]
sub_RhesusGeneCount <- sub_RhesusGeneCount[apply(sub_RhesusGeneCount, 1, function(x){sum(x>1)>3}), ]
sub_colData <- Rhesus_colData[colnames(sub_RhesusGeneCount),]

Rhesus_GeneObject_STR <- DESeqDataSetFromMatrix(sub_RhesusGeneCount, 
                                                sub_colData, 
                                                ~ID+SR)
Rhesus_GeneObject_STR <- DESeq(Rhesus_GeneObject_STR, 
                               test = "LRT",
                               reduced = ~ID,
                               parallel = T, 
                               BPPARAM=MulticoreParam(40))
save(Rhesus_GeneObject_STR, 
     file="revise0530/gene_count_res/Rhesus_GeneObject_STR.RData")


### dif genes
library(dendsort)
library(RColorBrewer)
load("revise0530/gene_count_res/Rhesus_GeneObject_STR.RData")

###
object <- Rhesus_GeneObject_STR
region_hyper_genes <- list()
for(r1 in levels(object$SR)){
  hyper_genes <- list()
  for(r2 in levels(object$SR)[-which(levels(object$SR)==r1)]){
    res <- results(object, contrast = c("SR", r1, r2))
    res <- res[complete.cases(res),]
    res <- res[order(res$padj),]
    res <- res[res$padj<0.05 & res$log2FoldChange > 1,]
    hyper_genes[[r2]] <- rownames(res)
  }
  hyper_genes <- Reduce(intersect, hyper_genes)
  region_hyper_genes[[r1]] <- hyper_genes
}


genes <- melt(region_hyper_genes)
colnames(genes) <- c("gene_id", "region")
genes$gene_id <- as.character(genes$gene_id)
save(genes, file="revise0530/gene_count_res/STR_genes.RData")



##### hyper genes expression
load("revise0530/gene_count_res/STR_genes.RData")
Rhesus_norCounts <- counts(Rhesus_GeneObject_STR, normalized=TRUE)
dim(Rhesus_norCounts) # 18858    24
df <- Rhesus_norCounts
dfheat <- df[genes$gene_id,]
dim(dfheat) # 741  24


###
colInfo <- Rhesus_colData[colnames(df),c("SR"), drop=F]
colInfo <- colInfo[order(colInfo$SR),,drop=F]
head(colInfo)
str(colInfo)

ann_colors <- list(SR = setNames(brewer.pal(9, "Pastel1")[1:3], unique(colInfo$SR)))



## scale function
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

# scale_df <- scale_mat(df, "row")
scale_df <- scale_mat(dfheat, "row")
scale_df[scale_df>4] <- 3


filename= paste("revise0530/gene_count_res/LRT_log2FC_STR_profile.pdf", sep="")
pheatmap(
  # dfheat,
  scale_df[, rownames(colInfo)],
  col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
  breaks = c(-3,seq(-2,2,length=20),3),
  legend_breaks = c(-2,-1,1,2),
  filename = filename,
  # scale="row",
  # main = title,
  # width = 10,
  # height = 10,
  # clustering_method = "complete",
  border_color=NA,
  # fontsize_col = 0.5,
  # annotation_row = geneCluster,
  annotation_col = colInfo,
  annotation_colors = ann_colors,
  gaps_row = cumsum(table(genes$region)),
  annotation_legend = T,
  show_rownames = F,
  show_colnames = F,
  cluster_rows = F,
  cluster_cols = F)




######################### function 
library(gProfileR)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/STR_genes.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


gprofiler_res_list <- list()
for(r in unique(genes$region)){
  gprofiler_res <- gprofiler(
    query=genes$gene_id[genes$region==r], 
    organism = "mmulatta")
  gprofiler_res <- gprofiler_res[order(gprofiler_res[,3]),]
  gprofiler_res$gene_name <- gprofiler_res$intersection
  if(dim(gprofiler_res)[1]==0) next
  else{for(i in 1:dim(gprofiler_res)[1]){
    x <- gprofiler_res$gene_name[i]
    gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],"gene_name"], collapse=",")
  }
    gprofiler_res_list[[r]] <- gprofiler_res
    filename=paste("revise0530/gene_count_res/STR_LRT_", r,"_function.txt", sep="")
    write.table(gprofiler_res, filename,
                quote = F, sep="\t", row.names = F, col.names = T)
  }
}



######################### function plot
library(ggplot2)
## ACb
cluster_res <- gprofiler_res_list[["ACb"]]
cluster_res <- cluster_res[cluster_res$domain %in% c("BP", "keg"), c(3,10,12)]
rownames(cluster_res) <- NULL
ACb <- cluster_res[c(1:5),]
ACb <- data.frame(ACb, region="ACb")
## "CN"
cluster_res <- gprofiler_res_list[["CN"]]
cluster_res <- cluster_res[cluster_res$domain %in% c("BP", "keg"), c(3,10,12)]
rownames(cluster_res) <- NULL
CN <- cluster_res[c(1:5),]
CN <- data.frame(CN, region="CN")
## "PUT"
cluster_res <- gprofiler_res_list[["PUT"]]
cluster_res <- cluster_res[cluster_res$domain %in% c("BP", "keg"), c(3,10,12)]
rownames(cluster_res) <- NULL
PUT <- cluster_res[c(1:5),]
PUT <- data.frame(PUT, region="PUT")

##
pp <- rbind(ACb, CN, PUT)
pp <- pp[order(pp$p.value, decreasing = T),]
pp$term.name <- factor(pp$term.name, levels = unique(pp$term.name))

ggplot(pp, aes(x=term.name, y=-log10(p.value))) +
  geom_bar(stat="identity", fill="#C1A5CE") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.7, color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(color = "black"),
        # strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white", color="white", size=1)) +
  # scale_colour_manual(values = c(male=brewer.pal(9,"Blues")[8], female=brewer.pal(9,"Reds")[8])) +
  # scale_fill_manual(values = c(male=brewer.pal(9,"Blues")[8], female=brewer.pal(9,"Reds")[8])) +
  facet_grid(region~., scales="free", space="free_y") +
  coord_flip() +
  xlab("") +
  ggtitle("")

ggsave("revise0530/gene_count_res/STR_LRT_function_bar.pdf",
       width = 8, height = 5)



############
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/STR_genes.RData")
load("revise0530/gene_count_res/Rhesus_GeneObject_STR.RData")

###
object <- Rhesus_GeneObject_STR
res <- results(object)

dfheat <- dfheat[, rownames(colInfo)]
mean <- sapply(unique(colInfo$SR), function(x){rowMeans(dfheat[,colInfo$SR==x])})
head(mean)
mean <- data.frame(gene_name=gtf_ensembl_gene[rownames(mean), "gene_name"],
                   genes, mean)
mean <- data.frame(mean, res[rownames(mean), 
                             c("log2FoldChange", "pvalue", "padj")])

write.table(mean, "revise0530/gene_count_res/STR_LRT_gene_expression.txt", 
            row.names = F, col.names = T, quote = F, sep="\t")




#################################################################
#                         CB LRT test                           #
#################################################################
library(DESeq2)
library(pheatmap)
library(BiocParallel)
library(reshape2)
library(ConsensusClusterPlus)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/gene_count_res/RhesusGeneCount.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


sub_RhesusGeneCount <- RhesusGeneCount[, rownames(Rhesus_colData)[Rhesus_colData$SR_merge=="CB"]]
sub_RhesusGeneCount <- sub_RhesusGeneCount[apply(sub_RhesusGeneCount, 1, function(x){sum(x>1)>3}), ]
sub_colData <- Rhesus_colData[colnames(sub_RhesusGeneCount),]

Rhesus_GeneObject_CB <- DESeqDataSetFromMatrix(sub_RhesusGeneCount, 
                                                sub_colData, 
                                                ~ID+SR)
Rhesus_GeneObject_CB <- DESeq(Rhesus_GeneObject_CB, 
                               test = "LRT",
                               reduced = ~ID,
                               parallel = T, 
                               BPPARAM=MulticoreParam(40))
save(Rhesus_GeneObject_CB, 
     file="revise0530/gene_count_res/Rhesus_GeneObject_CB.RData")


### dif genes
library(dendsort)
library(RColorBrewer)
load("revise0530/gene_count_res/Rhesus_GeneObject_CB.RData")

object <- Rhesus_GeneObject_CB
region_hyper_genes <- list()
for(r1 in levels(object$SR)){
  hyper_genes <- list()
  for(r2 in levels(object$SR)[-which(levels(object$SR)==r1)]){
    res <- results(object, contrast = c("SR", r1, r2))
    res <- res[complete.cases(res),]
    res <- res[order(res$padj),]
    res <- res[res$padj<0.05 & res$log2FoldChange > 1,]
    hyper_genes[[r2]] <- rownames(res)
  }
  hyper_genes <- Reduce(intersect, hyper_genes)
  region_hyper_genes[[r1]] <- hyper_genes
}


genes <- melt(region_hyper_genes)
colnames(genes) <- c("gene_id", "region")
genes$gene_id <- as.character(genes$gene_id)
save(genes, file="revise0530/gene_count_res/CB_genes.RData")



##### hyper genes expression
load("revise0530/gene_count_res/CB_genes.RData")
Rhesus_norCounts <- counts(Rhesus_GeneObject_CB, normalized=TRUE)
dim(Rhesus_norCounts) # 18858    16
df <- Rhesus_norCounts
dfheat <- df[genes$gene_id,]


###
colInfo <- Rhesus_colData[colnames(df),c("SR"), drop=F]
colInfo <- colInfo[order(colInfo$SR),,drop=F]
head(colInfo)
str(colInfo)

ann_colors <- list(SR = setNames(brewer.pal(9, "Pastel1")[1:2], unique(colInfo$SR)))



## scale function
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

# scale_df <- scale_mat(df, "row")
scale_df <- scale_mat(dfheat, "row")
scale_df[scale_df>4] <- 3

gene_name <- gtf_ensembl_gene[rownames(scale_df), "gene_name"]
gene_name[which(is.na(gene_name))] <- rownames(scale_df)[which(is.na(gene_name))]
rownames(scale_df) <- gene_name

filename= paste("revise0530/gene_count_res/LRT_log2FC_CB_profile.pdf", sep="")
pheatmap(
  # dfheat,
  scale_df[, rownames(colInfo)],
  col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
  breaks = c(-3,seq(-2,2,length=20),3),
  legend_breaks = c(-2,-1,1,2),
  filename = filename,
  # scale="row",
  # main = title,
  # width = 10,
  # height = 10,
  # clustering_method = "complete",
  border_color=NA,
  # fontsize_col = 0.5,
  # annotation_row = geneCluster,
  annotation_col = colInfo,
  annotation_colors = ann_colors,
  gaps_row = cumsum(table(genes$region)),
  annotation_legend = T,
  # show_rownames = F,
  show_colnames = F,
  cluster_rows = F,
  cluster_cols = F)




######################### function 
library(gProfileR)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/CB_genes.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


gprofiler_res_list <- list()
for(r in unique(genes$region)){
  gprofiler_res <- gprofiler(
    query=genes$gene_id[genes$region==r], 
    organism = "mmulatta")
  gprofiler_res <- gprofiler_res[order(gprofiler_res[,3]),]
  gprofiler_res$gene_name <- gprofiler_res$intersection
  if(dim(gprofiler_res)[1]==0) next
  else{for(i in 1:dim(gprofiler_res)[1]){
    x <- gprofiler_res$gene_name[i]
    gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],"gene_name"], collapse=",")
  }
    gprofiler_res_list[[r]] <- gprofiler_res
    filename=paste("revise0530/gene_count_res/CB_LRT_", r,"_function.txt", sep="")
    write.table(gprofiler_res, filename,
                quote = F, sep="\t", row.names = F, col.names = T)
  }
}




############
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/CB_genes.RData")
load("revise0530/gene_count_res/Rhesus_GeneObject_CB.RData")

###
object <- Rhesus_GeneObject_CB
res <- results(object)

dfheat <- dfheat[, rownames(colInfo)]
mean <- sapply(unique(colInfo$SR), function(x){rowMeans(dfheat[,colInfo$SR==x])})
head(mean)
mean <- data.frame(gene_name=gtf_ensembl_gene[rownames(mean), "gene_name"],
                   genes, mean)
mean <- data.frame(mean, res[rownames(mean), 
                             c("log2FoldChange", "pvalue", "padj")])

write.table(mean, "revise0530/gene_count_res/CB_LRT_gene_expression.txt", 
            row.names = F, col.names = T, quote = F, sep="\t")



#########################################################
#                    Age specific genes                 #
#########################################################
###### run DEseq2 and results
library(DESeq2)
library(BiocParallel)
library(gtools)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/RhesusGeneCount.RData")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
str(Rhesus_colData)


colData <- colData(Rhesus_GeneObject)
Age_res <- list()
Age_genes <- list()
for(r in unique(colData$SR_merge)){
  sub_colData <- colData[colData$SR_merge == r,]
  df <- RhesusGeneCount[, rownames(sub_colData)]
  df <- df[apply(df, 1, function(x){sum(x>1)>3}), ]
  if(length(unique(sub_colData$SR))>1){
    Object <- DESeqDataSetFromMatrix(df,
                                     sub_colData,
                                     ~Age_Stage+Sex+SR)
    Object <- DESeq(Object,
                    test = "LRT",
                    reduced = ~Sex+SR,
                    parallel = T,
                    BPPARAM=MulticoreParam(10))
  }else{
    Object <- DESeqDataSetFromMatrix(df,
                                     sub_colData,
                                     ~Age_Stage+Sex)
    Object <- DESeq(Object,
                    test = "LRT",
                    reduced = ~Sex,
                    parallel = T,
                    BPPARAM=MulticoreParam(10))
  }
  res <- results(Object,
                 contrast = c("Age_Stage", c("Mid", "Young")))
  res <- res[complete.cases(res),]
  hyper <- res[res$padj<0.05 & res$log2FoldChange>1, ]
  hypo <- res[res$padj<0.05 & res$log2FoldChange < -1, ]
  Age_res[[r]] <- Object
  Age_genes[[r]] <- list(Mid=hyper, Young=hypo)
}

save(Age_res, Age_genes, 
     file = "revise0530/gene_count_res/Age_res.RData")




########### plot the count of differential genes for Age_Stage
library(reshape2)
library(ggplot2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Age_res.RData")

Age_genes_InEachSR <- lapply(Age_genes, function(x){list(Mid=rownames(x[["Mid"]]), Young=rownames(x[["Young"]]))})
Age_genes_count <- melt(Age_genes_InEachSR)
colnames(Age_genes_count) <- c("gene", "Age", "SR")
Age_genes_count1 <- table(Age_genes_count$SR, Age_genes_count$Age)
Age_genes_count1 <- melt(Age_genes_count1)
colnames(Age_genes_count1) <- c("SR", "Age", "counts")
head(Age_genes_count1)

Age_genes_count2 <- table(Age_genes_count$SR)
Age_genes_count2 <- Age_genes_count2[order(Age_genes_count2, decreasing=T)]
head(Age_genes_count2)
Age_genes_count1$SR <- factor(Age_genes_count1$SR, levels=(names(Age_genes_count2)))


### plot
pp <- Age_genes_count1
ggplot(pp,aes(x=SR, y=counts, fill=Age)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = c(Mid=brewer.pal(9,"Blues")[8], Young=brewer.pal(9,"Reds")[8])) +
  ggtitle("the counts of DEG for age")

ggsave("revise0530/gene_count_res/Age_DEG_counts_bar_InSR.pdf",
       width=6, height = 5)




#####################################################################
#                  Age specific genes pheatmap                      #
#####################################################################
library(DESeq2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Age_res.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
# Rhesus_norCounts_ex1 <- counts(Rhesus_GeneObject_ex1, normalized=TRUE)


#################
## scale function
#################
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

#####
for(region in c("cortex", "HIP", "PIT", "MO", "STR", "CB", "PA")){
  Object <- Age_res[[region]]
  vsd <- varianceStabilizingTransformation(Object, blind=FALSE)
  vsd_data <- assay(vsd)
  ####
  rowInfo <-  melt(list(Mid=rownames(Age_genes[[region]][["Mid"]]), 
                        Young=rownames(Age_genes[[region]][["Young"]])))
  colnames(rowInfo) <- c("gene", "Age")
  rownames(rowInfo) <- rowInfo$gene
  head(rowInfo)
  ####
  colData <- colData(Object)
  colInfo <- data.frame(colData[, c("Sex", "Age_Stage")])
  head(colInfo)
  str(colInfo)
  
  ann_colors <- list(Age_Stage = anno_colors$Age_Stage,
                      Sex =anno_colors$Sex)
  
  dfheat <- vsd_data[rownames(rowInfo), rownames(colInfo)]
  scale_df <- scale_mat(dfheat, "row")
  scale_df[scale_df>4] <- 3
  filename = paste("revise0530/gene_count_res/Age_DEG_profile_in_", region, ".pdf", sep="")
  pheatmap(
    dfheat,
    col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
    breaks = c(-3,seq(-2,2,length=20),3),
    # legend_breaks = c(-2,-1,1,2),
    filename = filename,
    main = region,
    scale="row",
    # main = main,
    width = 10,
    height = 10,
    # clustering_method = "single",
    border_color=NA,
    # fontsize = 0.7,
    annotation_col = colInfo,
    annotation_colors = ann_colors,
    gaps_row = cumsum(table(rowInfo$Age)[c("Mid", "Young")]),
    annotation_legend = T,
    show_rownames = F,
    show_colnames = T,
    cluster_rows = F 
  )
}






################################################################
#     Age specific genes (permutations null distribution)      #
################################################################
## 同样适用于Sex
library(DESeq2)
library(reshape2)
library(BiocParallel)
library(gtools)
library(coin)
library(lmPerm)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/RhesusGeneCount.RData")
str(Rhesus_colData)


Permutaion_object_in_eachSR <- list()
Permutaion_res_in_eachSR <- list()
for(r in unique(Rhesus_colData$SR_merge)[16:16]){
  colData <- Rhesus_colData[Rhesus_colData$SR_merge == r, ]
  df <- RhesusGeneCount[, rownames(colData)]
  df <- df[apply(df, 1, function(x){sum(x>1)>3}), ]
  ###
  Permutaion_object_in_eachSR[[r]] <- list()
  Permutaion_res_in_eachSR[[r]] <- list()
  Age_test <- colData$Age_Stage
  Sex_test <- colData$Sex
  if(length(unique(colData$SR))>1){
    for(i in 1:20){
      colData$Age_test <- sample(Age_test)
      colData$Sex_test <- sample(Sex_test)
      #####
      Object <- DESeqDataSetFromMatrix(df, colData, ~Age_test+Sex_test+SR)
      Object <- DESeq(Object, 
                      test = "LRT",
                      reduced = ~Sex_test+SR,
                      parallel = T,
                      BPPARAM=MulticoreParam(30))
      res <- results(Object)[, "pvalue",drop=F]
      Permutaion_object_in_eachSR[[r]][i] <- Object
      Permutaion_res_in_eachSR[[r]][[i]] <- res
     }
    }else{
      for(i in 1:20){
        colData$Age_test <- sample(Age_test)
        colData$Sex_test <- sample(Sex_test)
        #####
        Object <- DESeqDataSetFromMatrix(df, colData, ~Age_test+Sex_test)
        Object <- DESeq(Object, 
                        test = "LRT",
                        reduced = ~Sex_test,
                        parallel = T,
                        BPPARAM=MulticoreParam(30))
        res <- results(Object)[, "pvalue",drop=F]
        Permutaion_object_in_eachSR[[r]][i] <- Object
        Permutaion_res_in_eachSR[[r]][[i]] <- res
      }
  }
    save(Permutaion_object_in_eachSR, Permutaion_res_in_eachSR,
       file="revise0530/gene_count_res/Permutaion20_res.RData")
}


save(Permutaion_object_in_eachSR, Permutaion_res_in_eachSR,
     file="revise0530/gene_count_res/Permutaion20_res.RData")




################################################################
#         Age pvalue adjusted using permutations pvalue        #
################################################################
library(fdrci)
library(DESeq2)
library(reshape2)
library(ggplot2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Permutaion20_res.RData")
load("revise0530/gene_count_res/Age_res.RData")


fdr_od_res <- list()
fdrTbl_res <- list()
pvalue_res <- list()
for(r in names(Permutaion_res_in_eachSR)){
  perm_list <- Permutaion_res_in_eachSR[[r]]
  perm_list <- lapply(perm_list, as.data.frame)
  observe <- Age_res[[r]]
  observe_pvalue <- results(observe)[, "pvalue"]
  pp <- data.frame(observe_pvalue, perm_list)
  pp <- pp[complete.cases(pp),]
  pp <- melt(pp)
  pp$variable <- as.character(pp$variable)
  pp$variable[pp$variable != "observe_pvalue"] <- "random_pvalue"
  pvalue_res[[r]] <- data.frame(pp, region=r)
  ######
  test <- fdrTbl(observe_pvalue, perm_list, "pvalue", length(observe_pvalue), 1, 50)
  test2 <- fdr_od(observe_pvalue, perm_list, "pvalue", length(observe_pvalue), 
                  thres=0.05, cl=.95)
  fdr_od_res[[r]] <- test2
  fdrTbl_res[[r]] <- test
}

save(pvalue_res, fdr_od_res, fdrTbl_res, 
     file="revise0530/gene_count_res/Age_fdrci_res.RData")


##### pvalue distribution
sample_pvalue_res <- lapply(pvalue_res, 
                            function(x){x[c(1:sum(x$variable=="observe_pvalue"), sample(which(x$variable=="random_pvalue"), sum(x$variable=="observe_pvalue"))), ]})
pp <- Reduce(rbind, sample_pvalue_res)

ggplot(pp, aes(x= value, fill= variable, color=variable)) +
  geom_histogram(position="identity", alpha=0.3) +
  # geom_line(stat = "density") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        # legend.position = c(0.8,0.8),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  facet_wrap(~region, ncol=4, scales = "free") +
  scale_color_manual(values=c("#F25F5C", "#247BA0")) +
  scale_fill_manual(values=c("#F25F5C", "#247BA0")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue for Age") 

# ggsave("revise0530/gene_count_res/Age_permutation_pvalue_density.pdf",
#        height = 7, width = 11)

ggsave("revise0530/gene_count_res/Age_permutation_pvalue_hist.pdf",
       height = 7, width = 11)


####
Age_permutation_genes <- list()
for(r in names(fdrTbl_res)){
  # pdf(paste("revise0530/gene_count_res/", r, "_FDRplot.pdf", sep=""))
  # FDRplot(fdrTbl_res[[r]],0,3,"My FDR Plot",lpos = "bottomleft")
  # dev.off()
  object <- Age_res[[r]]
  res <- results(object, contrast=c("Age_Stage", "Mid", "Young"))
  res <- res[complete.cases(res),]
  res$logP <- -log10(res$pvalue)
  res <- res[order(res$logP, decreasing = T),]
  ####
  fdrres <- fdrTbl_res[[r]]
  fdrres <- fdrres[complete.cases(fdrres),]
  fdrres <- fdrres[fdrres$fdr < 0.1,]
  if(dim(fdrres)[1]==0) next
  else{
    logP <- fdrres[1, "threshold"]
    ##
    hyper <- res[res$logP>logP & res$log2FoldChange>1, ]
    hypo <- res[res$logP>logP & res$log2FoldChange < -1, ]
    Age_permutation_genes[[r]] <- list(Mid=hyper, Young=hypo)
  }
 }

save(Age_permutation_genes, 
     file="revise0530/gene_count_res/Age_permutation_genes.RData")


for(r in names(Age_permutation_genes)){
  res <- Age_permutation_genes[[r]]
  for(state in c("Mid", "Young")){
    if(dim(res[[state]])[1]==0) next
    a <- data.frame(hyper=state, res[[state]])
    gene_name <- gtf_ensembl_gene[rownames(a), "gene_name"]
    a <- data.frame(gene_id=rownames(a), gene_name=gene_name, a)
    write.table(a, paste("revise0530/gene_count_res/Age_", r, "_", state, "_hyper_genes.txt", sep=""),
                row.names = F, col.names = T, quote = F, sep="\t")
  }
}




########### plot the count of differential genes for Age_Stage
library(reshape2)
library(ggplot2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Age_permutation_genes.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


Age_genes_InEachSR <- lapply(Age_permutation_genes, function(x){list(Mid=rownames(x[["Mid"]]), Young=rownames(x[["Young"]]))})
Age_genes_count <- melt(Age_genes_InEachSR)
colnames(Age_genes_count) <- c("gene", "Age", "SR")
Age_genes_count1 <- table(Age_genes_count$SR, Age_genes_count$Age)
Age_genes_count1 <- melt(Age_genes_count1)
colnames(Age_genes_count1) <- c("SR", "Age", "counts")
head(Age_genes_count1)

Age_genes_count2 <- table(Age_genes_count$SR)
Age_genes_count2 <- Age_genes_count2[order(Age_genes_count2, decreasing=T)]
head(Age_genes_count2)
Age_genes_count1$SR <- factor(Age_genes_count1$SR, levels=(names(Age_genes_count2)))


### plot
pp <- Age_genes_count1
ggplot(pp,aes(x=SR, y=counts, fill=Age)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = c(Mid=brewer.pal(9,"Blues")[8], Young=brewer.pal(9,"Reds")[8])) +
  ggtitle("the counts of DEG for age (after permutation)")

ggsave("revise0530/gene_count_res/Age_DEG_counts_bar_InSR_permutation_adjust.pdf",
       width=6, height = 5)




#####################################################################
#         Age specific genes pheatmap (after permutations)          #
#####################################################################
library(DESeq2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/gene_count_res/Age_permutation_genes.RData")
# Rhesus_norCounts_ex1 <- counts(Rhesus_GeneObject_ex1, normalized=TRUE)


#################
## scale function
#################
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

#####
for(region in c("cortex", "HIP", "PIT","STR", "CB")){
  Object <- Age_res[[region]]
  vsd <- varianceStabilizingTransformation(Object, blind=FALSE)
  vsd_data <- assay(vsd)
  ####
  rowInfo <-  melt(list(Mid=rownames(Age_permutation_genes[[region]][["Mid"]]), 
                        Young=rownames(Age_permutation_genes[[region]][["Young"]])))
  colnames(rowInfo) <- c("gene", "Age")
  rownames(rowInfo) <- rowInfo$gene
  head(rowInfo)
  ####
  colData <- colData(Object)
  colInfo <- data.frame(colData[, c("Sex", "Age_Stage")])
  head(colInfo)
  str(colInfo)
  
  ann_colors <- list(Age_Stage = anno_colors$Age_Stage,
                     Sex =anno_colors$Sex)
  
  dfheat <- vsd_data[rownames(rowInfo), rownames(colInfo)]
  scale_df <- scale_mat(dfheat, "row")
  scale_df[scale_df>4] <- 3
  filename = paste("revise0530/gene_count_res/Age_DEG_profile_in_", region, "_after_permutation.pdf", sep="")
  pheatmap(
    dfheat,
    col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
    breaks = c(-3,seq(-2,2,length=20),3),
    # legend_breaks = c(-2,-1,1,2),
    filename = filename,
    main = region,
    scale="row",
    # main = main,
    width = 10,
    height = 10,
    # clustering_method = "single",
    border_color=NA,
    # fontsize = 0.7,
    annotation_col = colInfo,
    annotation_colors = ann_colors,
    gaps_row = cumsum(table(rowInfo$Age)[c("Mid", "Young")]),
    annotation_legend = T,
    show_rownames = F,
    show_colnames = T,
    cluster_rows = F 
  )
}



#############################################################
#      Age specific genes fuction (after permutations)      #
#############################################################
library(gProfileR)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Age_permutation_genes.RData")


gprofiler_res_list <- list()
for(region in c("cortex", "PIT", "HIP", "STR", "CB")){
  for(age in c("Mid", "Young")){
    genes <- rownames(Age_permutation_genes[[region]][[age]])
    gprofiler_res <- gprofiler(query=genes, organism = "mmulatta")
    gprofiler_res$gene_name <- gprofiler_res$intersection
    if(dim(gprofiler_res)[1]==0) next
    else{for(i in 1:dim(gprofiler_res)[1]){
      x <- gprofiler_res$gene_name[i]
      gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],"gene_name"], collapse=",")
    }
      gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),]
      filename=paste("revise0530/gene_count_res/", region,"_", age, "_hyper_function_permutation.txt", sep="")
      write.table(gprofiler_res, filename,
                  quote = F, sep="\t", row.names = F, col.names = T)
      gprofiler_res_list[[region]][[age]] <- gprofiler_res
    }
  }
}

save(gprofiler_res_list,
     file="revise0530/gene_count_res/Age_permutation_function_gprofiler_res_list.RData")


#####
library(ggplot2)
load("revise0530/gene_count_res/Age_permutation_function_gprofiler_res_list.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")


### cortex
region <- "cortex"
gprofiler_res <- gprofiler_res_list[[region]][["Mid"]]
gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),c(3, 10,12)]
rownames(gprofiler_res) <- NULL
cortex_Mid <- gprofiler_res[c(1,4:7),]
cortex_Mid <- data.frame(cortex_Mid, Age="Mid", region="cortex")

gprofiler_res <- gprofiler_res_list[[region]][["Young"]]
gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),c(3, 10,12)]
rownames(gprofiler_res) <- NULL
cortex_Young <- gprofiler_res[1:4,]
cortex_Young <- data.frame(cortex_Young, Age="Young", region="cortex")



### PIT
# region <- "PIT"
# gprofiler_res <- gprofiler_res_list[[region]][["Mid"]]
# gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),c(3, 10,12)]
# rownames(gprofiler_res) <- NULL
# PIT_Mid <- gprofiler_res[c(1:5),]
# PIT_Mid <- data.frame(PIT_Mid, Age="Mid", region="PIT")
# 
# gprofiler_res <- gprofiler_res_list[[region]][["Young"]]
# gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),c(3, 10,12)]
# rownames(gprofiler_res) <- NULL
# PIT_Young <- gprofiler_res[c(1:3),]
# PIT_Young <- data.frame(PIT_Young, Age="Young", region="PIT")


### HIP
region <- "HIP"
gprofiler_res <- gprofiler_res_list[[region]][["Mid"]]
gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),c(3, 10,12)]
rownames(gprofiler_res) <- NULL
HIP_Mid <- gprofiler_res[c(1,2,3,4,6),]
HIP_Mid <- data.frame(HIP_Mid, Age="Mid", region="HIP")

gprofiler_res <- gprofiler_res_list[[region]][["Young"]]
gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),c(3, 10,12)]
rownames(gprofiler_res) <- NULL
HIP_Young <- gprofiler_res[c(1:5),]
HIP_Young <- data.frame(HIP_Young, Age="Young", region="HIP")



### STR
region <- "STR"
gprofiler_res <- gprofiler_res_list[[region]][["Mid"]]
gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),c(3, 10,12)]
rownames(gprofiler_res) <- NULL
STR_Mid <- gprofiler_res[c(1),]
STR_Mid <- data.frame(STR_Mid, Age="Mid", region="STR")

gprofiler_res <- gprofiler_res_list[[region]][["Young"]]
gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),c(3, 10,12)]
rownames(gprofiler_res) <- NULL
STR_Young <- gprofiler_res[c(2:4),]
STR_Young <- data.frame(STR_Young, Age="Young", region="STR")


pp <- rbind(
  # cortex_Mid, cortex_Young
  # MO_Mid, MO_Young
  HIP_Mid, HIP_Young
  # STR_Mid, STR_Young
)
order <- pp$term.name[order(pp$p.value, decreasing=T)]
pp$term.name <- factor(pp$term.name, levels = order)

ggplot(pp, aes(x=term.name, y=-log10(p.value), fill=Age)) +
  geom_bar(stat="identity", width = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        # panel.border = element_rect(size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        # legend.position = c(0.8,0.8),
        strip.text = element_text(color="black", size = 15),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  # geom_text(aes(label=symbol), vjust=0.5, hjust=0,size=4, colour="black") +
  # annotate("text", x=pp$NAME[2], y=-1.5,size=3,
  #          label="nominal P values: \n*P < 0.05 \n**P < 0.005 \n***P < 0.0005 \nNS: not significant") +
  # scale_x_discrete(limits=c(pp$term.name))+
  coord_flip() +
  xlab("") +
  facet_grid(Age~region, scales = "free_y", space = "free_x") + 
  scale_fill_manual(values = anno_colors$Age)  +
  guides(fill=F)


ggsave("revise0530/gene_count_res/HIP_Age_function_plot.pdf", 
       height = 5, width = 8.5)



##############################################################################
#                     Age specific : term example genes                      #
##############################################################################
library(ggplot2)
library(DESeq2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Age_permutation_genes.RData")  
load("revise0530/gene_count_res/Age_permutation_function_gprofiler_res_list.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")


vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)

###   
region="cortex"
Mid_term_name_list = list(term1="MHC protein complex",
                          term2="Antigen processing and presentation"
                          # term3="CCR chemokine receptor binding"
                          )
Young_term_name_list <- list(term1="receptor ligand activity",
                             term2="receptor regulator activity"
                             # term3="regulation of signaling receptor activity"
                             )
# region="HIP"
# Mid_term_name_list = list(term1="integral component of membrane",
#                           term2="signaling receptor activity")
# Young_term_name_list <- list(term1="Cell adhesion molecules (CAMs)",
#                              term2="Graft-versus-host disease")

## Mid
Mid <- gprofiler_res_list[[region]][["Mid"]]
Mid_res_list <- list()
for(term in names(Mid_term_name_list)){
  Mid_term_name <- Mid_term_name_list[[term]]
  Mid_intersection <- strsplit(Mid$intersection[Mid$term.name %in% Mid_term_name], ",")[[1]]
  Mid_res <- data.frame(Age_permutation_genes[[region]][["Mid"]])[Mid_intersection,]
  Mid_res <- Mid_res[order(Mid_res$padj),]
  Mid_res$gene_name <- gtf_ensembl_gene[rownames(Mid_res), "gene_name"]
  Mid_res <- Mid_res[complete.cases(Mid_res$gene_name),]
  if(term=="term2"){Mid_res <- Mid_res[-which(rownames(Mid_res) %in% rownames(Mid_res_list[[1]])[1:3]),]}
  # if(term=="term3"){Mid_res <- Mid_res[-which(rownames(Mid_res) %in% rownames(Mid_res_list[[1]])[1:3] | rownames(Mid_res) %in% rownames(Mid_res_list[[2]])[1:3]),]}
  Mid_res_list[[Mid_term_name]] <- Mid_res
}

Mid_top_genes <- lapply(Mid_res_list, function(x)(rownames(x)[1:3]))
Mid_top_genes <- melt(Mid_top_genes)
colnames(Mid_top_genes) <- c("gene", "term_name")
Mid_top_genes <- data.frame(Mid_top_genes, Age="Mid")
Mid_top_genes <- Mid_top_genes[complete.cases(Mid_top_genes$gene),]
dim(Mid_top_genes) # 9 3



## Young
Young <- gprofiler_res_list[[region]][["Young"]]
Young_res_list <- list()
for(term in names(Young_term_name_list)){
  Young_term_name <- Young_term_name_list[[term]]
  Young_intersection <- strsplit(Young$intersection[Young$term.name == Young_term_name], ",")[[1]]
  Young_res <- data.frame(Age_res_InEachSR[[region]])[Young_intersection,]
  Young_res <- Young_res[order(Young_res$padj),]
  Young_res$gene_name <- gtf_ensembl_gene[rownames(Young_res), "gene_name"]
  Young_res <- Young_res[complete.cases(Young_res$gene_name),]
  if(term=="term2"){Young_res <- Young_res[-which(rownames(Young_res) %in% rownames(Young_res_list[[1]])[1:3]),]}
  # if(term=="term3"){Young_res <- Young_res[-which(rownames(Young_res) %in% rownames(Young_res_list[[1]])[1:3] | rownames(Young_res) %in% rownames(Young_res_list[[2]])[1:3]),]}
  Young_res_list[[Young_term_name]] <- Young_res
}

Young_top_genes <- lapply(Young_res_list, function(x)(rownames(x)[1:3]))
Young_top_genes <- melt(Young_top_genes)
colnames(Young_top_genes) <- c("gene", "term_name")
Young_top_genes <- data.frame(Young_top_genes, Age="Young")
Young_top_genes <- Young_top_genes[complete.cases(Young_top_genes$gene),]
dim(Young_top_genes) # 9 3


##
genes <- rbind(Mid_top_genes, Young_top_genes)
rownames(genes) <- genes$gene
head(genes)
colData <- Rhesus_colData[colnames(vsd_data),]
colData <- colData[colData$SR_merge == region, ]
df <- vsd_data[rownames(genes), rownames(colData)]
dim(df)

pp <- melt(df)
colnames(pp) <- c("gene", "sample_id", "expression")
pp$Age <- colData[as.character(pp$sample_id), "Age_stage"]
pp$gene_name <- gtf_ensembl_gene[as.character(pp$gene), "gene_name"]
pp$term_name <- genes[as.character(pp$gene), "term_name"]
pp$term_name <- factor(pp$term_name, levels=unique(pp$term_name))

head(pp)
str(pp)


ggplot(pp, aes(x=gene_name, y=expression, fill=Age, colour=Age)) +
  geom_violin(alpha=0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.7, color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(color = "black"),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white", size=1)) +
  scale_colour_manual(values = c(Mid=brewer.pal(9,"Blues")[8], Young=brewer.pal(9,"Reds")[8])) +
  scale_fill_manual(values = c(Mid=brewer.pal(9,"Blues")[8], Young=brewer.pal(9,"Reds")[8])) +
  facet_wrap(~term_name, ncol=3, scales="free_x") +
  # annotate("text", x = 4, y = 8, label = paste("Mid: ", Mid_term_name, "\nYoung: ",Young_term_name, sep="")) +
  xlab("") +
  ylab("expression(VST)") +
  ggtitle(region)


ggsave(paste("matrix0420/gene_count_res/Age_", region, "_term_top3_genes_violin.pdf", sep=""),
       width=8, height = 5)




##############################################################################
#                                Age marker genes                            #
##############################################################################
library(ggplot2)
library(DESeq2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Age_res.RData")
load("revise0530/gene_count_res/Age_permutation_genes.RData")


vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)

genes <- c("CD74", "DRA", "DRB1")
gene_ids <- gtf_ensembl_gene$gene_id[gtf_ensembl_gene$gene_name %in% genes] 
data <- vsd_data[gene_ids,]
rownames(data) <- genes
pp <- data.frame(t(data), Rhesus_colData[colnames(data),c("Age_Stage", "SR_merge")])
pp <- pp[pp$SR_merge %in% c("cortex", "PIT", "HIP", "STR", "CB"),]
p <- melt(pp)
colnames(p)[3:4] <- c("gene", "expression")

ggplot(p, aes(x=SR_merge, y=expression, colour=Age_Stage, fill=Age_Stage)) +
  # geom_violin()+
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.7, color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(color = "black"),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white", size=1)) +
  scale_colour_manual(values = c(Mid=brewer.pal(9,"Blues")[8], Young=brewer.pal(9,"Reds")[8])) +
  scale_fill_manual(values = c(Mid=brewer.pal(9,"Blues")[8], Young=brewer.pal(9,"Reds")[8])) +
  facet_wrap(~gene, ncol=3, scales="free_y") +
  scale_x_discrete(limits=c("PIT", "STR", "CB", "HIP", "cortex"))+
  xlab("") +
  guides(fill=F, colour=F)

ggsave("revise0530/gene_count_res/Age_marker_genes_expression_boxplot.pdf", width = 7, height = 3)


genes <- c("CD74", "DRA", "DRB1")
pvalue_list <- list()
for(gene_name in genes){
  pvalue_list[[gene_name]] <- list()
  for(r in c("PIT", "STR", "CB", "HIP", "cortex")){
    res <- Age_permutation_genes[[r]]
    res <- Reduce(rbind, res)
    gene_id <- gtf_ensembl_gene$gene_id[gtf_ensembl_gene$gene_name %in% gene_name] 
    if(gene_id %in% rownames(res)){
      pvalue_list[[gene_name]][[r]] <- res[gene_id, "pvalue"]
    }

  }
}

  
  
  
##################################################################
#                      Age selected genes                        #
##################################################################
library(ggplot2)
library(DESeq2)
library(RColorBrewer)
library(reshape2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Age_res.RData")
load("revise0530/gene_count_res/Age_permutation_genes.RData")


vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)

selected_genes_list <- list(PIT_Mid_Hyper = c("CDKN2A","RXFP1","TLE4"),
                            PIT_Young_Hyper = c("GH1", "PRLR", "PLAGL1"),
                            cortex_Mid_hyper=c("HSD11B1", "APLNR", "IDO1"),
                            cortex_Young_hyper=c("TERT", "Srpx2", "MMP3"),
                            STR_Mid_hyper=c("GLI1", "SLC22A2", "IL33"),
                            STR_Young_hyper=c("NPTX2", "MCHR1", "NR4A1"),
                            HIP_Mid_Hyper = c("CD44","GALR1","TLR3"),
                            HIP_Young_Hyper = c("DRD2","DRD4","TLR9"))


selected_genes_expr_list <- list()
for(state in names(selected_genes_list)){
  region = strsplit(state, "_")[[1]][1]
  genes <- selected_genes_list[[state]]
  gene_ids <- gtf_ensembl_gene$gene_id[gtf_ensembl_gene$gene_name %in% genes] 
  samples <- rownames(Rhesus_colData)[Rhesus_colData$SR_merge %in% region]
  data <- vsd_data[gene_ids,samples]
  rownames(data) <- genes
  pp <- data.frame(t(data), Rhesus_colData[colnames(data),c("Age_Stage", "SR_merge")])
  p <- melt(pp)
  colnames(p)[3:4] <- c("gene", "expression")
  p <- data.frame(p, state=state)
  selected_genes_expr_list[[state]] <- p
}

p <- Reduce(rbind, selected_genes_expr_list)



ggplot(p, aes(x=gene, y=expression, colour=Age_Stage, fill=Age_Stage)) +
  # geom_violin()+
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.7, color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(color = "black"),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white", size=1)) +
  scale_colour_manual(values = c(Mid=brewer.pal(9,"Blues")[8], Young=brewer.pal(9,"Reds")[8])) +
  scale_fill_manual(values = c(Mid=brewer.pal(9,"Blues")[8], Young=brewer.pal(9,"Reds")[8])) +
  facet_wrap(~state, ncol=4, scales="free") +
  xlab("") +
  guides(fill=F, colour=F)

ggsave("revise0530/gene_count_res/Age_selected_hyper_genes_expression_boxplot.pdf", 
       width = 10, height = 6)




###########################################################
#                    Sex specific genes                   #
###########################################################
###### run DEseq2 and results
library(DESeq2)
library(BiocParallel)
library(gtools)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/RhesusGeneCount.RData")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
str(Rhesus_colData)


colData <- colData(Rhesus_GeneObject)
Sex_res <- list()
Sex_genes <- list()
for(r in unique(colData$SR_merge)){
  sub_colData <- colData[colData$SR_merge == r,]
  df <- RhesusGeneCount[, rownames(sub_colData)]
  df <- df[apply(df, 1, function(x){sum(x>1)>3}), ]
  str(Rhesus_colData)
  genes <- gtf_ensembl_gene[-which(gtf_ensembl_gene$seqnames %in% c("X", "Y")),]$gene_id
  df <- df[rownames(df) %in% genes, ]
  ####
  if(length(unique(sub_colData$SR))>1){
    Object <- DESeqDataSetFromMatrix(df, 
                                     sub_colData, 
                                     ~Age_Stage+Sex+SR)
    Object <- DESeq(Object, 
                    test = "LRT",
                    reduced = ~Age_Stage+SR,
                    parallel = T, 
                    BPPARAM=MulticoreParam(10))
  }else{
    Object <- DESeqDataSetFromMatrix(df, 
                                     sub_colData, 
                                     ~Age_Stage+Sex)
    Object <- DESeq(Object, 
                    test = "LRT",
                    reduced = ~Age_Stage,
                    parallel = T, 
                    BPPARAM=MulticoreParam(10))
  }
  res <- results(Object,
                 contrast = c("Sex", c("male", "female")))
  res <- res[complete.cases(res),]
  hyper <- res[res$padj<0.05 & res$log2FoldChange>1, ]
  hypo <- res[res$padj<0.05 & res$log2FoldChange < -1, ]
  Sex_res[[r]] <- Object
  Sex_genes[[r]] <- list(male=hyper, female=hypo)
}


save(Sex_res, Sex_genes, 
     file = "revise0530/gene_count_res/Sex_res.RData")



########### plot the count of differential genes for Sex
library(reshape2)
library(ggplot2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Sex_res.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")


Sex_genes_InEachSR <- lapply(Sex_genes, function(x){list(male=rownames(x[["male"]]), female=rownames(x[["female"]]))})
Sex_genes_count <- melt(Sex_genes_InEachSR)
colnames(Sex_genes_count) <- c("gene", "Sex", "SR")
Sex_genes_count1 <- table(Sex_genes_count$SR, Sex_genes_count$Sex)
Sex_genes_count1 <- melt(Sex_genes_count1)
colnames(Sex_genes_count1) <- c("SR", "Sex", "counts")
head(Sex_genes_count1)

Sex_genes_count2 <- table(Sex_genes_count$SR)
Sex_genes_count2 <- Sex_genes_count2[order(Sex_genes_count2, decreasing=T)]
head(Sex_genes_count2)
Sex_genes_count1$SR <- factor(Sex_genes_count1$SR, levels=(names(Sex_genes_count2)))


### plot
pp <- Sex_genes_count1
ggplot(pp,aes(x=SR, y=counts, fill=Sex)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = anno_colors$Sex) +
  ggtitle("the counts of DEG for Sex")

ggsave("revise0530/gene_count_res/Sex_DEG_counts_bar_InSR.pdf",
       width=6, height = 5)




################################################################
#               Sex specific genes pheatmap                    #
################################################################
library(DESeq2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Sex_res.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
# Rhesus_norCounts <- counts(Rhesus_GeneObject, normalized=TRUE)


#################
## scale function
#################
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

#####
for(region in c("cortex","PIT")){
  Object <- Sex_res[[region]]
  vsd <- varianceStabilizingTransformation(Object, blind=FALSE)
  vsd_data <- assay(vsd)
  ####
  rowInfo <-  melt(list(male=rownames(Sex_genes[[region]][["male"]]), 
                        female=rownames(Sex_genes[[region]][["female"]])))
  colnames(rowInfo) <- c("gene", "Sex")
  rownames(rowInfo) <- rowInfo$gene
  head(rowInfo)
  ####
  colData <- colData(Object)
  colInfo <- data.frame(colData[, c("Sex", "Age_Stage")])
  head(colInfo)
  str(colInfo)
  
  ann_colors <- list(Age_Stage = anno_colors$Age_Stage,
                     Sex =anno_colors$Sex)
  
  dfheat <- vsd_data[rownames(rowInfo), rownames(colInfo)]
  scale_df <- scale_mat(dfheat, "row")
  scale_df[scale_df>4] <- 3
  filename = paste("revise0530/gene_count_res/Sex_DEG_profile_in_", region, ".pdf", sep="")
  pheatmap(
    dfheat,
    col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
    breaks = c(-3,seq(-2,2,length=20),3),
    # legend_breaks = c(-2,-1,1,2),
    filename = filename,
    main = region,
    scale="row",
    # main = main,
    width = 10,
    height = 10,
    # clustering_method = "single",
    border_color=NA,
    # fontsize = 0.7,
    annotation_col = colInfo,
    annotation_colors = ann_colors,
    gaps_row = cumsum(table(rowInfo$Sex)[c("male", "female")]),
    annotation_legend = T,
    show_rownames = F,
    show_colnames = T,
    cluster_rows = F 
  )
}



################################################################
#     Sex specific genes (permutations null distribution)      #
################################################################
library(DESeq2)
library(reshape2)
library(BiocParallel)
library(gtools)
library(coin)
library(lmPerm)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/RhesusGeneCount.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
str(Rhesus_colData)



Permutaion_object_in_eachSR <- list()
Permutaion_res_in_eachSR <- list()
for(r in unique(Rhesus_colData$SR_merge)[5:16]){
  colData <- Rhesus_colData[Rhesus_colData$SR_merge == r, ]
  df <- RhesusGeneCount[, rownames(colData)]
  df <- df[apply(df, 1, function(x){sum(x>1)>3}), ]
  genes <- gtf_ensembl_gene[-which(gtf_ensembl_gene$seqnames %in% c("X", "Y")),]$gene_id
  df <- df[rownames(df) %in% genes, ]
  ###
  Permutaion_object_in_eachSR[[r]] <- list()
  Permutaion_res_in_eachSR[[r]] <- list()
  Age_test <- colData$Age_Stage
  Sex_test <- colData$Sex
  i <- 1
  while(i <= 20){
    colData$Age_test <- sample(Age_test)
    colData$Sex_test <- sample(Sex_test)
    len <- length(unique(paste(colData$Age_test, colData$Sex_test)))
    if(len!=2){
      if(length(unique(colData$SR))>1){
        Object <- DESeqDataSetFromMatrix(df, colData, ~Age_test+Sex_test+SR)
        Object <- DESeq(Object, 
                        test = "LRT",
                        reduced = ~Age_test+SR,
                        parallel = T,
                        BPPARAM=MulticoreParam(30))
        res <- results(Object)[, "pvalue",drop=F]
        Permutaion_object_in_eachSR[[r]][i] <- Object
        Permutaion_res_in_eachSR[[r]][[i]] <- res
      }else{
        #####
        Object <- DESeqDataSetFromMatrix(df, colData, ~Age_test+Sex_test)
        Object <- DESeq(Object, 
                        test = "LRT",
                        reduced = ~Age_test,
                        parallel = T,
                        BPPARAM=MulticoreParam(30))
        res <- results(Object)[, "pvalue",drop=F]
        Permutaion_object_in_eachSR[[r]][i] <- Object
        Permutaion_res_in_eachSR[[r]][[i]] <- res
      }
      i <- i + 1
      save(Permutaion_object_in_eachSR, Permutaion_res_in_eachSR,
           file="revise0530/gene_count_res/Permutaion20_res_for_Sex.RData")
    }else{
      i <- i
    }
  }
}

save(Permutaion_object_in_eachSR, Permutaion_res_in_eachSR,
     file="revise0530/gene_count_res/Permutaion20_res_for_Sex.RData")




#########################################################
#          Sex specific genes (after permutations)      #
#########################################################

##########
library(fdrci)
library(DESeq2)
library(reshape2)
library(ggplot2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Permutaion20_res_for_Sex.RData")
load("revise0530/gene_count_res/Sex_res.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


fdr_od_res <- list()
fdrTbl_res <- list()
pvalue_res <- list()
for(r in names(Permutaion_res_in_eachSR)){
  observe <- Sex_res[[r]]
  observe_pvalue <- results(observe)[, "pvalue", drop=F]
  ##
  perm_list <- Permutaion_res_in_eachSR[[r]]
  perm_list <- lapply(perm_list, function(x){x[rownames(x) %in% rownames(observe_pvalue),,drop=F]})
  ##
  pp <- data.frame(observe_pvalue, perm_list)
  pp <- pp[complete.cases(pp),]
  pp <- melt(pp)
  pp$variable <- as.character(pp$variable)
  pp$variable[pp$variable == "pvalue"] <- "observe_pvalue"
  pp$variable[pp$variable != "observe_pvalue"] <- "random_pvalue"
  pvalue_res[[r]] <- data.frame(pp, region=r)
  ######
  test <- fdrTbl(observe_pvalue[,"pvalue"], perm_list, "pvalue", length(observe_pvalue), 1, 50)
  test2 <- fdr_od(observe_pvalue[,"pvalue"], perm_list, "pvalue", length(observe_pvalue),
                  thres=0.05, cl=.95)
  fdr_od_res[[r]] <- test2
  fdrTbl_res[[r]] <- test
}


save(fdr_od_res, fdrTbl_res, pvalue_res, 
     file = "revise0530/gene_count_res/Sex_fdrci_res.RData")


#####
sample_pvalue_res <- lapply(pvalue_res, 
                            function(x){x[c(1:sum(x$variable=="observe_pvalue"), sample(which(x$variable=="random_pvalue"), sum(x$variable=="observe_pvalue"))), ]})
pp <- Reduce(rbind, sample_pvalue_res)

ggplot(pp, aes(x= value, fill= variable, color=variable)) +
  geom_histogram(position="identity", alpha=0.5) +
  # geom_line(stat = "density") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        # legend.position = c(0.8,0.8),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  facet_wrap(~region, ncol=4, scales = "free") +
  scale_color_manual(values=c("#1B9E77", "#D95F02")) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue for Sex") 


ggsave("revise0530/gene_count_res/Sex_permutation_pvalue_hist.pdf",
       height = 7, width = 11)


##########
load("revise0530/gene_count_res/Sex_res.RData")
Sex_permutation_genes <- list()
for(r in names(fdrTbl_res)){
  # pdf(paste("revise0530/gene_count_res/", r, "_FDRplot.pdf", sep=""))
  # FDRplot(fdrTbl_res[[r]],0,3,"My FDR Plot",lpos = "bottomleft")
  # dev.off()
  object <- Sex_res[[r]]
  res <- results(object, contrast=c("Sex", "male", "female"))
  res <- res[complete.cases(res),]
  res$logP <- -log10(res$pvalue)
  res <- res[order(res$logP, decreasing = T),]
  ####
  fdrres <- fdrTbl_res[[r]]
  fdrres <- fdrres[complete.cases(fdrres),]
  fdrres <- fdrres[fdrres$fdr < 0.1,]
  if(dim(fdrres)[1]==0) next
  else{
    logP <- fdrres[1, "threshold"]
    ##
    hyper <- res[res$logP>logP & res$log2FoldChange>1, ]
    hypo <- res[res$logP>logP & res$log2FoldChange < -1, ]
    Sex_permutation_genes[[r]] <- list(male=hyper, female=hypo)
  }
}

save(Sex_permutation_genes, 
     file="revise0530/gene_count_res/Sex_permutation_genes.RData")


for(r in names(Sex_permutation_genes)){
  res <- Sex_permutation_genes[[r]]
  for(state in c("male", "female")){
    if(dim(res[[state]])[1]==0) next
    a <- data.frame(hyper=state, res[[state]])
    gene_name <- gtf_ensembl_gene[rownames(a), "gene_name"]
    a <- data.frame(gene_id=rownames(a), gene_name=gene_name, a)
    write.table(a, paste("revise0530/gene_count_res/Sex_", r, "_", state, "_hyper_genes.txt", sep=""),
                row.names = F, col.names = T, quote = F, sep="\t")
  }
}


########### plot the count of differential genes for Sex
library(reshape2)
library(ggplot2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Sex_res.RData")
load("revise0530/gene_count_res/Sex_permutation_genes.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")


Sex_genes_InEachSR <- lapply(Sex_permutation_genes, function(x){list(male=rownames(x[["male"]]), female=rownames(x[["female"]]))})
Sex_genes_count <- melt(Sex_genes_InEachSR)
colnames(Sex_genes_count) <- c("gene", "Sex", "SR")
Sex_genes_count1 <- table(Sex_genes_count$SR, Sex_genes_count$Sex)
Sex_genes_count1 <- melt(Sex_genes_count1)
colnames(Sex_genes_count1) <- c("SR", "Sex", "counts")
head(Sex_genes_count1)

Sex_genes_count2 <- table(Sex_genes_count$SR)
Sex_genes_count2 <- Sex_genes_count2[order(Sex_genes_count2, decreasing=T)]
head(Sex_genes_count2)
Sex_genes_count1$SR <- factor(Sex_genes_count1$SR, levels=(names(Sex_genes_count2)))


### plot
pp <- Sex_genes_count1
ggplot(pp,aes(x=SR, y=counts, fill=Sex)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  scale_fill_manual(values = anno_colors$Sex) +
  ggtitle("the counts of DEG for Sex (after permutation)")

ggsave("revise0530/gene_count_res/Sex_DEG_counts_bar_InSR_permutation_adjust.pdf",
       width=6, height = 5)




#####################################################################
#                  Sex specific genes pheatmap                      #
#####################################################################
library(DESeq2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/gene_count_res/Sex_permutation_genes.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
# Rhesus_norCounts_ex1 <- counts(Rhesus_GeneObject_ex1, normalized=TRUE)


#################
## scale function
#################
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}


#####
for(region in c("cortex")){
  Object <- Sex_res[[region]]
  vsd <- varianceStabilizingTransformation(Object, blind=FALSE)
  vsd_data <- assay(vsd)
  ####
  # rowInfo <-  melt(list(male=rownames(Sex_permutation_genes[[region]][["male"]]), 
  #                       female=rownames(Sex_permutation_genes[[region]][["female"]])))
  rowInfo <-  melt(list(male=rownames(Sex_permutation_genes[[region]][["male"]])[1:10],
                        female=rownames(Sex_permutation_genes[[region]][["female"]])[1:10]))
  
  colnames(rowInfo) <- c("gene", "Sex")
  rownames(rowInfo) <- as.character(rowInfo$gene)
  rowInfo$gene_name <- gtf_ensembl_gene[as.character(rowInfo$gene), "gene_name"]
  rowInfo$gene_name[which(is.na(rowInfo$gene_name))] <- as.character(rowInfo$gene)[which(is.na(rowInfo$gene_name))]
  head(rowInfo)
  ####
  colData <- colData(Object)
  colInfo <- data.frame(colData[, c("Sex", "Age_Stage")])
  head(colInfo)
  str(colInfo)
  
  ann_colors <- list(Age_Stage = anno_colors$Age_Stage,
                     Sex =anno_colors$Sex)
  
  dfheat <- vsd_data[rownames(rowInfo), rownames(colInfo)]
  rownames(dfheat) <- rowInfo$gene_name
  scale_df <- scale_mat(dfheat, "row")
  scale_df[scale_df>4] <- 3
  filename = paste("revise0530/gene_count_res/Sex_DEG_profile_in_", region, "_after_permutation.pdf", sep="")
  pheatmap(
    dfheat,
    col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
    breaks = c(-3,seq(-2,2,length=20),3),
    # legend_breaks = c(-2,-1,1,2),
    filename = filename,
    main = region,
    scale="row",
    # main = main,
    width = 10,
    height = 10,
    # clustering_method = "single",
    border_color=NA,
    # fontsize = 0.7,
    annotation_col = colInfo,
    annotation_colors = ann_colors,
    gaps_row = cumsum(table(rowInfo$Sex)[c("male", "female")]),
    annotation_legend = T,
    show_rownames = T,
    show_colnames = F,
    cluster_rows = F 
  )
}




#################################################################
#                 Sex specific genes fuction                    #
#################################################################
library(gProfileR)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Sex_permutation_genes.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")

Sex_function_gprofiler_res_list <- list()
for(region in c("cortex", "HIP")){
  for(Sex in c("female", "male")){
    genes <- rownames(Sex_permutation_genes[[region]][[Sex]])
    gprofiler_res <- gprofiler(query=genes, organism = "mmulatta")
    gprofiler_res$gene_name <- gprofiler_res$intersection
    if(dim(gprofiler_res)[1]==0) next
    else{for(i in 1:dim(gprofiler_res)[1]){
      x <- gprofiler_res$gene_name[i]
      gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],5], collapse=",")
    }
      filename=paste("revise0530/gene_count_res/", region,"_", Sex, "_hyper_function_permutation.txt", sep="")
      write.table(gprofiler_res, filename,
                  quote = F, sep="\t", row.names = F, col.names = T)
      Sex_function_gprofiler_res_list[[region]][[Sex]] <- gprofiler_res
    }
  }
}

save(Sex_function_gprofiler_res_list, 
     file = "revise0530/gene_count_res/Sex_permutation_function_gprofiler_res_list.RData")



##################################################################
#                        Sex top genes                           #
##################################################################
library(ggplot2)
library(DESeq2)
library(RColorBrewer)
library(reshape2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Sex_res.RData")
load("revise0530/gene_count_res/Sex_permutation_genes.RData")


vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)


top_genes_expr_list <- list()
for(region in c("cortex")){
  res_list <- list()
  for(a in c("male", "female")){
    res <- Sex_permutation_genes[[region]][[a]]
    res$gene_name <- gtf_ensembl_gene[rownames(res), "gene_name"]
    res <- res[complete.cases(res), ]
    res <- res[1:10, ]
    samples <- rownames(Rhesus_colData)[Rhesus_colData$SR_merge %in% region]
    data <- vsd_data[rownames(res),samples]
    rownames(data) <- res$gene_name
    pp <- data.frame(t(data), Rhesus_colData[colnames(data),c("Sex", "SR_merge")])
    p <- melt(pp)
    colnames(p)[3:4] <- c("gene", "expression")
    p <- data.frame(p, hyper=a)
    res_list[[a]] <- p
  }
  top_genes_expr_list[[region]] <- Reduce(rbind, res_list)
}


p <- Reduce(rbind, top_genes_expr_list)

ggplot(p, aes(x=gene, y=expression, colour=Sex, fill=Sex)) +
  # geom_violin()+
  geom_boxplot(alpha=0.7)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.7, color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(color = "black"),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white", size=1)) +
  scale_colour_manual(values = c(male="#F8CE7A", female="#BF645D")) +
  scale_fill_manual(values = c(male="#F8CE7A", female="#BF645D")) +
  facet_wrap(~hyper, nrow=1, scales = "free_x") +
  xlab("") +
  ylab("Expression (VST)") +
  guides(fill=F, colour=F)

ggsave("revise0530/gene_count_res/Sex_top_dif_genes_expression_boxplot.pdf", 
       width = 12, height = 3)


## ND3
vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)
gene_name <- "ND3"
gene_id <- gtf_ensembl_gene[gtf_ensembl_gene$gene_name %in% gene_name,]
colData <- Rhesus_colData[Rhesus_colData$SR_merge %in% c("cortex","HIP","STR","CB"),]
df <- vsd_data[rownames(gene_id),rownames(colData)]

pp <- data.frame(expression=df)
pp <- data.frame(pp, colData[rownames(pp), c("SR_merge", "Sex")])
ggplot(pp, aes(x=SR_merge, y=expression, colour=Sex, fill=Sex)) +
  # geom_violin()+
  geom_boxplot(alpha=0.7)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.7, color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(color = "black"),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white", size=1)) +
  scale_colour_manual(values = c(male="#F8CE7A", female="#BF645D")) +
  scale_fill_manual(values = c(male="#F8CE7A", female="#BF645D")) +
  # facet_wrap(~hyper, nrow=2, scales = "free_x") +
  xlab("Region") +
  ylab("Expression (VST)") 
  # guides(fill=F, colour=F)

ggsave("revise0530/gene_count_res/Sex_ND3_expression_boxplot.pdf", 
       width = 6, height = 5)

## CRH
vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)
gene_name <- "CRH"
gene_id <- gtf_ensembl_gene[gtf_ensembl_gene$gene_name %in% gene_name,]
colData <- Rhesus_colData
df <- vsd_data[rownames(gene_id),rownames(colData)]

pp <- data.frame(expression=df)
pp <- data.frame(pp, colData[rownames(pp), c("SR_merge", "Sex")])
ggplot(pp, aes(x=SR_merge, y=expression, colour=Sex, fill=Sex)) +
  # geom_violin()+
  geom_boxplot(alpha=0.7)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.7, color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(color = "black"),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white", size=1)) +
  scale_colour_manual(values = c(male="#F8CE7A", female="#BF645D")) +
  scale_fill_manual(values = c(male="#F8CE7A", female="#BF645D")) +
  # facet_wrap(~hyper, nrow=2, scales = "free_x") +
  xlab("Region") +
  ylab("Expression (VST)") 
# guides(fill=F, colour=F)

ggsave("revise0530/gene_count_res/Sex_CRH_expression_boxplot.pdf", 
       width = 8, height = 5)


Sex_genes_InEachSR <- lapply(Sex_permutation_genes, function(x){list(male=rownames(x[["male"]]), female=rownames(x[["female"]]))})
Sex_genes_count <- melt(Sex_genes_InEachSR)
Sex_genes_count[Sex_genes_count$value %in% gene_id,]



#################################################################
#                  FoldChange (Age and Sex)                    #
#################################################################
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Sex_res.RData")
load("revise0530/gene_count_res/Age_res.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Sex_permutation_genes.RData")
load("revise0530/gene_count_res/Age_permutation_genes.RData")



for(r in names(Age_res)){
  a <- results(Age_res[[r]], contrast = c("Age_Stage", "Mid", "Young"))
  b <- results(Sex_res[[r]], contrast = c("Sex", "male", "female"))
  inter <- intersect(rownames(a), rownames(b))
  f <- data.frame(Age=a[inter,c("log2FoldChange", "pvalue")], 
                  Sex=b[inter,c("log2FoldChange", "pvalue")])
  f$difGenes <- "NO"
  f$difGenes[rownames(f) %in% rownames(Reduce(rbind, Age_permutation_genes[[r]]))] <- "Age"
  f$difGenes[rownames(f) %in% rownames(Reduce(rbind, Sex_permutation_genes[[r]]))] <- "Sex"
  # f <- f[abs(f$Age.log2FoldChange)>=1 | abs(f$Sex.log2FoldChange)>=1, ]
  ###
  ggplot(f, aes(x=abs(Age.log2FoldChange), 
                y=abs(Sex.log2FoldChange), color=difGenes)) +
    geom_point(size=0.9) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                      colour = "black", size = 11),
          axis.text.y  = element_text(colour = "black"),
          strip.background = element_rect(fill="white", colour = "white", size=1)) +
    scale_color_manual(values = c(Age="#F25F5C", NO="#8F908F", Sex="#FFE066")) +
    guides(color = guide_legend(override.aes = list(size=2))) +
    ggtitle(r)
    
  
  ggsave(paste("revise0530/gene_count_res/", r, "_logFC_point.pdf", sep=""),
         width = 5, height = 4.5)
}

             


#################################################################
#                  coefficients (Age and Sex)                   #
#################################################################
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Sex_res.RData")
load("revise0530/gene_count_res/Age_res.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Sex_permutation_genes.RData")
load("revise0530/gene_count_res/Age_permutation_genes.RData")


##
Age_table <- list()
for(r in names(Age_permutation_genes)){
  object <- Age_res[[r]]
  Mid <- Age_permutation_genes[[r]][["Mid"]][, c(2, 5)]
  Young <- Age_permutation_genes[[r]][["Young"]][, c(2, 5)]
  coef <- coef(object)[, c("Age_Stage_Young_vs_Mid", "Sex_male_vs_female")]
  if(dim(Mid)[1]>0){Mid <- data.frame(region=r, hyper="Mid", Mid, coef[rownames(Mid),, drop=F])}
  if(dim(Young)[1]>0){Young <- data.frame(region=r, hyper="Young", Young, coef[rownames(Young),, drop=F])}
  Age_table[[r]] <- rbind(Mid, Young)
}


Age_table <- Reduce(rbind, Age_table)
Age_table <- data.frame(gene=rownames(Age_table), gene_name=gtf_ensembl_gene[rownames(Age_table), "gene_name"], Age_table)
write.table(Age_table, "revise0530/gene_count_res/Age_table.txt",
            row.names = F, col.names = T, sep="\t", quote = F)


##
Sex_table <- list()
for(r in names(Sex_permutation_genes)){
  object <- Sex_res[[r]]
  male <- Sex_permutation_genes[[r]][["male"]][, c(2, 5)]
  female <- Sex_permutation_genes[[r]][["female"]][, c(2, 5)]
  coef <- coef(object)[, c("Age_Stage_Young_vs_Mid", "Sex_male_vs_female")]
  if(dim(male)[1]>0){male <- data.frame(region=r, hyper="male", male, coef[rownames(male),, drop=F])}
  if(dim(female)[1]>0){female <- data.frame(region=r, hyper="female", female, coef[rownames(female),, drop=F])}
  Sex_table[[r]] <- rbind(male, female)
}


Sex_table <- Reduce(rbind, Sex_table)
Sex_table <- data.frame(gene=rownames(Sex_table), gene_name=gtf_ensembl_gene[rownames(Sex_table), "gene_name"], Sex_table)
write.table(Sex_table, "revise0530/gene_count_res/Sex_table.txt",
            row.names = F, col.names = T, sep="\t", quote = F)
















