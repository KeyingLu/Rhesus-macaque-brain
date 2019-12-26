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
                           parallel = T, 
                           BPPARAM=MulticoreParam(40))
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
vsd_data_save <- data.frame(GeneID=rownames(vsd_data), 
                            vsd_data)
write.table(vsd_data_save, "SourceData/SourceData_Fig.1e.txt",
            row.names = F, col.names = T, quote = F, sep="\t")


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


Rhesus_norlCount <- counts(Rhesus_GeneObject, normalized=T)
save(Rhesus_norlCount, file="revise0530/gene_count_res/Rhesus_norlCount.RData")






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

