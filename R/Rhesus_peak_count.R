library(DESeq2)
library(cqn)
library(ggplot2)
library(RColorBrewer)
library(dendextend)
library(gtools)
library(BiocParallel)
library(data.table)
library(Rtsne)



setwd("/mnt/data2/Rhesus_brain/ATAC/matrix")
dir.create("peak_top10000_res")


##############################################################################
#                           loading Rhesus                                   #
##############################################################################
RhesusPeakCount <- read.table("/mnt/data2/Rhesus_brain/ATAC/liu_results/peak_without_summits/peak_top_10k/without_ca1/matrix523/fragment_final.matrix", header=T, row.names = 1)
dim(RhesusPeakCount) # [1] 26442    33
seq <- sapply(colnames(RhesusPeakCount), function(x){strsplit(x, "_")[[1]][3]})
RhesusPeakCount <- RhesusPeakCount[, -which(seq=="46")]
dim(RhesusPeakCount) # [1] 26442    26
chroms <- sapply(rownames(RhesusPeakCount), function(x){strsplit(x, ":")[[1]][1]})
RhesusPeakCount <- RhesusPeakCount[-which(chroms %in% c("X", "Y")),]
dim(RhesusPeakCount) # 25842    26
RhesusPeakCount <- RhesusPeakCount[apply(RhesusPeakCount, 1, function(x){sum(x>1)>3}), ]
dim(RhesusPeakCount) # 25842    26
save(RhesusPeakCount, file="peak_top10000_res/RhesusPeakCount.RData")


# load("peak_top10000_res/RhesusPeakCount.RData")
# write.table(RhesusPeakCount, "peak_top10000_res/RhesusPeakCount.txt", row.names = T, col.names = T, sep="\t", quote = F)



##########################################################
#                        colData                        #
##########################################################
load("/mnt/data2/Rhesus_brain/revise0530/gene_count_res/Rhesus_colData.RData")
Rhesus_peak_colData <- Rhesus_colData[colnames(RhesusPeakCount),]
save(Rhesus_peak_colData, file="peak_top10000_res/Rhesus_peak_colData.RData")



##########################################################
#                        DESeq2                          #
##########################################################
## run DESeq2 on R console
library(DESeq2)
library(BiocParallel)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain/ATAC/matrix")
load("peak_top10000_res/RhesusPeakCount.RData")
load("peak_top10000_res/Rhesus_peak_colData.RData")
str(Rhesus_peak_colData)
Rhesus_peak_colData$BR <- factor(Rhesus_peak_colData$BR)

Rhesus_PeakObject <- DESeqDataSetFromMatrix(RhesusPeakCount, 
                                            Rhesus_peak_colData, 
                                            ~ID+SR)
Rhesus_PeakObject <- DESeq(Rhesus_PeakObject,
                           test = "LRT",
                           reduced = ~ SR,
                           parallel = T, 
                           BPPARAM=MulticoreParam(10))
save(Rhesus_PeakObject, 
     file="peak_top10000_res/Rhesus_DESeq2_object.RData")



##########################################################
#                        t-SNE                          #
##########################################################
library(Rtsne)
library(DESeq2)
library(ggplot2)
setwd("/mnt/data2/Rhesus_brain/ATAC/matrix")
load("peak_top10000_res/Rhesus_peak_colData.RData")
load("peak_top10000_res/Rhesus_DESeq2_object.RData")
load("/mnt/data2/Rhesus_brain/revise0530/gene_count_res/Rhesus_anno_colors.RData")

vsd <- varianceStabilizingTransformation(Rhesus_PeakObject, blind=FALSE)
vsd_data <- assay(vsd)
dim(vsd_data) # 25842    26


for(seed in 1:50){
  seed=29
  set.seed(seed)
  perplexity = 2
  tsne <- Rtsne::Rtsne(t(vsd_data), perplexity=2)
  rownames(tsne$Y) <- colnames(vsd_data)
  colnames(tsne$Y) <- c("tsne1", "tsne2")
  
  
  p <- data.frame(tsne$Y, Rhesus_peak_colData[rownames(tsne$Y),])
  head(p)
  
  title = paste("Rhesus peak t-SNE \n seed:", seed, " perplexity:", perplexity, sep="")
  
  ggplot(p, aes(x=tsne1, y=tsne2, colour=ID, shape=SR)) +
    geom_point(size=4) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black", size=1),
          plot.title = element_text(hjust = 0.5, color="black"),
          axis.title = element_text(color="black"),
          axis.text = element_text(color="black", size = 10),
          legend.title = element_text(color="black")) +
    scale_colour_manual(values=anno_colors$ID) +
    # scale_colour_manual(values=colorRampPalette(brewer.pal(12,"Paired"))(52)) +
    # geom_text(aes(label=region2), vjust=-1.2, size=1, colour="black") +
    xlab("t-SNE1") +
    ylab("t-SNE2") +
    ggtitle(title)
  
  
  filename = paste("peak_top10000_res/Rhesus_Peak_t-SNE_ID_seed",seed,"_per", perplexity," .pdf", sep="")
  ggsave(filename,width = 6, height = 5)
}



################################################################
#         Likelihood ratio test for ID (in cortex)             #
################################################################
library(DESeq2)
library(ConsensusClusterPlus)
setwd("/mnt/data2/Rhesus_brain/ATAC/matrix")
load("peak_top10000_res/Rhesus_peak_colData.RData")
load("peak_top10000_res/Rhesus_DESeq2_object.RData")
load("/mnt/data2/Rhesus_brain/revise0530/gene_count_res/Rhesus_anno_colors.RData")


##############
#### dif peaks
##############
res <- results(Rhesus_PeakObject)
res <- data.frame(res)
res <- res[order(res$padj),]
res <- res[res$padj<0.05,]
dim(res) # 9780    6
peaks <- rownames(res)
length(peaks)
peaks_LRT_res <- res
save(peaks_LRT_res, file="peak_top10000_res/peaks_LRT_res.RData")


##################
##### peak cluster
##################
vsd <- varianceStabilizingTransformation(Rhesus_PeakObject, blind=FALSE)
vsd_data <- assay(vsd)

Rhesus_norCounts <- counts(Rhesus_PeakObject, normalized=TRUE)
dim(Rhesus_norCounts) # 25842    26
df <- Rhesus_norCounts
df <- vsd_data

top_df <- df[peaks,]
ttop_df <- t(top_df)
ttop_df = scale(ttop_df)
title="peak_top10000_res/ConsensusCluster_DifPeaks_0609"
DIfPeaks_Results = ConsensusClusterPlus(ttop_df,maxK=10,reps=100,pItem=0.8,pFeature=1,
                                        title=title,clusterAlg="km",
                                        innerLinkage = "ward.D2",
                                        finalLinkage = "ward.D2",
                                        seed = 1,
                                        distance="euclidean",
                                        plot="png")
# cal <- calcICL(DIfPeaks_Results)
# cal$clusterConsensus
save(DIfPeaks_Results, 
     file="peak_top10000_res/ConsensusCluster_DifPeaks.RData")


##############
#### heatmap
##############
library(RColorBrewer)
library(pheatmap)
library(dendsort)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix/peak_top10000_res/Rhesus_peak_colData.RData")
load("ATAC/matrix/peak_top10000_res/Rhesus_DESeq2_object.RData")
load("ATAC/matrix/peak_top10000_res/ConsensusCluster_DifPeaks.RData")
load("/mnt/data2/Rhesus_brain/revise0530/gene_count_res/Rhesus_anno_colors.RData")



k=4
res <- DIfPeaks_Results[[k]]
peakCluster <- data.frame(peakCluster=factor(sort(res$consensusClass)))


colInfo <- Rhesus_peak_colData[colnames(vsd_data),c("ID", "Age_Stage", "Sex", "SR")]
head(colInfo)
str(colInfo)

ann_colors <- list(Age_Stage = anno_colors$Age_Stage,
                   SR = setNames(brewer.pal(5, "Set3"), unique(colInfo$SR)),
                   Sex = anno_colors$Sex,
                   ID = anno_colors$ID)

peak_order <- rownames(peakCluster)
dfheat <- vsd_data[peak_order, ]


method1 = c("euclidean", "maximum", "manhattan", "canberra", "minkowski")
method2 = c("ward.D", "ward.D2", "single",  "complete", "average")
method3 = c("average","min")
for(i in method1){
  for(j in method2){
    for(z in method3){
      i = "maximum"
      j = "complete"
      z = "average"
      filename= paste("ATAC/matrix/peak_top10000_res/dif_peak_profile_", i, "_", j, "_",z,".pdf", sep="")
      title = paste("distance:", i, " hclust:",j, " type:",z, sep = "")
      pheatmap(
        dfheat,
        col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
        breaks = c(-3,seq(-2,2,length=20),3),
        legend_breaks = c(-2,-1,1,2),
        filename = filename,
        scale="row",
        main = title,
        width = 10,
        height = 10,
        # clustering_method = "complete",
        border_color=NA,
        # fontsize_col = 0.5,
        annotation_col = colInfo,
        annotation_colors = ann_colors,
        gaps_row = cumsum(table(peakCluster$peakCluster)),
        annotation_legend = T,
        show_rownames = F,
        show_colnames = F,
        cluster_rows = F,
        cluster_cols = as.hclust(dendsort(hclust(dist(t(dfheat), method=i), method = j), type=z))
      )
    }
  }
}



####
# library(pheatmap)
# library(dendsort)
# load("ATAC/matrix/ConsensusCluster_DifPeaks.RData")
# k=4
# res <- DIfPeaks_Results[[k]]
# peakCluster <- data.frame(peakCluster=factor(sort(res$consensusClass)))
# 
# correlation <- cor(t(top_df))
# colInfo <- peakCluster
# colInfo$peakCluster <- paste("peakCluster", colInfo$peakCluster, sep="")
# head(colInfo)
# 
# ann_colors <- list(peakCluster=c(peakCluster1="#3385BB",
#                                  peakCluster2="#FE8D19",
#                                  peakCluster3="#6DBD57",
#                                  peakCluster4="#F25F5C"))
# 
# 
# pheatmap(correlation,
#          col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
#          # filename = "ATAC/matrix/dif_peak_correlation_heatmap.pdf",
#          filename = "ATAC/matrix/dif_peak_correlation_heatmap.png",
#          annotation_col = colInfo,
#          annotation_colors = ann_colors,
#          show_rownames = F,
#          show_colnames = F,
#          width=6.5,
#          height=5
# )
# 
# ### help
# min(correlation)
# test <- data.frame(a=c(-0.9680889,-0.9680889), b=c(1,1))
# pheatmap(test,
#          col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
#          # filename = "ATAC/matrix/dif_peak_correlation_heatmap.pdf",
#          filename = "ATAC/matrix/dif_peak_correlation_heatmap_legend_help.pdf",
#          show_rownames = F,
#          show_colnames = F,
# )



################################################################
#                     peakCluster function                     #
################################################################ 
library(gProfileR)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("ATAC/matrix/peak_top10000_res/ConsensusCluster_DifPeaks.RData")
load("ATAC/matrix/peak_top10000_res/dif_peakAnno.RData")

##
k=4
res <- DIfPeaks_Results[[k]]
peakCluster <- data.frame(peakCluster=factor(sort(res$consensusClass)))


peakCluster_function_list <- list()
for(cluster in 1:k){
  peaks <- rownames(peakCluster)[peakCluster$peakCluster==cluster]
  genes <- peakAnno[peaks,]$geneId
  gprofiler_res <- gprofiler(query=genes, organism = "mmulatta")
  gprofiler_res <- gprofiler_res[order(gprofiler_res[,3]),]
  gprofiler_res$gene_name <- gprofiler_res$intersection
  if(dim(gprofiler_res)[1]==0) next
  else{for(i in 1:dim(gprofiler_res)[1]){
    x <- gprofiler_res$gene_name[i]
    gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],6], collapse=",")
  }
    filename=paste("ATAC/matrix/peak_top10000_res/LRT_peakCluster", cluster,"_function.txt", sep="")
    write.table(gprofiler_res, filename,
                quote = F, sep="\t", row.names = F, col.names = T)
  }
  peakCluster_function_list[[cluster]] <- gprofiler_res
}

#### cluster1
cluster <- 1
gprofiler_res <- peakCluster_function_list[[cluster]]
gprofiler_res <- gprofiler_res[gprofiler_res$domain=="BP",]
gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),]
cluster1 <- data.frame(gprofiler_res[1:15,], cluster="peakCluster1")

#### cluster1
cluster <- 3
gprofiler_res <- peakCluster_function_list[[cluster]]
gprofiler_res <- gprofiler_res[gprofiler_res$domain=="BP",]
gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),]
cluster3 <- data.frame(gprofiler_res[1:15,], cluster="peakCluster3")

##
gprofiler_res <- rbind(cluster1, cluster3)
order <- gprofiler_res$term.name[order(gprofiler_res$p.value, decreasing=T)]
gprofiler_res$term.name <- factor(gprofiler_res$term.name, levels=order)
ggplot(gprofiler_res, aes(x=term.name, y=-log10(p.value))) +
  geom_point(aes(size=overlap.size/query.size, colour=-log10(p.value))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.7, face="bold"),
        axis.title = element_text(face="bold"),
        axis.text.x  = element_text(face="bold", angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 10),
        axis.text.y  = element_text(face="bold", colour = "black",size=10),
        panel.border = element_rect(size = 1.5),
        strip.text = element_text(size=10, colour = "black", face="bold", angle=90),
        strip.background = element_rect(fill="white", colour = "white", size=1),legend.title = element_text(face="bold")) +
  coord_flip() +
  scale_color_gradientn(colours = brewer.pal(9,"Reds")[3:9])  +
  facet_grid(cluster~., scales = "free") +
  ylab("-log10") +
  xlab("") +
  labs(size='GeneRatio') 

# filename=paste("anova_peakCluster", cluster,"_function_top20.pdf", sep="")
filename=paste("anova_peakCluster_function_top.pdf", sep="")
ggsave(filename , width = 9, height = 6)




################################################################
#                       LRT  for genes                         #
################################################################
# RNA samples in accordance with ATAC samples
library(DESeq2)
library(BiocParallel)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/RhesusGeneCount.RData")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("ATAC/matrix/peak_top10000_res/Rhesus_peak_colData.RData")


genes <- gtf_ensembl_gene[-which(gtf_ensembl_gene$seqnames %in% c("X", "Y")),]$gene_id
length(genes)
df <- RhesusGeneCount[rownames(RhesusGeneCount) %in% genes, rownames(Rhesus_peak_colData)]
dim(df) # 22700    26
colData <- Rhesus_peak_colData[colnames(df),]
###
Object <- DESeqDataSetFromMatrix(df, colData, ~ID+SR)
Object <- DESeq(Object, 
                test = "LRT",
                reduced = ~SR,
                parallel = T,
                BPPARAM=MulticoreParam(30))
save(Object, file="ATAC/matrix/peak_top10000_res/Rhesus_RNA_Object_26.RData")

res <- results(Object)
res <- res[complete.cases(res),]
res <- res[order(res$padj),]
res <- res[res$padj<0.05,]
dim(res) #  6673    6
head(res)
gene_LRT_res <- res
save(gene_LRT_res, file="ATAC/matrix/peak_top10000_res/gene_LRT_res.RData")



##################################
#### gene差异比较 (RNA samples in accordance with ATAC samples)
###################################
### Rhesus6和7和8 ??? Rhesus10和11和3和4
# load("matrix0420/gene_count_res/Rhesus_DESeq2_object.RData")
# load("matrix0420/gene_count_res/gtf_ensembl_gene.RData")
# load("ATAC/matrix/peak_top10000_res/Rhesus_peak_colData.RData")
# genes <- gtf_ensembl_gene[-which(gtf_ensembl_gene$seqnames %in% c("X", "Y")),]$gene_id
# length(genes)
# vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
# vsd_data <- assay(vsd)
# df <- vsd_data[genes, rownames(Rhesus_peak_colData)]
# dim(df) # 29361    26
# colData <- Rhesus_peak_colData[colnames(df),]
# testInfo <- colData$ID
# testInfo[which(testInfo %in% c("Rhesus6", "Rhesus7", "Rhesus8"))] <- "Rhesus_a"
# testInfo[which(testInfo %in% c("Rhesus10", "Rhesus11", "Rhesus3", "Rhesus4"))] <- "Rhesus_b"
# pvalue <- apply(df, 1, function(x){d <- data.frame(x,testInfo);summary(aov(x~testInfo, data=d))[[1]][1,5]})
# padj <- p.adjust(pvalue, method = "fdr", n=length(pvalue))
# mean <- sapply(unique(testInfo), function(x){rowMeans(df[, testInfo==x])})
# head(mean)
# log2FC <- log2(mean[, "Rhesus_a"]/mean[, "Rhesus_b"])
# anova_res <- data.frame(pvalue, padj, mean, log2FC)
# anova_res <- anova_res[order(anova_res$padj),]
# anova_res <- anova_res[anova_res$padj<0.05,]
# dim(anova_res) #  13  2
# head(anova_res)
# anova_res$state <- anova_res$log2FC
# anova_res$state[which(anova_res$log2FC>0)] <- "Rhesus_a"
# anova_res$state[which(anova_res$log2FC<0)] <- "Rhesus_b"
# gene_anova_res_for_two_group <- anova_res
# gene_anova_res_for_two_group$gene_name <- gtf_ensembl_gene[rownames(gene_anova_res_for_two_group),"gene_name"]
# gene_anova_res_for_two_group[order(gene_anova_res_for_two_group$log2FC, decreasing = T),]




###########################################################
#                       ChIPseeker                        #
###########################################################
library(ChIPseeker)
library(GenomicFeatures)
library(DESeq2)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix/peak_top10000_res/Rhesus_DESeq2_object.RData")

txdb <- makeTxDbFromGFF("/mnt/data1/Ref/Rhesus_macaque_8.0.1/Macaca_mulatta.Mmul_8.0.1.91.chr.gtf",
                      format = "gtf",
                      dataSource = "ensemblgenomes",
                      organism = "Macaca mulatta")

promoter <- getPromoters(TxDb=txdb, upstream=2500, downstream=2500)


### all peaks 
vsd <- varianceStabilizingTransformation(Rhesus_PeakObject, blind=FALSE)
vsd_data <- assay(vsd)
all_peaks <- rownames(vsd_data)
all_peaks <- data.frame(chrom=sapply(all_peaks, function(x){strsplit(x, ":")[[1]][1]}),
                    start=sapply(all_peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1]}),
                    end=sapply(all_peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2]}))
save(all_peaks, file="ATAC/matrix/peak_top10000_res/all_peaks.RData")
all_peaks_save <- data.frame(paste("chr", all_peaks[,1], sep=""), all_peaks[,2:3], rownames(all_peaks))
write.table(all_peaks_save, "ATAC/matrix/peak_top10000_res/all_peaks.bed",
            row.names = F, col.names = F, sep="\t", quote = F)
write.table(all_peaks_save, "ATAC/liftover_res/all_peaks.bed",
            row.names = F, col.names = F, sep="\t", quote = F)
all_peaks <- GRanges(all_peaks)
all_peakAnno  <- annotatePeak(all_peaks, tssRegion=c(-2500, 2500), TxDb = txdb, level = "gene")
### dif peaks 
load("TAC/matrix/peak_top10000_res/peaks_LRT_res.RData")
peaks <- rownames(peaks_LRT_res)
peaks <- data.frame(chrom=sapply(peaks, function(x){strsplit(x, ":")[[1]][1]}),
                    start=sapply(peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1]}),
                    end=sapply(peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2]}))
peaks <- GRanges(peaks)
peakAnno  <- annotatePeak(peaks, tssRegion=c(-2500, 2500), TxDb = txdb, level = "gene")
### cluster peaks
load("ATAC/matrix/peak_top10000_res/ConsensusCluster_DifPeaks.RData")
k=4
res <- DIfPeaks_Results[[k]]
peakCluster <- data.frame(peakCluster=factor(sort(res$consensusClass)))
peakClusterAnno_list <- list()
for(i in 1:k){
  peaks <- rownames(peakCluster)[peakCluster$peakCluster==i]
  peaks <- data.frame(chrom=sapply(peaks, function(x){strsplit(x, ":")[[1]][1]}),
                      start=sapply(peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1]}),
                      end=sapply(peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2]}))
  peaks <- GRanges(peaks)
  peakClusterAnno  <- annotatePeak(peaks, tssRegion=c(-2500, 2500), TxDb = txdb, level = "gene")
  cluster <- paste("peakCluster", i, sep="")
  peakClusterAnno_list[[cluster]] <- peakClusterAnno
}


##########
peakAnno_list <- peakClusterAnno_list
peakAnno_list[["dif_peak"]] <- peakAnno
peakAnno_list[["all_peak"]] <- all_peakAnno
save(peakAnno_list, file="ATAC/matrix/peak_top10000_res/peakAnno_list.RData")


#
pdf("ATAC/matrix/peak_top10000_res/Distribution_of_peak_relative_to_TSS.pdf", width=5, height = 3)
plotDistToTSS(peakAnno_list, 
              title="Distribution of peak relative to TSS",
              ylab="Genomic Region (%) (5'->3')")
dev.off()
#
pdf("ATAC/matrix/peak_top10000_res/Feature_Distribution_of_peak.pdf", width=7, height = 3)
plotAnnoBar(peakAnno_list)
dev.off()


######################################
## chisq.test for Feature Distribution
######################################
library(gtools)
load("ATAC/matrix/peak_top10000_res/peakAnno_list.RData")
names(peakAnno_list)

comb <- combinations(length(names(peakAnno_list)), 2, names(peakAnno_list))
chisq_test_res_list <- list()
for(i in 1:dim(comb)[1]){
  peak_type1 <- comb[i, 1]
  peak_type2 <- comb[i, 2]
  comb_type <- paste(peak_type1, peak_type2, sep = "-")
  anno1 <- peakAnno_list[[peak_type1]]@anno$annotation
  anno2 <- peakAnno_list[[peak_type2]]@anno$annotation
  for(feature in c("Promoter", "Distal")){
    testInfo1 <- sapply(anno1, function(x){strsplit(x, " ")[[1]][1]})
    testInfo2 <- sapply(anno2, function(x){strsplit(x, " ")[[1]][1]})
    testInfo1[which(testInfo1!=feature)] <- "Other"
    testInfo2[which(testInfo2!=feature)] <- "Other"
    x <- matrix(c(as.vector(table(testInfo1)), as.vector(table(testInfo2))), ncol = 2)
    pvalue <- chisq.test(x)$p.value # 2.417162e-85(promoter) 7.399233e-64(Distal intergenic)
    chisq_test_res_list[[comb_type]][[feature]] <- pvalue
  }
}

t(data.frame(chisq_test_res_list))





#################################### function
library(gProfileR)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
library(ggplot2)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix/peak_top10000_res/peakAnno_list.RData")
load("ATAC/matrix/peak_top10000_res/peaks_LRT_res.RData")

peakAnno <- data.frame(peakAnno_list[['dif_peak']])
rownames(peakAnno) <- rownames(peaks_LRT_res)
save(peakAnno, file="ATAC/matrix/peak_top10000_res/dif_peakAnno.RData")
# 
# gene_peak_2.5kb <- peakAnno[abs(peakAnno$distanceToTSS)<2500,]
# gprofiler_res <- gprofiler(query=gene_peak_2.5kb$geneId,organism = "mmulatta")
# View(gprofiler_res)
# 
# 
# order <- gprofiler_res$term.name[order(gprofiler_res$p.value, decreasing=T)]
# gprofiler_res$term.name <- factor(gprofiler_res$term.name, levels = order)
# 
# ggplot(gprofiler_res, aes(x=term.name, y=-log10(p.value))) +
#   geom_point(aes(size=overlap.size/query.size, colour=-log10(p.value))) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         plot.title = element_text(hjust = 0.7, face="bold"),
#         axis.title = element_text(face="bold"),
#         axis.text.x  = element_text(face="bold", angle = 0, hjust=0.5, vjust=0.5, 
#                                     colour = "black", size = 10),
#         axis.text.y  = element_text(face="bold", colour = "black",size=10),
#         panel.border = element_rect(size = 1.5),
#         strip.background = element_rect(fill="white", colour = "black", size=1),legend.title = element_text(face="bold")) +
#   coord_flip() +
#   scale_color_gradientn(colours = brewer.pal(9,"Reds")[3:9])  +
#   ylab("-log10") +
#   xlab("") +
#   labs(size='GeneRatio') 
# 
# ggsave("different_peak_2.5kb_gene_function.pdf", width = 8, height = 5)

                                                                                      


##############################################################################
#                       dif-peak gene intersect dif-gene                     #
##############################################################################
###
# library(venn)
# setwd("/mnt/data1/Rhesus_brain/ATAC/matrix")
# load("peak_top10000_res/gene_anova_res.RData")
# 
# intersect(gene_peak_2.5kb$geneId, rownames(gene_anova_res))
# ####
# pdf("dif_peak_gene_2.5kb_intersect_dif_gene_venn.pdf")
# venn(list(gene_peak_2.5kb$geneId, rownames(gene_anova_res)),
#      snames = c("ATAC-seq", "RNA-seq"),
#      col=c("#3F60AC", "#8B2052"),
#      zcolor = c("#3F60AC", "#8B2052"),
#      cexil = 1.5,	# Character expansion for the intersection labels
#      cexsn = 1.5 #Character expansion for the set names,
# )
# dev.off()
# 
# pdf("all_dif_peak_gene_intersect_dif_gene_venn.pdf")
# venn(list(peakAnno$geneId, rownames(gene_anova_res)),
#      snames = c("ATAC-seq", "RNA-seq"),
#      col=c("#3F60AC", "#8B2052"),
#      zcolor = c("#3F60AC", "#8B2052"),
#      cexil = 1.5,	# Character expansion for the intersection labels
#      cexsn = 1.5 #Character expansion for the set names,
# )
# dev.off()
# 
# ##
# library(gProfileR)
# library(ggplot2)
# 
# genes <- intersect(gene_peak_2.5kb$geneId, rownames(gene_anova_res))
# gprofiler_res <- gprofiler(query=genes, organism = "mmulatta")
# View(gprofiler_res)



##############################################################################
#                   不同的distance, peak与基因(nearest)的cor                 #
##############################################################################
# library(hexbin)
# setwd("/mnt/data2/Rhesus_brain")
# load("ATAC/matrix/peak_top10000_res/dif_peakAnno.RData")
# load("peak_top10000_res/Rhesus_DESeq2_object.RData")
# load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
# 
# 
# vsd <- varianceStabilizingTransformation(Rhesus_PeakObject, blind=FALSE)
# peak_vsd_data <- assay(vsd)
# vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
# genes <- gtf_ensembl_gene[-which(gtf_ensembl_gene$seqnames %in% c("X", "Y")),]$gene_id
# length(genes)
# gene_vsd_data <- assay(vsd)
# gene_vsd_data <- gene_vsd_data[genes,colnames(peak_vsd_data)]
# 
# difPeakAnno <- peakAnno[,c("geneId", "distanceToTSS")]
# difPeakAnno <- data.frame(peak=rownames(difPeakAnno),difPeakAnno, stringsAsFactors = F)
# str(difPeakAnno)
# head(difPeakAnno)
# 
# # correlation <- apply(difPeakAnno, 1, function(x){cor(peak_vsd_data[x[1],], gene_vsd_data[x[2],])})
# correlation <- apply(difPeakAnno, 1, function(x){cor(rank(peak_vsd_data[x[1],]), rank(gene_vsd_data[x[2],]))})
# pp <- data.frame(difPeakAnno, correlation)
# head(pp)
# # test <- pp[pp$geneId %in% rownames(gene_anova_res),]
# # pp <- test
# 
# ggplot(pp, aes(x=distanceToTSS, y=correlation)) +
#   # geom_point() +
#   stat_binhex() +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         plot.title = element_text(hjust = 0.5, face="bold"),
#         axis.title = element_text(face="bold", size = 16),
#         axis.text.x  = element_text(face="bold", angle = 0, hjust=0.5, vjust=0.5, 
#                                     colour = "black", size = 15),
#         axis.text.y  = element_text(face="bold", colour = "black",size=15),
#         panel.border = element_rect(size = 1.5),
#         strip.background = element_rect(fill="white", colour = "black", size=1),
#         legend.title = element_text(face="bold")) +
#   scale_fill_gradient(low="lightblue", high="red") +
#   # scale_x_continuous(limits = c(-10000, 10000)) +
#   xlab("Genomic Region (5'->3')") +
#   ggtitle("the correlation of peak intensity and gene expression")
#   
#   
# ggsave("the_correlation_of_peak_intensity_and_gene_expression.pdf",
#        height = 5, width = 5.8)




##############################################################################
#                             peak与基因的cor                                #
##############################################################################
library(reshape2)
library(DESeq2)
library(WGCNA)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix/peak_top10000_res/dif_peakAnno.RData")
load("ATAC/matrix/peak_top10000_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("ATAC/matrix/peak_top10000_res/peaks_LRT_res.RData")
load("ATAC/matrix/peak_top10000_res/gene_LRT_res.RData")



##
peak_norCounts <- counts(Rhesus_PeakObject, normalized=TRUE)
gene_norCounts <- counts(Rhesus_GeneObject, normalized=TRUE)
genes <- gtf_ensembl_gene[-which(gtf_ensembl_gene$seqnames %in% c("X", "Y")),]$gene_id
length(genes) # 29361
gene_norCounts <- gene_norCounts[rownames(gene_norCounts) %in% genes,colnames(peak_norCounts)]
##
vsd <- varianceStabilizingTransformation(Rhesus_PeakObject, blind=FALSE)
peak_vsd_data <- assay(vsd)
vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
gene_vsd_data <- assay(vsd)
genes <- gtf_ensembl_gene[-which(gtf_ensembl_gene$seqnames %in% c("X", "Y")),]$gene_id
length(genes) # 29361
gene_vsd_data <- gene_vsd_data[rownames(gene_vsd_data) %in% genes, colnames(peak_vsd_data)]
##
# peak_df <- peak_vsd_data
# gene_df <- gene_vsd_data
peak_df <- peak_norCounts
gene_df <- gene_norCounts
cor(peak_df["6:179675817-179678513",], gene_df["ENSMMUG00000019858",])
# 0.6232772
peakCluster["6:179675817-179678513",] #3

########
all_peaks <- rownames(peak_df)
all_peaks <- data.frame(chrom=sapply(all_peaks, function(x){strsplit(x, ":")[[1]][1]}),
                        start=sapply(all_peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1]}),
                        end=sapply(all_peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2]}),
                        stringsAsFactors = F)
rownames(all_peaks) <- rownames(peak_df)


##### peak split chrom
peak_chrom_split_list <- list()
for(chrom in unique(all_peaks$chrom)){
  peak_chrom_split_list[[chrom]] <- rownames(all_peaks)[all_peaks$chrom == chrom]
}
##### gene split chrom
gene_chrom_split_list <- list()
for(chrom in unique(all_peaks$chrom)){
  gene_chrom_split_list[[chrom]] <- rownames(gene_df)[as.character(gtf_ensembl_gene[rownames(gene_df),"seqnames"]) == chrom]
}


### cor with dif chrom
peak_gene_cor_list <- list()
peak_gene_cor0.7_list <- list()
peak_gene_FDR_list <- list()
for(chrom in names(peak_chrom_split_list)){
  peak <- peak_chrom_split_list[[chrom]]
  gene <- gene_chrom_split_list[[chrom]]
  cor_value <- cor(t(peak_df[peak,]), t(gene_df[gene,]))
  corP = corPvalueStudent(cor_value, ncol(peak_df))
  dim(cor_value)
  peak_gene_cor <- melt(cor_value)
  peak_gene_p <- melt(corP)
  colnames(peak_gene_cor) <- c("peak", "gene", "correlation")
  colnames(peak_gene_p) <- c("peak", "gene", "pvalue")
  peak_gene_p$FDR <- p.adjust(peak_gene_p$pvalue, method = "fdr", n=length(peak_gene_p$pvalue))
  peak_gene_cor <- data.frame(peak_gene_cor, peak_gene_p[, c("pvalue", "FDR")])
  head(peak_gene_cor)
  peak_gene_cor <- peak_gene_cor[complete.cases(peak_gene_cor), ]
  peak_gene_cor <- data.frame(peak_gene_cor, all_peaks[as.character(peak_gene_cor$peak),])
  peak_gene_cor <- data.frame(peak_gene_cor,
                              gtf_ensembl_gene[as.character(peak_gene_cor$gene), c(1,2,3,4)])
  colnames(peak_gene_cor)[6:11] <- c("peak_chr", "peak_start", "peak_end", "gene_chr", "gene_start", "gene_end")
  peak_gene_cor$peak_start <- as.numeric(peak_gene_cor$peak_start)
  peak_gene_cor$peak_end <- as.numeric(peak_gene_cor$peak_end)
  peak_gene_cor$strand <- as.character(peak_gene_cor$strand)
  peak_gene_cor_list[[chrom]] <- peak_gene_cor
  peak_gene_cor0.7 <- peak_gene_cor[which(abs(peak_gene_cor$correlation)>=0.70),]
  dim(peak_gene_cor0.7)
  peak_gene_cor0.7_list[[chrom]] <- peak_gene_cor0.7
  peak_gene_FDR <- peak_gene_cor[which(peak_gene_cor$FDR<0.1),]
  peak_gene_FDR_list[[chrom]] <- peak_gene_FDR
}

lapply(peak_gene_cor0.7_list, function(x){dim(x)[1]})

save(peak_gene_cor_list, peak_gene_cor0.7_list, peak_gene_FDR_list,
     file="ATAC/matrix/peak_top10000_res/peak_gene_cor_list.RData")


###
load("ATAC/matrix/peak_top10000_res/peak_gene_cor_list.RData")

for(chorm in names(peak_gene_cor0.7_list)){
  peak_gene_cor0.7 <- peak_gene_cor0.7_list[[chorm]]
  distance <- c()
  if(dim(peak_gene_cor0.7)[1]==0) next
  for(i in 1:dim(peak_gene_cor0.7)[1]){
    x <- peak_gene_cor0.7[i,]
    if(x[12]=="+"){d <- unlist(c(x[7]-x[10],x[8]-x[10])); d <- unique(d[abs(d) == min(abs(d))])[1]}
    if(x[12]=="-"){d <- unlist(c(x[11]-x[7],x[11]-x[8])); d <- unique(d[abs(d) == min(abs(d))])[1]}
    if(length(d)>1){print(d)}
    distance <- c(distance, d)
  }
  peak_gene_cor0.7 <- data.frame(peak_gene_cor0.7, distance=distance)
  peak_gene_cor0.7_list[[chorm]] <- peak_gene_cor0.7
}


peak_gene_cor0.7 <- Reduce(rbind, peak_gene_cor0.7_list)
dim(peak_gene_cor0.7) # 10071    13
save(peak_gene_cor0.7, 
     file="ATAC/matrix/peak_top10000_res/peak_gene_cor0.7.RData")


#########################################
library(ggplot2)
library(dplyr)
library(ChIPseeker)
library(GenomicFeatures)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix/peak_top10000_res/peak_gene_cor0.7.RData")
load("ATAC/matrix/peak_top10000_res/all_peaks.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


peak_gene_cor_list[["6"]]$gene_name <- gtf_ensembl_gene[as.character(peak_gene_cor_list[["6"]]$gene), "gene_name"]
test <- peak_gene_cor_list[["6"]][peak_gene_cor_list[["6"]]$gene_name %in% "BTNL3",]
test[abs(test$correlation)>0.7,]


peak_gene_cor0.7$gene_name <- gtf_ensembl_gene[as.character(peak_gene_cor0.7$gene), "gene_name"]
peak_gene_cor0.7[peak_gene_cor0.7$gene_name %in% "BTNL3",]
peak_gene_cor0.7[peak_gene_cor0.7$gene_name %in% "IGKC",]

# length(unique(peak_gene_cor0.7$peak)) # 5148
# peak_gene_cor0.7_highest <- list()
# for(i in as.character(unique(peak_gene_cor0.7$peak))){
#   sub <- peak_gene_cor0.7[peak_gene_cor0.7$peak %in% i,]
#   if(dim(sub)[1]>1){
#     order <- abs(sub$correlation) == max(abs(sub$correlation))
#     sub <- sub[order,]
#   }
#   peak_gene_cor0.7_highest[[i]] <- sub
# }
# 
# pp <- Reduce(rbind, peak_gene_cor0.7_highest)
# peak_gene_cor0.7 <- pp



##############
### random function
##############
random_function <- function(peak_gene_pair){
  chr_num <- table(peak_gene_pair$peak_chr)
  chr_res_list <- list()
  for(chr in names(chr_num)){
    n <- chr_num[chr]
    chr_peaks <- all_peaks[all_peaks$chrom==chr,]
    # chr_peaks <- chr_peaks[sample(rownames(chr_peaks), n,  replace = T),]
    chr_peaks <- chr_peaks[sample(rownames(chr_peaks), n),]
    colnames(chr_peaks) <- c("peak_chr", "peak_start", "peak_end")
    chr_gtf <- gtf_ensembl_gene[gtf_ensembl_gene$seqnames==chr,c(1:5)]
    chr_gtf <- chr_gtf[sample(rownames(chr_gtf), n, replace = T), ]
    colnames(chr_gtf) <- c("gene_chr", "gene_start", "gene_end", "strand", "gene")
    chr_res <- cbind(chr_peaks, chr_gtf)
    chr_res$strand <- as.character(chr_res$strand)
    chr_res$peak_start <- as.numeric(chr_res$peak_start)
    chr_res$peak_end <- as.numeric(chr_res$peak_end)
    distance <- c()
    for(i in 1:dim(chr_res)[1]){
      x <- chr_res[i,]
      if(x[7]=="+"){d <- unlist(c(x[2]-x[5],x[3]-x[5])); d <- unique(d[abs(d) == min(abs(d))])[1]}
      if(x[7]=="-"){d <- unlist(c(x[5]-x[2],x[5]-x[3])); d <- unique(d[abs(d) == min(abs(d))])[1]}
      # if(length(d)>1){print(d)}
      distance <- c(distance, d)
    }
    chr_res <- data.frame(chr_res, distance=distance)
    chr_res_list[[chr]] <- chr_res
  }
  distance_res <- Reduce(rbind, chr_res_list)
  random_peaks_distance <- data.frame(distance=distance_res$distance, peak_type="random_peak")
  return(random_peaks_distance)
}




####################################  peak_gene_cor0.7
random_peaks_distance <- random_function(peak_gene_cor0.7)
peak_gene_cor0.7_distance <- data.frame(distance=peak_gene_cor0.7$distance,
                                        peak_type="peak_gene_cor0.7")
pp <- rbind(peak_gene_cor0.7_distance, random_peaks_distance)
pp <- data.frame(peak_gene_cor0.7_distance=peak_gene_cor0.7_distance$distance,
                 random_peaks_distance=random_peaks_distance$distance)
head(pp)
ggplot(pp, aes(x=peak_gene_cor0.7_distance)) +
  # geom_histogram(aes(x=random_peaks_distance), binwidth=10000000,boundary = 0.5, fill="#AAAAAA", color="#AAAAAA", alpha=0.5) +
  # geom_histogram(binwidth=10000000,boundary = 0.5,fill="#6F8DC8", color="#6F8DC8", alpha=0.5) +
  # geom_histogram(aes(x=random_peaks_distance), binwidth=10000,boundary = 0.5, fill="#AAAAAA", color="#AAAAAA", alpha=0.5) +
  # geom_histogram(binwidth=10000,boundary = 0.5,fill="#6F8DC8", color="#6F8DC8", alpha=0.5) +
  # geom_histogram(binwidth=10000000, color="black", fill="#6F8DC8",boundary = 0.5) +
  geom_density(stat = "density", size=0.5, color="#247BA0", fill="#247BA0") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1, color = "black"),
        axis.text = element_text(size = 15, colour = "black"),
        plot.title = element_text(hjust = 0.5,size = 15),
        axis.title.x =element_text(size=15), 
        axis.title.y=element_text(size=15))+
  ggtitle("Peak-gene distance with cor0.7 ") +
  # scale_x_continuous(limits = c(-1000000, 1000000), labels = c("-1Mb", "-500kb", "TSS", "500kb", "1Mb")) +
  # scale_x_continuous(labels = c("-200Mb", "-100Mb", "TSS", "100Mb", "200Mb")) +
  # scale_y_log10() +
  xlab("Genomic Region (5'->3')")
  
ggsave("ATAC/matrix/peak_top10000_res/Peak-gene_distance_with_cor0.7.pdf", width = 5.5, height = 5)



####################
#### distance < 50Mb
####################
ggplot(pp, aes(x=peak_gene_cor0.7_distance)) +
  geom_histogram(aes(x=random_peaks_distance), binwidth=1000000,boundary = 0.5, fill="#AAAAAA", color="#AAAAAA", alpha=0.5) +
  geom_histogram(binwidth=1000000,boundary = 0.5,fill="#6F8DC8", color="#6F8DC8", alpha=0.5) +
  # geom_histogram(binwidth=10000000, color="black", fill="#6F8DC8",boundary = 0.5) +
  # geom_line(stat = "density", size=0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1, color = "black"),
        axis.text = element_text(size = 15, colour = "black"),
        plot.title = element_text(hjust = 0.5,size = 15),
        axis.title.x =element_text(size=15), 
        axis.title.y=element_text(size=15))+
  ggtitle("Peak-gene distance with correlation 0.7 ") +
  scale_x_continuous(limits = c(-50000000, 50000000), labels = c("-50Mb", "-25Mb", "TSS", "25Mb", "50Mb")) +
  xlab("Genomic Region (5'->3')")

ggsave("ATAC/matrix/peak_top10000_res/Peak-gene_distance_with_correlation_0.7_distance_50Mb_and_random_peak.pdf", width = 5.5, height = 5)


# 
# ####################
# #### distance < 10Mb
# ####################
# ggplot(pp, aes(x=peak_gene_cor0.7_distance)) +
#   geom_histogram(aes(x=random_peaks_distance), binwidth=500000,boundary = 0.5, fill="#AAAAAA", color="#AAAAAA", alpha=0.5) +
#   geom_histogram(binwidth=500000,boundary = 0.5,fill="#6F8DC8", color="#6F8DC8", alpha=0.5) +
#   # geom_histogram(binwidth=10000000, color="black", fill="#6F8DC8",boundary = 0.5) +
#   # geom_line(stat = "density", size=0.5) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         panel.border = element_rect(size=1, color = "black"),
#         axis.text = element_text(size = 15, colour = "black"),
#         plot.title = element_text(hjust = 0.5,size = 15),
#         axis.title.x =element_text(size=15), 
#         axis.title.y=element_text(size=15))+
#   ggtitle("Peak-gene distance with correlation 0.7 ") +
#   scale_x_continuous(limits = c(-10000000, 10000000), labels = c("-10Mb", "-5Mb", "TSS", "5Mb", "10Mb")) +
#   xlab("Genomic Region (5'->3')")
# 
# 
# ggsave("ATAC/matrix/peak_top10000_res/Peak-gene_distance_with_correlation_0.7_distance_10Mb.pdf", width = 5.5, height = 5)
# 
# 
# 
# ###################
# #### distance < 1Mb
# ###################
# ggplot(pp, aes(x=peak_gene_cor0.7_distance)) +
#   geom_histogram(aes(x=random_peaks_distance), binwidth=50000,boundary = 0.5, fill="#AAAAAA", color="#AAAAAA", alpha=0.5) +
#   geom_histogram(binwidth=50000,boundary = 0.5,fill="#6F8DC8", color="#6F8DC8", alpha=0.5) +
#   # geom_histogram(binwidth=10000000, color="black", fill="#6F8DC8",boundary = 0.5) +
#   # geom_line(stat = "density", size=0.5) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         panel.border = element_rect(size=1, color = "black"),
#         axis.text = element_text(size = 15, colour = "black"),
#         plot.title = element_text(hjust = 0.5,size = 15),
#         axis.title.x =element_text(size=15), 
#         axis.title.y=element_text(size=15))+
#   ggtitle("Peak-gene distance with correlation 0.7 ") +
#   scale_x_continuous(limits = c(-1000000, 1000000), labels = c("-1Mb", "-500kb", "TSS", "500kb", "1Mb")) +
#   xlab("Genomic Region (5'->3')")
# 
# 
# ggsave("ATAC/matrix/peak_top10000_res/Peak-gene_distance_with_correlation_0.7_distance_1Mb_and_random_peak.pdf", width = 5.5, height = 5)
# 
# 
# 
# 
# ###################
# #### distance < 0.5Mb
# ###################
# ggplot(pp, aes(x=peak_gene_cor0.7_distance)) +
#   geom_histogram(aes(x=random_peaks_distance), binwidth=5000,boundary = 0.5, fill="#AAAAAA", color="#AAAAAA", alpha=0.5) +
#   geom_histogram(binwidth=5000,boundary = 0.5,fill="#6F8DC8", color="#6F8DC8", alpha=0.5) +
#   # geom_histogram(binwidth=10000000, color="black", fill="#6F8DC8",boundary = 0.5) +
#   # geom_line(stat = "density", size=0.5) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         panel.border = element_rect(size=1, color = "black"),
#         axis.text = element_text(size = 15, colour = "black"),
#         plot.title = element_text(hjust = 0.5,size = 15),
#         axis.title.x =element_text(size=15), 
#         axis.title.y=element_text(size=15))+
#   ggtitle("Peak-gene distance with correlation 0.7 ") +
#   scale_x_continuous(limits = c(-500000, 500000), labels = c("-0.5Mb", "-250kb", "TSS", "250kb", "0.5Mb")) +
#   xlab("Genomic Region (5'->3')")
# 
# 
# ggsave("ATAC/matrix/peak_top10000_res/Peak-gene_distance_with_correlation_0.7_distance_0.5Mb_and_random_peak.pdf", width = 5.5, height = 5)
# 
# 
# 
# #########
# peak_bar <- table(table(as.character(peak_gene_cor0.7$peak)))
# peak_bar <- data.frame(peak_bar)
# colnames(peak_bar) <- c("time", "count")
# 
# ggplot(peak_bar, aes(x=time, y=count)) +
#   geom_bar(stat="identity", color="black", fill="white") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         panel.border = element_rect(size=1, color = "black"),
#         axis.text = element_text(size = 15, colour = "black"),
#         plot.title = element_text(hjust = 0.5, face="bold",size = 15),
#         axis.title = element_text(face="bold"),
#         axis.title.x =element_text(size=15), 
#         axis.title.y=element_text(size=15),
#         legend.title = element_text(face="bold"))+
#   xlab("the times of the same peak occurs")+
#   ylab("peak counts of the same times") +
#   ggtitle("")
# 
# ggsave("peak_counts_of_the_same_times.pdf", width = 6, height = 5)
# 
# ## help
# ggplot(peak_bar[6:dim(peak_bar)[1],], aes(x=time, y=count)) +
#   geom_bar(stat="identity", color="black", fill="white") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         panel.border = element_rect(size=1, color = "black"),
#         axis.text = element_text(size = 15, colour = "black"),
#         plot.title = element_text(hjust = 0.5, face="bold",size = 15),
#         axis.title = element_text(face="bold"),
#         axis.title.x =element_text(size=15), 
#         axis.title.y=element_text(size=15),
#         legend.title = element_text(face="bold"))+
#   xlab("")+
#   ylab("") +
#   ggtitle("")
# 
# ggsave("peak_counts_of_the_same_times_help.pdf", width = 5, height = 4)
# 
# 
# 
# 
# #########
# gene_bar <- table(table(as.character(peak_gene_cor0.7$gene)))
# gene_bar <- data.frame(gene_bar)
# colnames(gene_bar) <- c("time", "count")
# head(gene_bar)
# 
# ggplot(gene_bar, aes(x=time, y=count)) +
#   geom_bar(stat="identity", color="black", fill="white") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         panel.border = element_rect(size=1, color = "black"),
#         axis.text = element_text(size = 15, colour = "black"),
#         plot.title = element_text(hjust = 0.5, face="bold",size = 15),
#         axis.title = element_text(face="bold"),
#         axis.title.x =element_text(size=15), 
#         axis.title.y=element_text(size=15),
#         legend.title = element_text(face="bold"))+
#   xlab("the times of the same gene occurs")+
#   ylab("gene counts of the same times") +
#   ggtitle("")
# 
# ggsave("gene_counts_of_the_same_times.pdf", width = 12, height = 5)
# 
# ## help
# ggplot(gene_bar[7:dim(gene_bar)[1],], aes(x=time, y=count)) +
#   geom_bar(stat="identity", color="black", fill="white") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         panel.border = element_rect(size=1, color = "black"),
#         axis.text = element_text(size = 15, colour = "black"),
#         plot.title = element_text(hjust = 0.5, face="bold",size = 15),
#         axis.title = element_text(face="bold"),
#         axis.title.x =element_text(size=15), 
#         axis.title.y=element_text(size=15),
#         legend.title = element_text(face="bold"))+
#   xlab("")+
#   ylab("") +
#   ggtitle("")
# 
# ggsave("gene_counts_of_the_same_times_help.pdf", width = 12, height = 4)
# 
# 

############################################################
#                         biotype                          #
############################################################
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix/peak_top10000_res/peak_gene_cor0.7.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("ATAC/matrix/peak_top10000_res/peaks_LRT_res.RData")


colors <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(13), 
                   unique(gtf_ensembl_gene$gene_biotype)[2:14])

#######################
###### peak_gene_FDR
#######################

## protein_coding and noncoding
sub <- peak_gene_cor0.7[peak_gene_cor0.7$peak %in% rownames(peaks_LRT_res), ]
gene_biotype <- gtf_ensembl_gene[as.character(sub$gene),"gene_biotype"]
gene_biotype[which(gene_biotype!="protein_coding")] <- "noncoding"
pp <- data.frame(table(gene_biotype))

pp = pp[order(pp$Freq, decreasing = TRUE),] 
myLabel = as.vector(pp$gene_biotype)   ## 转成向量，否则图例的标签可能与实际顺序不一???
myLabel = paste(myLabel, "(", round(pp$Freq / sum(pp$Freq) * 100, 2), "%)", sep = "")   ## ??? round() 对结果保留两位小???

ggplot(pp, aes(x="", y=Freq, fill=gene_biotype)) +
  geom_bar(stat = "identity", width = 1, alpha=0.6) +
  coord_polar(theta = "y", start = 1.8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.7),
        # axis.text.x  = element_text(face="bold", angle = 0, hjust=0.5, vjust=0.5, 
        #                             colour = "black", size = 10),
        axis.text.y  = element_text(colour = "black",size=10),
        panel.border = element_rect(size = 1.5),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size=10, colour = "black", angle=90),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  scale_fill_manual(values = c(protein_coding="#3F60AC", noncoding="#8B2052"), breaks =pp$gene_biotype, labels = myLabel) +
  labs(x = "", y = "", title = "") +
  ggtitle("peak-gene cor0.7 : gene_biotype ")
  

ggsave("ATAC/matrix/peak_top10000_res/difPeak_gene_cor0.7_gene_biotype_percent.pdf",
       width = 6, height = 6)


##### noncoding detail
sub <- peak_gene_cor0.7[peak_gene_cor0.7$peak %in% rownames(peaks_LRT_res), ]
gene_biotype <- gtf_ensembl_gene[as.character(sub$gene),"gene_biotype"]
gene_biotype <- gene_biotype[which(gene_biotype!="protein_coding")]
pp <- data.frame(table(gene_biotype))

pp = pp[order(pp$Freq, decreasing = TRUE),] 
myLabel = as.vector(pp$gene_biotype)   ## 转成向量，否则图例的标签可能与实际顺序不一???
myLabel = paste(myLabel, "(", round(pp$Freq / sum(pp$Freq) * 100, 2), "%)", sep = "")   ## ??? round() 对结果保留两位小???

ggplot(pp, aes(x="", y=Freq, fill=gene_biotype)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.7, face="bold"),
        axis.title = element_text(face="bold"),
        # axis.text.x  = element_text(face="bold", angle = 0, hjust=0.5, vjust=0.5, 
        #                             colour = "black", size = 10),
        axis.text.y  = element_text(face="bold", colour = "black",size=10),
        panel.border = element_rect(size = 1.5),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size=10, colour = "black", face="bold", angle=90),
        strip.background = element_rect(fill="white", colour = "white", size=1),legend.title = element_text(face="bold")) +
  scale_fill_manual(values = colors, breaks =pp$gene_biotype, labels = myLabel) +
  labs(x = "", y = "", title = "") +
  ggtitle("peak-gene cor0.7: gene_biotype ")


ggsave("ATAC/matrix/peak_top10000_res/difPeak_gene_cor0.7_gene_biotype_detail_percent.pdf", 
       width = 6, height = 6)


##################
######   all genes
##################

## protein_coding and noncoding
gene_biotype <- gtf_ensembl_gene[-which(gtf_ensembl_gene$seqnames %in% c("M", "X", "Y")),"gene_biotype"]
gene_biotype[which(gene_biotype!="protein_coding")] <- "noncoding"
pp <- data.frame(table(gene_biotype))
pp
pp = pp[order(pp$Freq, decreasing = TRUE),] 
myLabel = as.vector(pp$gene_biotype)   ## 转成向量，否则图例的标签可能与实际顺序不一???
myLabel = paste(myLabel, "(", round(pp$Freq / sum(pp$Freq) * 100, 2), "%)", sep = "")   ## ??? round() 对结果保留两位小???

ggplot(pp, aes(x="", y=Freq, fill=gene_biotype)) +
  geom_bar(stat = "identity", width = 1, alpha=0.6) +
  coord_polar(theta = "y", start = 2.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.7, face="bold"),
        axis.title = element_text(face="bold"),
        # axis.text.x  = element_text(face="bold", angle = 0, hjust=0.5, vjust=0.5, 
        #                             colour = "black", size = 10),
        axis.text.y  = element_text(face="bold", colour = "black",size=10),
        panel.border = element_rect(size = 1.5),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size=10, colour = "black", face="bold", angle=90),
        strip.background = element_rect(fill="white", colour = "white", size=1),legend.title = element_text(face="bold")) +
  scale_fill_manual(values = c(protein_coding="#3F60AC", noncoding="#8B2052"), breaks =pp$gene_biotype, labels = myLabel) +
  labs(x = "", y = "", title = "") +
  ggtitle("all gene gene_biotype ")


ggsave("ATAC/matrix/peak_top10000_res/all_gene_biotype_percent.pdf", width = 6, height = 6)

##  noncoding detail
gene_biotype <- gtf_ensembl_gene[-which(gtf_ensembl_gene$seqnames %in% c("M", "X", "Y")),"gene_biotype"]
gene_biotype <- gene_biotype[which(gene_biotype!="protein_coding")]
pp <- data.frame(table(gene_biotype))
pp
pp = pp[order(pp$Freq, decreasing = TRUE),] 
myLabel = as.vector(pp$gene_biotype)   ## 转成向量，否则图例的标签可能与实际顺序不一???
myLabel = paste(myLabel, "(", round(pp$Freq / sum(pp$Freq) * 100, 2), "%)", sep = "")   ## ??? round() 对结果保留两位小???

ggplot(pp, aes(x="", y=Freq, fill=gene_biotype)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.7, face="bold"),
        axis.title = element_text(face="bold"),
        # axis.text.x  = element_text(face="bold", angle = 0, hjust=0.5, vjust=0.5, 
        #                             colour = "black", size = 10),
        axis.text.y  = element_text(face="bold", colour = "black",size=10),
        panel.border = element_rect(size = 1.5),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size=10, colour = "black", face="bold", angle=90),
        strip.background = element_rect(fill="white", colour = "white", size=1),legend.title = element_text(face="bold")) +
  scale_fill_manual(values = colors, breaks =pp$gene_biotype, labels = myLabel) +
  labs(x = "", y = "", title = "") +
  ggtitle("all gene gene_biotype ")


ggsave("ATAC/matrix/peak_top10000_res/all_gene_biotype_detail_percent.pdf", width = 6, height = 6)



## chisq.test
biotype_all <- gtf_ensembl_gene[-which(gtf_ensembl_gene$seqnames %in% c("M", "X", "Y")),"gene_biotype"]
biotype_all[which(biotype_all!="protein_coding")] <- "noncoding"
biotype_cor <- gtf_ensembl_gene[as.character(peak_gene_FDR$gene),"gene_biotype"]
biotype_cor[which(biotype_cor!="protein_coding")] <- "noncoding"

x <- matrix(c(as.vector(table(biotype_all)), as.vector(table(biotype_cor))), ncol = 2)
chisq.test(x)$p.value



##############################################################################
#                       dif peak and gene (peak-gene pair)                   #
##############################################################################
library(gProfileR)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
library(ggplot2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix/peak_top10000_res/peak_gene_cor0.7.RData")
load("ATAC/matrix/peak_top10000_res/ConsensusCluster_DifPeaks.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("ATAC/matrix/peak_top10000_res/dif_peakAnno.RData")

###
k=4
res <- DIfPeaks_Results[[k]]
geneCluster <- data.frame(geneCluster=factor(sort(res$consensusClass)))

##
sub_peak_gene_cor0.7 <- peak_gene_cor0.7[as.character(peak_gene_cor0.7$peak) %in% rownames(geneCluster), ]
distance <- head(sort(abs(sub_peak_gene_cor0.7$distance)))
sub_peak_gene_cor0.7[abs(sub_peak_gene_cor0.7$distance)==distance[1],]


##
peakCluster_function_list <- list()
cluster_gene_list <- list()
for(cluster in 1:k){
  peaks <- rownames(geneCluster)[geneCluster$geneCluster==cluster]
  sub_peak_gene_cor0.7 <- peak_gene_cor0.7[as.character(peak_gene_cor0.7$peak) %in% peaks, ]
  genes <- sub_peak_gene_cor0.7[sub_peak_gene_cor0.7$correlation>0, "gene"]
  # genes <- sub_peak_gene_cor0.7$gene
  genes <- unique(as.character(genes))
  length(genes)
  cluster_gene_list[[cluster]] <- genes
  gprofiler_res <- gprofiler(query=genes, organism = "mmulatta")
  gprofiler_res <- gprofiler_res[order(gprofiler_res[,3]),]
  gprofiler_res$gene_name <- gprofiler_res$intersection
  if(dim(gprofiler_res)[1]==0) next
  else{for(i in 1:dim(gprofiler_res)[1]){
    x <- gprofiler_res$gene_name[i]
    gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],"gene_name"], collapse=",")
  }
    filename=paste("ATAC/matrix/peak_top10000_res/LRT_peakCluster", cluster,"_gene_function_with_cor0.7.txt", sep="")
    write.table(gprofiler_res, filename,
                quote = F, sep="\t", row.names = F, col.names = T)
  }
  peakCluster_function_list[[cluster]] <- gprofiler_res
}


### cluster1
cluster <- 1
gprofiler_res <- peakCluster_function_list[[cluster]]
View(gprofiler_res)
gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),]
cluster1 <- data.frame(gprofiler_res, cluster="peakCluster1")

#### cluster2
# cluster <- 2
# gprofiler_res <- peakCluster_function_list[[cluster]]
# View(gprofiler_res)
# gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),]
# cluster2 <- data.frame(gprofiler_res, cluster="peakCluster2")

#### cluster3
cluster <- 3
gprofiler_res <- peakCluster_function_list[[cluster]]
gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),]
cluster3 <- data.frame(gprofiler_res, cluster="peakCluster3")


#### cluster4
cluster <- 4
gprofiler_res <- peakCluster_function_list[[cluster]]
gprofiler_res <- gprofiler_res[-which(gprofiler_res$domain %in% "hp"), ]
gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),]
cluster4 <- data.frame(gprofiler_res, cluster="peakCluster4")


##
gprofiler_res <- rbind(cluster1, cluster3, cluster4)
order <- gprofiler_res$term.name[order(gprofiler_res$p.value, decreasing=T)]
gprofiler_res$term.name <- factor(gprofiler_res$term.name, levels=order)
ggplot(gprofiler_res, aes(x=term.name, y=-log10(p.value))) +
  geom_bar(stat = "identity", fill="#BEA8D0", width = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.7),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 10),
        axis.text.y  = element_text(colour = "black",size=10),
        panel.border = element_rect(size = 1.5),
        strip.text = element_text(size=10, colour = "black", angle=90),
        strip.background = element_rect(fill="white", colour = "white", size=1),legend.title = element_text(face="bold")) +
  coord_flip() +
  scale_color_gradientn(colours = brewer.pal(9,"Reds")[3:9])  +
  facet_grid(cluster~., scales = "free_y",space = "free_y") +
  ylab("-log10(pvalue)") +
  xlab("") +
  labs(size='GeneRatio') 

# filename=paste("anova_peakCluster", cluster,"_function_top20.pdf", sep="")
filename=paste("ATAC/matrix/peak_top10000_res/LRT_peakCluster_cor_gene_function.pdf", sep="")
ggsave(filename, width = 5, height = 4)


### peak_gene_FDR 中所有差异peak对应的基因的功能
peaks <- rownames(geneCluster)[geneCluster$geneCluster==cluster]
genes <- peak_gene_cor0.7$gene[as.character(peak_gene_cor0.7$peak) %in% rownames(geneCluster)]
genes <- unique(as.character(genes))
length(genes) # 304
gprofiler_res <- gprofiler(query=genes, organism = "mmulatta")
View(gprofiler_res[, c("p.value", "term.name")])




####### cluster1(p) and cluster3(n)
peaks_1 <- rownames(geneCluster)[geneCluster$geneCluster==1]
genes_1 <- peak_gene_cor0.7$gene[(as.character(peak_gene_cor0.7$peak) %in% peaks_1) & peak_gene_cor0.7$correlation>0.7]
genes_1 <- unique(as.character(genes_1))

peaks_3 <- rownames(geneCluster)[geneCluster$geneCluster==3]
genes_3 <- peak_gene_cor0.7$gene[(as.character(peak_gene_cor0.7$peak) %in% peaks_3) & peak_gene_cor0.7$correlation< -0.7]
genes_3 <- unique(as.character(genes_3))

p_n <-intersect(genes_1, genes_3)
gtf_ensembl_gene[p_n, 6]


####### cluster1(n) and cluster3(p)
peaks_1 <- rownames(geneCluster)[geneCluster$geneCluster==1]
genes_1 <- peak_gene_cor0.7$gene[(as.character(peak_gene_cor0.7$peak) %in% peaks_1) & peak_gene_cor0.7$correlation< -0.7]
genes_1 <- unique(as.character(genes_1))

peaks_3 <- rownames(geneCluster)[geneCluster$geneCluster==3]
genes_3 <- peak_gene_cor0.7$gene[(as.character(peak_gene_cor0.7$peak) %in% peaks_3) & peak_gene_cor0.7$correlation> 0.7]
genes_3 <- unique(as.character(genes_3))

n_p <- intersect(genes_1, genes_3)
gtf_ensembl_gene[n_p, 6]





##############################################################################
#                      peak and gene module correlation                      #
##############################################################################
library(WGCNA)
library(reshape2)
library(DESeq2)
setwd("/mnt/xdlab1/home/looking/brain_project")
load("ATAC/matrix/peak_anova_res.RData")
load("matrix0420/gene_count_res/WGCNA_Rhesus_datExpr_exclude_MB_ex1.RData")
load("matrix0420/gene_count_res/WGCNA_Rhesus_net_modules_exclude_MB_ex1.RData")
load("ATAC/matrix/peak_top10000_res/Rhesus_DESeq2_object.RData")
load("matrix0420/gene_count_res/gtf_ensembl_gene.RData")

moduleLabels = net$colors
MEs0 = moduleEigengenes(datExpr, moduleLabels)$eigengenes # MEs = net$MEs
MEsWW = orderMEs(MEs0)
dim(MEsWW) # 407 56
rownames(MEsWW) <- rownames(datExpr)


vsd <- varianceStabilizingTransformation(Rhesus_PeakObject, blind=FALSE)
peak_vsd_data <- assay(vsd)

MEsWW_sub <- MEsWW[colnames(peak_vsd_data),]


cor <- cor(MEsWW_sub, t(peak_vsd_data))
dim(cor)
cor <- t(cor)
cor <- melt(cor)
head(cor)

cor_sub <- cor[abs(cor$value)>0.7,]
dif_peak <- rownames(peak_anova_res)
head(cor_sub)
sum(dif_peak %in% cor_sub$Var1)
sum(cor_sub$Var1 %in% dif_peak)


a <- hist(cor$value)















#######################
## peak length distribution


### all peaks 
vsd <- varianceStabilizingTransformation(Rhesus_PeakObject, blind=FALSE)
vsd_data <- assay(vsd)
all_peaks <- rownames(vsd_data)
all_peaks <- data.frame(chrom=sapply(all_peaks, function(x){strsplit(x, ":")[[1]][1]}),
                        start=sapply(all_peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1]}),
                        end=sapply(all_peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2]}),
                        stringsAsFactors = F)
all_peaks$length <- as.numeric(all_peaks$end)-as.numeric(all_peaks$start)
View(all_peaks)


ggplot(all_peaks, aes(x=length)) +
  geom_density() +
  scale_x_log10()





