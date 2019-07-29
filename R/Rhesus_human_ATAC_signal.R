

########################################################
#                      colData                         #
########################################################
setwd("/mnt/data2/Rhesus_brain/ATAC")
human_ATAC_colData1 <- read.table("liftover_res/GR_human_ATAC_SraRunTable.txt", header=T, sep="\t", stringsAsFactors = F)
human_ATAC_colData2 <- read.table("GSE96949/Supplemental_Table_S3.txt", header=T, sep="\t", stringsAsFactors = F)[, c(2,3)]
human_ATAC_colData3 <- read.table("liftover_res/genome_average_signal.txt", header = F, sep="\t", stringsAsFactors = F)
human_ATAC_colData3$Sample_Name <- sapply(human_ATAC_colData3$V1, function(x){strsplit(x, "_")[[1]][1]})
human_ATAC_colData3$ID <- paste("human", sapply(human_ATAC_colData3$V1, function(x){strsplit(x, "_")[[1]][2]}), sep="")
rownames(human_ATAC_colData3) <- human_ATAC_colData3$Sample_Name
human_ATAC_colData2 <- human_ATAC_colData2[order(human_ATAC_colData2$Brain.region..Abbreviation.),]
human_ATAC_colData2 <- human_ATAC_colData2[!duplicated(human_ATAC_colData2$Brain.region..Abbreviation.),]
rownames(human_ATAC_colData2) <- human_ATAC_colData2$Brain.region..Name.
human_ATAC_colData1$region <- human_ATAC_colData2[human_ATAC_colData1$source_name, 1]
rownames(human_ATAC_colData1) <- human_ATAC_colData1$Run
human_ATAC_colData <- human_ATAC_colData1[, c("gender", "region", "age_of_death", "cell_type", "source_name", "Sample_Name")]
human_ATAC_colData$ID <- human_ATAC_colData3[human_ATAC_colData$Sample_Name, "ID"]
human_ATAC_colData <- human_ATAC_colData[human_ATAC_colData$region %in% c("ACC", "DLPFC", "INS", "ITC", "OFC", "PMC", "PVC", "STC", "VLPFC"),]


save(human_ATAC_colData, file="GSE96949/matrix/human_ATAC_colData.RData")

human_ATAC_colData <- data.frame(SRA=rownames(human_ATAC_colData), human_ATAC_colData)
write.table(human_ATAC_colData, "GSE96949/matrix/human_ATAC_colData.txt",
            col.names = T, row.names = F, quote = F, sep="\t")


########################################################
#                        rowData                       #
########################################################
setwd("/mnt/data2/Rhesus_brain")

rowData <- read.table("ATAC/liftover_res/hg38_RheMac8/all_peaks_sorted_orthology_hg38.bed")
colnames(rowData) <- c("human_chr", "human_start", "human_end", "Rhesus_peak")
rowData$human_chr <- gsub("chr", "", rowData$human_chr)
rowData$Rhesus_peak <- gsub("chr", "", rowData$Rhesus_peak)
names <- paste(rowData$human_chr, ":", rowData$human_start, "-", rowData$human_end, sep="")
rownames(rowData) <- names

save(rowData, file="GSE96949/matrix/rowData.Rdata")



########################################################
#               load human ATAC peak                   #
########################################################
setwd("/mnt/data2/Rhesus_brain/ATAC")
load("GSE96949/matrix/human_ATAC_colData.RData")


human_peak_binary <- read.table("GSE96949/matrix/human_conserve_region_peak.matrix",
                                header=T, row.names = 1)
human_peak_binary <- human_peak_binary[,rownames(human_ATAC_colData)]


region_sharing <- apply(human_peak_binary, 1, function(x){sum(x>0)})
head(region_sharing)
pp <- data.frame(table(region_sharing))
pp$percent <- pp$Freq/dim(human_peak_binary)[1]
head(pp)

percent <- c()
percent <- setNames(pp$percent[1], "0")
for(i in seq(1, 88, 8)){
  j <- i + 7
  res <- pp
  res$region_sharing <- as.integer(res$region_sharing)-1
  res <- res[res$region_sharing>=i & res$region_sharing<=j,]
  name = paste(i, j, sep="-")
  percent_sum <- sum(res$percent)
  names(percent_sum) <- name
  percent <- c(percent, percent_sum)
}

percent <- data.frame(percent)
percent$percent <- percent$percent*100
percent$region_sharing <- rownames(percent)
percent$region_sharing <- factor(percent$region_sharing, levels = percent$region_sharing)

ggplot(percent, aes(x=region_sharing, y=percent)) +
  geom_bar(stat = "identity", fill="#247BA0") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.7, color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(color = "black")) 

ggsave("GSE96949/matrix/human_conserve_region_sharing_bar.pdf",
       width = 6, height = 5)


#####
# load("matrix/peak_top10000_res/peaks_LRT_res.RData")
# 
# 
# region_sharing_tmp <- region_sharing
# names(region_sharing_tmp) <- rowData[names(region_sharing), "Rhesus_peak"]
# 
# region_sharing_tmp <- data.frame(region_sharing_tmp)
# colnames(region_sharing_tmp) <- c("region_sharing")
# region_sharing_tmp$difPeak <- "NO"
# region_sharing_tmp$difPeak[rownames(region_sharing_tmp) %in% rownames(peaks_LRT_res)] <- "YES"
# 
# 
# pp <- table(region_sharing_tmp$region_sharing, region_sharing_tmp$difPeak)
# pp <- data.frame(pp)
# colnames(pp)[1:2] <- c("region_sharing", "difPeak")
# pp$percent <- pp$Freq/dim(human_peak_binary)[1]
# head(pp)
# 
# 
# a <- pp[pp$difPeak=="YES",]
# a <- pp[pp$difPeak=="NO",]
# percent <- c()
# percent <- setNames(a$percent[1], "0")
# for(i in seq(1, 115, 23)){
#   j <- i + 22
#   res <- a
#   res$region_sharing <- as.integer(res$region_sharing)-1
#   res <- res[res$region_sharing>=i & res$region_sharing<=j,]
#   name = paste(i, j, sep="-")
#   percent_sum <- sum(res$percent)
#   names(percent_sum) <- name
#   percent <- c(percent, percent_sum)
# }
# NO_percent <- percent
# YES_percent <- percent
# 
# percent <- rbind(data.frame(percent=NO_percent, difPeak="NO"),
#                  data.frame(percent=YES_percent, difPeak="YES"))
# percent$percent <- percent$percent*100
# percent$region_sharing <- rownames(percent)
# percent$region_sharing[7:12] <- percent$region_sharing[1:6]
# 
# 
# ggplot(percent, aes(x=region_sharing, y=percent, fill=difPeak)) +
#   geom_bar(stat = "identity") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.7, color = "black"),
#         axis.text.y = element_text(color = "black"),
#         text = element_text(color = "black")) 
# 
# 
# 
# ggsave("GSE96949/matrix/human_conserve_region_sharing_bar2.pdf",
#        width = 6.3, height = 5)
 


########################################################
#               load human ATAC count                  #
########################################################
library(DESeq2)
library(BiocParallel)
setwd("/mnt/data2/Rhesus_brain/ATAC")
load("GSE96949/matrix/human_ATAC_colData.RData")

###
human_peak_binary <- read.table("GSE96949/matrix/human_conserve_region_peak.matrix",
                                header=T, row.names = 1)
human_peak_binary <- human_peak_binary[,rownames(human_ATAC_colData)]
region_sharing <- apply(human_peak_binary, 1, function(x){sum(x>0)})
###
human_ATAC_count <- read.table("GSE96949/matrix/human_conserve_region_count.matrix",
                         header=T, row.names = 1)
human_ATAC_count <- human_ATAC_count[names(region_sharing)[region_sharing>0], rownames(human_ATAC_colData)]
dim(human_ATAC_count) # 7842   88



### Neuronal
count <- human_ATAC_count[, rownames(human_ATAC_colData)[human_ATAC_colData$cell_type=="Neuronal"]]
count <- count[apply(count, 1, function(x){sum(x>1)>3}),]
dim(count) # 7819   44
human_ATAC_Neuronal <- DESeqDataSetFromMatrix(count, 
                                            human_ATAC_colData[human_ATAC_colData$cell_type=="Neuronal",], 
                                            ~region+ID)
human_ATAC_Neuronal <- DESeq(human_ATAC_Neuronal, 
                             test="LRT",
                             reduced = ~region,
                           parallel = T, 
                           BPPARAM=MulticoreParam(40))
save(human_ATAC_Neuronal, 
     file="GSE96949/matrix/human_ATAC_Neuronal.RData")


### Non-Neuronal 
count <- human_ATAC_count[, rownames(human_ATAC_colData)[human_ATAC_colData$cell_type=="Non-Neuronal"]]
count <- count[apply(count, 1, function(x){sum(x>1)>3}),]
dim(count) # 7834   44
human_ATAC_Non_Neuronal <- DESeqDataSetFromMatrix(count, 
                                              human_ATAC_colData[human_ATAC_colData$cell_type=="Non-Neuronal",], 
                                              ~region+ID)
human_ATAC_Non_Neuronal <- DESeq(human_ATAC_Non_Neuronal, 
                                 test="LRT",
                                 reduced = ~region,
                             parallel = T, BPPARAM=MulticoreParam(40))
save(human_ATAC_Non_Neuronal, 
     file="GSE96949/matrix/human_ATAC_Non_Neuronal.RData")




########################################################
#                     Neuronal                         #
########################################################
library(DESeq2)
library(VennDiagram)
setwd("/mnt/data2/Rhesus_brain/ATAC")
load("GSE96949/matrix/human_ATAC_colData.RData")
load("GSE96949/matrix/human_ATAC_Neuronal.RData")
load("GSE96949/matrix/rowData.Rdata")
load("matrix/peak_top10000_res/peaks_LRT_res.RData")



res <- results(human_ATAC_Neuronal)
res <- res[complete.cases(res),]
res <- res[order(res$padj), ]
Neuronal_res <- res[res$padj<0.05,] 
dim(Neuronal_res) # 3398    6

inter <- intersect(rownames(peaks_LRT_res), 
                   rowData[rownames(Neuronal_res), "Rhesus_peak"])
length(inter) # 923

length(inter)/dim(Neuronal_res)[1] # 0.2716304



##
area1 = length(rownames(peaks_LRT_res))
area2 = length(rowData[rownames(Neuronal_res), "Rhesus_peak"])
cross.area = length(inter)

pdf("GSE96949/matrix/Neuronal_dif_peaks_human_Rhesus_venn.pdf", width=5, height = 5)
draw.pairwise.venn(area1, area2, cross.area,
                   category = c("Rhesus_dif_peaks", "Human_dif_peaks"),
                   fill=c("#3F60AC", "#8B2052"),
                   col=c("#3F60AC", "#8B2052"),
                   alpha=0.5,
                   print.mode = "raw",
                   cat.pos = c(210,150),
                   lwd = c(0,0),
                   scaled = T,
                   cat.cex=2,
                   main= "Neuronal",
                   cat.dist = c(0.05,0.05),
                   margin=0.05,
                   ind = TRUE,
                   cex = 1.5)
dev.off()



########################################################
#                     Non_Neuronal                     #
########################################################
library(DESeq2)
library(VennDiagram)
setwd("/mnt/data2/Rhesus_brain/ATAC")
load("GSE96949/matrix/human_ATAC_colData.RData")
load("GSE96949/matrix/human_ATAC_Non_Neuronal.RData")
load("GSE96949/matrix/rowData.Rdata")
load("matrix/peak_top10000_res/peaks_LRT_res.RData")



res <- results(human_ATAC_Non_Neuronal)
res <- res[complete.cases(res),]
res <- res[order(res$padj), ]
Non_Neuronal_res <- res[res$padj<0.05,] 
dim(Non_Neuronal_res) # 4798    6

inter <- intersect(rownames(peaks_LRT_res), 
                   rowData[rownames(Non_Neuronal_res), "Rhesus_peak"])
length(inter) # 1357

length(inter)/dim(Non_Neuronal_res)[1] # 0.2828262

##
area1 = length(rownames(peaks_LRT_res))
area2 = length(rowData[rownames(Neuronal_res), "Rhesus_peak"])
cross.area = length(inter)

pdf("GSE96949/matrix/Non_Neuronal_dif_peaks_human_Rhesus_venn.pdf", width=5, height = 5)
draw.pairwise.venn(area1, area2, cross.area,
                   category = c("Rhesus_dif_peaks", "Human_dif_peaks"),
                   fill=c("#3F60AC", "#8B2052"),
                   col=c("#3F60AC", "#8B2052"),
                   alpha=0.5,
                   print.mode = "raw",
                   cat.pos = c(210,150),
                   lwd = c(0,0),
                   scaled = T,
                   cat.cex=2,
                   main= "Non-Neuronal",
                   cat.dist = c(0.05,0.05),
                   margin=0.05,
                   ind = TRUE,
                   cex = 1.5)
dev.off()


########################################################
#                        t-SNE                         #
########################################################
library(Rtsne)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
setwd("/mnt/data2/Rhesus_brain/ATAC")
load("GSE96949/matrix/human_ATAC_Object.RData")
load("GSE96949/matrix/human_ATAC_colData.RData")


vsd <- varianceStabilizingTransformation(human_ATAC_Object, blind=FALSE)
vsd_data <- assay(vsd)


# PCA-based t-SNE
seed = 1:50
per = c(5, 10, 15, 20, 25, 30, 35)
for(i in seed){
  # for(perplexity in per){
  i=2
  set.seed(i)
  perplexity = 3
  tsne <- Rtsne::Rtsne(t(vsd_data), initial_dims=50,
                       max_iter = 5000,perplexity=perplexity)
  rownames(tsne$Y) <- colnames(vsd_data)
  colnames(tsne$Y) <- c("tsne1", "tsne2")
  colData <- data.frame(colData(human_ATAC_Object))
  p <- data.frame(tsne$Y, colData[rownames(tsne$Y),])
  head(p)
  title = paste("t-SNE : seed=", i, " perplexity=", perplexity, sep="")
  
 
  ggplot(p, aes(x=tsne1, y=tsne2, colour=cell_type)) +
  # ggplot(p, aes(x=tsne1, y=tsne2, colour=region)) +
  # ggplot(p, aes(x=tsne1, y=tsne2, colour=ID)) +
    geom_point(size=1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    # scale_colour_manual(values=anno_colors$ID) +
    # scale_colour_manual(values=anno_colors$SR) +
    scale_colour_manual(values=brewer.pal(8,"Dark2")[c(1,2,3,4,6)]) +
    # scale_colour_manual(values=colorRampPalette(brewer.pal(12,"Paired"))(14)) +
    # geom_text(aes(label=rownames(p)), vjust=-1.2, size=0.5, colour="black") +
    # geom_text(aes(label=SR), vjust=-1.2, size=1, colour="black") +
    ggtitle(title) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    xlab("t-SNE1") +
    ylab("t-SNE2")
  
  # 
  # filename = paste("GSE96949/matrix/human_ATAC_t-SNE_region_seed", i, "_per", perplexity, ".pdf", sep="")
  # ggsave(filename,width = 5.6, height = 5)
  # filename = paste("GSE96949/matrix/human_ATAC_t-SNE_ID_seed", i,  "_per", perplexity, ".pdf", sep="")
  # ggsave(filename,width = 6, height = 5)
  filename = paste("GSE96949/matrix/human_ATAC_t-SNE_cell_type_seed", i, "_per", perplexity, ".pdf", sep="")
  ggsave(filename,width = 6, height = 5)
  # }
}









##################################################################
############
library(pheatmap)
library(RColorBrewer)

normal_matrix_N <- normal_matrix[, colData$cell_type=="N"]
normal_matrix_N <- normal_matrix_N[apply(normal_matrix_N, 1, function(x){mean(x)>0}),]
dim(normal_matrix_N)


pheatmap(
  normal_matrix_N,
  col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
  # breaks = c(-4,seq(-3,3,length=20),4),
  # legend_breaks = c(-3,-1,1,3),
  filename = "test_N.pdf",
  scale="row",
  main = "Neuronal",
  border_color=NA,
  annotation_col = colData[colData$cell_type=="N", c("patient", "region", "gender")],
  # annotation_colors = ann_colors,
  # gaps_row = cumsum(table(geneCluster$geneCluster)),
  annotation_legend = T,
  show_rownames = F,
  show_colnames = F,
  # cluster_rows = F
  )




#######  neuronal
cortex_region <- c("VLPFC", "ITC", "STC", "PMC", "DLPFC", "OFC", "ACC", "PVC")

normal_matrix_N <- normal_matrix[, colData$region %in% cortex_region & colData$cell_type == "N"]
normal_matrix_N <- normal_matrix_N[apply(normal_matrix_N, 1, function(x){mean(x)>0}),]
dim(normal_matrix_N)


pheatmap(
  normal_matrix_N,
  col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
  # breaks = c(-4,seq(-3,3,length=20),4),
  # legend_breaks = c(-3,-1,1,3),
  filename = "test_N_cortex.pdf",
  scale="row",
  main = "Neuronal",
  border_color=NA,
  annotation_col = colData[colData$cell_type=="N", c("patient", "region", "gender")],
  # annotation_colors = ann_colors,
  # gaps_row = cumsum(table(geneCluster$geneCluster)),
  annotation_legend = T,
  show_rownames = F,
  show_colnames = F,
  # cluster_rows = F
)






########################################################
library(ggplot2)

seed = 65
set.seed(seed)
perplexity = 2
tsne <- Rtsne::Rtsne(t(normal_matrix_N), perplexity=2)
rownames(tsne$Y) <- colnames(normal_matrix_N)
colnames(tsne$Y) <- c("tsne1", "tsne2")


p <- data.frame(tsne$Y, colData[rownames(tsne$Y),])
head(p)

title = paste("human ATAC t-SNE \n seed:", seed, " perplexity:", perplexity, sep="")

ggplot(p, aes(x=tsne1, y=tsne2, colour=patient)) +
  geom_point(size=4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", size=1),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.title = element_text(face="bold"),
        axis.text = element_text(face="bold", size = 10),
        legend.title = element_text(face="bold")) +
  # scale_colour_manual(values=anno_colors$ID) +
  # scale_colour_manual(values=colorRampPalette(brewer.pal(12,"Paired"))(52)) +
  # geom_text(aes(label=region2), vjust=-1.2, size=1, colour="black") +
  xlab("t-SNE1") +
  ylab("t-SNE2") +
  ggtitle(title)


filename = paste("human_ATAC_t-SNE_ID_seed",seed,"_per", perplexity," .pdf", sep="")
ggsave(filename,width = 6, height = 5)


#####################################
# anova
df <- normal_matrix_N
testInfo <- colData[colnames(normal_matrix_N),]$patient
pvalue <- apply(df, 1, function(x){d <- data.frame(x,testInfo);summary(aov(x~testInfo, data=d))[[1]][1,5]})
padj <- p.adjust(pvalue, method = "fdr", n=length(pvalue))
anova_res <- data.frame(pvalue, padj)
anova_res <- anova_res[order(anova_res$padj),]
anova_res <- anova_res[anova_res$padj<0.05,]
dim(anova_res) #  56  2
head(anova_res)
N_anova_res <- anova_res


###
pheatmap(
  normal_matrix_N[rownames(peak_anova_res),],
  col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
  # breaks = c(-4,seq(-3,3,length=20),4),
  # legend_breaks = c(-3,-1,1,3),
  filename = "test_N_cortex_pheatmap.pdf",
  scale="row",
  main = "Neuronal",
  border_color=NA,
  annotation_col = colData[colData$cell_type=="N", c("patient", "gender")],
  # annotation_colors = ann_colors,
  # gaps_row = cumsum(table(geneCluster$geneCluster)),
  annotation_legend = T,
  show_rownames = F,
  show_colnames = F,
  # cluster_rows = F
)








#######  non neuronal
cortex_region <- c("VLPFC", "ITC", "STC", "PMC", "DLPFC", "OFC", "ACC", "PVC")

normal_matrix_G <- normal_matrix[, colData$region %in% cortex_region & colData$cell_type == "G"]
normal_matrix_G <- normal_matrix_G[apply(normal_matrix_G, 1, function(x){mean(x)>0}),]
dim(normal_matrix_G)





########################################################
library(ggplot2)
df <- normal_matrix_G
seed = 65
set.seed(seed)
perplexity = 2
tsne <- Rtsne::Rtsne(t(df), perplexity=2)
rownames(tsne$Y) <- colnames(df)
colnames(tsne$Y) <- c("tsne1", "tsne2")


p <- data.frame(tsne$Y, colData[rownames(tsne$Y),])
head(p)

title = paste("human ATAC t-SNE \n seed:", seed, " perplexity:", perplexity, sep="")

ggplot(p, aes(x=tsne1, y=tsne2, colour=patient)) +
  geom_point(size=4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", size=1),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.title = element_text(face="bold"),
        axis.text = element_text(face="bold", size = 10),
        legend.title = element_text(face="bold")) +
  # scale_colour_manual(values=anno_colors$ID) +
  # scale_colour_manual(values=colorRampPalette(brewer.pal(12,"Paired"))(52)) +
  # geom_text(aes(label=region2), vjust=-1.2, size=1, colour="black") +
  xlab("t-SNE1") +
  ylab("t-SNE2") +
  ggtitle(title)


filename = paste("human_ATAC_t-SNE_ID_seed",seed,"_per", perplexity," .pdf", sep="")
ggsave(filename,width = 6, height = 5)


#####################################
# anova
df <- normal_matrix_G
testInfo <- colData[colnames(df),]$patient
pvalue <- apply(df, 1, function(x){d <- data.frame(x,testInfo);summary(aov(x~testInfo, data=d))[[1]][1,5]})
padj <- p.adjust(pvalue, method = "fdr", n=length(pvalue))
anova_res <- data.frame(pvalue, padj)
anova_res <- anova_res[order(anova_res$padj),]
anova_res <- anova_res[anova_res$padj<0.05,]
dim(anova_res) #  1905    2
head(anova_res)
G_anova_res <- anova_res


###
pheatmap(
  normal_matrix_N[rownames(peak_anova_res),],
  col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
  # breaks = c(-4,seq(-3,3,length=20),4),
  # legend_breaks = c(-3,-1,1,3),
  filename = "test_N_cortex_pheatmap.pdf",
  scale="row",
  main = "NON-Neuronal",
  border_color=NA,
  annotation_col = colData[colData$cell_type=="N", c("patient", "gender")],
  # annotation_colors = ann_colors,
  # gaps_row = cumsum(table(geneCluster$geneCluster)),
  annotation_legend = T,
  show_rownames = F,
  show_colnames = F,
  # cluster_rows = F
)



#################################################################
#                       peak distrubution                       #
#################################################################
library(ChIPseeker)
library(GenomicFeatures)
library(DESeq2)
setwd("/mnt/data2/Rhesus_brain/ATAC/liftover_res")
load("ATAC/matrix/peak_top10000_res/Rhesus_DESeq2_object.RData")

txdb <- makeTxDbFromGFF("D:/????????????/????/Macaca_mulatta.Mmul_8.0.1.91.chr.gtf.gz",
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
save(all_peaks, file="ATAC/matrix/all_peaks.RData")
write.table(all_peaks, "ATAC/matrix/all_peaks.bed",
            row.names = F, col.names = F, sep="\t", quote = F)
all_peaks <- GRanges(all_peaks)
all_peakAnno  <- annotatePeak(all_peaks, tssRegion=c(-2500, 2500), TxDb = txdb, level = "gene")
### dif peaks 
load("ATAC/matrix/peak_anova_res.RData")
peaks <- rownames(peak_anova_res)
peaks <- data.frame(chrom=sapply(peaks, function(x){strsplit(x, ":")[[1]][1]}),
                    start=sapply(peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1]}),
                    end=sapply(peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2]}))
peaks <- GRanges(peaks)
peakAnno  <- annotatePeak(peaks, tssRegion=c(-2500, 2500), TxDb = txdb, level = "gene")
### cluster peaks
load("ATAC/matrix/ConsensusCluster_DifPeaks.RData")
k=3
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
save(peakAnno_list, file="ATAC/matrix/peakAnno_list.RData")


#
pdf("Distribution_of_peak_relative_to_TSS.pdf", width=5, height = 3)
plotDistToTSS(peakAnno_list, 
              title="Distribution of peak relative to TSS",
              ylab="Genomic Region (%) (5'->3')")
dev.off()
#
pdf("Feature_Distribution_of_peak.pdf", width=7, height = 3)
plotAnnoBar(peakAnno_list)
dev.off()


######################################
## chisq.test for Feature Distribution
######################################
library(gtools)
load("ATAC/matrix/peakAnno_list.RData")
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















