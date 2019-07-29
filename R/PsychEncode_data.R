library(ggplot2)



##############################################################################
#                            loading file                                    #
##############################################################################
setwd("/mnt/data2/Rhesus_brain")


RPKM <- read.table("human/PsychENCODE/nhp_development_RPKM_rmTechRep.txt", 
                   header = T, sep= "\t", row.names = 1)
expr_count <- read.table("human/PsychENCODE/nhp_development_count_rmTechRep.txt", 
                         header = T, sep= "\t", row.names = 1)

colData <- read.table("human/PsychENCODE/mRNA-seq_QC.txt", 
                      sep="\t", header=T, row.names = 1, stringsAsFactors = F)
colData <- colData[, c("Species", "Brain", "Region", "NCXRegion", "Age", "Days", "Sequencing.site")]
colData <- colData[colData$Species=="Human",]
# colData <- colData[(colData$Species=="Macaque" & colData$Age %in% c("4Y", "7Y",  "5Y",  "11Y", "3.5Y")) |
#                      (colData$Species=="Human" & colData$Age %in% c("11 Y", "13 Y", "15 Y", "19 Y", "21 Y", "23 Y", "30 Y", "36 Y", "37 Y", "40 Y")),]

type = gsub("[0-9]|[.]|[ ]+", "", colData$Age)
colData <- colData[type %in% c("Y", "M"), ]
unique(colData$Age)

re_Brain <- names(table(colData$Brain))[table(colData$Brain)>=16]
colData <- colData[colData$Brain %in% re_Brain, ]

table(colData$Brain, colData$Age)
table(colData$NCXRegion, colData$Brain)
# colData <- colData[colData$Age %in% c("11 Y", "13 Y", "15 Y", "19 Y", "21 Y", "23 Y", "30 Y", "36 Y", "37 Y", "40 Y"),]
colData$NCXRegion[colData$NCXRegion=="NCX"] <- "cortex"
colData$NCXRegion[colData$NCXRegion=="MD"] <- "THA"
colData <- colData[colData$Sequencing.site=="YALE",]

dim(colData) # 80  7

colnames(expr_count) <- gsub("[.]", "_", colnames(expr_count))
rownames(colData) <- gsub("[.]", "_", rownames(colData))
colnames(RPKM) <- gsub("[.]", "_", colnames(RPKM))
expr_count <- expr_count[,rownames(colData)]
RPKM <- RPKM[,rownames(colData)]
human_colData <- colData

save(expr_count, file="human/PsychENCODE/expr_count.RData")
save(human_colData, file="human/PsychENCODE/human_colData.RData")
save(RPKM, file="human/PsychENCODE/RPKM.RData")






#######################################################
#                human t-SNE with RPKM                #
#######################################################
library(Rtsne)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
setwd("/mnt/data2/Rhesus_brain")
load("human/PsychENCODE/RPKM.RData")
load("human/PsychENCODE/human_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")


human_logRPKM <- log2(RPKM + 1)
dim(human_logRPKM) # 27932   80

df <- human_logRPKM
# PCA-based t-SNE
seed = 11:30
for(i in seed){
  i=19
  set.seed(i)
  perplexity = 4
  tsne <- Rtsne::Rtsne(t(df), initial_dims=50, 
                       max_iter = 5000,perplexity=perplexity)
  rownames(tsne$Y) <- colnames(df)
  colnames(tsne$Y) <- c("tsne1", "tsne2")
  p <- data.frame(tsne$Y, human_colData[rownames(tsne$Y),])
  str(p)
  head(p)
  title = paste("t-SNE : seed=", i, " perplexity=", perplexity, sep="")
  
  ## NCXRegion, Brain
  ggplot(p, aes(x=tsne1, y=tsne2, colour=Brain)) +
    geom_point(size=0.3) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face="bold"),
          axis.title = element_text(face="bold"),
          axis.text = element_text(face="bold")) +
    # scale_colour_manual(values=c(anno_colors$SR, anno_colors$merge)[unique(p$NCXRegion)]) +
    scale_colour_manual(values=brewer.pal(12,"Paired")) +
    # scale_colour_manual(values=brewer.pal(12,"Paired")[c(2,4,6,8,10,12)]) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    ggtitle(title)
  
  # filename = paste("human/PsychENCODE/PsychENCODE_region_t-SNE_seed", i, "_per", perplexity, ".pdf", sep="")
  filename = paste("human/PsychENCODE/PsychENCODE_ID_t-SNE_seed", i, "_per", perplexity, ".pdf", sep="")

  ggsave(filename, width = 5, height = 4)
}















