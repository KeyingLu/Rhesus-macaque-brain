library(DESeq2)



###################################################
#                     human                       #
###################################################
library(gtools)
library(DESeq2)
library(BiocParallel)
setwd("/mnt/data2/Rhesus_brain")
load("human/PsychENCODE/RPKM.RData")
load("human/PsychENCODE/human_colData.RData")
load("human/PsychENCODE/expr_count.RData")

human_logRPKM <- log2(RPKM + 1)
dim(human_logRPKM) # 27932    80
human_count <- expr_count
dim(human_count)
head(human_colData)

rownames(human_count) <- sapply(rownames(human_count), function(x){strsplit(x, "[.]")[[1]][1]})
write.table(human_colData, "human/PsychENCODE/human_sample_information.txt",
            col.names = T, row.names = F, quote=F, sep="\t")


#####################################################
#                       Rhesus                      #
#####################################################
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
colData <- Rhesus_colData

cortex_mapping_list <- list(A1C="STG", DFC="MFG", IPC=c("AMG&PMG", "ALG"), 
                            ITC="ITG", M1C=c("SPrCG", "IPrCG"), MFC="ACG",
                            OFC="MOG", S1C=c("SPCG", "IPCG"), STC="STG",
                            V1C=c("aSOG", "pSOG"), VFC="IFG")
cortex_mapping <- melt(cortex_mapping_list)
colnames(cortex_mapping) <- c("Rhesus", "human")

colData <- colData[colData$SR %in% c(as.character(cortex_mapping$Rhesus), "AMY", "CBC", "THA") | colData$SR_merge %in% c("HIP", "STR"), ]


Rhesus_count <- read.table("revise0530/Rhesus_count.matrix", header=T, row.names=1)
sub_Rhesus_count <- Rhesus_count[, rownames(colData)]
dim(sub_Rhesus_count) # 30807   184



############################################################
#                      orth count                          #
############################################################
orth_genes <- read.table("ANNOTATION/humanGRCH38_Rhesus_orth_genes_one2one.txt", stringsAsFactors = F)
colnames(orth_genes) <- c("Human_gene_id", "Rhesus_gene_id", "Homologous_type")
rownames(orth_genes) <- orth_genes$Human_gene_id
head(orth_genes)
dim(orth_genes) # 19810     3


sub_Rhesus <- sub_Rhesus_count
orth_Rhesus <- sub_Rhesus[rownames(sub_Rhesus) %in% orth_genes$Rhesus_gene_id,]
orth_human <- human_count[rownames(human_count) %in% orth_genes$Human_gene_id,]
rownames(orth_human) <- orth_genes[rownames(orth_human),]$Rhesus_gene_id
genes <- intersect(rownames(orth_Rhesus), rownames(orth_human))
length(genes) # 15202


orth_Counts <- data.frame(orth_Rhesus[genes,], orth_human[genes,])
dim(orth_Counts) # 15202   264


remain = apply(orth_Counts, 1, function(x){sum(x>1)>=3})
orth_Counts <- orth_Counts[remain,]
dim(orth_Counts)  # 14528   264
save(orth_Counts, file= "human/PsychENCODE/orth_Counts.RData")




#################################################################
#                      orth_colData                             #
#################################################################
load("human/PsychENCODE/human_colData.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")


##
cortex_mapping_list <- list(A1C="STG", DFC="MFG", IPC=c("AMG&PMG", "ALG"), 
                            ITC="ITG", M1C=c("SPrCG", "IPrCG"), MFC="ACG",
                            OFC="MOG", S1C=c("SPCG", "IPCG"), STC="STG",
                            V1C=c("aSOG", "pSOG"), VFC="IFG")
cortex_mapping <- melt(cortex_mapping_list)
colnames(cortex_mapping) <- c("Rhesus", "human")
cortex_mapping <- cortex_mapping[-which(cortex_mapping$Rhesus=="STG"),]
rownames(cortex_mapping) <- as.character(cortex_mapping$Rhesus)

##
colData <- Rhesus_colData[colnames(orth_Rhesus), ]
colData$SR[colData$SR_merge %in% c("HIP", "STR")] <- colData$SR_merge[colData$SR_merge %in% c("HIP", "STR")]
colData_1 <- colData[, c("SR", "ID")]
colnames(colData_1) <- c("region", "ID")
colData_1 <- data.frame(colData_1, Species="Rhesus")
colData_1$region[colData_1$region %in% cortex_mapping$Rhesus] <- cortex_mapping[colData_1$region[colData_1$region %in% cortex_mapping$Rhesus], "human"]


colData <- human_colData[colnames(orth_human), ]
colData_2 <- colData[, c("Region", "Brain")]
colnames(colData_2) <- c("region", "ID")
colData_2 <- data.frame(colData_2, Species="human")
colData_2$region[colData_2$region %in% c("A1C", "STC")] <- "STG"
colData_2$region[colData_2$region %in% c("MD")] <- "THA"

orth_colData <- rbind(colData_1, colData_2)
orth_colData$Region <- orth_colData$region
orth_colData$Region[!(orth_colData$region %in% c("STR", "HIP", "AMY", "THA", "CBC"))] <- "cortex"

save(orth_colData, file= "human/PsychENCODE/orth_colData.RData")




#################################################################
#                            DEseq2                             #
#################################################################
library(DESeq2)
library(BiocParallel)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("human/PsychENCODE/orth_Counts.RData")
load("human/PsychENCODE/orth_colData.RData")

dim(orth_Counts) #14528   264
orth_colData$SpeciesRegion <- paste(orth_colData$region, orth_colData$Species, sep="_")
Rhesus_human_Object <- DESeqDataSetFromMatrix(orth_Counts, 
                                              orth_colData, 
                                              ~SpeciesRegion)
Rhesus_human_Object <- DESeq(Rhesus_human_Object, 
                              parallel = T, BPPARAM=MulticoreParam(20))
save(Rhesus_human_Object, 
     file="human/PsychENCODE/Rhesus_human_Object.RData")




#######################################################
#               t-SNE using DEseq2 results            #
#######################################################
library(Rtsne)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(limma)
setwd("/mnt/data2/Rhesus_brain")
load("human/PsychENCODE/orth_colData.RData")
load("human/PsychENCODE/Rhesus_human_Object.RData")

vsd <- varianceStabilizingTransformation(Rhesus_human_Object, blind=FALSE)
vsd_data <- limma::removeBatchEffect(assay(vsd), vsd$Species)
df <- vsd_data
dim(df) # 14528   264

vsd_data_save <- data.frame(GeneID=rownames(vsd_data), 
                            vsd_data)
write.table(vsd_data_save, "SourceData/SourceData_Fig.2a.txt",
            row.names = F, col.names = T, quote = F, sep="\t")


seed = 1
for(i in seed){
  i=42
  set.seed(i)
  perplexity = 4
  tsne <- Rtsne::Rtsne(t(df), initial_dims=50, 
                       max_iter = 5000,perplexity=perplexity)
  rownames(tsne$Y) <- colnames(df)
  colnames(tsne$Y) <- c("tsne1", "tsne2")
  p <- data.frame(tsne$Y, orth_colData[rownames(tsne$Y),])
  p$Region <- p$region
  p$Region[!(p$region %in% c("STR", "HIP", "AMY", "THA", "CBC"))] <- "cortex"
  head(p)
  title = paste("t-SNE : seed=", i, " perplexity=", perplexity, sep="")
  
  ggplot(p, aes(x=tsne1, y=tsne2, colour=Region, shape=Species)) +
    geom_point(size=1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(colour="black"),
          axis.text = element_text(colour="black"),
          legend.text = element_text(colour="black"),
          legend.title = element_text(colour="black")) +
    scale_colour_manual(values=c(anno_colors$SR, anno_colors$merge)[unique(p$Region)]) +
    scale_shape_manual(values = c(19,3))+
    # scale_colour_manual(values=anno_colors$Region_cortex) +
    # scale_colour_manual(values=colorRampPalette(brewer.pal(12,"Paired"))(length(unique(p$region)))) +
    # scale_colour_manual(values=colorRampPalette(brewer.pal(12,"Paired"))(length(unique(p$ID)))) +
    ggtitle(title) +
    xlab("t-SNE1") +
    ylab("t-SNE2") +
    guides(colour = guide_legend(override.aes = list(size=2)))
  
  
  filename = paste("human/PsychENCODE/human_Rhesus_Gene_t-SNE_Region_seed", i, "_per", perplexity, ".pdf", sep="")
  ggsave(filename,width = 5.8, height = 5)
  # filename = paste("human/PsychENCODE/human_Gene_t-SNE_ID_seed", i, ".pdf", sep="")
  # ggsave(filename,width = 5.8, height = 5)
}


p <- data.frame(Sample_ID=rownames(p), p)
write.table(p, "SourceData/Fig.2a.txt",
            col.names = T, row.names = F, quote=F, sep="\t")



################################################################
#         species specific genes in region (wilcox test)       #
################################################################
library(gclus)
library(limma)
setwd("/mnt/data2/Rhesus_brain")
load("human/PsychENCODE/orth_colData.RData")
load("human/PsychENCODE/Rhesus_human_Object.RData")

vsd <- varianceStabilizingTransformation(Rhesus_human_Object, blind=FALSE)
norCounts <- counts(Rhesus_human_Object, normalized=T)
norCounts <- limma::removeBatchEffect(norCounts, vsd$Species)
dim(norCounts) # 14528   264


human_Rhesus_inEachRegion_genes_list <- list()
human_Rhesus_inEachRegion_res <- list()
for(region in unique(orth_colData$region)){
  df <- norCounts[, orth_colData$region %in% region]
  colData <- orth_colData[colnames(df),]
  testInfo <- as.character(colData$Species)
  pvalue <- apply(df, 1, function(x){wilcox.test(x ~ testInfo)$p.value})
  padj_fdr <- p.adjust(pvalue, method = "fdr", n=length(pvalue))
  padj_BH <- p.adjust(pvalue, method = "BH", n=length(pvalue))
  mean <-  sapply(unique(testInfo), function(x){rowMeans(df[,testInfo==x])})
  p <- data.frame(pvalue, padj=padj_fdr)
  log2FoldChange <- log2(mean[,"human"]/mean[,"Rhesus"])
  res <- data.frame(p, log2FoldChange)
  res <- res[order(res$padj),]
  # hyper <- rownames(res[which(res$padj < 0.05 & res$log2FoldChange > 0.5),])
  # hypo <- rownames(res[which(res$padj < 0.05 & res$log2FoldChange < -0.5),])
  hyper <- rownames(res[which(res$padj < 0.05 & res$log2FoldChange > 1),])
  hypo <- rownames(res[which(res$padj < 0.05 & res$log2FoldChange < -1),])
  human_Rhesus_inEachRegion_genes_list[[region]] <- list(human=hyper, Rhesus=hypo)
  human_Rhesus_inEachRegion_res[[region]] <- res
}

save(human_Rhesus_inEachRegion_genes_list, 
     human_Rhesus_inEachRegion_res,
     file="human/PsychENCODE/human_Rhesus_inEachRegion_dif_genes.RData")


####
library(reshape2)
load("human/PsychENCODE/human_Rhesus_inEachRegion_dif_genes.RData")

pp <- melt(human_Rhesus_inEachRegion_genes_list)
pp <- melt(table(pp$L1, pp$L2))
colnames(pp) <-c("region", "Species", "count")
pp$count[which(pp$Species=="Rhesus")] <- -pp$count[which(pp$Species=="Rhesus")]
pp <- pp[order(pp$count), ]
pp$region <- as.character(pp$region)
pp$region <- factor(pp$region, levels = unique(pp$region))


ggplot(pp, aes(x=region, y=count, fill=Species)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, color = "black"),
        axis.text.y = element_text(color = "black"),
        # strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white")) +
  scale_fill_manual(values =  c(Rhesus=brewer.pal(9,"Blues")[6], human=brewer.pal(9,"Reds")[6])) 


ggsave("human/PsychENCODE/species_specific_genes_bar.pdf", 
       width = 8, height = 5)



###
lapply(human_Rhesus_inEachRegion_res, dim)
pp <- melt(human_Rhesus_inEachRegion_genes_list)
pp <- melt(table(pp$L1, pp$L2))
colnames(pp) <-c("region", "Species", "count")
pp$VS <- as.character(pp$Species)
pp$VS[pp$Species=="human"] <- "H>M"
pp$VS[pp$Species=="Rhesus"] <- "H<M"

region_sum <- sapply(levels(pp$region), function(x){sum(pp$count[pp$region==x])})
H_M <- 14528 - region_sum
pp <- rbind(pp, data.frame(region=names(H_M), Species="-", count=H_M, VS="H=M"))
pp <- pp[order(pp$count), ]
pp$region <- as.character(pp$region)
pp$region <- factor(pp$region, levels = unique(pp$region))
pp$size <- pp$count
pp$size[pp$count<=100] <- "0-100"
pp$size[pp$count>100 & pp$count<=200] <- "100-200"
pp$size[pp$count>200 & pp$count<=300] <- "200-300"
pp$size[pp$count>300 & pp$count<=600] <- "300-600"
pp$size[pp$count>600 & pp$count<=700] <- "600-700"
pp$size[pp$count>700 & pp$count<=900] <- "700-900"
pp$size[pp$count>900 & pp$count<=1100] <- "900-1100"
pp$size[pp$count>1100 & pp$count<=15000] <- "1100-15000"
pp$size <- factor(pp$size, levels=unique(pp$size))


ggplot(pp, aes(x=region, y=VS, size=size, color=VS)) +
  geom_point() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.background = element_rect(fill="white")) +
  scale_size_discrete(breaks=c("0-100","100-200","200-300","300-600",
                               "600-700","700-900", "900-1100", 
                               "1100-15000"),
                      labels=c("0-100","100-200","200-300","300-600",
                               "600-700","700-900", "900-1100", 
                               "1100-15000")) +
  scale_color_manual(values =c("H=M"="#BFBEBE", "H>M"="#E8192C", "H<M"="#4E83B3")) +
  scale_y_discrete(limits=c("H<M", "H>M", "H=M")) 


ggsave("human/PsychENCODE/species_specific_genes_point.pdf", 
       width = 8, height = 5)




ggplot(pp, aes(x=cluster, y=motif, size=logQ)) + 
  geom_point(colour="#7D53A3", alpha=0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        plot.title = element_text(hjust = 0.5, colour="black"),
        axis.title = element_text(size=16, colour="black"),
        axis.text = element_text(size=13, colour="black")) +
  scale_size_continuous(breaks = c(1,100, 200, 300,400),
                        # range=c(0.5,3.5),
                        labels=c("<2","100","200","300", ">400")) +
  labs(size="-logQ") +
  scale_y_discrete(limits=rev(levels(pp$motif))) +
  xlab("") +
  ylab("")

ggsave("ATAC/matrix/peak_top10000_res/cluster_motif_Qvalue_least_200.pdf",
       width=6, height = 10)



################################################################
#            write species specific genes in region            #
################################################################
setwd("/mnt/data2/Rhesus_brain")
load("human/PsychENCODE/human_Rhesus_inEachRegion_dif_genes.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


for(r in names(human_Rhesus_inEachRegion_res)){
  human <- human_Rhesus_inEachRegion_genes_list[[r]][["human"]]
  Rhesus <- human_Rhesus_inEachRegion_genes_list[[r]][["Rhesus"]]
  res <- human_Rhesus_inEachRegion_res[[r]]
  human <- res[human,]
  human <- data.frame(rownames(human), 
                      gene_name= gtf_ensembl_gene[rownames(human), "gene_name"],
                      human)
  Rhesus <- res[Rhesus,]
  Rhesus <- data.frame(rownames(Rhesus), 
                      gene_name= gtf_ensembl_gene[rownames(Rhesus), "gene_name"],
                      Rhesus)
  write.table(human, paste("human/PsychENCODE/human_", r, "_hyper_genes.txt"),
              row.names = F, col.names = T, sep="\t", quote = F)
  write.table(Rhesus, paste("human/PsychENCODE/Rhesus_", r, "_hyper_genes.txt"),
              row.names = F, col.names = T, sep="\t", quote = F)
}


################################################################
#        write species specific genes in region to a txt       #
################################################################
setwd("/mnt/data2/Rhesus_brain")
load("human/PsychENCODE/human_Rhesus_inEachRegion_dif_genes.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


write_list <- list()
for(r in names(human_Rhesus_inEachRegion_res)){
  human <- human_Rhesus_inEachRegion_genes_list[[r]][["human"]]
  Rhesus <- human_Rhesus_inEachRegion_genes_list[[r]][["Rhesus"]]
  res <- human_Rhesus_inEachRegion_res[[r]]
  human <- res[human,]
  human <- data.frame(gene_id=rownames(human), 
                      gene_name=gtf_ensembl_gene[rownames(human), "gene_name"],
                      human, Region=r,hyper_species="Human")
  Rhesus <- res[Rhesus,]
  Rhesus <- data.frame(gene_id=rownames(Rhesus), 
                       gene_name= gtf_ensembl_gene[rownames(Rhesus), "gene_name"],
                       Rhesus, Region=r, hyper_species="Rhesus")
  write_list[[r]] <- rbind(human, Rhesus)
}

pp <- Reduce(rbind, write_list)
write.table(pp, "human/PsychENCODE/species_specific_genes_in_each_region.txt",
      row.names = F, col.names = T, sep="\t", quote = F)




###################################################################
#      AMY/HIP vs Cortex in species, respectively (distance)      #
###################################################################
library("lmPerm")
library("coin")
library("gtools")
library(reshape2)
library(venn)
library(ggplot2)
setwd("/mnt/data2/Rhesus_brain")
load("human/PsychENCODE/orth_colData.RData")
load("human/PsychENCODE/Rhesus_human_Object.RData")

vsd <- varianceStabilizingTransformation(Rhesus_human_Object, blind=FALSE)
vsd_data <- limma::removeBatchEffect(assay(vsd), vsd$Species)
df <- vsd_data
dim(df) # 14528   264


### distance
distance_res <- list()
combs <- combinations(length(unique(orth_colData$Region)), 2, unique(orth_colData$Region))

for(i in 1:dim(combs)[1]){
  region1 <- combs[i, 1]
  region2 <- combs[i, 2]
  for(species in c("human", "Rhesus")){
    sub_colData <- orth_colData[orth_colData$Region %in% c(region1, region2) & orth_colData$Species == species, ]
    df <- vsd_data[, rownames(sub_colData)]
    comb <- paste(region1, region2, species, sep="_")
    vs <- paste(region1, region2, sep="_")
    a <- df[,sub_colData$Region==region1]
    b <- df[, sub_colData$Region==region2]
    colnames(a) <- paste(region1, colnames(a), sep="_")
    colnames(b) <- paste(region2, colnames(b), sep="_")
    distance <- dist(t(cbind(a, b)))
    distance <- melt(as.matrix(distance))
    n <- sapply(as.character(distance[,1]), function(x){strsplit(x, "_")[[1]][1]})
    m <- sapply(as.character(distance[,2]), function(x){strsplit(x, "_")[[1]][1]})
    p1 <- sapply(as.character(distance[,1]), function(x){strsplit(x, "_")[[1]][2]})
    p2 <- sapply(as.character(distance[,2]), function(x){strsplit(x, "_")[[1]][2]})
    distance <- distance[n!=m & p1==p2, ]
    distance_res[[comb]] <- data.frame(distance, region=vs, species=species)
  }
}


pp <- Reduce(rbind, distance_res)


sapply(levels(pp$region), function(x){wilcox.test(value~species, d=pp[pp$region==x,])$p.value})

pp <- pp[pp$region %in% c("AMY_cortex", "AMY_CBC", "cortex_HIP", "CBC_HIP"), 
        c("value", "region", "species")]

ggplot(pp, aes(y=value, x=species, fill=species, color=species)) +
  geom_boxplot(alpha=0.5) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.5), 
             color="black", size=0.2, alpha=0.5) +
  # geom_jitter(size=0.1, color="#3F60AC", width = 0.2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(hjust = 0.5, vjust=0.5, color = "black"),
        axis.text.y = element_text(color = "black", size=10),
        strip.text = element_text(size=9.2),
        strip.background = element_rect(fill="white", color = "white")) +
  scale_color_manual(values=c(brewer.pal(12, "Paired")[2], brewer.pal(12, "Paired")[6])) +
  scale_fill_manual(values=c(brewer.pal(12, "Paired")[2], brewer.pal(12, "Paired")[6])) +
  facet_wrap(~region, ncol=5,scales = "free_x")+
  ylab("Distance")


ggsave("human/PsychENCODE/Distance_AMY_or_HIP_vs_cortex_CBC_in_same_species1226.pdf",
       width =10, height = 3)




##############################################################################
#                                 correlation                                #
##############################################################################
library("lmPerm")
library("coin")
library("gtools")
library(reshape2)
library(venn)
library(ggplot2)
setwd("/mnt/data2/Rhesus_brain")
load("human/PsychENCODE/orth_colData.RData")
load("human/PsychENCODE/Rhesus_human_Object.RData")

vsd <- varianceStabilizingTransformation(Rhesus_human_Object, blind=FALSE)
vsd_data <- limma::removeBatchEffect(assay(vsd), vsd$Species)
df <- vsd_data
dim(df) # 14528   264


for(top in c(1000, 2000, 5000)){
  top=1000
  df <- vsd_data
  mean <- rowMeans(df)
  df <- df[order(mean, decreasing=T)[1:top],]
  pp <- cor(df)
  pp <- melt(pp)
  pp$species1 <- orth_colData[pp$Var1,]$Species
  pp$species2 <- orth_colData[pp$Var2,]$Species
  pp <- pp[pp$species1!=pp$species2,]
  pp$region1 <- orth_colData[pp$Var1,]$Region
  pp$region2 <- orth_colData[pp$Var2,]$Region
  pp <- pp[1:(dim(pp)[1]/2),]
  head(pp)
  pp$region1 <- paste(pp$species1, pp$region1, sep="_")
  pp$region2 <- paste(pp$species2, pp$region2, sep="_")
  
  ggplot(pp, aes(y=value, x="")) +
    geom_violin(color="black") +
    geom_jitter(size=0.1, color="#4292C6", width = 0.2) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face="bold"),
          axis.title.y = element_text(size = 12),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = "black"),
          axis.text.y = element_text(color = "black", size=10),
          text = element_text(face="bold"),
          strip.text = element_text(size=9.2),
          strip.background = element_rect(fill="white", color = "white")) +
    facet_grid(region1~region2, scales = "free_x")+
    xlab("") +
    ylab("Correlation")
  
  filename = paste("human/PsychENCODE/Rhesus_human_cor_violin_ExprTop", top, ".pdf")
  ggsave(filename, height = 6.4, width = 6.6)
}




