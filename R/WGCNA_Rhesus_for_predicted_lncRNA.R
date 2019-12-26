library(ggplot2)
library(WGCNA)
library(pheatmap)
require(RColorBrewer)
library(DESeq2)
options(stringsAsFactors = FALSE)


##############################################################################
#                          loading Rhesus DESeq object                       #
##############################################################################
#=============================================================================
#  Code chunk 1
#=============================================================================
getwd()
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/novel_lncRNA_trans.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_DESeq2_object.RData")
Rhesus_norCounts <- counts(Rhesus_TransObject,normalized=TRUE)
dim(Rhesus_norCounts) # 51406   408
#=====================================================================================
#  Code chunk 2
#=====================================================================================
datExpr0=as.data.frame(t(Rhesus_norCounts))
datExpr0=datExpr0[, colnames(datExpr0) %in% rownames(novel_lncRNA_trans)]
dim(datExpr0) # 408 2708
#=====================================================================================
#  Code chunk 3
#=====================================================================================
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK # TRUE
#=====================================================================================
#  Code chunk 4
#=====================================================================================
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr0)  ## 408 2708
#=====================================================================================
#  Code chunk 5
#=====================================================================================
sampleTree = hclust(dist(datExpr0), method = "complete");  # complete\average\ward.D2
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(20,5)
# pdf(file = "revise0530/gene_count_res/WGCNA_Rhesus_sampleClustering.pdf", 
#     width = 20, height = 5);
par(mar = c(0,3,2,0))
plot(sampleTree, main = "Rhesus Sample clustering to detect outliers", sub="", 
     xlab="", cex = 0.3, cex.axis = 1.5, cex.main =1.5)
# abline(h = 3000000, col = "red") # Plot a line to show the cut
# dev.off()

# # Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = 3000000, minSize = 5)
# table(clust)
# # clust 1 contains the samples we want to keep.
# keepSamples = (clust==1)
# datExpr = datExpr0[keepSamples, ]
datExpr = datExpr0
nGenes = ncol(datExpr) # 2708
nSamples = nrow(datExpr) # 408
save(datExpr, 
     file="revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA.RData")
#=====================================================================================
#  Code chunk 7
#=====================================================================================
load("revise0530/gene_count_res/Rhesus_colData.RData") # Rhesus_colData
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA.RData")
datTraits <- Rhesus_colData[rownames(datExpr),]
head(datTraits)
dim(datTraits) # [1] 408   9
names(datTraits)
collectGarbage()
#=====================================================================================
#  Code chunk 8
#=====================================================================================
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "complete")
traitColors = labels2colors(datTraits) 
# Plot the sample dendrogram and the colors underneath.
# pdf(file = "revise0530/gene_count_res/WGCNA_Rhesus_sampleClustering_2.pdf", 
#     width = 20, height = 5);
plotDendroAndColors(sampleTree2, traitColors,
                    cex.dendroLabels = 0.5,
                    dendroLabels = datTraits$region2,
                    groupLabels = names(datTraits), 
                    main = "Rhesus Sample dendrogram and trait heatmap")
# dev.off()




#############################################################################################
#           Automatic, one-step network construction and module detection                   #
#############################################################################################
#=====================================================================================
#  Code chunk 1
#=====================================================================================
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()
#=====================================================================================
#  Code chunk 2
#=====================================================================================
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
pdf("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_lncRNA_soft_threshold.pdf",
    width = 9, height = 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#=====================================================================================
#  Code chunk 3
#=====================================================================================
library(WGCNA)
setwd("/mnt/xdlab1/home/looking/brain_project")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA.RData")
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = 3,
                       TOMType = "unsigned",
                       # minModuleSize = 10,
                       # minModuleSize = 20,
                       minModuleSize = 50,
                       reassignThreshold = 0, 
                       # mergeCutHeight = 0.15,
                       # deepSplit=4,
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE
)
save(net, file="revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_net_modules_lncRNA.RData")

######### TOM
enableWGCNAThreads()
TOM = TOMsimilarityFromExpr(datExpr, power = 3);
save(TOM, file="revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_TOM_lncRNA.RData")

####
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_net_modules_lncRNA.RData")
table(net$colors) ##    7 modules
# net$dendrograms[[1]]


#=====================================================================================
#  Code chunk 4
#=====================================================================================
# library(WGCNA)
# setwd("/mnt/xdlab1/home/looking/brain_project")
# load("revise0530/gene_count_res/WGCNA_Rhesus_net_modules_exclude_MB.RData")
# # open a graphics window
# sizeGrWindow(12, 9)
# # Convert labels to colors for plotting
# mergedColors = labels2colors(net$colors)
# # Plot the dendrogram and the module colors underneath
# pdf("revise0530/gene_count_res/WGCNA_47module_cluster_dendrogram.pdf", width=12, height=9)
# plotDendroAndColors(net$dendrograms[[1]], 
#                     mergedColors[net$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = F, hang = 0.02,
#                     autoColorHeight=F,colorHeight = 0.09,
#                     addGuide = TRUE, guideHang = 0.02)
# # dend1 <- net$colors[net$blockGenes[[1]]]
# # names(dend1) <- mergedColors[net$blockGenes[[1]]]
# # legend("topright", legend = unique(dend1),
# #        fill = unique(names(dend1)))
# dev.off()


#=====================================================================================
#  Code chunk 5
#=====================================================================================
library(WGCNA)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_net_modules_lncRNA.RData")


moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs = net$MEs; 
dim(MEs) # 408   8
geneTree = net$dendrograms[[1]]


############################################################################
#                      cor between modules and traits                      #
############################################################################


#############################
#### SR_merge
#############################
load("revise0530/gene_count_res/Rhesus_colData.RData")
datTraits <- Rhesus_colData[rownames(datExpr),]

SR_merge_traits <- matrix(0, ncol=length(unique(datTraits$SR_merge)), 
                          nrow=nrow(datTraits))
rownames(SR_merge_traits) <- rownames(datTraits)
colnames(SR_merge_traits) <- unique(datTraits$SR_merge)
for(a in unique(datTraits$SR_merge)){
  info <- as.character(datTraits$SR_merge)
  info[info==a] <- 1
  info[info!=1] <- 0
  info <- as.integer(info)
  SR_merge_traits[,a] <- info
}
dim(SR_merge_traits) # 408  16

####
moduleColorsWW = moduleColors
MEs0 = moduleEigengenes(datExpr, moduleLabels)$eigengenes # MEs = net$MEs
dim(MEs0) # [1]  408   8
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, SR_merge_traits, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nrow(datExpr))
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
save(modTraitCor, modTraitP, 
     file="revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_modTraitCor_for_SR_merge_lncRNA.RData")


##### plot cor
library(pheatmap)
library(RColorBrewer)
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_modTraitCor_for_SR_merge_lncRNA.RData")
pheatmap(
  modTraitCor,
  main="Module-Region_trait relationships",
  filename="revise0530/transcripts_count_res_with_TACO_gtf/Net_modules_and_SR_merge_lncRNA_cor_pheatmap.pdf",
  col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(40),
  breaks = seq(-1,1,length=40),
  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
  # cellwidth = 9,
  # cellheight = 9.5,
  fontsize=12,
  width = 10, 
  height = 10,
  cluster_rows = T,
  cluster_cols = T,
  show_rownames = T,
  show_colnames = T,
  border_color = NA
)

modTraitCor_help <- modTraitCor
sum(modTraitCor>=0.80 & modTraitP<0.01) 

p <- data.frame(gene=rownames(modTraitCor), modTraitCor)
write.table(p, "SourceData/Fig.S6a.txt",
            col.names = T, row.names = F, quote=F, sep="\t")




########### 保存 module trans_id
library(WGCNA)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_net_modules_lncRNA.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_modTraitCor_for_SR_merge_lncRNA.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_TOM_lncRNA.RData")
dim(datExpr) # 248 2782
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
moduleLabels2 <- paste("Lnc-MA", net$colors, sep="")

rownames(modTraitCor) <- gsub("ME", "Lnc-MA", rownames(modTraitCor))
region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)
probes = names(datExpr)



for(module in region_modules){
  region <- names(region_modules)[which(region_modules==module)]
  inModule = is.finite(match(moduleLabels2, module));
  modProbes = probes[inModule];
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  kIM = intramodularConnectivity(modTOM, moduleColors[inModule], scaleByMax = TRUE)
  kIM <- kIM[order(kIM$kWithin, decreasing=T),1:2]
  kIM$module <- module
  kIM <- data.frame(gene_id=rownames(kIM), kIM, individual=region)
  filename=paste("revise0530/transcripts_count_res_with_TACO_gtf/module_", module, "_", region, "_info.txt", sep = "")
  write.table(kIM, filename, row.names = F, col.names = T, quote = F, sep = "\t")
}



#### bar plots showing the expression of MEs across brain
library(reshape2)
library(wesanderson)
library(reshape2)
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_modTraitCor_for_SR_merge_lncRNA.RData")
region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)

pp <- data.frame(MEsWW[, region_modules], SR=datTraits$SR_merge)
pp <- melt(pp)
head(pp)


colors <- c(ME3="#cfbddc",
            ME2="#247ba0",
            ME7="#74a089", 
            ME4="#f25f5c",
            ME5="#ffe066")

ggplot(pp, aes(x=SR, y=value, fill=variable, color=variable)) +
  geom_boxplot(outlier.size = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        text = element_text(face="bold"),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white")) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  facet_wrap(~variable,ncol=2, scales="free_y") +
  stat_summary(geom = "point", shape=23, fill="white",size=1,fun.y="median") +
  guides(fill=F, color=F) +
  xlab("") +
  ylab("Module eigengene expression")

ggsave("revise0530/transcripts_count_res_with_TACO_gtf/Module_eigengene_expression_lncRNA.pdf", 
       width = 4, height = 4.5)


















