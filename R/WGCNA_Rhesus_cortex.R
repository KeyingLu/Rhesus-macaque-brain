library(ggplot2)
library(WGCNA)
library(pheatmap)
require(RColorBrewer)
library(DESeq2)
options(stringsAsFactors = FALSE)


################################################################
#               loading Rhesus DESeq object                    #
################################################################
#===============================================================
#  Code chunk 1
#===============================================================
getwd()
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_GeneObject_cortex.RData") ## Rhesus_object
Rhesus_Gene_norCounts <- counts(Rhesus_GeneObject_cortex,normalized=TRUE)
dim(Rhesus_Gene_norCounts) # [1] 21306   248
#===============================================================
#  Code chunk 2
#===============================================================
datExpr0=t(Rhesus_Gene_norCounts)
dim(datExpr0) # 248 21306
#===============================================================
#  Code chunk 3
#===============================================================
gsg = goodSamplesGenes(datExpr0, verbose = 3)
# Flagging genes and samples with too many missing values...
# ..step 1
gsg$allOK
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
dim(datExpr0)  ##  248 21306

#=====================================================================================
#  Code chunk 5
#=====================================================================================
sampleTree = hclust(dist(datExpr0), method = "complete");  # complete\average\ward.D2
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(15,5)
# pdf(file = "revise0530/gene_count_res/WGCNA_Rhesus_sampleClustering.pdf", 
#     width = 20, height = 5);
par(mar = c(0,3,2,0))
plot(sampleTree, main = "Rhesus Sample clustering to detect outliers", sub="", 
     xlab="", cex = 0.3, cex.axis = 1.5, cex.main =1.5)
# abline(h = 3000000, col = "red") # Plot a line to show the cut
dev.off()

# # Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = 3000000, minSize = 5)
# table(clust)
# # clust 1 contains the samples we want to keep.
# keepSamples = (clust==1)
# datExpr = datExpr0[keepSamples, ]
datExpr = datExpr0
nGenes = ncol(datExpr) # 21306
nSamples = nrow(datExpr) # 248
save(datExpr, file="revise0530/gene_count_res/WGCNA_Rhesus_datExpr_cortex.RData")
#=====================================================================================
#  Code chunk 7
#=====================================================================================
load("revise0530/gene_count_res/Rhesus_colData.RData") # Rhesus_colData
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_cortex.RData")
datTraits <- Rhesus_colData[rownames(datExpr),]
head(datTraits)
dim(datTraits) # 248   9
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
dev.off()




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
pdf("revise0530/gene_count_res/WGCNA_Rhesus_soft_threshold_for_cortex.pdf",
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
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_cortex.RData")
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = 9,
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
save(net, file="revise0530/gene_count_res/WGCNA_Rhesus_net_modules_cortex.RData")

table(net$colors) ## 32
# net$dendrograms[[1]]
TOM = TOMsimilarityFromExpr(datExpr, power = 9);
save(TOM, file="revise0530/gene_count_res/WGCNA_Rhesus_TOM_cortex.RData")


#=====================================================================================
#  Code chunk 4
#=====================================================================================
library(WGCNA)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/WGCNA_Rhesus_net_modules_cortex.RData")
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
pdf("revise0530/gene_count_res/WGCNA_47module_cluster_dendrogram.pdf", width=12, height=9)
plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = F, hang = 0.02,
                    autoColorHeight=F,colorHeight = 0.09,
                    addGuide = TRUE, guideHang = 0.02)
# dend1 <- net$colors[net$blockGenes[[1]]]
# names(dend1) <- mergedColors[net$blockGenes[[1]]]
# legend("topright", legend = unique(dend1),
#        fill = unique(names(dend1)))
dev.off()






###################################################################
#                 cor between modules and traits                  #
###################################################################
library(WGCNA)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/WGCNA_Rhesus_net_modules_cortex.RData")

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs = net$MEs; 
dim(MEs) # 248  33
geneTree = net$dendrograms[[1]]


##############
#### lobe(in cortex)
##############
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_cortex.RData")
datTraits <- Rhesus_colData[rownames(datExpr),]
Info <- datTraits$BR

ID_cortex_traits <- matrix(0, ncol=length(unique(Info)), nrow=length(Info))
rownames(ID_cortex_traits) <- rownames(datTraits)
colnames(ID_cortex_traits) <- unique(unique(Info))
for(a in unique(Info)){
  info <- Info
  info[info==a] <- 1
  info[info!=1] <- 0
  info <- as.integer(info)
  ID_cortex_traits[,a] <- info
}
dim(ID_cortex_traits) # 248   5

####
moduleColorsWW = moduleColors
MEs0 = moduleEigengenes(datExpr, moduleLabels)$eigengenes # MEs = net$MEs
dim(MEs0) # 248  32
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, ID_cortex_traits, use = "p", method = "pearson")
modTraitP = corPvalueStudent(modTraitCor, nrow(datExpr))
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
save(modTraitCor, modTraitP, file="revise0530/gene_count_res/WGCNA_modTraitCor_for_lobe_cortex.RData")


##### plot cor
library(pheatmap)
library(RColorBrewer)
load("revise0530/gene_count_res/WGCNA_modTraitCor_for_lobe_cortex.RData")
pheatmap(
  t(modTraitCor),
  main="Module-lobe(cortex)_trait relationships",
  filename="revise0530/gene_count_res/Net_modules_and_lobe_cortex_traits_cor_pheatmap.pdf",
  col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(40),
  breaks = seq(-1,1,length=40),
  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
  # cellwidth = 9,
  # cellheight = 9.5,
  # fontsize_row = 6,
  width = 10, 
  height = 10,
  cluster_rows = T,
  cluster_cols = T,
  show_rownames = T,
  show_colnames = T,
  border_color = NA
)

### plot cor help
sum(modTraitCor>=0.80 & modTraitP<0.01) 
sum(modTraitCor <= -0.80 & modTraitP<0.01)
modTraitCor_help <- modTraitCor

p <- data.frame(Module=rownames(modTraitCor), modTraitCor)
write.table(p, "SourceData/Fig.3d.txt",
            col.names = T, row.names = F, quote=F, sep="\t")


#### bar plots showing the expression of MEs across brain
library(reshape2)
library(wesanderson)
region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)

pp <- data.frame(MEsWW[, region_modules, drop=F], BR=datTraits$BR)
pp <- melt(pp)
head(pp)
colnames(pp)[2:3] <- c("Module", "expression")
pp$Module <- as.character(pp$Module)
p <- matrix(0, ncol=4, nrow=length(unique(pp$BR))*length(region_modules))
colnames(p) <- c("BR", "Module", "mean", "se")

j=1
for(r in unique(pp$BR)){
  for(i in unique(pp$Module)){
    mean <- mean(pp$expression[pp$BR==r & pp$Module==i])
    se <- sd(pp$expression[pp$BR==r & pp$Module==i])
    p[j, ] <- c(r, i, mean, se)
    j = j+1
  }
}

p <- data.frame(p)
p$mean <- as.numeric(as.character(p$mean))
p$se <- as.numeric(as.character(p$se))
head(p)
str(p)
ggplot(p, aes(x=BR, y=mean, fill=Module)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = c(wes_palette("Royal2")[c(1,3,5)], 
                               wes_palette("GrandBudapest1"))) +
  facet_wrap(~Module,ncol=1, scales="free_y") +
  geom_errorbar((aes(ymin=mean-se, ymax=mean+se))) +
  guides(fill=F) +
  xlab("") +
  ylab("Module eigengene expression")

ggsave("revise0530/gene_count_res/Module_eigengene_expression_cortex_lobe.pdf", 
       width = 6, height = 5)


### signature genes expression
library(reshape2)
library(ggplot2)
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_cortex.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
genes <- c("PVALB", "OSTN", "FSTL1")
gene_info <- gtf_ensembl_gene[gtf_ensembl_gene$gene_name %in% genes, ]
df <- datExpr[,gene_info$gene_id]
colnames(df) <- gene_info$gene_name
head(df)

pp <- melt(t(df))
colnames(pp) <- c("gene", "sample_id", "expression")
pp$BR <- Rhesus_colData[as.character(pp$sample_id),"BR"]
pp$BR <- factor(pp$BR, levels = c("occipital_lobe", "parietal_lobe", "Limbic_cortex", "temporal_lobe", "frontal_lobe"))
head(pp)
str(pp)


ggplot(pp, aes(x=BR, y=log2(expression+1),fill=BR, colour=BR)) +
  # geom_violin(alpha=0.1) +
  # geom_jitter(width = 0.1, height = 0.3, size=0.01) +
  geom_boxplot(alpha=0.1) + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        line = element_line(colour = "black"),
        rect = element_rect(colour = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white", colour = "white")) +
  scale_fill_manual(values = anno_colors$BR)+
  scale_colour_manual(values = anno_colors$BR)+
  facet_wrap(~gene, nrow = 1) +
  # coord_flip() +
  # guides(fill=F, colour=F) +
  # guides(fill=guide_legend(title = NULL)) +
  xlab("") 


# ggsave("revise0530/gene_count_res/WGCNA_for_lobe_signature_genes_expression_violin.pdf", 
#        width = 5, height = 3)
ggsave("revise0530/gene_count_res/WGCNA_for_lobe_signature_genes_expression_boxplot.pdf", 
       width = 7, height = 3)
       
###  
p <- pp
colnames(p)[4] <- c("lobe")
write.table(p, "SourceData/Fig.3e.txt",
            col.names = T, row.names = F, quote=F, sep="\t")





############################################################################
#                                network                                   #
############################################################################
###导出网络到Cytoscape###
# Recalculate topological overlap if needed
enableWGCNAThreads()
library(WGCNA)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_cortex.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_TOM_cortex.RData")
load("revise0530/gene_count_res/WGCNA_modTraitCor_for_ID_cortex.RData")
load("revise0530/gene_count_res//WGCNA_Rhesus_net_modules_cortex.RData")
load("revise0530/gene_count_res/WGCNA_modTraitCor_for_lobe_cortex.RData")

# TOM = TOMsimilarityFromExpr(datExpr, power = 9);
# save(TOM, file="revise0530/gene_count_res/WGCNA_Rhesus_TOM_cortex.RData")


#####  hub gene
dim(datExpr) # 248 27694
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
moduleLabels2 <- paste("ME", net$colors, sep="")
hub_res <- chooseTopHubInEachModule(datExpr,
                                    colorh=moduleLabels2,
                                    power=9, 
                                    type = "signed")
hub_res <- data.frame(hub_gene=hub_res)
hub_res$gene_name <- gtf_ensembl_gene[hub_res$hub_gene,2]
head(hub_res)


# Select modules，选择需要导出的模块颜色
region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)


# Select module probes选择模块探测
probes = names(datExpr)
inModule = is.finite(match(moduleLabels2, region_modules));
modProbes = probes[inModule];
length(modProbes)

mean_expression <- sapply(unique(datTraits$ID), function(x) {colMeans(datExpr[which(datTraits$ID==x),])})
mean_expression <- sapply(unique(datTraits$BR), function(x) {colMeans(datExpr[which(datTraits$BR==x),])})

for(module in region_modules){
  region <- names(region_modules)[which(region_modules==module)]
  inModule = is.finite(match(moduleLabels2, module));
  modProbes = probes[inModule];
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  kIM = intramodularConnectivity(modTOM, moduleColors[inModule], scaleByMax = TRUE)
  kIM <- kIM[order(kIM$kWithin, decreasing=T),1:2]
  kIM <- data.frame(gene_name=gtf_ensembl_gene[rownames(kIM),5], kIM)
  kIM$module <- module
  kIM <- data.frame(kIM, mean_expression[rownames(kIM),])
  kIM <- data.frame(gene_id=rownames(kIM), kIM)
  filename=paste("revise0530/gene_count_res/module_", module, "_", region, "_info.txt", sep = "")
  write.table(kIM, filename, row.names = F, col.names = T, quote = F, sep = "\t")
}


# # For modular membership:
# MM = as.data.frame(cor(datExpr, MEs, use = "p"));

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("revise0530/gene_count_res/Rhesus-CytoscapeInput-edges-", paste(region,region_modules, sep="-"), ".txt", sep=""),
                               nodeFile = paste("revise0530/gene_count_res/Rhesus-CytoscapeInput-nodes-", paste(region,region_modules, sep="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);

# # For intramodular connectivity:
kIM = intramodularConnectivity(modTOM, moduleColors[inModule], scaleByMax = TRUE)
edge <- cyt$edgeData
dim(cyt$nodeData)
degree <- table(c(as.character(edge$fromNode),as.character(edge$toNode)))
degree <- data.frame(degree)
colnames(degree) <- c("nodeName", "degree")
eg = bitr(degree$nodeName, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mmu.eg.db")
degree$AltName <- as.character(degree$nodeName)
for(i in degree$AltName){if(i %in% eg$ENSEMBL){degree$AltName[degree$AltName==i] <- eg$SYMBOL[eg$ENSEMBL==i]}}
degree$KTotal <- kIM[as.character(degree$nodeName), "kTotal"]
degree <- degree[order(degree$KTotal, decreasing=T),]
head(degree)
write.table(degree,
            paste("revise0530/gene_count_res/Rhesus-CytoscapeInput-degree-", paste(region,region_modules, sep="-"), ".txt", sep=""),
            row.names = F,
            col.names = T,
            quote = F,
            sep="\t"
)





#############################################################
#                        function                           #
#############################################################
library(gProfileR)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/WGCNA_Rhesus_net_modules_cortex.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_cortex.RData")
load("revise0530/gene_count_res/WGCNA_modTraitCor_for_ID_cortex.RData")
load("revise0530/gene_count_res//WGCNA_Rhesus_net_modules_cortex.RData")
load("revise0530/gene_count_res/WGCNA_modTraitCor_for_lobe_cortex.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_TOM_cortex.RData")


moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
moduleLabels2 <- paste("ME", net$colors, sep="")
probes <- colnames(datExpr)

region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)


for(module in region_modules){
  region <- names(region_modules)[which(region_modules==module)]
  inModule = is.finite(match(moduleLabels2, module));
  modProbes = probes[inModule];
  gprofiler_res <- gprofiler(query=modProbes, organism = "mmulatta")
  gprofiler_res$gene_name <- gprofiler_res$intersection
  if(dim(gprofiler_res)[1]==0) next
  else{for(i in 1:dim(gprofiler_res)[1]){
    x <- gprofiler_res$gene_name[i]
    gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],5], collapse=",")
  }
    filename=paste("revise0530/gene_count_res/WGCNA_", module,"_", region, "_function.txt", sep="")
    write.table(gprofiler_res, filename,
                quote = F, sep="\t", row.names = F, col.names = T)
  }
}



############################################################################
#                          保存 module trans_id                            #
############################################################################
library(WGCNA)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_cortex.RData")
load("revise0530/gene_count_res//WGCNA_Rhesus_net_modules_cortex.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_TOM_cortex.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")

############
load("revise0530/gene_count_res/WGCNA_modTraitCor_for_ID_cortex.RData")
dim(datExpr) # 248 21306
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
moduleLabels2 <- paste("ME", net$colors, sep="")

region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)
probes = colnames(datExpr)

for(module in region_modules){
  region <- names(region_modules)[which(region_modules==module)]
  inModule = is.finite(match(moduleLabels2, module));
  modProbes = probes[inModule];
  modProbes <- data.frame(gene_id=modProbes, 
                          gene_name=gtf_ensembl_gene[modProbes, "gene_name"])
  write.table(modProbes, paste("revise0530/gene_count_res/", module, "_", region, "_gene_id.txt", sep=""),
              quote = F, sep="\t", row.names = F, col.names = F)
}

################
load("revise0530/gene_count_res/WGCNA_modTraitCor_for_lobe_cortex.RData")
dim(datExpr) # 248 21306
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
moduleLabels2 <- paste("ME", net$colors, sep="")

region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)
probes = colnames(datExpr)

for(module in region_modules){
  region <- names(region_modules)[which(region_modules==module)]
  inModule = is.finite(match(moduleLabels2, module));
  modProbes = probes[inModule];
  modProbes <- data.frame(gene_id=modProbes, 
                          gene_name=gtf_ensembl_gene[modProbes, "gene_name"])
  write.table(modProbes, paste("revise0530/gene_count_res/", module, "_", region, "_gene_id.txt", sep=""),
              quote = F, sep="\t", row.names = F, col.names = F)
}



for(module in region_modules){
  region <- names(region_modules)[which(region_modules==module)]
  inModule = is.finite(match(moduleLabels2, module));
  modProbes = probes[inModule];
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  kIM = intramodularConnectivity(modTOM, moduleColors[inModule], scaleByMax = TRUE)
  kIM <- kIM[order(kIM$kWithin, decreasing=T),1:2]
  kIM <- data.frame(gene_name=gtf_ensembl_gene[rownames(kIM),"gene_name"], kIM)
  kIM$module <- module
  kIM <- data.frame(gene_id=rownames(kIM), kIM)
  filename=paste("revise0530/gene_count_res/module_", module, "_", region, "_gene_id_info.txt", sep = "")
  write.table(kIM, filename, row.names = F, col.names = T, quote = F, sep = "\t")
}

