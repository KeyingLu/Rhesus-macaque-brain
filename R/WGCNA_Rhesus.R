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
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData") 

Rhesus_Gene_norCounts <- counts(Rhesus_GeneObject,normalized=TRUE)
dim(Rhesus_Gene_norCounts) # 23651   408
#=====================================================================================
#  Code chunk 2
#=====================================================================================
datExpr0=as.data.frame(t(Rhesus_Gene_norCounts))
dim(datExpr0) # 408 23651
#=====================================================================================
#  Code chunk 3
#=====================================================================================
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
dim(datExpr0)  ## 407 23651
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
nGenes = ncol(datExpr) # 
nSamples = nrow(datExpr) # 408
save(datExpr, 
     file="revise0530/gene_count_res/WGCNA_Rhesus_datExpr_exclude_MB.RData")
#=====================================================================================
#  Code chunk 7
#=====================================================================================
load("revise0530/gene_count_res/Rhesus_colData.RData") # Rhesus_colData
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_exclude_MB.RData")
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
pdf("revise0530/gene_count_res/WGCNA_Rhesus_soft_threshold.pdf",
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
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_exclude_MB.RData")
cor <- WGCNA::cor # aviod a conflict between the WGCNA and the other packages
net = blockwiseModules(datExpr, power = 6,
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
save(net, file="revise0530/gene_count_res/WGCNA_Rhesus_net_modules_exclude_MB.RData")

######### TOM
enableWGCNAThreads()
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
save(TOM, file="revise0530/gene_count_res/WGCNA_Rhesus_TOM_exclude_MB.RData")

####
load("revise0530/gene_count_res/WGCNA_Rhesus_net_modules_exclude_MB.RData")
table(net$colors) ##  59 modules
# net$dendrograms[[1]]


#=====================================================================================
#  Code chunk 4
#=====================================================================================
# library(WGCNA)
# setwd("/mnt/data2/Rhesus_brain")
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
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_exclude_MB.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_net_modules_exclude_MB.RData")


moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs = net$MEs; 
dim(MEs) # 408  60
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
dim(MEs0) # [1] 408  60
MEsWW = orderMEs(MEs0)
colnames(MEsWW) <- gsub("E", "", colnames(MEsWW))
modTraitCor = cor(MEsWW, SR_merge_traits, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nrow(datExpr))
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
save(modTraitCor, modTraitP, 
     file="revise0530/gene_count_res/WGCNA_modTraitCor_for_SR_merge_exclude_MB.RData")


##### plot cor
library(pheatmap)
library(RColorBrewer)
pheatmap(
  modTraitCor,
  main="Module-Region_trait relationships",
  filename="revise0530/gene_count_res/Net_modules_and_SR_merge_exMB_traits_cor_pheatmap.pdf",
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
modTraitCor_help <- modTraitCor
sum(modTraitCor>=0.80 & modTraitP<0.01)  # 11
region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)



#### bar plots showing the expression of MEs across brain
library(ggplot2)
library(reshape2)
library(wesanderson) # names(wesanderson)
load("revise0530/gene_count_res/WGCNA_modTraitCor_for_SR_merge_exclude_MB.RData")
region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)

pp <- data.frame(MEsWW[, region_modules], SR=datTraits$SR_merge)
pp <- melt(pp)
colnames(pp)[2:3] <- c("Module", "expression")
head(pp)

p <- matrix(0, ncol=4, nrow=length(unique(pp$SR))*length(region_modules))
colnames(p) <- c("SR", "Module", "mean", "se")
j=1
for(r in unique(pp$SR)){
  for(i in unique(pp$Module)){
    mean <- mean(pp$expression[pp$SR==r & pp$Module==i])
    se <- sd(pp$expression[pp$SR==r & pp$Module==i])
    p[j, ] <- c(r, i, mean, se)
    j = j+1
  }
}


p <- data.frame(p)
p$mean <- as.numeric(as.character(p$mean))
p$se <- as.numeric(as.character(p$se))
head(p)
str(p)
ggplot(p, aes(x=SR, y=mean, fill=Module)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        text = element_text(face="bold"),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white")) +
  scale_fill_manual(values = c(wes_palette("Royal2"), 
                               wes_palette("GrandBudapest1"),
                               wes_palette("Chevalier1")[3:4])) +
  facet_wrap(~Module, ncol=4, scales="free_y") +
  geom_errorbar((aes(ymin=mean-se, ymax=mean+se))) +
  guides(fill=F) +
  xlab("") +
  ylab("Module eigengene expression")

ggsave("revise0530/gene_count_res/Module_eigengene_expression_exclude_MB.pdf", 
       width = 12, height = 8)


#### bar plots showing the expression of MEs across brain (all)
library(reshape2)
library(wesanderson)
pp <- data.frame(MEsWW, SR=datTraits$SR_merge)
pp <- melt(pp)
head(pp)
colnames(pp)[2:3] <- c("Module", "expression")
p <- matrix(0, ncol=4, nrow=length(levels(pp$SR))*50)
colnames(p) <- c("SR", "Module", "mean", "se")
j=1
for(r in unique(pp$SR)){
  for(i in unique(pp$Module)){
    mean <- mean(pp$expression[pp$SR==r & pp$Module==i])
    se <- sd(pp$expression[pp$SR==r & pp$Module==i])
    p[j, ] <- c(r, i, mean, se)
    j = j+1
  }
}

p <- data.frame(p)
p$mean <- as.numeric(as.character(p$mean))
p$se <- as.numeric(as.character(p$se))
head(p)
str(p)

ggplot(p, aes(x=SR, y=mean)) +
  geom_bar(stat="identity", fill="grey") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        text = element_text(face="bold", size=6),
        strip.text = element_text(size=rel(1.5)),
        strip.background = element_rect(fill="white",colour = "white")) +
  facet_wrap(~Module,ncol=8, scales="free_y") +
  geom_errorbar((aes(ymin=mean-se, ymax=mean+se))) +
  guides(fill=F) +
  xlab("") +
  ylab("Module eigengene expression")

ggsave("revise0530/gene_count_res/ALL_Module_eigengene_expression_exclude_MB.pdf", width = 16, height = 12)



##########################################################
#                      network                           #
##########################################################
enableWGCNAThreads()
library(WGCNA)
setwd("/mnt/data2/Rhesus_brain/")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_exclude_MB.RData")
load("revise0530/gene_count_res/WGCNA_modTraitCor_for_SR_merge_exclude_MB.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_net_modules_exclude_MB.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData") # Rhesus_colData
load("revise0530/gene_count_res/WGCNA_Rhesus_TOM_exclude_MB.RData")
datTraits <- Rhesus_colData[rownames(datExpr),]

#####  hub gene
dim(datExpr) # 408 23651
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
moduleLabels2 <- paste("M", net$colors, sep="")
# hub_res <- chooseTopHubInEachModule(datExpr,
#                                     colorh=moduleLabels2,
#                                     power=6, 
#                                     type = "signed")
# hub_res <- data.frame(hub_gene=hub_res)
# hub_res$gene_name <- gtf_ensembl_gene[hub_res$hub_gene,"gene_name"]
# head(hub_res)



region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)
probes = names(datExpr)

mean_expression <- sapply(unique(datTraits$SR_merge), function(x) {colMeans(datExpr[which(datTraits$SR_merge==x),])})

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
  kIM <- data.frame(kIM, mean_expression[rownames(kIM),])
  kIM <- data.frame(gene_id=rownames(kIM), kIM)
  filename=paste("revise0530/gene_count_res/module_", module, "_", region, "_info.txt", sep = "")
  write.table(kIM, filename, row.names = F, col.names = T, quote = F, sep = "\t")
}



# Select modules，选择需要导出的模块颜色
probes = names(datExpr)
region = "OB" # HTHA  CB PIT OB
module <- rownames(modTraitCor)[modTraitCor[,region]>=0.8 & modTraitP[,region]<0.01]
module
inModule = is.finite(match(moduleLabels2, module));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
dim(modTOM)

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("revise0530/gene_count_res/Rhesus-CytoscapeInput-edges-", paste(region,paste(module,collapse = "-"), sep="-"), ".txt", sep=""),
                               nodeFile = paste("revise0530/gene_count_res/Rhesus-CytoscapeInput-nodes-", paste(region,paste(module,collapse = "-"), sep="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.25,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);
edge <- cyt$edgeData
dim(cyt$nodeData)
degree <- table(c(as.character(edge$fromNode),as.character(edge$toNode)))
degree <- data.frame(degree)
colnames(degree) <- c("nodeName", "degree")

kIM = intramodularConnectivity(modTOM, moduleColors[inModule], scaleByMax = TRUE)
kIM <- kIM[order(kIM$kWithin, decreasing=T),1:2]
kIM <- data.frame(gene_name=gtf_ensembl_gene[rownames(kIM),"gene_name"], kIM)
degree <- data.frame(degree, kIM[as.character(degree$nodeName), c("kTotal", "gene_name")])
degree <- degree[order(degree$kTotal, decreasing=T),]
degree$gene_name <- as.character(degree$gene_name)
degree$gene_name[which(is.na(degree$gene_name))] <- as.character(degree$nodeName)[which(is.na(degree$gene_name))]
head(degree)

write.table(degree,
            paste("revise0530/gene_count_res/Rhesus-CytoscapeInput-degree-", paste(region,paste(module,collapse = "-"), sep="-"), ".txt", sep=""),
            row.names = F,
            col.names = T,
            quote = F,
            sep="\t"
)



###################################################
#                     function                    #
###################################################
library(gProfileR)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
# set_base_url("http://biit.cs.ut.ee/gprofiler")
library(WGCNA)
setwd("/mnt/data2/Rhesus_brain/")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_exclude_MB.RData")
load("revise0530/gene_count_res/WGCNA_modTraitCor_for_SR_merge_exclude_MB.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_net_modules_exclude_MB.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_TOM_exclude_MB.RData")


moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
moduleLabels2 <- paste("M", net$colors, sep="")
probes <- colnames(datExpr)

region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)


module_function_list <- list()
for(module in region_modules){
  region <- names(region_modules)[which(region_modules==module)]
  inModule = is.finite(match(moduleLabels2, module));
  modProbes = probes[inModule];
  gprofiler_res <- gprofiler(query=modProbes, organism = "mmulatta")
  gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),]
  gprofiler_res$gene_name <- gprofiler_res$intersection
  if(dim(gprofiler_res)[1]==0) next
  else{for(i in 1:dim(gprofiler_res)[1]){
    x <- gprofiler_res$gene_name[i]
    gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],"gene_name"], collapse=",")
  }
    module_function_list[[module]] <- gprofiler_res
    filename=paste("revise0530/gene_count_res/WGCNA_", module,"_", region, "_function.txt", sep="")
    write.table(gprofiler_res, filename,
                quote = F, sep="\t", row.names = F, col.names = T)
  }
}

save(module_function_list, file="revise0530/gene_count_res/WGCNA_module_funciton.RData")



library(ggplot2)
load("revise0530/gene_count_res/WGCNA_module_funciton.RData")
#### term  plot
names(module_function_list)
M7_res <- module_function_list[["M7"]][, c("p.value", "domain", "term.name")]
rownames(M7_res) <- NULL
# View(M7_res)
M7_res_sub <- M7_res[1, ]
M7_res_sub <- data.frame(M7_res_sub, module="M7")


M28_res <- module_function_list[["M28"]][, c("p.value", "domain", "term.name")]
rownames(M28_res) <- NULL
M28_res_sub <- M28_res[c(1:5), ]
M28_res_sub <- data.frame(M28_res_sub, module="M28")


M45_res <- module_function_list[["M45"]][, c("p.value", "domain", "term.name")]
rownames(M45_res) <- NULL
M45_res_sub <- M45_res[c(1:5), ]
M45_res_sub <- data.frame(M45_res_sub, module="M45")


M2_res <- module_function_list[["M2"]][, c("p.value", "domain", "term.name")]
rownames(M2_res) <- NULL
M2_res_sub <- M2_res[c(1:5), ]
M2_res_sub <- data.frame(M2_res_sub, module="M2")


M36_res <- module_function_list[["M36"]][, c("p.value", "domain", "term.name")]
rownames(M36_res) <- NULL
M36_res_sub <- M36_res[c(2:5, 7), ]
M36_res_sub <- data.frame(M36_res_sub, module="M36")


M11_res <- module_function_list[["M11"]][, c("p.value", "domain", "term.name")]
rownames(M11_res) <- NULL
M11_res_sub <- M11_res[c(3,9:12), ]
M11_res_sub <- data.frame(M11_res_sub, module="M11")


M6_res <- module_function_list[["M6"]][, c("p.value", "domain", "term.name")]
rownames(M6_res) <- NULL
M6_res_sub <- M6_res[c(4, 5, 6, 7,12), ]
M6_res_sub <- data.frame(M6_res_sub, module="M6")


M25_res <- module_function_list[["M25"]][, c("p.value", "domain", "term.name")]
rownames(M25_res) <- NULL
M25_res_sub <- M25_res[c(1:5), ]
M25_res_sub <- data.frame(M25_res_sub, module="M25")




pp <- rbind(M7_res_sub,M28_res_sub,M45_res_sub,
            M2_res_sub,M36_res_sub,
            M11_res_sub,M6_res_sub,M25_res_sub)

pp$term.name <- factor(pp$term.name, levels=c(rev(unique(pp$term.name))))


ggplot(pp, aes(x=term.name, y=-log10(p.value))) +
  geom_bar(stat="identity", fill="#C1AFD5", width = 0.8) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        # legend.position = c(0.8,0.8),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  # geom_text(aes(label=symbol), vjust=0.5, hjust=0,size=4, colour="black") +
  # annotate("text", x=pp$NAME[2], y=-1.5,size=3,
  #          label="nominal P values: \n*P < 0.05 \n**P < 0.005 \n***P < 0.0005 \nNS: not significant") +
  # scale_x_discrete(limits=c(pp$term.name))+
  coord_flip() +
  xlab("") +
  facet_grid(module~., scales = "free", space = "free")



ggsave("revise0530/gene_count_res/WGCNA_for_region_module_main_function_bar.pdf",
       height = 10, width = 6)
##




############################################################################
#                         selected genes expression                        #
############################################################################
library(ggplot2)
library(reshape2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain/")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_exclude_MB.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/module_genes_order_list.RData")


genes <- c("CDK5R1", "CDK5", "LMTK2", "SP8", "SP9", "TBX21", 
           "CBLN1", "CBLN3", "FAT2", "AVP", "OXT", "GNRH", "FEZF1",
           "GH1", "GHRHR", "CGA", "POU1F1", "ADCY5", "PPP1R1B",
           "CA12", "COCH", "PDE1B")


# res <- module_genes_order_list[["M9_cortex1"]]
res <- module_genes_order_list[["M10_cortex2"]]

res <- res[complete.cases(res), ]
genes <- as.character(res$gene_name[1:10])

genes[which(genes %in% gtf_ensembl_gene$gene_name)]
genes <- genes[which(genes %in% gtf_ensembl_gene$gene_name)]
gene_info <- gtf_ensembl_gene[gtf_ensembl_gene$gene_name %in% genes, ]


# vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
# vsd_data <- assay(vsd)
pp <- t(datExpr)[gene_info$gene_id,]
# pp <- vsd_data[gene_info$gene_id,]
rownames(pp) <- gene_info$gene_name
pp <- pp[genes,]
pp <- melt(pp)
colnames(pp) <- c("gene", "sample_id", "expression")
head(pp)
pp$SR <- Rhesus_colData[as.character(pp$sample_id), ]$SR_merge
# pp$gene <- as.character(pp$gene)
# pp$SR <- as.character(pp$SR)
# sums <- aggregate(expression~gene, data=pp, sum)
# sums <- sums[order(sums$expression, decreasing = T),]
# pp$gene <- factor(pp$gene, levels = sums$gene)
# means <- aggregate(expression~SR,data=pp, mean)
# means <- means[order(means$expression, decreasing = T),]
# pp$SR <- factor(pp$SR, levels = means$SR)
# 

order = rev(c("cortex", "OB" , "CB" , "HTHA", "PIT", "STR","HIP",
              "AMY", "PA", "CC", "THA", "MO", "PON", "VTA",
              "SN" , "GP"))

ggplot(pp, aes(x=SR, y=expression, fill=SR, color=SR)) +
  geom_violin() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill="white", color = "white")) +
  facet_grid(.~gene, scales="free_x") +
  coord_flip() +
  guides(fill=F, color=F) +
  xlab("") +
  # ylab("expression(VST)") +
  # scale_x_discrete(limits=rev(names())) +
  scale_x_discrete(limits=order) +
  scale_fill_manual(values = setNames(colorRampPalette(brewer.pal(12,"Paired"))(16), order)) +
  scale_color_manual(values = setNames(colorRampPalette(brewer.pal(12,"Paired"))(16), order)) +
  ggtitle("")


ggsave("revise0530/gene_count_res/WGCNA_module_selected_genes_violin_M10_top10.pdf",
       width = 8, height = 4)



############################################################################
#                           hub genes (top)                                #
############################################################################

enableWGCNAThreads()
library(DESeq2)
library(WGCNA)
library(pheatmap)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain/")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_exclude_MB.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_TOM_exclude_MB.RData")
load("revise0530/gene_count_res/WGCNA_modTraitCor_for_SR_merge_exclude_MB.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_net_modules_exclude_MB.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")


vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)

#####  hub gene
dim(datExpr) # 408 23651
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
moduleLabels2 <- paste("M", net$colors, sep="")

region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)
probes = names(datExpr)

module_genes_order_list <- list()
for(module in region_modules){
  region <- names(region_modules)[which(region_modules==module)]
  inModule = is.finite(match(moduleLabels2, module));
  modProbes = probes[inModule];
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  kIM = intramodularConnectivity(modTOM, moduleColors[inModule], scaleByMax = TRUE)
  kIM <- kIM[order(kIM$kWithin, decreasing=T),1:2]
  kIM <- data.frame(gene_name=gtf_ensembl_gene[rownames(kIM),"gene_name"], kIM)
  names <- paste(module, region, sep="_")
  module_genes_order_list[[names]] <- kIM
}

save(module_genes_order_list, 
     file="revise0530/gene_count_res/module_genes_order_list.RData")


##
names(module_genes_order_list)

##
M7_CB <- module_genes_order_list[["M7_CB"]]
M7_CB <- M7_CB[complete.cases(M7_CB),]
M7_CB <- data.frame(M7_CB, module="M7_CB")
M7_CB <- M7_CB[1:8,]
M7_CB
##
M28_THA <- module_genes_order_list[["M28_THA"]]
M28_THA <- M28_THA[complete.cases(M28_THA),]
M28_THA <- data.frame(M28_THA, module="M28_THA")
M28_THA <- M28_THA[c(1:8),]
M28_THA
##
M45_HTHA <- module_genes_order_list[["M45_HTHA"]]
M45_HTHA <- M45_HTHA[complete.cases(M45_HTHA),]
M45_HTHA <- data.frame(M45_HTHA, module="M45_HTHA")
M45_HTHA <- M45_HTHA[c(1:8),]
M45_HTHA
##
M2_PIT <- module_genes_order_list[["M2_PIT"]]
M2_PIT <- M2_PIT[complete.cases(M2_PIT),]
M2_PIT <- data.frame(M2_PIT, module="M2_PIT")
M2_PIT <- M2_PIT[c(1:8),]
M2_PIT
##
# M23_PON <- module_genes_order_list[["M23_PON"]]
# M23_PON <- M23_PON[complete.cases(M23_PON),]
# M23_PON <- data.frame(M23_PON, module="M23_PON")
# M23_PON <- M23_PON[c(1:8),]
# M23_PON
##
M36_SN <- module_genes_order_list[["M36_SN"]]
M36_SN <- M36_SN[complete.cases(M36_SN),]
M36_SN <- data.frame(M36_SN, module="M36_SN")
M36_SN <- M36_SN[c(1:8),]
M36_SN
##
M11_OB <- module_genes_order_list[["M11_OB"]]
M11_OB <- M11_OB[complete.cases(M11_OB),]
M11_OB <- data.frame(M11_OB, module="M11_OB")
M11_OB <- M11_OB[c(1:8),]
M11_OB
##
M6_STR <- module_genes_order_list[["M6_STR"]]
M6_STR <- M6_STR[complete.cases(M6_STR),]
M6_STR <- data.frame(M6_STR, module="M6_STR")
M6_STR <- M6_STR[c(1:8),]
M6_STR
##
M25_HIP <- module_genes_order_list[["M25_HIP"]]
M25_HIP <- M25_HIP[complete.cases(M25_HIP),]
M25_HIP <- data.frame(M25_HIP, module="M25_HIP")
M25_HIP <- M25_HIP[c(1:8),]
M25_HIP
##
# ME6_cortex <- module_genes_order_list[["ME6_cortex1"]]
# ME6_cortex <- ME6_cortex[complete.cases(ME6_cortex),]
# ME6_cortex <- data.frame(ME6_cortex, module="ME6_cortex")
# ME6_cortex <- ME6_cortex[c(1:8),]
# ME6_cortex
# ##
# ME10_cortex <- module_genes_order_list[["ME10_cortex2"]]
# ME10_cortex <- ME10_cortex[complete.cases(ME10_cortex),]
# ME10_cortex <- data.frame(ME10_cortex, module="ME10_cortex")
# ME10_cortex <- ME10_cortex[c(1:8),]
# ME10_cortex


##
top_genes <- rbind(M7_CB, M28_THA, M45_HTHA, M2_PIT, M36_SN,
                   M11_OB, M6_STR, M25_HIP)
region <- unique(gsub("M[0-9]+_", "", top_genes$module))
colData <- Rhesus_colData[rownames(datExpr),]
samples <- c()
for(r in region){s <- rownames(colData)[colData$SR_merge == r];samples <- c(samples, s)}
df <- t(datExpr)[rownames(top_genes), samples]
rownames(df) <- top_genes$gene_name



dfheat <- vsd_data[rownames(top_genes), samples]
rownames(dfheat) <- top_genes$gene_name

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


pheatmap(
  scale_df,
  col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
  breaks = c(-4,seq(-3,3,length=20),4),
  # legend_breaks = c(-2,-1,1,2),
  filename = "revise0530/gene_count_res/WGCNA_modules_top_genes_heatmap.pdf",
  # scale="row",
  # main = title,
  width = 10,
  height = 10,
  # clustering_method = "complete",
  border_color=NA,
  gaps_row = cumsum(table(top_genes$module)),
  # fontsize_col = 0.5,
  # annotation_row = geneCluster,
  # annotation_col = colInfo,
  # annotation_colors = ann_colors,
  # annotation_legend = T,
  show_rownames = T,
  show_colnames = F,
  cluster_rows = F,
  cluster_cols= F
)





############################################################################
#                          保存 module trans_id                            #
############################################################################
library(WGCNA)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_exclude_MB.RData")
load("revise0530/gene_count_res/WGCNA_modTraitCor_for_SR_merge_exclude_MB.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_net_modules_exclude_MB.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")


dim(datExpr) # 408 28568
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
moduleLabels2 <- paste("M", net$colors, sep="")

region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)
probes = names(datExpr)

for(module in region_modules){
  region <- names(region_modules)[which(region_modules==module)]
  inModule = is.finite(match(moduleLabels2, module));
  modProbes = probes[inModule];
  modProbes <- data.frame(gene_id=modProbes, 
                          gene_name=gtf_ensembl_gene[modProbes, "gene_name"])
  write.table(modProbes, paste("revise0530/gene_count_res/", module, "_", region, "_gene_id.txt", sep=""),
              quote = F, sep="\t", row.names = F, col.names = F)
  
}



