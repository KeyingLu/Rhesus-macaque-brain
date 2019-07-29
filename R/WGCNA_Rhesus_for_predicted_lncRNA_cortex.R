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
datExpr0=as.data.frame(t(Rhesus_norCounts[, rownames(Rhesus_colData)[Rhesus_colData$SR_cortex == "cortex"]]))
datExpr0=datExpr0[,colnames(datExpr0) %in% rownames(novel_lncRNA_trans)]
dim(datExpr0) # 248 2708
#=====================================================================================
#  Code chunk 3
#=====================================================================================
gsg = goodSamplesGenes(datExpr0, verbose = 3)
# Flagging genes and samples with too many missing values...
# ..step 1
# ..Excluding 63 genes from the calculation due to too many missing samples or zero variance.
# ..step 2
gsg$allOK # FALSE
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
dim(datExpr0)  ##  248 2678
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
nGenes = ncol(datExpr) # 2678
nSamples = nrow(datExpr) # 248
save(datExpr, 
     file="revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA_cortex.RData")
#=====================================================================================
#  Code chunk 7
#=====================================================================================
load("revise0530/gene_count_res/Rhesus_colData.RData") # Rhesus_colData
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA_cortex.RData")
datTraits <- Rhesus_colData[rownames(datExpr),]
head(datTraits)
dim(datTraits) # [1] 248   9
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
pdf("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_lncRNA_cortex_soft_threshold.pdf",
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
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA_cortex.RData")
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = 7,
                       TOMType = "unsigned",
                       # minModuleSize = 10,
                       minModuleSize = 20,
                       # minModuleSize = 50,
                       reassignThreshold = 0, 
                       # mergeCutHeight = 0.15,
                       # deepSplit=4,
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE
)
save(net, file="revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_net_modules_lncRNA_cortex.RData")

######### TOM
enableWGCNAThreads()
TOM = TOMsimilarityFromExpr(datExpr, power = 7);
save(TOM, file="revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_TOM_lncRNA_cortex.RData")

####
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_net_modules_lncRNA_cortex.RData")
table(net$colors) ##   11 modules
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
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA_cortex.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_net_modules_lncRNA_cortex.RData")


moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs = net$MEs; 
dim(MEs) # 248  12
geneTree = net$dendrograms[[1]]



############################################################################
#                      cor between modules and traits                      #
############################################################################


##############
#### lobe(in cortex)
##############
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA_cortex.RData")
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
dim(MEs0) # 
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, ID_cortex_traits, use = "p", method = "pearson")
modTraitP = corPvalueStudent(modTraitCor, nrow(datExpr))
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
save(modTraitCor, modTraitP, 
     file="revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_modTraitCor_for_lobe_lncRNA_cortex.RData")

##### plot cor
library(pheatmap)
library(RColorBrewer)
pheatmap(
  t(modTraitCor),
  main="Module-lobe(cortex)_trait relationships",
  # filename="revise0530/transcripts_count_res_with_TACO_gtf/Net_modules_and_lobe_cortex_traits_cor_pheatmap.pdf",
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






################################
#### ID(in cortex)
################################
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA_cortex.RData")
datTraits <- Rhesus_colData[rownames(datExpr),]
Info <- datTraits$ID

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
dim(ID_cortex_traits) # 248   8

####
moduleColorsWW = moduleColors
MEs0 = moduleEigengenes(datExpr, moduleLabels)$eigengenes # MEs = net$MEs
dim(MEs0) # 
MEsWW = orderMEs(MEs0)
modTraitCor = cor(MEsWW, ID_cortex_traits, use = "p", method = "pearson")
modTraitP = corPvalueStudent(modTraitCor, nrow(datExpr))
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
save(modTraitCor, modTraitP, 
     file="revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_modTraitCor_for_ID_lncRNA_cortex.RData")

##### plot cor
library(pheatmap)
library(RColorBrewer)
pheatmap(
  modTraitCor,
  main="Module-ID(cortex)_trait relationships",
  filename="revise0530/transcripts_count_res_with_TACO_gtf/Net_modules_and_ID_cortex_traits_cor_pheatmap_lncRNA.pdf",
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
sum(modTraitCor>=0.80 & modTraitP<0.01) # 7
sum(modTraitCor <= -0.80 & modTraitP<0.01)
modTraitCor_help <- modTraitCor


#### bar plots showing the expression of MEs across brain
library(reshape2)
library(wesanderson)
region_modules <- apply(modTraitCor, 2, function(x){rownames(modTraitCor)[which(x>=0.8)]})
region_modules <- unlist(region_modules)

pp <- data.frame(MEsWW[, region_modules, drop=F], ID=datTraits$ID)
pp <- melt(pp)
head(pp)
colnames(pp)[2:3] <- c("Module", "expression")
pp$Module <- as.character(pp$Module)
p <- matrix(0, ncol=4, nrow=8*length(region_modules))
colnames(p) <- c("ID", "Module", "mean", "se")

j=1
for(r in unique(pp$ID)){
  for(i in unique(pp$Module)){
    mean <- mean(pp$expression[pp$ID==r & pp$Module==i])
    se <- sd(pp$expression[pp$ID==r & pp$Module==i])
    p[j, ] <- c(r, i, mean, se)
    j = j+1
  }
}

p <- data.frame(p)
p$mean <- as.numeric(as.character(p$mean))
p$se <- as.numeric(as.character(p$se))
head(p)
str(p)
ggplot(p, aes(x=ID, y=mean, fill=Module)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = c(wes_palette("Royal2")[c(1,3,5)], wes_palette("GrandBudapest1"))) +
  facet_wrap(~Module,ncol=3, scales="free_y") +
  geom_errorbar((aes(ymin=mean-se, ymax=mean+se))) +
  guides(fill=F) +
  xlab("") +
  ylab("Module eigengene expression")

ggsave("revise0530/transcripts_count_res_with_TACO_gtf/Module_eigengene_expression_ID_lncRNA_cortex.pdf", 
       width = 6, height = 4)



########### 保存 module trans_id
library(WGCNA)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA_cortex.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_net_modules_lncRNA_cortex.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_modTraitCor_for_ID_lncRNA_cortex.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_TOM_lncRNA_cortex.RData")
dim(datExpr) # 248 2782
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
moduleLabels2 <- paste("Lnc-MB", net$colors, sep="")

rownames(modTraitCor) <- gsub("ME", "Lnc-MB", rownames(modTraitCor))
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



############################################################################
#           correlation between lncRNA module and RNA                      #
############################################################################
library(WGCNA)
library(gProfileR)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_datExpr_lncRNA_cortex.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/WGCNA_Rhesus_net_modules_lncRNA_cortex.RData")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


##
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs = net$MEs; 
dim(MEs) # 248  12
rownames(MEs) <- rownames(datExpr)

##
norCounts <- counts(Rhesus_GeneObject, normalized=T)
sub_norCounts <- norCounts[, rownames(Rhesus_colData)[Rhesus_colData$SR_merge=="cortex"]]

cor_res <- cor(t(sub_norCounts), MEs)
cor_res <- cor_res[, c("ME10", "ME3", "ME4", "ME5", "ME6", "ME7", "ME8")] 
  
##
pp <- melt(cor_res)
pp <- pp[complete.cases(pp),]
colnames(pp) <- c("gene", "module", "correlation")
pp_p <- pp[pp$correlation>0.8,]
dim(pp_p)
pp_n <- pp[pp$correlation< -0.8,]
dim(pp_n)

pp_pn <- pp[abs(pp$correlation) > 0.8,]  
dim(pp_pn)
  
  
lncRNA_module_function_list <- list()
for(module in levels(pp_pn$module)){
  gprofiler_res <- gprofiler(
    query=as.character(pp_pn[pp_pn$module==module, "gene"]), 
    organism = "mmulatta")
  gprofiler_res <- gprofiler_res[order(gprofiler_res[,3]),]
  gprofiler_res$gene_name <- gprofiler_res$intersection
  if(dim(gprofiler_res)[1]==0) next
  else{for(i in 1:dim(gprofiler_res)[1]){
    x <- gprofiler_res$gene_name[i]
    gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],"gene_name"], collapse=",")
  }
    # re <- paste("cluster", cluster, sep="")
    lncRNA_module_function_list[[module]] <- gprofiler_res
    # filename=paste("revise0530/gene_count_res/LRT_ID_geneCluster", cluster,"_function.txt", sep="")
    # write.table(gprofiler_res, filename,
    #             quote = F, sep="\t", row.names = F, col.names = T)
  }
}


lapply(lncRNA_module_function_list, function(x){x[,"term.name"]})


