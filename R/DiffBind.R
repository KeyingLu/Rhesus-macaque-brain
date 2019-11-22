# /mnt/data1/Tools/R-3.6.0/bin/R
library(DiffBind)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix1020/score5_res/Rhesus_peak_colData.RData")



###########################################################
#                       sampleSheet                       #
###########################################################
sampleSheet <- Rhesus_peak_colData[, c("ID", "SR", "SR_merge", "Age_Stage", "Sex")]
sampleSheet <- sampleSheet[order(sampleSheet$SR),]
sampleSheet$Tissue <- sampleSheet$SR
sampleSheet$Condition <- sampleSheet$SR_merge
sampleSheet$SampleID <- gsub("X", "", rownames(sampleSheet))
sampleSheet$bamReads <- paste("ATAC/bam/", sampleSheet$SampleID, "_selected.bam", sep="")
# sampleSheet$Peaks <- paste("ATAC/peak/", sampleSheet$SampleID, "_normalScore_peaks.bed", sep="")  # PeakCaller = "raw"
sampleSheet$Peaks <- paste("ATAC/Genrich_res/", sampleSheet$SR, "_a200.narrowPeak", sep="")  # PeakCaller = "narrowPeak"
# sampleSheet$PeakCaller <- "raw"
sampleSheet$PeakCaller <- "narrow"
write.csv(sampleSheet, "ATAC/matrix1020/score5_res/sampleSheet.csv")


######
sub_sampleSheet <- sampleSheet[sampleSheet$SR %in% c("CA1", "PCG"), ]
write.csv(sub_sampleSheet, "ATAC/matrix1020/score5_res/CA1_and_PCG_sampleSheet.csv")






###########################################################
#                        33 samples                       #
###########################################################
library(ChIPseeker)
library(GenomicFeatures)
setwd("/mnt/data2/Rhesus_brain")
txdb <- makeTxDbFromGFF("/mnt/data1/Ref/Rhesus_macaque_8.0.1/Macaca_mulatta.Mmul_8.0.1.91.chr.gtf",
                        format = "gtf",
                        dataSource = "ensemblgenomes",
                        organism = "Macaca mulatta")

promoter <- getPromoters(TxDb=txdb, upstream=2500, downstream=2500)


#####
RegionPeak <- dba(sampleSheet="ATAC/matrix1020/score5_res/sampleSheet.csv")
RegionPeak <- dba.count(RegionPeak ,bParallel=T)
lapply(RegionPeak$peaks, dim)
save(RegionPeak, file="ATAC/matrix1020/score5_res/RegionPeak_33.RData")

###
called <- RegionPeak$called
colnames(called) <- colnames(RegionPeak$binding)[4:36]
allPeaks <- RegionPeak$peaks[[1]][,1:3]
rownames(allPeaks) <- paste(allPeaks[,1], ":", allPeaks[,2], "-", allPeaks[,3],sep="")
rownames(called) <- paste(allPeaks[,1], ":", allPeaks[,2], "-", allPeaks[,3],sep="")
colData <- RegionPeak$samples
rownames(colData) <- colData$SampleID

region_occur <- called[,cumsum(table(colData$SR))]
colnames(region_occur) <- names(cumsum(table(colData$SR)))
save(region_occur, file="ATAC/matrix1020/score5_res/region_occur.RData")


num <- apply(region_occur,1,function(x){sum(x>=1)})

sharePeaksAnno <- list()
for(i in 1:6){
  sub_peaks <- allPeaks[num==i,]
  sub_peaks <- GRanges(sub_peaks)
  sub_peaks  <- annotatePeak(sub_peaks, tssRegion=c(-2500, 2500), TxDb = txdb, level = "gene")
  ##
  sharePeaksAnno[[paste("reproducible",i,sep="_")]] <- sub_peaks
}
all_peakAnno <-  annotatePeak(GRanges(allPeaks), tssRegion=c(-2500, 2500), TxDb = txdb, level = "gene")
sharePeaksAnno[["all_peak"]] <- all_peakAnno
save(sharePeaksAnno, file="ATAC/matrix1020/score5_res/sharePeaksAnno.RData")


pdf("ATAC/matrix1020/score5_res/Feature_Distribution_of_peak_for_reproducible_region.pdf", 
    width=7, height = 10)
plotAnnoBar(sharePeaksAnno)
dev.off()


library(gtools)
load("ATAC/matrix1020/score5_res/sharePeaksAnno.RData")
names(sharePeaksAnno)

comb <- combinations(length(names(sharePeaksAnno)), 2, names(sharePeaksAnno))
chisq_test_res_list <- list()
for(i in 1:6){
  peak_type1 <- comb[i, 1]
  peak_type2 <- comb[i, 2]
  comb_type <- paste(peak_type1, peak_type2, sep = "-")
  anno1 <- sharePeaksAnno[[peak_type1]]@anno$annotation
  anno2 <- sharePeaksAnno[[peak_type2]]@anno$annotation
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


anno <- data.frame(sharePeaksAnno[["reproducible_6"]]@anno)
# feature <- sapply(anno$annotation, function(x){strsplit(x, " ")[[1]][1]})
feature <- anno$annotation
sum(feature %in% "Promoter (<=1kb)")/length(feature)



##########
# cor(ATAC and RNA) 
#########
library(reshape2)
library(DESeq2)
library(WGCNA)
library(DiffBind)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData") # Rhesus_GeneObject
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("ATAC/matrix1020/score5_res/region_occur.RData")
load("ATAC/matrix1020/score5_res/RegionPeak_33.RData")
load("ATAC/matrix1020/score5_res/sharePeaksAnno.RData")


vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)

all_peak <- RegionPeak$peaks[[1]][, 1:3]
rownames(all_peak) <- paste(all_peak[,1], ":", all_peak[,2], "-", all_peak[,3], sep="")
colnames(all_peak) <- c("chrom", "start", "end") 
binding <- data.frame(RegionPeak$binding)
rownames(binding) <- rownames(all_peak)
binding <- binding[,4:dim(binding)[2]]

##
peak_df <- log2(binding+1)
gene_df <- vsd_data[, colnames(binding)]


########
num <- apply(region_occur,1,function(x){sum(x>=1)})

sharePeaks_cor <- list()
for(i in 1:6){
  sub_peaks <- all_peak[num==i,]
  ##### peak split chrom
  peak_chrom_split_list <- list()
  for(chrom in unique(all_peak$chrom)){
    peak_chrom_split_list[[chrom]] <- rownames(sub_peaks)[sub_peaks$chrom == chrom]
  }
  ##### gene split chrom
  gene_chrom_split_list <- list()
  for(chrom in unique(all_peak$chrom)){
    gene_chrom_split_list[[chrom]] <- rownames(gene_df)[as.character(gtf_ensembl_gene[rownames(gene_df),"seqnames"]) == chrom]
  }
  ### cor with dif chrom
  peak_gene_cor_list <- list()
  # peak_gene_cor0.7_list <- list()
  # peak_gene_FDR_list <- list()
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
    peak_gene_cor <- data.frame(peak_gene_cor, sub_peaks[as.character(peak_gene_cor$peak),])
    peak_gene_cor <- data.frame(peak_gene_cor,
                                gtf_ensembl_gene[as.character(peak_gene_cor$gene), c(1,2,3,4)])
    colnames(peak_gene_cor)[6:11] <- c("peak_chr", "peak_start", "peak_end", "gene_chr", "gene_start", "gene_end")
    peak_gene_cor$peak_start <- as.numeric(peak_gene_cor$peak_start)
    peak_gene_cor$peak_end <- as.numeric(peak_gene_cor$peak_end)
    peak_gene_cor$strand <- as.character(peak_gene_cor$strand)
    peak_gene_cor_list[[chrom]] <- peak_gene_cor
    # peak_gene_cor0.7 <- peak_gene_cor[which(abs(peak_gene_cor$correlation)>=0.70),]
    # dim(peak_gene_cor0.7)
    # peak_gene_cor0.7_list[[chrom]] <- peak_gene_cor0.7
    # peak_gene_FDR <- peak_gene_cor[which(peak_gene_cor$FDR<0.1),]
    # peak_gene_FDR_list[[chrom]] <- peak_gene_FDR
  }
  ##
  sharePeaks_cor[[paste("reproducible",i,sep="_")]] <- peak_gene_cor_list
  
}


save(sharePeaks_cor,
     file="ATAC/matrix1020/score5_res/sharePeaks_cor.RData")


###
load("ATAC/matrix1020/score5_res/sharePeaks_cor.RData")


sharePeaks_cor_FDR <- list()
for(name in names(sharePeaks_cor)){
  rep_cor <- sharePeaks_cor[[name]]
  for(chorm in names(rep_cor)){
    peak_gene <- rep_cor[[chorm]]
    peak_gene_FDR <- peak_gene[abs(peak_gene$correlation)>0.7 & peak_gene$FDR<0.1, ]
    distance <- c()
    if(dim(peak_gene_FDR)[1]==0) next
    for(i in 1:dim(peak_gene_FDR)[1]){
      x <- peak_gene_FDR[i,]
      if(x[12]=="+"){d <- unlist(c(x[7]-x[10],x[8]-x[10]));
      d <- unique(d[abs(d) == min(abs(d))])[1]}
      if(x[12]=="-"){d <- unlist(c(x[11]-x[7],x[11]-x[8]));
      d <- unique(d[abs(d) == min(abs(d))])[1]}
      if(length(d)>1){print(d)}
      distance <- c(distance, d)
    }
    peak_gene_FDR <- data.frame(peak_gene_FDR, distance=distance)
    sharePeaks_cor_FDR[[name]][[chorm]] <- peak_gene_FDR
  }
}





#########################################
library(gProfileR)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
library(ggplot2)
library(dplyr)
library(ChIPseeker)
library(GenomicFeatures)
library(DiffBind)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix1020/score5_res/peak_gene_FDR.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("ATAC/matrix1020/score5_res/CA1_PCG_gene_res.RData")


sharePeaks_cor_FDR <- list()
for(name in names(sharePeaks_cor)){
  rep_cor <- sharePeaks_cor[[name]]
  for(chorm in names(rep_cor)){
    peak_gene <- rep_cor[[chorm]]
    peak_gene_FDR <- peak_gene[abs(peak_gene$correlation)>0.7 & peak_gene$FDR<0.1, ]
    distance <- c()
    if(dim(peak_gene_FDR)[1]==0) next
        sharePeaks_cor_FDR[[name]][[chorm]] <- peak_gene_FDR
  }
}


#######
load("ATAC/matrix1020/score5_res/sharePeaksAnno.RData")
function_list <- list()
for(name in names(sharePeaks_cor_FDR)){
  # rep6 <- sharePeaksAnno[[name]]
  # rep6 <- data.frame(rep6@anno)
  # rep6_promoter <- rep6[abs(rep6$distanceToTSS)<3000,]
  # genes_tmp <- unique(rep6_promoter[, "geneId"])
  ##
  sub <- Reduce(rbind, sharePeaks_cor_FDR[[name]])
  genes <- as.character(unique(sub[sub$correlation>0.7, "gene"]))
  # genes <- intersect(genes_tmp, genes)
  # gtf_ensembl_gene[genes, "gene_name"]
  gprofiler_res <- gprofiler(query=genes, organism = "mmulatta")
  gprofiler_res <- gprofiler_res[order(gprofiler_res[,3]),]
  gprofiler_res$gene_name <- gprofiler_res$intersection
  if(dim(gprofiler_res)[1]>0){
    for(i in 1:dim(gprofiler_res)[1]){
      x <- gprofiler_res$gene_name[i]
      gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],6], collapse=",")
    }
  }
  function_list[[name]]  <- gprofiler_res[, c("p.value", "domain", "term.name")]
}

rep6 <- sharePeaksAnno[[name]]
rep6 <- data.frame(rep6@anno)
rep6_promoter <- rep6[abs(rep6$distanceToTSS)<3000,]
genes_tmp <- unique(rep6_promoter[, "geneId"])
gprofiler_res <- gprofiler(query=genes_tmp, organism = "mmulatta")




###############################################################
#                         t-SNE                               #
###############################################################
library(Rtsne)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix1020/score5_res/RegionPeak.RData")
load("ATAC/matrix1020/score5_res/Rhesus_peak_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")

binding <- RegionPeak$binding[,4:36]
variance <- rowVars(binding)
sub <- binding[order(variance, decreasing = T)[1:1000],]

# PCA-based t-SNE
seed = 1:50
per = c(5, 10, 15, 20, 25, 30, 35)
for(i in seed){
  # for(perplexity in per){
  i=46
  set.seed(i)
  perplexity = 3
  tsne <- Rtsne::Rtsne(t(sub), initial_dims=50,
                       max_iter = 5000,perplexity=perplexity)
  rownames(tsne$Y) <- colnames(sub)
  colnames(tsne$Y) <- c("tsne1", "tsne2")
  p <- data.frame(tsne$Y, Rhesus_peak_colData[paste("X", rownames(tsne$Y),sep=""),])
  head(p)
  title = paste("t-SNE : seed=", i, " perplexity=", perplexity, sep="")
  
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
  
  
  filename = paste("ATAC/matrix1020/score5_res/Rhesus_Peak_t-SNE_top1000_seed",i,"_per", perplexity," .pdf", sep="")
  ggsave(filename,width = 6, height = 5)
  # }
}





###########################################################
#                       SR compare                        #
###########################################################
# /mnt/data1/Tools/R-3.6.0/bin/R
library(gtools)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix1020/score5_res/RegionPeak_33.RData")
sampleSheet <- read.csv("ATAC/matrix1020/score5_res/sampleSheet.csv",
                        stringsAsFactors = F)


######
combs <- combinations(length(unique(sampleSheet$SR)), 2, unique(sampleSheet$SR))
dba.report_list <- list()
for(i in 1:dim(combs)[1]){
  region1 <- combs[i, 1]
  region2 <- combs[i, 2]
  # sub_sampleSheet <- sampleSheet[sampleSheet$SR %in% c(region1, region2), ]
  # write.csv(sub_sampleSheet, "ATAC/matrix1020/score5_res/sub_sampleSheet.csv")
  # RegionPeak <- dba(sampleSheet="ATAC/matrix1020/score5_res/sub_sampleSheet.csv")
  # RegionPeak <- dba.count(RegionPeak ,bParallel=T, score=DBA_SCORE_TMM_READS_FULL)
  ###
  test <- dba.contrast(RegionPeak, 
                       # categories=DBA_CONDITION,
                       group1=RegionPeak$masks[[region1]],
                       group2=RegionPeak$masks[[region2]],
                       name1=region1,
                       name2=region2
                             )
  test <- dba.analyze(test)
  # dba.show(test, bContrast=T)
  test <- dba.report(test)
  # test2 <- test[abs(test$Fold)>0.8 & test$FDR<0.1, ]
  # dim(data.frame(test2))
  dba.report_list[[paste(region1, region2, sep="_")]] <- test
}

save(dba.report_list, file="ATAC/matrix1020/score5_res/SR_compare_dba.report_list.RData")


diff_num <- list()
for(vs in names(dba.report_list)){
  test <- dba.report_list[[vs]]
  test2 <- test[abs(test$Fold)>1 & test$FDR<0.1, ]
  diff_num[[vs]] <- dim(data.frame(test2))[1]
}

melt(diff_num)

###########################################################
#                       CA1 vs PCG                        #
###########################################################
RegionPeak <- dba(sampleSheet="ATAC/matrix1020/score5_res/CA1_and_PCG_sampleSheet.csv")
RegionPeak <- dba.count(RegionPeak ,bParallel=T, score=DBA_SCORE_TMM_READS_FULL)
RegionPeak <- dba.contrast(RegionPeak, categories=DBA_CONDITION)
RegionPeak <- dba.analyze(RegionPeak)
RegionPeak.DB <- dba.report(RegionPeak)


dba.show(RegionPeak,bContrast=TRUE)
RegionPeak$config$AnalysisMethod
pdf("ATAC/matrix1020/score5_res/MA_plot_using_DB.pdf")
dba.plotMA(RegionPeak_DB, method=RegionPeak_DB$config$AnalysisMethod,
           th=0.1,fold=1)
dev.off()


RegionPeak_DB <- RegionPeak
save(RegionPeak_DB, 
     file="ATAC/matrix1020/score5_res/RegionPeak_DB.RData")


#################
# EnrichedHeatmap
#################
library(DESeq2)
library(DiffBind)
library(EnrichedHeatmap)
library(rtracklayer)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix1020/score5_res/Rhesus_peak_colData.RData")
load("ATAC/matrix1020/score5_res/RegionPeak_DB.RData")


RegionPeak.DB <- dba.report(RegionPeak_DB)
df <- data.frame(RegionPeak.DB)
rownames(df) <- paste(df[,1], ":", df[,2], "-", df[,3], sep="")
CA1_hyper <- rownames(df)[df$Fold>1 & df$FDR<0.1]
CA1_hypo <- rownames(df)[df$Fold< -1 & df$FDR<0.1]
peaks <- c(CA1_hyper, CA1_hypo)
order <- c(setNames(rep("hyper", length(CA1_hyper)), CA1_hyper), 
           setNames(rep("hypo", length(CA1_hypo)), CA1_hypo))


###############
tmp.extension <- 3000
targets <- GRanges(df[peaks,1:3])
targets_extended <- resize(targets, fix = "center", width = tmp.extension*2)
###
matrix_mean <- function(r){
  samples <- rownames(Rhesus_peak_colData)[Rhesus_peak_colData$SR %in% r]
  mat_list <- NULL
  for(sample_id in samples){
    bw_file = paste("ATAC/Genrich_res/bw/", gsub("X", "", sample_id), ".bw", sep="")
    tmp.bigwig <- rtracklayer::import(bw_file, 
                                      format = "BigWig", 
                                      selection = BigWigSelection(targets_extended))
    normMatrix <- normalizeToMatrix(signal = tmp.bigwig, 
                                    target = resize(targets, fix = "center", width = 1), 
                                    background = 0, 
                                    keep = c(0, 0.99),      ## minimal value to the 99th percentile
                                    target_ratio = 0,
                                    mean_mode = "w0",       ## see ?EnrichedHeatmap on other options
                                    value_column = "score", ## = the name of the 4th column of the bigwig
                                    extend = tmp.extension)
    mat_list[[sample_id]] <- normMatrix
  }
  ##
  mat_mean = getSignalsFromList(mat_list)
  return(mat_mean)
}

#########
enrich_func <- function(r, mat_mean, value){
  # col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21)
  col_fun = circlize::colorRamp2(c(0, value), c("#6BAED6", "#E31A1C"))
  ## heatmap function:
  enrHtmp <- EnrichedHeatmap( mat = mat_mean, 
                              pos_line = FALSE, ## no dashed lines around the start
                              border = FALSE,   ## no box around heatmap
                              col = col_fun,    ## color gradients from above
                              column_title = r, 
                              column_title_gp = gpar(fontsize = 15, fontfamily = "sans"),
                              # row_order = ,
                              split = order,
                              use_raster = TRUE, raster_quality = 10, 
                              rect_gp = gpar(col = "transparent"), 
                              heatmap_legend_param = list(legend_direction = "vertical", 
                                                          title = "legend"),
                              top_annotation = HeatmapAnnotation(enriched = anno_enriched(gp = gpar(col = c("#1F78B4", "#E31A1C"), lty = 1, lwd=2),  
                                                                                          ylim=c(0.03,0.18)
                                                                                          ))
  ) 
  return(enrHtmp)
}

###########
CA1_mat_mean <- matrix_mean("CA1")
PCG_mat_mean <- matrix_mean("PCG")
other_mat_mean <- matrix_mean(c("SPL", "pSOG", "MFG", "ITG"))
# value <- max(max(CA1_mat_mean), max(PCG_mat_mean))
value <- max(max(CA1_mat_mean), max(PCG_mat_mean), max(other_mat_mean))
CA1_enrHtmp <- enrich_func("CA1", CA1_mat_mean, value)
PCG_enrHtmp <- enrich_func("PCG", PCG_mat_mean, value)
other_enrHtmp <- enrich_func(c("other cortex"), other_mat_mean, value)



##########
pdf(paste("ATAC/matrix1020/score5_res/EnrichedHeatmap_CA1_vs_PCG.pdf", sep=""), 
    width = 5, height = 6)
draw(CA1_enrHtmp + PCG_enrHtmp + other_enrHtmp,                            
     heatmap_legend_side = "right",     
     annotation_legend_side = "right",
     padding = unit(c(4, 4, 4, 4), "mm") ## some padding to avoid labels beyond plot borders
)
dev.off() 



##############
# ChIPseeker 
##############
library(ChIPseeker)
library(GenomicFeatures)
library(DESeq2)
library(GenomicRanges)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix1020/score5_res/RegionPeak_DB.RData")



txdb <- makeTxDbFromGFF("/mnt/data1/Ref/Rhesus_macaque_8.0.1/Macaca_mulatta.Mmul_8.0.1.91.chr.gtf",
                        # "TACO_minExpr_5.5/assembly.refcomp.gtf",
                        format = "gtf",
                        dataSource = "ensemblgenomes",
                        organism = "Macaca mulatta")

promoter <- getPromoters(TxDb=txdb, upstream=2500, downstream=2500)


#############
### all peaks 
RegionPeak_DB
all_peaks <- RegionPeak_DB$peaks[[1]][, 1:3]
save(all_peaks, file="ATAC/matrix1020/score5_res/all_peaks.RData")


all_peakAnno  <- annotatePeak(GRanges(all_peaks), tssRegion=c(-2500, 2500), TxDb = txdb, level = "gene")
save(all_peakAnno, file="ATAC/matrix1020/score5_res/all_peakAnno.RData")


############
### dif peaks 
RegionPeak.DB <- dba.report(RegionPeak_DB)
df <- data.frame(RegionPeak.DB)
rownames(df) <- paste(df[,1], ":", df[,2], "-", df[,3], sep="")
CA1_hyper <- df[df$Fold>1 & df$FDR<0.1,]
CA1_hypo <- df[df$Fold< -1 & df$FDR<0.1,]
peaks <- rbind(CA1_hyper, CA1_hypo)


peakAnno  <- annotatePeak(GRanges(peaks), tssRegion=c(-2500, 2500), TxDb = txdb, level = "gene")
HIPpeakAnno  <- annotatePeak(GRanges(CA1_hyper), tssRegion=c(-2500, 2500), TxDb = txdb, level = "gene")
CortexpeakAnno  <- annotatePeak(GRanges(CA1_hypo), tssRegion=c(-2500, 2500), TxDb = txdb, level = "gene")


##########
peakAnno_list <- list()
peakAnno_list[["HIP_peaks"]] <- HIPpeakAnno
peakAnno_list[["Cortex_peaks"]] <- CortexpeakAnno
peakAnno_list[["dif_peaks"]] <- peakAnno
peakAnno_list[["all_peaks"]] <- all_peakAnno
save(peakAnno_list, file="ATAC/matrix1020/score5_res/peakAnno_list.RData")


#
pdf("ATAC/matrix1020/score5_res/Feature_Distribution_of_peak.pdf", width=7, height = 3)
plotAnnoBar(peakAnno_list)
dev.off()


######################################
## chisq.test for Feature Distribution
######################################
library(gtools)
load("ATAC/matrix1020/score5_res/peakAnno_list.RData")
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







###########################################################
#             correlation  ATAC and RNA                   #
###########################################################
library(reshape2)
library(DESeq2)
library(WGCNA)
library(DiffBind)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData") # Rhesus_GeneObject
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("ATAC/matrix1020/score5_res/RegionPeak_DB.RData")

res <- results(Rhesus_GeneObject, contrast = c("SR", "CA1", "PCG"))
# res <- res[complete.cases(res),]
CA1_PCG_gene_res <- res
save(CA1_PCG_gene_res, file="ATAC/matrix1020/score5_res/CA1_PCG_gene_res.RData")





vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)

all_peak <- RegionPeak_DB$peaks[[1]][, 1:3]
rownames(all_peak) <- paste(all_peak[,1], ":", all_peak[,2], "-", all_peak[,3], sep="")
colnames(all_peak) <- c("chrom", "start", "end") 
binding <- data.frame(RegionPeak_DB$binding)
rownames(binding) <- rownames(all_peak)
binding <- binding[,4:dim(binding)[2]]

##
peak_df <- log2(binding+1)
gene_df <- vsd_data[, colnames(binding)]


########
RegionPeak.DB <- dba.report(RegionPeak_DB)
df <- data.frame(RegionPeak.DB)
rownames(df) <- paste(df[,1], ":", df[,2], "-", df[,3], sep="")
DifPeaks <- df[abs(df$Fold)>1 & df$FDR<0.1, 1:3]
colnames(DifPeaks)[1] <- "chrom"


##### peak split chrom
peak_chrom_split_list <- list()
for(chrom in unique(all_peak$chrom)){
  peak_chrom_split_list[[chrom]] <- rownames(DifPeaks)[DifPeaks$chrom == chrom]
}
##### gene split chrom
gene_chrom_split_list <- list()
for(chrom in unique(all_peak$chrom)){
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
  peak_gene_cor <- data.frame(peak_gene_cor, DifPeaks[as.character(peak_gene_cor$peak),])
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

save(peak_gene_cor0.7_list, peak_gene_FDR_list,
     file="ATAC/matrix1020/score5_res/DifPeak_gene_cor_list.RData")


###
load("ATAC/matrix1020/score5_res/DifPeak_gene_cor_list.RData")

# for(chorm in names(peak_gene_cor0.7_list)){
#   peak_gene_cor0.7 <- peak_gene_cor0.7_list[[chorm]]
#   distance <- c()
#   if(dim(peak_gene_cor0.7)[1]==0) next
#   for(i in 1:dim(peak_gene_cor0.7)[1]){
#     x <- peak_gene_cor0.7[i,]
#     if(x[12]=="+"){d <- unlist(c(x[7]-x[10],x[8]-x[10])); 
#     d <- unique(d[abs(d) == min(abs(d))])[1]}
#     if(x[12]=="-"){d <- unlist(c(x[11]-x[7],x[11]-x[8]));
#     d <- unique(d[abs(d) == min(abs(d))])[1]}
#     if(length(d)>1){print(d)}
#     distance <- c(distance, d)
#   }
#   peak_gene_cor0.7 <- data.frame(peak_gene_cor0.7, distance=distance)
#   peak_gene_cor0.7_list[[chorm]] <- peak_gene_cor0.7
# }

for(chorm in names(peak_gene_FDR_list)){
  peak_gene_FDR <- peak_gene_FDR_list[[chorm]]
  distance <- c()
  if(dim(peak_gene_FDR)[1]==0) next
  for(i in 1:dim(peak_gene_FDR)[1]){
    x <- peak_gene_FDR[i,]
    if(x[12]=="+"){d <- unlist(c(x[7]-x[10],x[8]-x[10]));
    d <- unique(d[abs(d) == min(abs(d))])[1]}
    if(x[12]=="-"){d <- unlist(c(x[11]-x[7],x[11]-x[8]));
    d <- unique(d[abs(d) == min(abs(d))])[1]}
    if(length(d)>1){print(d)}
    distance <- c(distance, d)
  }
  peak_gene_FDR <- data.frame(peak_gene_FDR, distance=distance)
  peak_gene_FDR_list[[chorm]] <- peak_gene_FDR
}

peak_gene_FDR <- Reduce(rbind, peak_gene_FDR_list)

# distance <- distance(GRanges(peak_gene_FDR[,c(6:8,12)]), GRanges(peak_gene_FDR[,9:12]), ignore.strand=F)
# peak_gene_FDR$distance <- distance
dim(peak_gene_FDR) # 24578    13

save(peak_gene_FDR, 
     file="ATAC/matrix1020/score5_res/peak_gene_FDR.RData")



#########################################
library(gProfileR)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
library(ggplot2)
library(dplyr)
library(ChIPseeker)
library(GenomicFeatures)
library(DiffBind)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix1020/score5_res/peak_gene_FDR.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("ATAC/matrix1020/score5_res/RegionPeak_DB.RData")
load("ATAC/matrix1020/score5_res/CA1_PCG_gene_res.RData")
res <- CA1_PCG_gene_res[complete.cases(CA1_PCG_gene_res),]
DifGenes <- res[res$padj<0.05 & abs(res$log2FoldChange)>1,]


RegionPeak.DB <- dba.report(RegionPeak_DB)
df <- data.frame(RegionPeak.DB)
rownames(df) <- paste(df[,1], ":", df[,2], "-", df[,3], sep="")
CA1_hyper <- df[df$Fold>1 & df$FDR<0.1,]
CA1_hypo <- df[df$Fold< -1 & df$FDR<0.1,]
# peaks <- rbind(CA1_hyper, CA1_hypo)


sub_list <- list()
sub_list[["CA1_hyper"]] <- peak_gene_FDR[peak_gene_FDR$peak %in% rownames(CA1_hyper), ]
sub_list[["CA1_hypo"]] <- peak_gene_FDR[peak_gene_FDR$peak %in% rownames(CA1_hypo), ]


for(name in names(sub_list)){
  sub <- sub_list[[name]]
  sub <- sub[order(sub$correlation, decreasing = T),]
  sub <- sub[sub$correlation>0.7,]
  sub$gene_name <- gtf_ensembl_gene[as.character(sub$gene), "gene_name"]
  sub$gene_FC <- CA1_PCG_gene_res[as.character(sub$gene), "log2FoldChange"]
  # table(is.na(CA1_PCG_gene_res[as.character(sub$gene), "log2FoldChange"]))
  sub$peak_FC <- df[as.character(sub$peak), "Fold"]
  sub <- sub[order(sub$peak),]
  write.table(sub, paste("ATAC/matrix1020/score5_res/", name, "_link_table.txt", sep=""),
              col.names=T, row.names = F, sep="\t", quote = F)
}




function_list <- list()
for(state in names(sub_list)){
  sub <- sub_list[[state]]
  genes <- as.character(unique(sub[sub$correlation>0.7, "gene"]))
  print(length(genes))
  # gprofiler_res <- gprofiler(query=genes, organism = "mmulatta")
  # gprofiler_res <- gprofiler_res[order(gprofiler_res[,3]),]
  # gprofiler_res$gene_name <- gprofiler_res$intersection
  # if(dim(gprofiler_res)[1]>0){
  #   for(i in 1:dim(gprofiler_res)[1]){
  #     x <- gprofiler_res$gene_name[i]
  #     gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],6], collapse=",")
  #   }
  # }
  # function_list[[state]]  <- gprofiler_res[gprofiler_res$domain == "BP", c("p.value", "domain", "term.name", "gene_name")]
}

for(name in names(function_list)){
  fuc <- function_list[[name]]
  write.table(fuc, paste("ATAC/matrix1020/score5_res/", name, "_function_.txt", sep=""), 
              col.names=T, row.names = F, sep="\t", quote = F)
}


####
pp <- function_list[["CA1_hyper"]]
order = pp$term.name[order(pp$p.value, decreasing=T)]
pp$term.name <- factor(pp$term.name, levels = order)

ggplot(pp[1:10,], aes(x=term.name, y=-log10(p.value))) +
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
  ggtitle("CA1_hyper") +
  xlab("") 


ggsave("ATAC/matrix1020/score5_res/CA1_hyper_function_bar.pdf",
       height = 7, width = 6)


###
pp <- function_list[["CA1_hypo"]]
order = pp$term.name[order(pp$p.value, decreasing=T)]
pp$term.name <- factor(pp$term.name, levels = order)

ggplot(pp[1:10,], aes(x=term.name, y=-log10(p.value))) +
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
  ggtitle("CA1_hypo") +
  xlab("") 


ggsave("ATAC/matrix1020/score5_res/CA1_hypo_function_bar.pdf",
       height = 7, width = 6)



################################################
#                 HOMER motif                  #                                      
################################################
# /mnt/data1/Tools/R-3.6.0/bin/R
library(DiffBind)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix1020/score5_res/RegionPeak_DB.RData")
### dif peaks 
RegionPeak.DB <- dba.report(RegionPeak_DB)
df <- data.frame(RegionPeak.DB)
rownames(df) <- paste(df[,1], ":", df[,2], "-", df[,3], sep="")
CA1_hyper <- df[df$Fold>1 & df$FDR<0.1,]
CA1_hypo <- df[df$Fold< -1 & df$FDR<0.1,]

DifPeaks <- list()
DifPeaks[["CA1_hyper"]] <- CA1_hyper[,1:3]
DifPeaks[["CA1_hypo"]] <- CA1_hypo[,1:3]


###
setwd("/mnt/data2/Rhesus_brain/ATAC/matrix1020/score5_res")
## all peaks
load("all_peaks.RData")
head(all_peaks)
all_peaks$peak_name <- paste(all_peaks[,1],":",all_peaks[,2],"-",all_peaks[,3],sep="")
filename = paste("all_peaks.bed", sep="")
write.table(all_peaks, filename, row.names = F, col.names = F, quote = F, sep="\t")
system(paste("sort -k1,1 -k2,2n all_peaks.bed -o all_peaks.sorted.bed", sep=""))
system(paste("/mnt/data1/Tools/bedtools2/bedtools2/bin/bedtools getfasta ",  
             "-fi /mnt/data1/Ref/Rhesus_macaque_8.0.1/Macaca_mulatta.Mmul_8.0.1.dna_sm.fa ", 
             "-bed all_peaks.sorted.bed ", 
             "-fo all_peaks.sorted.fasta", sep=""))

for(state in names(DifPeaks)){
  peaks <- DifPeaks[[state]]
  # for(i in 1:dim(peaks)[1]){
  #   start <- peaks[i,2]
  #   end <- peaks[i,3]
  #   mid <- round((start + end)/2,0)
  #   peaks[i,2] <- mid-250
  #   peaks[i,3] <- mid+250
  # }
  ###
  peaks <- data.frame(peaks, peak_name=rownames(peaks))
  head(peaks)
  filename = paste(state, "_peaks.bed", sep="")
  write.table(peaks, filename, row.names = F, col.names = F, quote = F, sep="\t")
  system(paste("sort -k1,1 -k2,2n ", state, "_peaks.bed -o ", state, "_peaks.sorted.bed", sep=""))
  system(paste("/mnt/data1/Tools/bedtools2/bedtools2/bin/bedtools getfasta ",  
               "-fi /mnt/data1/Ref/Rhesus_macaque_8.0.1/Macaca_mulatta.Mmul_8.0.1.dna_sm.fa ", 
               "-bed ",state,"_peaks.sorted.bed ", 
               "-fo ", state, "_peaks.sorted.fasta", sep=""))
  system(paste("/home/looking/software/HOMER/bin/findMotifs.pl ", 
         state, "_peaks.sorted.fasta rheMac8 ", state, "2/ -fastaBg all_peaks.sorted.fasta -fdr", sep=""))
}



################################################
#                 AME motif                  #                                      
################################################
setwd("/mnt/data2/Rhesus_brain/ATAC/matrix1020/score5_res")
for(state in names(DifPeaks)){
  peaks <- DifPeaks[[state]]
  # for(i in 1:dim(peaks)[1]){
  #   start <- peaks[i,2]
  #   end <- peaks[i,3]
  #   mid <- round((start + end)/2,0)
  #   peaks[i,2] <- mid-250
  #   peaks[i,3] <- mid+250
  # }
  ##
  peaks <- data.frame(peaks, peak_name=rownames(peaks))
  head(peaks)
  filename = paste(state, "_peaks.bed", sep="")
  write.table(peaks, filename, row.names = F, col.names = F, quote = F, sep="\t")
  system(paste("sort -k1,1 -k2,2n ", state, "_peaks.bed -o ", state, "_peaks.sorted.bed", sep=""))
  system(paste("/mnt/data1/Tools/bedtools2/bedtools2/bin/bedtools getfasta ",
               "-fi /mnt/data1/Ref/Rhesus_macaque_8.0.1/Macaca_mulatta.Mmul_8.0.1.dna_sm.fa ",
               "-bed ",state,"_peaks.sorted.bed ",
               "-fo ", state, "_peaks.sorted.fasta", sep=""))
  system(paste("/home/looking/software/meme-5.1.0/src/ame --normalise-affinity ",
              "--o AME_", state, "_normalise/ --control all_peaks.sorted.fasta ",
               state, "_peaks.sorted.fasta /home/looking/software/meme-5.1.0/motif_databases/EUKARYOTE/*.meme ",
               sep=""))
}

state <- "CA1_hyper"
system(paste("/home/looking/software/meme-5.1.0/src/centrimo ", 
             "--o centrimo_", state, "/ --neg CA1_hypo_peaks.sorted.fasta ",
             state, "_peaks.sorted.fasta /home/looking/software/meme-5.1.0/motif_databases/EUKARYOTE/*.meme ",
             sep=""))

state <- "CA1_hypo"
system(paste("/home/looking/software/meme-5.1.0/src/centrimo ", 
             "--o centrimo_", state, "/ --neg CA1_hyper_peaks.sorted.fasta ",
             state, "_peaks.sorted.fasta /home/looking/software/meme-5.1.0/motif_databases/EUKARYOTE/*.meme ",
             sep=""))


################################################
#                 GimmeMotifs                  #                                      
################################################
# /mnt/data1/Tools/R-3.6.0/bin/R
library(DiffBind)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix1020/score5_res/RegionPeak_DB.RData")
### dif peaks 
RegionPeak.DB <- dba.report(RegionPeak_DB)
df <- data.frame(RegionPeak.DB)
rownames(df) <- paste(df[,1], ":", df[,2], "-", df[,3], sep="")
CA1_hyper <- df[df$Fold>1 & df$FDR<0.1,]
CA1_hypo <- df[df$Fold< -1 & df$FDR<0.1,]

DifPeaks <- rbind(data.frame(loc=rownames(CA1_hyper), cluster="CA1_hyper"),
                  data.frame(loc=rownames(CA1_hypo), cluster="CA1_hypo"))
DifPeaks$loc <- paste("chr", DifPeaks$loc, sep = "")
write.table(DifPeaks, "ATAC/matrix1020/score5_res/DifPeaks.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")

setwd("/mnt/data2/Rhesus_brain/ATAC/matrix1020/score5_res")
cmd =paste("gimme maelstrom DifPeaks.txt ",
       "rheMac8 maelstrom.out/", sep="")


##
df <- read.table("maelstrom.out/final.out.csv")
m2f <- read.table("maelstrom.out/gimme.vertebrate.v5.0.motif2factors.txt", sep="\t", header = T)


df <- df[abs(df$CA1_hyper-df$CA1_hypo) >2 & (abs(df$CA1_hyper)>2 | abs(df$CA1_hypo)>2),]

factors <- c()
for(motif in rownames(df)){
  factor <- unique(toupper(as.character(m2f[m2f$Motif %in% motif, "Factor"])))
  factors <- c(factors, paste(factor, collapse = ","))
}

rownames(df) <- paste(rownames(df), ":", factors, sep="")

pheatmap(
  df,
  col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
  # breaks = c(-3,seq(-2,2,length=20),3),
  # legend_breaks = c(-2,-1,1,2),
  filename = "motif_enrichment.pdf",
  # scale="row",
  fontsize = 4,
  # clustering_method = "complete",
  border_color=NA,
  cluster_rows = T
)



################################################
#                   enhancer                   #                                      
################################################
# /mnt/data1/Tools/R-3.6.0/bin/R
library(DiffBind)
setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix1020/score5_res/RegionPeak_DB.RData")
load("ATAC/matrix1020/score5_res/all_peaks.RData")

### dif peaks 
RegionPeak.DB <- dba.report(RegionPeak_DB)
df <- data.frame(RegionPeak.DB)
rownames(df) <- paste(df[,1], ":", df[,2], "-", df[,3], sep="")
CA1_hyper <- df[df$Fold>1 & df$FDR<0.1,]
CA1_hypo <- df[df$Fold< -1 & df$FDR<0.1,]



Mac3 <- read.table("ATAC/liftover_res/rheMac3_enhancer.txt", 
                   sep="\t", header=T)
rownames(Mac3) <- paste(gsub("chr", "", Mac3[,1]), ":", Mac3[,2], "-", Mac3[,3], sep="")
##
Mac8_Mac3 <- read.table("ATAC/liftover_res/rheMac8_enhancer.bed",
                        sep="\t", header=F)
colnames(Mac8_Mac3) <- c("chrom", "start", "end", "Mac3")
Mac8_Mac3 <- data.frame(Mac8_Mac3, Mac3[Mac8_Mac3$Mac3, c("Cortex", "Subcortical.Structures")])


# 
# awk -F "chr" '{print $2,$3}' rheMac3_enhancer_orthology_rheMac8.bed > rheMac8_enhancer.bed
system(paste("/mnt/data1/Tools/bedtools2/bedtools2/bin/bedtools intersect ",
"-a ATAC/matrix1020/score5_res/CA1_hyper_peaks.sorted.bed ",
"-b ATAC/liftover_res/rheMac8_enhancer.bed ", 
"-wa -wb > ATAC/matrix1020/score5_res/CA1_hyper_peaks.enhancer.bed", sep=""))


CA1.hyper.inter <- read.table("ATAC/matrix1020/score5_res/CA1_hyper_peaks.enhancer.bed",
                              stringsAsFactors = F)
CA1.hyper.inter$Subcortical.Structures <- Mac3[CA1.hyper.inter[,8], "Subcortical.Structures"]
enhancer_peak <- unique(CA1.hyper.inter$V4[CA1.hyper.inter$Subcortical.Structures==1])
length(enhancer_peak)/dim(CA1_hyper)[1] # 0.238806
CA1_hyper_enhancer_peak <- enhancer_peak

pdf("ATAC/matrix1020/score5_res/CA1.hyper_intersect_enhancer_pie.pdf",
    width = 6, height = 5)
pp <- c(length(enhancer_peak), dim(CA1_hyper)[1]-length(enhancer_peak))
pie(pp,
    labels=paste(c("in Subcortical.Structures enhancer\n", "not in Subcortical.Structures enhancer\n"), "(", pp, ")", sep=""),
    border=F,
    clockwise = T,
    init.angle = -35,
    main = "CA1.hyper",
    col = c("#B33B49", "#5784B4"))
dev.off()


##
system(paste("/mnt/data1/Tools/bedtools2/bedtools2/bin/bedtools intersect ",
             "-a ATAC/matrix1020/score5_res/CA1_hypo_peaks.sorted.bed ",
             "-b ATAC/liftover_res/rheMac8_enhancer.bed ", 
             "-wa -wb > ATAC/matrix1020/score5_res/CA1_hypo_peaks.enhancer.bed", sep=""))


CA1.hypo.inter <- read.table("ATAC/matrix1020/score5_res/CA1_hypo_peaks.enhancer.bed",
                              stringsAsFactors = F)
CA1.hypo.inter$Cortex <- Mac3[CA1.hypo.inter[,8], "Cortex"]
enhancer_peak <- unique(CA1.hypo.inter$V4[CA1.hypo.inter$Cortex==1])
length(enhancer_peak)/dim(CA1_hypo)[1] # 0.6082192
CA1_hypo_enhancer_peak <- enhancer_peak

pdf("ATAC/matrix1020/score5_res/CA1.hypo_intersect_enhancer_pie.pdf",
    width = 6, height = 5)
pp <- c(length(enhancer_peak), dim(CA1_hypo)[1]-length(enhancer_peak))
pie(pp,
    labels=paste(c("in Cortex enhancer\n", "not in Cortex enhancer\n"), "(", pp, ")", sep=""),
    border=F,
    clockwise = T,
    init.angle = -35,
    main = "CA1.hypo",
    col = c("#B33B49", "#5784B4"))
dev.off()


##
system(paste("/mnt/data1/Tools/bedtools2/bedtools2/bin/bedtools intersect ",
             "-a ATAC/matrix1020/score5_res/all_peaks.sorted.bed ",
             "-b ATAC/liftover_res/rheMac8_enhancer.bed ", 
             "-wa -wb > ATAC/matrix1020/score5_res/all_peaks.enhancer.bed", sep=""))


all_peaks.inter <- read.table("ATAC/matrix1020/score5_res/all_peaks.enhancer.bed",
                             stringsAsFactors = F)
all_peaks.inter <- data.frame(all_peaks.inter, Mac3[all_peaks.inter[,8], c("Cortex", "Subcortical.Structures")])
Cortex_enhancer <- unique(all_peaks.inter$V4[all_peaks.inter$Cortex==1])
Subcortical_enhancer <- unique(all_peaks.inter$V4[all_peaks.inter$Subcortical.Structures==1])

length(Cortex_enhancer)/dim(all_peaks)[1] # 0.160428
length(Subcortical_enhancer)/dim(all_peaks)[1] # 0.1690076


pdf("ATAC/matrix1020/score5_res/all_peaks_intersect_enhancer_pie.pdf",
    width = 6, height = 5)
pp <- c(length(Subcortical_enhancer), dim(all_peaks)[1]-length(Subcortical_enhancer))
pie(pp,
    labels=paste(c("in Subcortical.Structures enhancer\n", "not in Subcortical.Structures enhancer\n"), "(", pp, ")", sep=""),
    border=F,
    clockwise = T,
    init.angle = -35,
    main = "all_peaks",
    col = c("#B33B49", "#5784B4"))

pp <- c(length(Cortex_enhancer), dim(all_peaks)[1]-length(Cortex_enhancer))
pie(pp,
    labels=paste(c("in Cortex enhancer\n", "not in Cortex enhancer\n"), "(", pp, ")", sep=""),
    border=F,
    clockwise = T,
    init.angle = -35,
    main = "all_peaks",
    col = c("#B33B49", "#5784B4"))
dev.off()






pp1 <- c(length(CA1_hyper_enhancer_peak), dim(CA1_hyper)[1]-length(CA1_hyper_enhancer_peak))
pp2 <- c(length(CA1_hypo_enhancer_peak), dim(CA1_hypo)[1]-length(CA1_hypo_enhancer_peak))
pp3 <- c(length(Subcortical_enhancer), dim(all_peaks)[1]-length(Subcortical_enhancer))
pp4 <- c(length(Cortex_enhancer), dim(all_peaks)[1]-length(Cortex_enhancer))


chisq.test(matrix(c(pp1,pp3), ncol=2))$p.value
chisq.test(matrix(c(pp2,pp4), ncol=2))$p.value











