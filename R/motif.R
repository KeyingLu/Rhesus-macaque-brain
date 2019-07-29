
#############################################################
#                           ATAC                            #
#############################################################

setwd("/mnt/data2/Rhesus_brain")
load("ATAC/matrix/peak_top10000_res/ConsensusCluster_DifPeaks.RData")

###
k=4
res <- DIfPeaks_Results[[k]]
geneCluster <- data.frame(geneCluster=factor(sort(res$consensusClass)))
for(cluster in 1:k){
  peaks <- rownames(geneCluster)[geneCluster$geneCluster==cluster]
  peaks <- data.frame(chrom=sapply(peaks, function(x){strsplit(x, ":")[[1]][1]}),
                      start=as.numeric(sapply(peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1]})),
                      end=as.numeric(sapply(peaks, function(x){strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2]})),
                      stringsAsFactors = F)
  # for(i in 1:dim(peaks)[1]){
  #   start <- peaks[i,2]
  #   end <- peaks[i,3]
  #   mid <- round((start + end)/2,0)
  #   peaks[i,2] <- mid-250
  #   peaks[i,3] <- mid+250
  # }
  peaks <- data.frame(peaks, peak_name=rownames(geneCluster)[geneCluster$geneCluster==cluster])
  filename = paste("ATAC/matrix/peak_top10000_res/dif_peak_cluster", cluster, "_0710.bed", sep="")
  write.table(peaks, filename, row.names = F, col.names = F, quote = F, sep="\t")
}


# sort -k1,1 -k2,2n dif_peak_cluster1_0710.bed -o dif_peak_cluster1_sorted.bed
# /mnt/data1/Tools/bedtools2/bedtools2/bin/bedtools getfasta
# -fi /mnt/data1/Ref/Rhesus_macaque_8.0.1/Macaca_mulatta.Mmul_8.0.1.dna_sm.fa
# -bed dif_peak_cluster1_sorted.bed
# -fo dif_peak_cluster1_sorted.fasta




################
k = 4
motif_res <- list()
for(i in 1:k){
  cluster <- paste("cluster", i, sep="")
  filename = paste("ATAC/matrix/peak_top10000_res/ame_", cluster, "_0710.tsv", sep="")
  if(file.exists(filename)){res <- read.table(filename, sep="\t", header=T)[,c(3:4, 6:7)]}
  if(file.exists(filename)==F){res <- data.frame("none", "none", 1, 1)}
  head(sort(res$adj_p.value))
  res[,3:4] <- -log10(res[,3:4])
  colnames(res) <- c("ID", "Alt_ID", "logP", "logQ")
  rownames(res) <- res$ID
  res$logQ[which(res$logQ>=400)] <- 400
  motif_res[[cluster]] <- res
}


############
motif_ID <- Reduce(union, lapply(motif_res, function(x){x[,1][x[,4]>250]}))

motif_res_modify <- list()
for(cluster in names(motif_res)){
  res <- motif_res[[cluster]][motif_ID,]
  rownames(res) <- motif_ID
  res <- res[,c(2,4)]
  res[res$logQ<2 | is.na(res$logQ), 2] <- 1
  motif_res_modify[[cluster]] <- res
}


test <- data.frame(motif_res_modify, stringsAsFactors = F)
test_name <- test[,c(1:k)*2-1]
Alt_ID <- list()
for(i in 1:dim(test_name)[1]){
  x <- as.character(test_name[i, ])
  index <- rownames(test_name)[i]
  x[which(x=="NA")] <- "1"
  if(all(x=="1")){Alt_ID[[index]] <- index}
  if(all(x!="1")){
    x <- test_name[i, ]
    Alt_ID[[index]] <- as.character(x[1,1])
    next
  }
  if(any(x!="1")){
    x <- test_name[i, ]
    Alt_ID[[index]] <- as.character(unique(x[-which(is.na(x))])[1,1])}

}

test_name <- matrix("a", ncol=1, nrow=length(Alt_ID))
rownames(test_name) <- names(Alt_ID)
colnames(test_name) <- "Alt_ID"
for(i in rownames(test_name)){test_name[i, 1] <- Alt_ID[[i]]}


test_qvale <- test[,c(1:k)*2]
rownames(test_qvale) <- test_name[rownames(test_qvale),]
colnames(test_qvale) <- gsub(".logQ", "", colnames(test_qvale))


##
library(reshape2)
library(ggplot2)
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
pp <- test_qvale
name <- sapply(rownames(pp), function(x){strsplit(x, "_")[[1]][1]})
index1 <- toupper(name) %in% gtf_ensembl_gene$gene_name
name <- sapply(rownames(pp), function(x){strsplit(x, ":")[[1]][1]})
index2 <- toupper(name) %in% gtf_ensembl_gene$gene_name
rownames(pp)[index1 | index2]
pp <- test_qvale[index1 | index2,]
# sub_gene <- gtf_ensembl_gene[gtf_ensembl_gene$gene_name %in% rownames(pp),]
# sub_norCounts <- Rhesus_norCounts[sub_gene$gene_id,-which()]
# rownames()


pp <- melt(t(pp))
colnames(pp) <- c("cluster", "motif", "logQ")

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
   




###
load("/mnt/xdlab1/home/looking/brain_project/matrix0420/gene_count_res/Rhesus_DESeq2_object.RData")

pp <- test_qvale
name1 <- sapply(rownames(pp), function(x){strsplit(x, "_")[[1]][1]})
name2 <- sapply(rownames(pp), function(x){strsplit(x, "::")[[1]][1]})
name3 <- sapply(rownames(pp), function(x){strsplit(x, "::")[[1]][2]})
name4 <- sapply(rownames(pp), function(x){strsplit(x, "::")[[1]][3]})
name <- Reduce(union, list(name1, name2, name3, name4))
name <- name[-which(is.na(name))]
index <- toupper(name) %in% gtf_ensembl_gene$gene_name[-which(gtf_ensembl_gene$seqnames %in% c("X", "Y", "MT"))]
name <- name[index]



norCounts <- counts(Rhesus_GeneObject, normalized=TRUE)
sub_gene <- gtf_ensembl_gene[gtf_ensembl_gene$gene_name %in% toupper(name),]
sub_norCounts <- norCounts[sub_gene$gene_id, colnames(peak_vsd_data)]
rownames(sub_norCounts) <- sub_gene$gene_name
dim(sub_norCounts)

pp <- melt(sub_norCounts)
colnames(pp) <- c("gene", "sample_id", "expression")
pp$ID <- Rhesus_colData[as.character(pp$sample_id),"ID"]
head(pp)

##

pp$gene <- factor(pp$gene, levels=rownames(sub_norCounts)[order(rowVars(sub_norCounts), decreasing = T)])


ggplot(pp, aes(x=ID, y=expression, fill=ID, colour=ID)) +
  geom_violin(alpha=0.8) +
  geom_jitter(size=0.5, width = 0.1, height = 0.8, colour="black") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        rect = element_rect(colour = "black"),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = -90, hjust = 0.5, vjust=0.5, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "black", size=1),
        text = element_text(face="bold", colour = "black")) +
  facet_grid(.~gene, scales = "free_x") +
  scale_fill_manual(values = setNames(c("#838B8B", "#713F27", "#9B3D4A",
                    "#CD8463","#AFC88C", "#789EC1","#6783B1"), unique(pp$ID))) +
  scale_colour_manual(values = setNames(c("#838B8B", "#713F27", "#9B3D4A",
                    "#CD8463","#AFC88C", "#789EC1","#6783B1"), unique(pp$ID))) +
  coord_flip() +
  guides(fill=F, colour=F) +
  xlab("") +
  ylab("Gene Expression")

ggsave("ATAC/matrix/motif_gene_expression.pdf",
       width = 15, height = 4)




#############################################################
#                           lncRNA                          #
#############################################################

setwd("/mnt/data2/Rhesus_brain")

################
lncRNA_module_name <- c("ME2_CB", "ME3_cortex", 
                        "ME4_STR", "ME5_PIT", "ME7_HIP")

motif_res <- list()
for(module in lncRNA_module_name){
  module_file <- paste("ame_", module, "_0721.tsv", sep="")
  filename = paste("revise0530/transcripts_count_res_with_TACO_gtf/", module_file, sep="")
  if(file.exists(filename)){res <- read.table(filename, sep="\t", header=T)[,c(3:4, 6:7)]; if(dim(res)[1]>10){res <- res[1:10,]}}
  # if(file.exists(filename)){res <- read.table(filename, sep="\t", header=T)[,c(3:4, 6:7)]}
  if(file.exists(filename)==F){res <- data.frame("none", "none", 1, 1)}
  head(sort(res$adj_p.value))
  res[,3:4] <- -log10(res[,3:4])
  colnames(res) <- c("ID", "Alt_ID", "logP", "logQ")
  rownames(res) <- res$ID
  motif_res[[module]] <- res
}


############
motif_ID <- Reduce(union, lapply(motif_res, function(x){x[,1][x[,4]>1]}))

motif_res_modify <- list()
for(cluster in names(motif_res)){
  res <- motif_res[[cluster]][motif_ID,]
  rownames(res) <- motif_ID
  res <- res[,c(2,4)]
  res[res$logQ<1 | is.na(res$logQ), 2] <- 1
  motif_res_modify[[cluster]] <- res
}


test <- data.frame(motif_res_modify, stringsAsFactors = F)
k = 5
test_name <- test[,c(1:k)*2-1]
Alt_ID <- list()
for(i in 1:dim(test_name)[1]){
  x <- as.character(test_name[i, ])
  xindex <- rownames(test_name)[i]
  x[which(is.na(x))] <- "NO"
  if(all(x=="NO")){Alt_ID[[index]] <- index}
  if(all(x!="NO")){
    x <- test_name[i, ]
    Alt_ID[[index]] <- as.character(x[1,1])
    next
  }
  if(any(x!="NO")){
    x <- test_name[i, ]
    Alt_ID[[xindex]] <- as.character(unique(x[-which(is.na(x))])[1,1])}
  
}

test_name <- matrix("a", ncol=1, nrow=length(Alt_ID))
rownames(test_name) <- names(Alt_ID)
colnames(test_name) <- "Alt_ID"
for(i in rownames(test_name)){test_name[i, 1] <- Alt_ID[[i]]}


test_qvale <- test[,c(1:k)*2]
rownames(test_qvale) <- test_name[rownames(test_qvale),]
colnames(test_qvale) <- gsub(".logQ", "", colnames(test_qvale))


##
library(reshape2)
library(ggplot2)
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
pp <- test_qvale
name <- sapply(rownames(pp), function(x){strsplit(x, "_")[[1]][1]})
name <- gsub("[)]", "", gsub("[(]", "", name))
rownames(pp) <- name
# index1 <- toupper(name) %in% gtf_ensembl_gene$gene_name
# name <- sapply(rownames(pp), function(x){strsplit(x, ":")[[1]][1]})
# index2 <- toupper(name) %in% gtf_ensembl_gene$gene_name
# rownames(pp)[index1 | index2]
# pp <- test_qvale[index1 | index2,]
# sub_gene <- gtf_ensembl_gene[gtf_ensembl_gene$gene_name %in% rownames(pp),]
# sub_norCounts <- Rhesus_norCounts[sub_gene$gene_id,-which()]
# rownames()


pp <- melt(t(pp))
colnames(pp) <- c("cluster", "motif", "logQ")

ggplot(pp, aes(x=cluster, y=motif, size=logQ)) + 
  geom_point(colour="#7D53A3", alpha=0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1),
        plot.title = element_text(hjust = 0.5, colour="black"),
        axis.title = element_text(size=16, colour="black"),
        axis.text.x = element_text(size=12, colour="black"),
        axis.text.y = element_text(size=7, colour="black")) +
  scale_size_continuous(breaks = c(1, 2, 3, 4, 5, 6),
                        # range=c(0.5,3.5),
                        labels=c("<1", "2", "3", "4", "5", "6")) +
  labs(size="-logQ") +
  scale_y_discrete(limits=rev(levels(pp$motif))) +
  xlab("") +
  ylab("")

ggsave("revise0530/transcripts_count_res_with_TACO_gtf/module_motif_Qvalue_least_1.pdf",
       width=6, height = 5.7)


###
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")

pp <- test_qvale
name <- sapply(rownames(pp), function(x){strsplit(x, "_")[[1]][1]})
name <- gsub("[)]", "", gsub("[(]", "", name))
name <- c(name, "PTBP2")

vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)
norCounts <- counts(Rhesus_GeneObject, normalized=TRUE)
sub_gene <- gtf_ensembl_gene[gtf_ensembl_gene$gene_name %in% toupper(name) | gtf_ensembl_gene$gene_id %in% name,]
sub_norCounts <- vsd_data[sub_gene$gene_id,]
gene_name = sub_gene$gene_name
gene_name[which(is.na(gene_name))] <- sub_gene$gene_id[which(is.na(gene_name))]
rownames(sub_norCounts) <- gene_name
dim(sub_norCounts)

pp <- melt(sub_norCounts)
colnames(pp) <- c("gene", "sample_id", "expression")
pp$SR <- Rhesus_colData[as.character(pp$sample_id),"SR_merge"]
head(pp)

##

pp$gene <- factor(pp$gene, levels=rownames(sub_norCounts)[order(rowVars(sub_norCounts), decreasing = T)])

region_order <- c("SN","HTHA","MO","PA", "AMY", "HIP", 
                  "PON",  "VTA", 
                  "GP","CC","cortex", "STR","THA", "CB","OB","PIT")
ggplot(pp, aes(x=SR, y=expression, fill=SR, color=SR)) +
  geom_boxplot(outlier.shape=16, outlier.size=0.5,lwd=0.2, alpha=0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black"),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color="black"),
        text = element_text(color="black"),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white")) +
  facet_grid(.~gene, scales = "free_x") +
  coord_flip() +
  guides(fill=F, color=F) +
  xlab("") +
  ylab("expression(VST)") +
  scale_x_discrete(limits=region_order) +
  scale_fill_manual(values = setNames(colorRampPalette(rev(brewer.pal(12,"Paired")))(16), region_order)) +
  scale_color_manual(values = setNames(colorRampPalette(rev(brewer.pal(12,"Paired")))(16), region_order)) +
  ggtitle("")



ggsave("revise0530/transcripts_count_res_with_TACO_gtf/motif_expression.pdf",
       width = 10, height = 4)








