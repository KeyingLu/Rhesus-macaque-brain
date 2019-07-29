
library(miscTools)

setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/gtf_ensembl_trans.RData")
load("human/GTEx_gene_count_res/human_gtf_ensembl_trans.RData")



# ##########   score
# score <- read.table("ANNOTATION/Rhesus_human_trans_similarity_score.txt")
# colnames(score) <- c("Rhesus_transcript_id", "Human_transcript_id", 
#                      "identity_percent", "alignment_length", "mismatchs",
#                      "gap_openings", "qurey_start", "qurey_end", "Subject_start",
#                      "Subject_end", "e_value", "bit_score")
# score <- insertCol(as.matrix(score), 2, 
#                    v=gtf_ensembl_trans[score$Rhesus_transcript_id,]$gene_id, 
#                    cName="Rhesus_gene_id")
# score <- insertCol(as.matrix(score), 4, 
#                    v=human_gtf_ensembl_trans[score[,"Human_transcript_id"],]$gene_id, 
#                    cName="Human_gene_id")
# score <- data.frame(score)
# score <- data.frame(link=paste(score$Rhesus_gene_id, score$Human_gene_id, sep = "_"), score)
# head(score)
# 
# 
# ######### orth
# orth_genes <- read.table("ANNOTATION/humanGRCH38_Rhesus_orth_genes_one2one.txt")
# colnames(orth_genes) <- c("Human_gene_id", "Rhesus_gene_id", "Homologous_type")
# rownames(orth_genes) <- paste(orth_genes$Rhesus_gene_id, orth_genes$Human_gene_id, sep="_")
# head(orth_genes)
# dim(orth_genes) # 19810     3
# 
# 
# orth_score <- score[score$link %in% rownames(orth_genes),]
# head(orth_score)
# dim(orth_score) # 587111     15
# 
# 
# hist(as.numeric(score$identity_percent))
# hist(as.numeric(score$bit_score))
# hist(as.numeric(orth_score$identity_percent))
# hist(as.numeric(orth_score$bit_score))



#######################################################################
#                             human gtf                               #            
#######################################################################
library(devtools)
library(rtracklayer)
human_gtf <- rtracklayer::import("/mnt/data1/Ref/human_hg38/GRCh38_ensembl/Homo_sapiens.GRCh38.91.chr.gtf")
human_gtf <- data.frame(human_gtf)
human_gtf_gene <- human_gtf[human_gtf$type=="gene",c(10,12)]
rownames(human_gtf_gene) <- human_gtf_gene$gene_id
dim(human_gtf_gene)
human_gtf_trans <- human_gtf[human_gtf$type=="transcript",c(10,12,15)]
rownames(human_gtf_trans) <- human_gtf_trans$transcript_id
save(human_gtf_gene, file="ANNOTATION/human_gtf_gene.RData")
save(human_gtf_trans, file="ANNOTATION/human_gtf_trans.RData")



#######################################################################
#                   human and denovo Rhesus ortholog                  #            
#######################################################################
library(reshape2)
setwd("/mnt/data2/Rhesus_brain")

denovo_orth <- read.table("TACO_minExpr_5.5/CPC2_CPAT_inparanoid_try2/human_denovo_coding_Rhesus_ortholog_one2one.txt",
                          header=F, stringsAsFactors = F)[c(2,3,5)]
colnames(denovo_orth) <- c("score", "human", "Rhesus")
human_protein <- sapply(denovo_orth$human, function(x){strsplit(x, "[|]")[[1]][3]})
denovo_orth$human_protein <- human_protein

## 先将蛋白id转成基因id，再利用GRCh38 gtf转化成为gene_symbol
load("ANNOTATION/human_gtf_gene.RData")
human_protein_gene_id <- read.table("TACO_minExpr_5.5/CPC2_CPAT_inparanoid_try2/human_protein_to_gene_id.txt", header=T, stringsAsFactors = F)
human_protein_gene_id_list <- sapply(unique(human_protein_gene_id$From), function(x){human_protein_gene_id$To[human_protein_gene_id$From %in% x]})
human_protein_gene_symbol_list <- lapply(human_protein_gene_id_list, function(x){paste(human_gtf_gene[x, "gene_name"][complete.cases(human_gtf_gene[x, "gene_name"])], collapse = ",")})

##
library(reshape2)
pp <- melt(human_protein_gene_symbol_list)
rownames(pp) <- pp$L1
denovo_orth$human_gene_name <- pp[denovo_orth$human_protein, "value"]
human_gene_name <- names(table(denovo_orth$human_gene_name))[table(denovo_orth$human_gene_name)==1]
denovo_orth <- denovo_orth[denovo_orth$human_gene_name %in% human_gene_name,]
head(denovo_orth)
dim(denovo_orth) # 7688    5



######
setwd("/mnt/data2/Rhesus_brain")
load("ANNOTATION/human_gtf_gene.RData")
load("ANNOTATION/human_gtf_trans.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")
denovo_orth$Rhesus_denovo_trans_id <- gsub(".p[0-9]+", "", as.character(denovo_orth$Rhesus))
denovo_orth <- data.frame(denovo_orth, TACO_gtf_trans[denovo_orth$Rhesus_denovo_trans_id,])
which(denovo_orth$human_gene_name == denovo_orth$ref_gene_name & denovo_orth$shared_splicing=="True")
View(denovo_orth)
head(denovo_orth)
denovo_orth <- denovo_orth[denovo_orth$score>100, ]
dim(denovo_orth) # 7685   38
save(denovo_orth, file="revise0530/de_novo_TACO_minExpr5.5_res/denovo_orth.RData")


###
load("revise0530/de_novo_TACO_minExpr5.5_res/denovo_orth.RData")
View(denovo_orth)

test <- denovo_orth[which(denovo_orth$shared_splicing=="True" & denovo_orth$category_relative_detail=="same_strand"), 
                    c("category_relative_detail", "transcript_id", "ref_gene_id", "ref_gene_name", "human_gene_name", "shared_splicing")]
test <- test[complete.cases(test),]
test <- test[-which(test$ref_gene_id==test$ref_gene_name),]
dim(test) # 3433    6
View(test)

#### Accuracy
length(which(test$ref_gene_name == test$human_gene_name))/dim(test)[1] # 3349




###############################################################
#              write novel&coding information                 #            
###############################################################
library(DESeq2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/denovo_orth.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_TransObject_LRT_for_SR_merge.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")


## coding
coding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="coding" & TACO_gtf_trans$CPAT_potential=="coding"),]
coding_trans <- coding_trans[which(coding_trans$shared_splicing!="True" | is.na(coding_trans$shared_splicing)),]
coding_trans <- coding_trans[, c("transcript_id", "category_relative_detail", "gene_id", "mammals91_score", "amniotes54_score", "CPC2_potential", "CPAT_potential")]
dim(coding_trans) # 27792    7

##
rownames(denovo_orth) <- denovo_orth$Rhesus_denovo_trans_id
coding_trans$inParanoid_score <- NA
coding_trans$Human_protein <- NA
coding_trans$predicted_gene_name <- NA
for(i in 1:dim(coding_trans)[1]){
  transcript_id = coding_trans[i, "transcript_id"]
  if(transcript_id %in% rownames(denovo_orth)){
    coding_trans[i, "inParanoid_score"] <- denovo_orth[transcript_id, "score"]
    coding_trans[i, "Human_protein"] <- denovo_orth[transcript_id, "human"]
    coding_trans[i, "predicted_gene_name"] <- as.character(denovo_orth[transcript_id, "human_gene_name"])
  }
}

write.table(coding_trans, 
            "revise0530/de_novo_TACO_minExpr5.5_res/de_novo_novel_coding_trans_information.txt",
            row.names = F, col.names = T, quote = F, sep="\t")





#######################################################################
#                  novel&coding&orth trans expression                 #            
#######################################################################
library(DESeq2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/denovo_orth.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_TransObject_LRT_for_SR_merge.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")


## coding
coding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="coding" & TACO_gtf_trans$CPAT_potential=="coding"),]
coding_trans <- coding_trans[which(coding_trans$shared_splicing!="True" | is.na(coding_trans$shared_splicing)),]
dim(coding_trans) # 27792    32

## coding&orth
rownames(denovo_orth) <- denovo_orth$Rhesus_denovo_trans_id
coding_orth_trans <- denovo_orth[rownames(denovo_orth) %in% coding_trans$transcript_id,]
coding_orth_trans <- coding_orth_trans[complete.cases(as.character(coding_orth_trans$human_gene_name)),]
dim(coding_orth_trans) # 3947   38
# View(coding_orth_trans)
length(unique(coding_orth_trans$gene_id)) # 3910
sum(is.na(coding_orth_trans$ref_gene_name)) # 222
length(unique(coding_orth_trans$human_gene_name)) # 3947

## coding&orth&intergenic
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/gtf_ensembl_trans.RData")
rownames(denovo_orth) <- denovo_orth$Rhesus_denovo_trans_id
coding_orth_trans <- denovo_orth[rownames(denovo_orth) %in% coding_trans$transcript_id,]
coding_orth_trans <- coding_orth_trans[complete.cases(coding_orth_trans$human_gene_name),]
intergenic_coding_orth_trans <- coding_orth_trans[coding_orth_trans$category_relative_detail=="intergenic",]
dim(intergenic_coding_orth_trans) # 222  38

index <- -which(intergenic_coding_orth_trans$human_gene_name %in% gtf_ensembl_gene$gene_name) # 8
intergenic_coding_orth_trans$human_gene_name[which(intergenic_coding_orth_trans$human_gene_name %in% gtf_ensembl_gene$gene_name)]
# View(intergenic_coding_orth_trans[which(intergenic_coding_orth_trans$human_gene_name %in% gtf_ensembl_gene$gene_name),])
refsame_trans <- intergenic_coding_orth_trans$transcript_id[which(intergenic_coding_orth_trans$human_gene_name %in% gtf_ensembl_gene$gene_name)]
write.table(refsame_trans, "revise0530/transcripts_count_res_with_TACO_gtf/refsame_trans_8.txt",
            row.names = F, col.names = F, sep="\t", quote=F)

View(gtf_ensembl_gene[which(gtf_ensembl_gene$gene_name %in% intergenic_coding_orth_trans$human_gene_name),])
gene_id <- gtf_ensembl_gene$gene_id[which(gtf_ensembl_gene$gene_name %in% intergenic_coding_orth_trans$human_gene_name)]
refsame_reftrans <- gtf_ensembl_trans$transcript_id[which(gtf_ensembl_trans$gene_id %in% gene_id)]
write.table(refsame_reftrans, "revise0530/transcripts_count_res_with_TACO_gtf/refsame_reftrans_9.txt",
            row.names = F, col.names = F, sep="\t", quote=F)


###
score <- read.table("revise0530/transcripts_count_res_with_TACO_gtf/refsame_reftrans_noveltrans_similarity_score_unique.txt",
                    header=F, stringsAsFactors = F)
score$TACO_gene_name <- intergenic_coding_orth_trans[score$V1,"human_gene_name"]
score$ref_gene_name <- gtf_ensembl_trans[score$V2, "gene_name"]
write.table(score[,c(1,2,3,13)], "revise0530/transcripts_count_res_with_TACO_gtf/refsame_reftrans_noveltrans_similarity_score_unique_gene_name.txt",
            row.names = F, col.names = F, quote = F, sep="\t")


###
factor(intergenic_coding_orth_trans$human_gene_name[index])
rownames(intergenic_coding_orth_trans)[intergenic_coding_orth_trans$human_gene_name=="CCKAR"]
intergenic_coding_orth_trans["TU39329", "gene_id"]
df["TU23195",colData$SR_merge=="PIT"]
write.table(intergenic_coding_orth_trans$human_gene_name, 
            "revise0530/transcripts_count_res_with_TACO_gtf/novel_intergenic_coding_trans.txt",
            row.names=F, col.names=F, quote=F, sep="\t")
View(intergenic_coding_orth_trans)


###
object <- Rhesus_TransObject
res <- results(object)
res <- res[complete.cases(res),]
res <- res[order(res$padj), ]
res <- res[res$padj<0.05,]
##
vsd <- varianceStabilizingTransformation(object, blind=FALSE)
vsd_data <- assay(vsd)
df <- vsd_data[rownames(res),]
dim(df) # 35640 408
df <- vsd_data[rownames(vsd_data) %in% rownames(intergenic_coding_orth_trans),]
df <- df[apply(df, 1, sd)!=0, ]
dim(df) # 222 408


genes <- c("ISL2","DLX2", "KCNG3",
           "DLX5",  "FABP7",  "BCL3",
           "B4GALT1", "ADRA1B", "GABRQ", "CCKAR","MRGPRE")
sub_df <- df
rownames(sub_df) <- intergenic_coding_orth_trans[rownames(sub_df), "human_gene_name"]
sub_df <- sub_df[genes,]




##########
library(reshape2)
library(ggplot2)
library(RColorBrewer)

pp <- melt(sub_df)
colnames(pp) <- c("gene", "sample_id", "expression")
head(pp)
pp$SR <- Rhesus_colData[as.character(pp$sample_id), ]$SR_merge

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
  facet_grid(.~gene) +
  coord_flip() +
  guides(fill=F, color=F) +
  xlab("") +
  ylab("expression(VST)") +
  scale_x_discrete(limits=region_order) +
  scale_fill_manual(values = setNames(colorRampPalette(rev(brewer.pal(12,"Paired")))(16), region_order)) +
  scale_color_manual(values = setNames(colorRampPalette(rev(brewer.pal(12,"Paired")))(16), region_order)) +
  ggtitle("")



ggsave("revise0530/transcripts_count_res_with_TACO_gtf/novel_trans_expr_example_violin.pdf",
       width = 10, height = 4)



##############
library(reshape2)
library(ggplot2)
library(RColorBrewer)
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_TransObject_LRT_for_SR_merge.RData")


object <- Rhesus_TransObject
res <- results(object)
res <- res[complete.cases(res),]
res <- res[order(res$padj), ]
res <- res[res$padj<0.05,]
##
vsd <- varianceStabilizingTransformation(object, blind=FALSE)
vsd_data <- assay(vsd)
norCounts <- counts(object, normalized=T)
# df <- norCounts[rownames(res),]
df <- vsd_data[rownames(res),]
dim(df) # 35640 408
df <- df[rownames(df) %in% rownames(intergenic_coding_orth_trans),]
df <- df[apply(df, 1, sd)!=0, ]
dim(df) # 222 408




genes <- c("ZNF816", "ENO1", "IGBP1", "VAPB", "MAD2L1", "MRPL40", "ARF1", "TOMM20")

sub_df <- df
rownames(sub_df) <- intergenic_coding_orth_trans[rownames(sub_df), "human_gene_name"]
sub_df <- sub_df[rownames(sub_df) %in% genes,]
rownames(sub_df) # "TOMM20" "ENO1"   "ARF1"   "VAPB"   "MRPL40" "MAD2L1" "IGBP1"



pp <- melt(sub_df)
colnames(pp) <- c("gene", "sample_id", "expression")
head(pp)
pp$SR <- Rhesus_colData[as.character(pp$sample_id), ]$SR_merge

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
  facet_grid(.~gene) +
  coord_flip() +
  guides(fill=F, color=F) +
  xlab("") +
  ylab("expression(VST)") +
  scale_x_discrete(limits=region_order) +
  scale_fill_manual(values = setNames(colorRampPalette(rev(brewer.pal(12,"Paired")))(16), region_order)) +
  scale_color_manual(values = setNames(colorRampPalette(rev(brewer.pal(12,"Paired")))(16), region_order)) +
  ggtitle("")


ggsave("revise0530/transcripts_count_res_with_TACO_gtf/novel_sameRef_difSite_expr_violin.pdf",
       width = 10, height = 4)





