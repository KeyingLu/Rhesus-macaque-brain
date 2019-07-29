
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
library(rtracklayer)
library(SummarizedExperiment)


setwd("/mnt/data2/Rhesus_brain")
dir.create("revise0530/de_novo_TACO_minExpr5.5_res")



###################################################################
#                          assembly  gtf                          #
###################################################################
library(devtools)
library(rtracklayer)
setwd("/mnt/data2/Rhesus_brain")

TACO_gtf <- rtracklayer::import("TACO_minExpr_5.5/assembly.refcomp.gtf")
TACO_gtf <- as.data.frame(TACO_gtf)
dim(TACO_gtf) # 526588     30
## conservation score
mammals91 <- read.table("Conservation_score/assembly_conservation_scores_91_mammals.tab", stringsAsFactors = F)
colnames(mammals91) <- c("name", "size", "covered", "sum", "mean0", "mean")
trans_id <- sapply(mammals91$name, function(x){strsplit(x, "[|]")[[1]][2]})
trans_id <- sapply(trans_id, function(x){strsplit(x, "[(]")[[1]][1]})
rownames(mammals91) <- trans_id
amniotes54 <- read.table("Conservation_score/assembly_conservation_scores_54_amniotes.tab", stringsAsFactors = F)
colnames(amniotes54) <- c("name", "size", "covered", "sum", "mean0", "mean")
trans_id <- sapply(amniotes54$name, function(x){strsplit(x, "[|]")[[1]][2]})
trans_id <- sapply(trans_id, function(x){strsplit(x, "[(]")[[1]][1]})
rownames(amniotes54) <- trans_id
## transcripts
TACO_gtf_trans <- TACO_gtf[which(TACO_gtf$type=="transcript"), ]
rownames(TACO_gtf_trans) <- TACO_gtf_trans$transcript_id
TACO_gtf_trans$mammals91_score <- mammals91[rownames(TACO_gtf_trans), "mean"]
TACO_gtf_trans$amniotes54_score <- amniotes54[rownames(TACO_gtf_trans), "mean"]
save(TACO_gtf_trans, file="revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")


###
unique(TACO_gtf_trans$category_relative_detail)
table(TACO_gtf_trans$category_relative_detail, TACO_gtf_trans$category_relative)
table(TACO_gtf_trans$category_relative_detail, TACO_gtf_trans$category)
table(TACO_gtf_trans$category_relative_detail, TACO_gtf_trans$annotation)
same_strand <- TACO_gtf_trans[TACO_gtf_trans$category_relative_detail == "same_strand",]
read_through <- TACO_gtf_trans[TACO_gtf_trans$category_relative_detail == "read_through",]
intergenic <- TACO_gtf_trans[TACO_gtf_trans$category_relative_detail == "intergenic",]
intronic_same_strand <- TACO_gtf_trans[TACO_gtf_trans$category_relative_detail == "intronic_same_strand",]
opp_strand <- TACO_gtf_trans[TACO_gtf_trans$category_relative_detail == "opp_strand",]
##
test <- TACO_gtf_trans[which(TACO_gtf_trans$shared_same_strand_bp == TACO_gtf_trans$ref_length | TACO_gtf_trans$shared_opp_strand_bp == TACO_gtf_trans$ref_length),]
test <- TACO_gtf_trans[which(TACO_gtf_trans$shared_same_strand_bp > 0 | TACO_gtf_trans$shared_opp_strand_bp > 0),]
View(read_through)
View(TACO_gtf[10000:10100,])
# TACO_gtf_trans <- TACO_gtf[which(TACO_gtf$type=="transcript"), 
#                            c(1,2,3,15, 18, 25, 26, 29, 23)]




###################################
#### category_relative_detail
library(ggplot2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")

pp <- data.frame(table(TACO_gtf_trans$category_relative_detail), stringsAsFactors = F)
colnames(pp) <- c("category_relative_detail", "count")
pp$category_relative_detail <- as.character(pp$category_relative_detail)
str(pp)
order <- unique(pp$category_relative_detail[order(pp$count, decreasing = T)])
pp$category_relative_detail <- factor(pp$category_relative_detail, levels=order)


### bar 
ggplot(pp, aes(x=category_relative_detail, y=count, fill=category_relative_detail)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", size = 13),
        axis.title = element_text(face="bold"),
        # axis.text.x  = element_text(face="bold", angle = 0, hjust=0.5, vjust=0.5,
        #                             colour = "black", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        axis.text.y  = element_text(face="bold", colour = "black",size=10),
        panel.border = element_rect(size = 1),
        strip.background = element_rect(fill="white", colour = "black", size=1),
        legend.position = c(0.7,0.7),
        legend.title = element_text(face="bold")) +
  scale_fill_manual(values = c("#803C20", "#B33B49", "#E38561", "#A4C98A",
                               "#639FC3", "#5784B4", "#442478", "#898989")) +
  geom_text(aes(label=count), vjust=-0.2) +
  xlab("") +
  ylab("") +
  ggtitle("")

ggsave("revise0530/de_novo_TACO_minExpr5.5_res/assembly_transcripts_category_relative_detail_bar.pdf",
       width = 5, height = 4.8)


### 饼环 
ggplot(pp, aes(x="", y=count, fill=category_relative_detail)) +
  geom_bar(stat="identity",width = 0.2) +
  coord_polar(theta = "y") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", size = 13),
        axis.title = element_text(face="bold"),
        # axis.text.x  = element_text(face="bold", angle = 0, hjust=0.5, vjust=0.5,
        #                             colour = "black", size = 12),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.text.y = element_blank(),
        axis.text.y  = element_text(face="bold", colour = "black",size=10),
        panel.border = element_rect(size = 1),
        strip.background = element_rect(fill="white", colour = "black", size=1),
        # legend.position = c(0.7,0.7),
        legend.title = element_text(face="bold")) +
  # scale_fill_manual(values = c("#803C20", "#B33B49", "#E38561", "#A4C98A", 
  #                              "#639FC3", "#5784B4", "#442478", "#898989")) +
  scale_fill_manual(values = brewer.pal(11, "RdYlGn")[]) +
  # geom_text(aes(label=count), vjust=-0.2) +
  xlab("") +
  ylab("") +
  ggtitle("")

ggsave("revise0530/de_novo_TACO_minExpr5.5_res/assembly_transcripts_category_relative_detail_Pie_Donut.pdf",
       width = 5, height = 4.8)



########################################################################
#                         same_strand detail                           #
########################################################################
same_strand <- TACO_gtf_trans[TACO_gtf_trans$category_relative_detail == "same_strand",]

pp <- data.frame(table(same_strand$shared_splicing))

pdf("revise0530/de_novo_TACO_minExpr5.5_res/assembly_transcripts_same_strand_detail_pie.pdf",
    width = 6, height = 5)
pie(pp$Freq,
    labels=paste(c("Partial Match\n", "Complete Macth\n"), "(", pp$Freq, ")", sep=""),
    border=F,
    clockwise = T,
    init.angle = -35,
    main = "same_strand",
    col = c("#B33B49", "#5784B4"))
dev.off()




########################################################################
#                 coding potential (CPC2 & CPAT)                       #
########################################################################
library(devtools)
library(rtracklayer)
library(VennDiagram)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")

### intersect trans (denovo and ref)
TACO_gtf_trans_refPart <- TACO_gtf_trans[which(TACO_gtf_trans$category_relative_detail=="same_strand" & TACO_gtf_trans$shared_splicing=="True"),]
# TACO_gtf_trans_refPart <- TACO_gtf_trans[which(TACO_gtf_trans$category_relative_detail=="same_strand" | TACO_gtf_trans$category_relative_detail=="read_through"),]
# TACO_gtf_trans_refPart <- TACO_gtf_trans[which(TACO_gtf_trans$category_relative_detail=="same_strand"),]
# TACO_gtf_trans_refPart <- TACO_gtf_trans[which(TACO_gtf_trans$shared_splicing=="True"),]
dim(TACO_gtf_trans_refPart) # 6714   30
unique(TACO_gtf_trans_refPart$ref_gene_type)
TACO_gtf_trans_refPart$ref_gene_type[TACO_gtf_trans_refPart$ref_gene_type=="protein_coding"] <- "coding"
TACO_gtf_trans_refPart$ref_gene_type[TACO_gtf_trans_refPart$ref_gene_type!="coding"] <- "noncoding"

area1 = dim(TACO_gtf_trans)[1]
area2 = 55057
cross.area = dim(TACO_gtf_trans_refPart)[1]

pdf("revise0530/de_novo_TACO_minExpr5.5_res/TACO_minExpr5.5_assembly_knowned_trans_inersect_venn.pdf", width=5, height = 5)
draw.pairwise.venn(area1, area2, cross.area,
                   category = c("assembly_trans", "ref_trans"),
                   fill=c("#3F60AC", "#8B2052"),
                   col=c("#3F60AC", "#8B2052"),
                   alpha=0.5,
                   print.mode = "raw",
                   cat.pos = c(210,150),
                   lwd = c(0,0),
                   scaled = T,
                   cat.cex=2,
                   main= "coding",
                   cat.dist = c(0.05,0.05),
                   margin=0.05,
                   ind = TRUE,
                   cex = 1.5)
dev.off()


################
### CPC2 results
################
CPC2_res <- read.table("TACO_minExpr_5.5/TACO_minExpr_5.5_gtf_CPC2_out.txt", header=F, sep="\t")
rownames(CPC2_res) <- CPC2_res[,1]
View(CPC2_res[1:10,])
CPC2_res_trans <- CPC2_res[,c("V1", "V8")]
colnames(CPC2_res_trans) <- c("transcript_id", "transcript_type")
head(CPC2_res_trans)
dim(CPC2_res_trans) # 53661     2


################
### CPAT results
################
CPAT_res <- read.table("TACO_minExpr_5.5/TACO_minExpr_5.5_gtf_CPAT_out.txt", header=T, row.names = 1)
CPAT_res$transcript_type <- CPAT_res$coding_prob
CPAT_res$transcript_type[CPAT_res$coding_prob > 0.6067] <- "coding"
CPAT_res$transcript_type[CPAT_res$coding_prob <= 0.6067] <- "noncoding"
CPAT_res_trans <- CPAT_res[, "transcript_type", drop=F]
table(CPAT_res_trans$transcript_type)


################
### combine ref trans
################
TACO_gtf_trans_refPart$CPC2_potential <- CPC2_res_trans[rownames(TACO_gtf_trans_refPart),2]
TACO_gtf_trans_refPart$CPAT_potential <- CPAT_res_trans[rownames(TACO_gtf_trans_refPart),1]
head(TACO_gtf_trans_refPart)
str(TACO_gtf_trans_refPart)
dim(TACO_gtf_trans_refPart)

sum(TACO_gtf_trans_refPart$CPC2_potential==TACO_gtf_trans_refPart$ref_gene_type)/dim(TACO_gtf_trans_refPart)[1]
sum(TACO_gtf_trans_refPart$CPAT_potential==TACO_gtf_trans_refPart$ref_gene_type)/dim(TACO_gtf_trans_refPart)[1]

#### CPC2 : plot Venn
library(VennDiagram)

### coding
predicted_coding <- rownames(TACO_gtf_trans_refPart)[TACO_gtf_trans_refPart$CPC2_potential=="coding"]
ref_coding <- rownames(TACO_gtf_trans_refPart)[TACO_gtf_trans_refPart$ref_gene_type=="coding"]

area1 = length(predicted_coding)
area2 = length(ref_coding)
cross.area = length(intersect(predicted_coding, ref_coding))

pdf("gtf_for_assembled/TACO_minExpr5.5_assembly_predicted_ref_inersect_for_coding_protein_CPC2.pdf",
    width = 5, height = 5)
draw.pairwise.venn(area1, area2, cross.area,
                   category = c("predicted", "ref"),
                   fill=c("#3F60AC", "#8B2052"),
                   col=c("#3F60AC", "#8B2052"),
                   alpha=0.5,
                   print.mode = "raw",
                   cat.pos = c(210,150),
                   lwd = c(0,0),
                   scaled = T,
                   cat.cex=2,
                   main= "coding",
                   cat.dist = c(0.05,0.05),
                   margin=0.05,
                   ind = TRUE,
                   cex = 1.5)
dev.off()

## nocoding
predicted_noncoding <- rownames(TACO_gtf_trans_refPart)[TACO_gtf_trans_refPart$CPC2_potential=="noncoding"]
ref_noncoding <- rownames(TACO_gtf_trans_refPart)[TACO_gtf_trans_refPart$ref_gene_type=="noncoding"]

area1 = length(predicted_noncoding)
area2 = length(ref_noncoding)
cross.area = length(intersect(predicted_noncoding, ref_noncoding))

pdf("gtf_for_assembled/TACO_minExpr5.5_assembly_predicted_ref_inersect_for_noncoding_CPC2.pdf",
    width = 5, height = 5)
draw.pairwise.venn(area1, area2, cross.area,
                   category = c("predicted", "ref"),
                   fill=c("#3F60AC", "#8B2052"),
                   col=c("#3F60AC", "#8B2052"),
                   alpha=0.5,
                   print.mode = "raw",
                   cat.pos = c(210,150),
                   lwd = c(0,0),
                   scaled = T,
                   cat.cex=2,
                   main= "coding",
                   cat.dist = c(0.05,0.05),
                   margin=0.05,
                   ind = TRUE,
                   cex = 1.5)
dev.off()


#### CPAT : plot Venn
library(VennDiagram)

### coding
predicted_coding <- rownames(TACO_gtf_trans_refPart)[TACO_gtf_trans_refPart$CPAT_potential=="coding"]
ref_coding <- rownames(TACO_gtf_trans_refPart)[TACO_gtf_trans_refPart$ref_gene_type=="coding"]

area1 = length(predicted_coding)
area2 = length(ref_coding)
cross.area = length(intersect(predicted_coding, ref_coding))

pdf("gtf_for_assembled/TACO_minExpr5.5_assembly_predicted_ref_inersect_for_coding_protein_CPAT.pdf",
    width = 5, height = 5)
draw.pairwise.venn(area1, area2, cross.area,
                   category = c("predicted", "ref"),
                   fill=c("#3F60AC", "#8B2052"),
                   col=c("#3F60AC", "#8B2052"),
                   alpha=0.5,
                   print.mode = "raw",
                   cat.pos = c(210,150),
                   lwd = c(0,0),
                   scaled = T,
                   cat.cex=2,
                   main= "coding",
                   cat.dist = c(0.05,0.05),
                   margin=0.05,
                   ind = TRUE,
                   cex = 1.5)
dev.off()

## nocoding
predicted_noncoding <- rownames(TACO_gtf_trans_refPart)[TACO_gtf_trans_refPart$CPAT_potential=="noncoding"]
ref_noncoding <- rownames(TACO_gtf_trans_refPart)[TACO_gtf_trans_refPart$ref_gene_type=="noncoding"]

area1 = length(predicted_noncoding)
area2 = length(ref_noncoding)
cross.area = length(intersect(predicted_noncoding, ref_noncoding))

pdf("gtf_for_assembled/TACO_minExpr5.5_assembly_predicted_ref_inersect_for_noncoding_CPAT.pdf",
    width = 5, height = 5)
draw.pairwise.venn(area1, area2, cross.area,
                   category = c("predicted", "ref"),
                   fill=c("#3F60AC", "#8B2052"),
                   col=c("#3F60AC", "#8B2052"),
                   alpha=0.5,
                   print.mode = "raw",
                   cat.pos = c(210,150),
                   lwd = c(0,0),
                   scaled = T,
                   cat.cex=2,
                   main= "coding",
                   cat.dist = c(0.05,0.05),
                   margin=0.05,
                   ind = TRUE,
                   cex = 1.5)
dev.off()



########################################################################
#                             CPC2 CPAT ref                            #
########################################################################
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")
TACO_gtf_trans$CPC2_potential <- CPC2_res_trans[rownames(TACO_gtf_trans),2]
TACO_gtf_trans$CPAT_potential <- CPAT_res_trans[rownames(TACO_gtf_trans),1]
TACO_gtf_trans[rownames(TACO_gtf_trans_refPart), ] <- TACO_gtf_trans_refPart[, colnames(TACO_gtf_trans)]
save(TACO_gtf_trans, file="revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")

##
library(venn)
pdf("gtf_for_assembled/TACO_minExpr_5.5/CPC2_CPAT_coding_ref.pdf")
venn(list(ref_noncoding=rownames(TACO_gtf_trans)[which(TACO_gtf_trans$ref_gene_type=="noncoding")],
          #  CPC2_noncoding=rownames(TACO_gtf_trans)[TACO_gtf_trans$CPC2_potential=="noncoding"],
          # CPAT_noncoding=rownames(TACO_gtf_trans)[TACO_gtf_trans$CPAT_potential=="noncoding"],
          CPC2_coding=rownames(TACO_gtf_trans)[TACO_gtf_trans$CPC2_potential=="coding"],
          CPAT_coding=rownames(TACO_gtf_trans)[TACO_gtf_trans$CPAT_potential=="coding"],
          ref_coding=rownames(TACO_gtf_trans)[TACO_gtf_trans$ref_gene_type=="coding"]),
     # size = 30,
     col=c("#778899","#3F60AC", "#8B2052","#778899"),
     zcolor = c("#778899","#3F60AC", "#8B2052","#778899"),
     cexil = 1.5,	# Character expansion for the intersection labels
     cexsn = 1.5 #Character expansion for the set names,
)
dev.off()



pdf("gtf_for_assembled/TACO_minExpr_5.5/CPC2_CPAT_coding_and_all_trans.pdf")
venn(list(# ref_noncoding=rownames(TACO_gtf_trans)[which(TACO_gtf_trans$ref_gene_type=="noncoding")],
  #  CPC2_noncoding=rownames(TACO_gtf_trans)[TACO_gtf_trans$CPC2_potential=="noncoding"],
  # CPAT_noncoding=rownames(TACO_gtf_trans)[TACO_gtf_trans$CPAT_potential=="noncoding"]
  All_trans=rownames(TACO_gtf_trans),
  CPC2_coding=rownames(TACO_gtf_trans)[TACO_gtf_trans$CPC2_potential=="coding"],
  CPAT_coding=rownames(TACO_gtf_trans)[TACO_gtf_trans$CPAT_potential=="coding"]
  # ref_coding=rownames(TACO_gtf_trans)[TACO_gtf_trans$ref_gene_type=="coding"]
),
# size = 30,
borders=F,
col=c("#778899","#3F60AC", "#8B2052"),
zcolor = c("#778899","#3F60AC", "#8B2052"),
cexil = 1.5,	# Character expansion for the intersection labels
cexsn = 1.5 #Character expansion for the set names,
)
dev.off()



########################################################################
#                           novel trans location                       #
########################################################################
# TACO_gtf_trans_novel <- TACO_gtf_trans[which(is.na(TACO_gtf_trans$ref_gene_id)),]
# TACO_gtf_trans_novel_CPAT_CPCP2 <- TACO_gtf_trans_novel[TACO_gtf_trans_novel$CPC2_potential==TACO_gtf_trans_novel$CPAT_potential,]
# 
# write.table(TACO_gtf_trans_novel_CPAT_CPCP2[,c("seqnames", "start", "end", "transcript_id", "CPC2_potential")],
#             "gtf_for_assembled/TACO_minExpr_5.5/TACO_gtf_trans_novel_CPAT_CPCP2.bed",
#             col.names = F, row.names = F, quote = F, sep="\t")
# 
# ### shell
# # /mnt/xdlab1/tools/bedtools2/bedtools2/bin/bedtools intersect 
# # -a TACO_gtf_trans_novel_CPAT_CPCP2.bed
# # -b /mnt/data2/Rhesus_brain/ANNOTATION/Rhesus_gene.bed
# # -wao > TACO_gtf_trans_novel_CPAT_CPCP2_intersect_ref_gene.bed
# 
# ## 
# data <- read.table("gtf_for_assembled/TACO_minExpr_5.5/TACO_gtf_trans_novel_CPAT_CPCP2_intersect_ref_gene.bed")
# head(data)
# 
# intergenic_trans <- data[,c(4,5)][data[,11]==0,]
# colnames(intergenic_trans) <- c("transcript_id", "biotype")
# intergenic_trans = intergenic_trans[!duplicated(intergenic_trans$transcript_id), ]
# dim(intergenic_trans)
# genic_trans <- data[,c(4,5)][data[,11]>0,]
# # genic_trans <- data[data[,11]>0,]
# colnames(genic_trans) <- c("transcript_id", "biotype")
# genic_trans = genic_trans[!duplicated(genic_trans$transcript_id), ]
# dim(genic_trans) # 824   2
# 
# 
# pp <- rbind(table(intergenic_trans$biotype),table(genic_trans$biotype))
# rownames(pp) <- c("intergenic", "overlap_ref_gene")
# 
# pdf("revise0530/de_novo_TACO_minExpr5.5_res/novel_coding_intergenic_genic.pdf",
#     width = 5, height = 5)
# p <- pp[,1]
# percent <- paste(round(p/sum(p)*100,2), "%", sep = "")
# pie(pp[,1],
#     labels=paste(names(p), "(",percent , ")",sep=""),
#     border=F,
#     clockwise = T,
#     init.angle = -15,
#     main = paste("coding", "(",sum(p), ")",sep=""),
#     col = c("#B0ABD4", "#FAB895"))
# dev.off()
# 
# 
# pdf("revise0530/de_novo_TACO_minExpr5.5_res/novel_noncoding_intergenic_genic.pdf",
#     width = 5, height = 5)
# p <- pp[,2]
# percent <- paste(round(p/sum(p)*100,2), "%", sep = "")
# pie(pp[,1],
#     labels=paste(names(p), "(",percent , ")",sep=""),
#     border=F,
#     clockwise = T,
#     init.angle = -15,
#     main = paste("noncoding", "(",sum(p), ")",sep=""),
#     col = c("#B0ABD4", "#FAB895"))
# dev.off()



#######################################################################
#         novel(coding/noncoding) trans conservation score            #            
#######################################################################
library(DESeq2)
library(BiocParallel)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")


## coding
coding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="coding" & TACO_gtf_trans$CPAT_potential=="coding"),]
coding_trans <- coding_trans[which(coding_trans$shared_splicing!="True" | is.na(coding_trans$shared_splicing)),]
dim(coding_trans) # 27792    34


test1 <- hist(coding_trans$mammals91_score)
test2 <- hist(coding_trans$amniotes54_score)
cor(coding_trans$mammals91_score, coding_trans$amniotes54_score)



## noncoding
noncoding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="noncoding" | TACO_gtf_trans$CPAT_potential=="noncoding"),]
noncoding_trans <- noncoding_trans[which(noncoding_trans$shared_splicing!="True" | is.na(noncoding_trans$shared_splicing)),]
dim(noncoding_trans) # 18822    34

pp <- rbind(data.frame(coding_trans, coding_type="coding"),
            data.frame(noncoding_trans, coding_type="noncoding"))

ggplot(pp, aes(x=mammals91_score, fill=coding_type)) +
  geom_density(stat = "density", size=0.5, alpha=0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1, color = "black"),
        axis.text = element_text(size = 15, colour = "black"),
        plot.title = element_text(hjust = 0.5,size = 15),
        axis.title.x =element_text(size=15), 
        axis.title.y=element_text(size=15))+
  scale_fill_manual(values = c("#F25F5C", "#247BA0")) +
  xlab("conservation score") +
  ggtitle("mammals91_score") 

ggsave("revise0530/de_novo_TACO_minExpr5.5_res/mammals91_score_novel_trans.pdf", 
       width = 5, height = 4)



ggplot(pp, aes(x=amniotes54_score, fill=coding_type)) +
  geom_density(stat = "density", size=0.5, alpha=0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1, color = "black"),
        axis.text = element_text(size = 15, colour = "black"),
        plot.title = element_text(hjust = 0.5,size = 15),
        axis.title.x =element_text(size=15), 
        axis.title.y=element_text(size=15))+
  scale_fill_manual(values = c("#F25F5C", "#247BA0")) +
  xlab("conservation score") +
  ggtitle("amniotes54_score ") 

ggsave("revise0530/de_novo_TACO_minExpr5.5_res/amniotes54_score_novel_trans.pdf", 
       width = 5, height = 4)



#######################################################################
#               novel(coding/noncoding) trans expression              #            
#######################################################################
library(DESeq2)
library(BiocParallel)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_DESeq2_object.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")
object <- Rhesus_TransObject

## coding
coding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="coding" & TACO_gtf_trans$CPAT_potential=="coding"),]
coding_trans <- coding_trans[which(coding_trans$shared_splicing!="True" | is.na(coding_trans$shared_splicing)),]
dim(coding_trans) # 27792    32

## noncoding
noncoding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="noncoding" | TACO_gtf_trans$CPAT_potential=="noncoding"),]
noncoding_trans <- noncoding_trans[which(noncoding_trans$shared_splicing!="True" | is.na(noncoding_trans$shared_splicing)),]
dim(noncoding_trans) # 18822    32

#### LRT for SR_merge
Rhesus_TransObject <- DESeq(Rhesus_TransObject, 
                            test = "LRT",
                            reduced = ~ID,
                            parallel = T, 
                            BPPARAM=MulticoreParam(40))
save(Rhesus_TransObject, 
     file="revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_TransObject_LRT_for_SR_merge.RData")




############
#   coding
############
res <- results(Rhesus_TransObject)
res <- res[rownames(res) %in% rownames(coding_trans), ]
dim(res) # 26766     6






#### gene cluster
library(ConsensusClusterPlus)
ttop_df <- t(top_df)
ttop_df = scale(ttop_df)
title="revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_coding_anovaTop5000_trans"
novel_coding_trans = ConsensusClusterPlus(ttop_df,maxK=10,reps=100,pItem=0.8,pFeature=1,
                                          title=title,clusterAlg="km",
                                          innerLinkage = "ward.D2",
                                          finalLinkage = "ward.D2",
                                          seed = 1,
                                          distance="euclidean",
                                          plot="png")
# cal <- calcICL(DIfGenes_Results_SR_merge)
# cal$clusterConsensus
save(novel_coding_trans,
     file="revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_coding_anovaTop5000_trans.RData")

###
load("revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_coding_anovaTop5000_trans.RData")
k=7
res <- novel_coding_trans[[k]]
geneCluster <- data.frame(geneCluster=factor(sort(res$consensusClass)))


colInfo <- Rhesus_colData[colnames(df),c("SR_merge", "ID")]
head(colInfo)
str(colInfo)

ann_colors <- list(SR_merge = c(anno_colors$SR,anno_colors$merge)[unique(colInfo$SR_merge)],
                   ID = anno_colors$ID)
dfheat <- df[rownames(geneCluster), ]
# dfheat <- df[rownames(anova_res), ]
# dfheat <- df


##
library(dendsort)
library(pheatmap)
library(RColorBrewer)
method1 = c("euclidean", "maximum", "manhattan", "canberra", "minkowski")
method2 = c("ward.D", "ward.D2", "single",  "complete", "average")
method3 = c("average","min")
for(i in method1){
  for(j in method2){
    for(z in method3){
      i = "canberra"
      j = "average"
      z = "min"
      filename= paste("revise0530/transcripts_count_res_with_TACO_gtf/AnovaTop5000_novel_coding_trans_profile_", i, "_", j, "_",z,".pdf", sep="")
      title = paste("novel&coding trans expression\ndistance:", i, " hclust:",j, " type:",z, sep = "")
      pheatmap(
        dfheat,
        col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
        breaks = c(-3,seq(-2,2,length=20),3),
        # legend_breaks = c(-3,-1,1,3),
        filename = filename,
        scale="row",
        main = title,
        # width = 10,
        # height = 10,
        fontsize=6,
        # clustering_method = "complete",
        border_color=NA,
        # fontsize_col = 0.5,
        annotation_col = colInfo,
        annotation_colors = ann_colors,
        gaps_row = cumsum(table(geneCluster$geneCluster)),
        annotation_legend = T,
        show_rownames = F,
        show_colnames = F,
        cluster_rows = F,
        cluster_cols = as.hclust(dendsort(hclust(dist(t(dfheat), method=i), method = j), type=z))
      )
    }
  }
}




#### expression 
mean <- sapply(unique(testInfo), function(x){rowMeans(top_df[,testInfo==x])})
mean <- data.frame(gene_name=gtf_ensembl_gene[rownames(mean), 2], kruskal_res[1:5000,], mean)
mean <- data.frame(gene_id=rownames(mean), mean)
mean <- data.frame(geneCluster, mean)
write.table(mean, "revise0530/gene_count_res/kruskal_for_region_gene_expression.txt", 
            row.names = F, col.names = T, quote = F, sep="\t")


###################################
#   noncoding
###################################
vsd <- varianceStabilizingTransformation(object, blind=FALSE)
vsd_data <- assay(vsd)
vsd_data_ex1 <- vsd_data[as.character(noncoding_trans$transcript_id), 
                         -which(colnames(vsd_data)=="Rhesus_6_36")]
df <- vsd_data_ex1
dim(df) #  18822   407
df <- df[apply(df, 1, sd)!=0, ]
dim(df) # 18782   407


#### anova
colData <- Rhesus_colData[colnames(df),]
testInfo <- colData$SR_merge
pvalue <- apply(df, 1, function(x){d <- data.frame(x,testInfo);summary(aov(x~testInfo, data=d))[[1]][1,5]})
padj <- p.adjust(pvalue, method = "fdr", n=length(pvalue))
anova_res <- data.frame(pvalue, padj)
anova_res <- anova_res[anova_res$padj<0.05,]
anova_res <- anova_res[order(anova_res$padj),]
dim(anova_res) #  12179     2
head(anova_res)
noncoding_anova_res <- anova_res

top_df <- df[rownames(anova_res)[1:5000],]

#### gene cluster
library(ConsensusClusterPlus)
ttop_df <- t(top_df)
ttop_df = scale(ttop_df)
title="revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_noncoding_anovaTop5000_trans"
novel_noncoding_trans = ConsensusClusterPlus(ttop_df,maxK=10,reps=100,pItem=0.8,pFeature=1,
                                             title=title,clusterAlg="km",
                                             innerLinkage = "ward.D2",
                                             finalLinkage = "ward.D2",
                                             seed = 1,
                                             distance="euclidean",
                                             plot="png")
# cal <- calcICL(DIfGenes_Results_SR_merge)
# cal$clusterConsensus
save(novel_noncoding_trans,
     file="revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_noncoding_anovaTop5000_trans.RData")

###
load("revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_noncoding_anovaTop5000_trans.RData")
k=8
res <- novel_noncoding_trans[[k]]
geneCluster <- data.frame(geneCluster=factor(sort(res$consensusClass)))


colInfo <- Rhesus_colData[colnames(df),c("SR_merge", "ID")]
head(colInfo)
str(colInfo)

ann_colors <- list(SR_merge = c(anno_colors$SR,anno_colors$merge)[unique(colInfo$SR_merge)],
                   ID = anno_colors$ID)
dfheat <- df[rownames(geneCluster), ]
# dfheat <- df[rownames(anova_res), ]
# dfheat <- df


##
library(dendsort)
library(pheatmap)
library(RColorBrewer)
method1 = c("euclidean", "maximum", "manhattan", "canberra", "minkowski")
method2 = c("ward.D", "ward.D2", "single",  "complete", "average")
method3 = c("average","min")
for(i in method1){
  for(j in method2){
    for(z in method3){
      i = "canberra"
      j = "ward.D"
      z = "min"
      filename= paste("revise0530/transcripts_count_res_with_TACO_gtf/AnovaTop5000_novel_noncoding_trans_profile_", i, "_", j, "_",z,".pdf", sep="")
      title = paste("novel&noncoding trans expression\ndistance:", i, " hclust:",j, " type:",z, sep = "")
      pheatmap(
        dfheat,
        col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
        breaks = c(-3,seq(-2,2,length=20),3),
        # legend_breaks = c(-3,-1,1,3),
        filename = filename,
        scale="row",
        main = title,
        # width = 10,
        # height = 10,
        fontsize=6,
        # clustering_method = "complete",
        border_color=NA,
        # fontsize_col = 0.5,
        annotation_col = colInfo,
        annotation_colors = ann_colors,
        gaps_row = cumsum(table(geneCluster$geneCluster)),
        annotation_legend = T,
        show_rownames = F,
        show_colnames = F,
        cluster_rows = F,
        cluster_cols = as.hclust(dendsort(hclust(dist(t(dfheat), method=i), method = j), type=z))
      )
    }
  }
}



#### expression 
mean <- sapply(unique(testInfo), function(x){rowMeans(top_df[,testInfo==x])})
mean <- data.frame(gene_name=gtf_ensembl_gene[rownames(mean), 2], kruskal_res[1:5000,], mean)
mean <- data.frame(gene_id=rownames(mean), mean)
mean <- data.frame(geneCluster, mean)
write.table(mean, "revise0530/gene_count_res/kruskal_for_region_gene_expression.txt", 
            row.names = F, col.names = T, quote = F, sep="\t")







#######################################################################
#         intergenic novel(coding/noncoding) trans expression         #            
#######################################################################

### anova (region)
library(DESeq2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_DESeq2_object.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")
object <- Rhesus_TransObject

## intergenic coding
coding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="coding" & TACO_gtf_trans$CPAT_potential=="coding"),]
coding_trans <- coding_trans[which(coding_trans$category_relative_detail=="intergenic"),]
dim(coding_trans) # 1701   32

## intergenic noncoding
noncoding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="noncoding" | TACO_gtf_trans$CPAT_potential=="noncoding"),]
noncoding_trans <- noncoding_trans[which(noncoding_trans$category_relative_detail=="intergenic"),]
dim(noncoding_trans) # 8002   32


####################
#  intergenic coding
####################
vsd <- varianceStabilizingTransformation(object, blind=FALSE)
vsd_data <- assay(vsd)
vsd_data_ex1 <- vsd_data[as.character(coding_trans$transcript_id), 
                         -which(colnames(vsd_data)=="Rhesus_6_36")]
df <- vsd_data_ex1
dim(df) # 1701  407
df <- df[apply(df, 1, sd)!=0, ]
dim(df) # 1699  407

#### anova
colData <- Rhesus_colData[colnames(df),]
testInfo <- colData$SR_merge
pvalue <- apply(df, 1, function(x){d <- data.frame(x,testInfo);summary(aov(x~testInfo, data=d))[[1]][1,5]})
padj <- p.adjust(pvalue, method = "fdr", n=length(pvalue))
anova_res <- data.frame(pvalue, padj)
anova_res <- anova_res[anova_res$padj<0.05,]
anova_res <- anova_res[order(anova_res$padj),]
dim(anova_res) # 1248    2
head(anova_res)
coding_anova_res <- anova_res

### mean
mean <- sapply(unique(testInfo), function(x){rowMeans(df[,testInfo==x])})
head(mean)
region_high_expression_list <- list()
for(region in colnames(mean)){
  fc <- mean[, region]/mean
  trans <- rownames(fc)[apply(fc, 1, function(x){all(x>=1) & any(x>1.5)})]
  trans <- intersect(rownames(coding_anova_res), trans)
  if(length(trans)==0){trans <- rownames(fc)[apply(fc, 1, function(x){all(x>=1) & any(x>1.4)})];
  trans <- intersect(rownames(coding_anova_res), trans)}
  if(length(trans)==0){trans <- rownames(fc)[apply(fc, 1, function(x){all(x>=1) & any(x>1.3)})];
  trans <- intersect(rownames(coding_anova_res), trans)}
  region_high_expression_list[[region]] <- trans
}
lengths(region_high_expression_list)



##########
library(ggplot2)
library(reshape2)
library(RColorBrewer)
trans <- unlist(region_high_expression_list)
sub_df <- df[trans,]
pp <- melt(sub_df)
colnames(pp) <- c("gene", "sample_id", "expression")
head(pp)
pp$SR <- colData[as.character(pp$sample_id), ]$SR_merge


ggplot(pp, aes(x=SR, y=expression, fill=SR, color=SR)) +
  geom_violin() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        text = element_text(face="bold"),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill="white", color = "white")) +
  facet_grid(.~gene, scales="free_x") +
  coord_flip() +
  guides(fill=F, color=F) +
  xlab("") +
  ylab("expression(VST)") +
  scale_x_discrete(limits=rev(names(region_high_expression_list))) +
  scale_fill_manual(values = colorRampPalette(rev(brewer.pal(11,"Paired")))(16)) +
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(11,"Paired")))(16)) +
  ggtitle("")


ggsave("revise0530/transcripts_count_res_with_TACO_gtf/novel_coding_trans_region_high_expr_violin.pdf",
       width = 16, height = 4)







#### gene cluster
library(ConsensusClusterPlus)
top_df <- df[rownames(anova_res),]
ttop_df <- t(top_df)
ttop_df = scale(ttop_df)
title="revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_intergenic_coding_anova_trans"
intergenic_coding_trans = ConsensusClusterPlus(ttop_df,maxK=10,reps=100,pItem=0.8,pFeature=1,
                                               title=title,clusterAlg="km",
                                               innerLinkage = "ward.D2",
                                               finalLinkage = "ward.D2",
                                               seed = 1,
                                               distance="euclidean",
                                               plot="png")
# cal <- calcICL(DIfGenes_Results_SR_merge)
# cal$clusterConsensus
save(intergenic_coding_trans,
     file="revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_intergenic_coding_anova_trans.RData")

###
load("revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_intergenic_coding_anova_trans.RData")
k=5
res <- intergenic_coding_trans[[k]]
geneCluster <- data.frame(geneCluster=factor(sort(res$consensusClass)))


colInfo <- Rhesus_colData[colnames(df),c("SR_merge", "ID")]
head(colInfo)
str(colInfo)

ann_colors <- list(SR_merge = c(anno_colors$SR,anno_colors$merge)[unique(colInfo$SR_merge)],
                   ID = anno_colors$ID)
dfheat <- df[rownames(geneCluster), ]
# dfheat <- df[rownames(anova_res), ]
# dfheat <- df


##
library(dendsort)
library(pheatmap)
library(RColorBrewer)
method1 = c("euclidean", "maximum", "manhattan", "canberra", "minkowski")
method2 = c("ward.D", "ward.D2", "single",  "complete", "average")
method3 = c("average","min")
for(i in method1){
  for(j in method2){
    for(z in method3){
      # i = "canberra"
      # j = "single"
      # z = "min"
      filename= paste("revise0530/transcripts_count_res_with_TACO_gtf/Anova_novel_intergenic_coding_trans_profile_", i, "_", j, "_",z,".pdf", sep="")
      title = paste("novel&coding&intergenic trans expression\ndistance:", i, " hclust:",j, " type:",z, sep = "")
      pheatmap(
        dfheat,
        col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
        breaks = c(-3,seq(-2,2,length=20),3),
        # legend_breaks = c(-3,-1,1,3),
        filename = filename,
        scale="row",
        main = title,
        # width = 10,
        # height = 10,
        fontsize=6,
        # clustering_method = "complete",
        border_color=NA,
        # fontsize_col = 0.5,
        annotation_col = colInfo,
        annotation_colors = ann_colors,
        gaps_row = cumsum(table(geneCluster$geneCluster)),
        annotation_legend = T,
        show_rownames = F,
        show_colnames = F,
        cluster_rows = F,
        cluster_cols = as.hclust(dendsort(hclust(dist(t(dfheat), method=i), method = j), type=z))
      )
    }
  }
}



###################################
#  intergenic noncoding
###################################
vsd <- varianceStabilizingTransformation(object, blind=FALSE)
vsd_data <- assay(vsd)
vsd_data_ex1 <- vsd_data[as.character(noncoding_trans$transcript_id), 
                         -which(colnames(vsd_data)=="Rhesus_6_36")]
df <- vsd_data_ex1
dim(df) #  8002  407
df <- df[apply(df, 1, sd)!=0, ]
dim(df) # 7999  407


#### anova
colData <- Rhesus_colData[colnames(df),]
testInfo <- colData$SR_merge
pvalue <- apply(df, 1, function(x){d <- data.frame(x,testInfo);summary(aov(x~testInfo, data=d))[[1]][1,5]})
padj <- p.adjust(pvalue, method = "fdr", n=length(pvalue))
anova_res <- data.frame(pvalue, padj)
anova_res <- anova_res[anova_res$padj<0.05,]
anova_res <- anova_res[order(anova_res$padj),]
dim(anova_res) #  5242    2
head(anova_res)
noncoding_anova_res <- anova_res


### mean
mean <- sapply(unique(testInfo), function(x){rowMeans(df[,testInfo==x])})
median <- sapply(unique(testInfo), function(x){rowMedians(df[, testInfo==x])})
rownames(median) <- rownames(df)
head(mean)
head(median)
region_high_expression_list <- list()
for(region in colnames(mean)){
  fc <- mean[, region]/mean
  trans <- rownames(fc)[apply(fc, 1, function(x){all(x>=1) & any(x>1.5)})]
  trans <- intersect(rownames(noncoding_anova_res), trans)
  if(length(trans)==0){trans <- rownames(fc)[apply(fc, 1, function(x){all(x>=1) & any(x>1.4)})];
  trans <- intersect(rownames(noncoding_anova_res), trans)}
  if(length(trans)==0){trans <- rownames(fc)[apply(fc, 1, function(x){all(x>=1) & any(x>1.3)})];
  trans <- intersect(rownames(noncoding_anova_res), trans)}
  region_high_expression_list[[region]] <- trans
}
lengths(region_high_expression_list)

# region_high_expression_list_median <- list()
# for(region in colnames(mean)){
#   fc <- median[, region]/median
#   trans <- rownames(fc)[apply(fc, 1, function(x){all(x>=1) & any(x>1.5)})]
#   trans <- intersect(rownames(noncoding_anova_res), trans)
#   if(length(trans)==0){trans <- rownames(fc)[apply(fc, 1, function(x){all(x>=1) & any(x>1.4)})];
#   trans <- intersect(rownames(noncoding_anova_res), trans)}
#   if(length(trans)==0){trans <- rownames(fc)[apply(fc, 1, function(x){all(x>=1) & any(x>1.3)})];
#   trans <- intersect(rownames(noncoding_anova_res), trans)}
#   region_high_expression_list_median[[region]] <- trans
# }
# lengths(region_high_expression_list_median)


##########
library(ggplot2)
library(reshape2)
library(RColorBrewer)
trans <- unlist(region_high_expression_list)
# trans <- unlist(region_high_expression_list_median)
sub_df <- df[trans,]
pp <- melt(sub_df)
colnames(pp) <- c("gene", "sample_id", "expression")
head(pp)
pp$SR <- colData[as.character(pp$sample_id), ]$SR_merge


ggplot(pp, aes(x=SR, y=expression, fill=SR, color=SR)) +
  geom_violin() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        text = element_text(face="bold"),
        strip.text = element_text(color = "black"),
        strip.background = element_rect(fill="white", color = "white")) +
  facet_grid(.~gene, scales="free_x") +
  coord_flip() +
  guides(fill=F, color=F) +
  xlab("") +
  ylab("expression(VST)") +
  scale_x_discrete(limits=rev(names(region_high_expression_list))) +
  scale_fill_manual(values = colorRampPalette(rev(brewer.pal(11,"Paired")))(16)) +
  scale_color_manual(values = colorRampPalette(rev(brewer.pal(11,"Paired")))(16)) +
  ggtitle("noncoding")


ggsave("revise0530/transcripts_count_res_with_TACO_gtf/novel_noncoding_trans_region_high_expr_violin.pdf",
       width = 16, height = 4)






#### gene cluster
library(ConsensusClusterPlus)
top_df <- df[rownames(anova_res),]
ttop_df <- t(top_df)
ttop_df = scale(ttop_df)
title="revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_intergenic_noncoding_anova_trans"
intergenic_noncoding_trans = ConsensusClusterPlus(ttop_df,maxK=10,reps=100,pItem=0.8,pFeature=1,
                                                  title=title,clusterAlg="km",
                                                  innerLinkage = "ward.D2",
                                                  finalLinkage = "ward.D2",
                                                  seed = 1,
                                                  distance="euclidean",
                                                  plot="png")
# cal <- calcICL(DIfGenes_Results_SR_merge)
# cal$clusterConsensus
save(intergenic_noncoding_trans,
     file="revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_intergenic_noncoding_anova_trans.RData")

###
load("revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_intergenic_noncoding_anova_trans.RData")
k=6
res <- intergenic_noncoding_trans[[k]]
geneCluster <- data.frame(geneCluster=factor(sort(res$consensusClass)))


colInfo <- Rhesus_colData[colnames(df),c("SR_merge", "ID")]
head(colInfo)
str(colInfo)

ann_colors <- list(SR_merge = c(anno_colors$SR,anno_colors$merge)[unique(colInfo$SR_merge)],
                   ID = anno_colors$ID)
dfheat <- df[rownames(geneCluster), ]
# dfheat <- df[rownames(anova_res), ]
# dfheat <- df


##
library(dendsort)
library(pheatmap)
library(RColorBrewer)
method1 = c("euclidean", "maximum", "manhattan", "canberra", "minkowski")
method2 = c("ward.D", "ward.D2", "single",  "complete", "average")
method3 = c("average","min")
for(i in method1){
  for(j in method2){
    for(z in method3){
      # i = "canberra"
      # j = "ward.D"
      # z = "min"
      filename= paste("revise0530/transcripts_count_res_with_TACO_gtf/Anova_novel_intergenic_noncoding_trans_profile_", i, "_", j, "_",z,".pdf", sep="")
      title = paste("novel&noncoding&intergenic trans expression\ndistance:", i, " hclust:",j, " type:",z, sep = "")
      pheatmap(
        dfheat,
        col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
        breaks = c(-3,seq(-2,2,length=20),3),
        # legend_breaks = c(-3,-1,1,3),
        filename = filename,
        scale="row",
        main = title,
        # width = 10,
        # height = 10,
        fontsize=6,
        # clustering_method = "complete",
        border_color=NA,
        # fontsize_col = 0.5,
        annotation_col = colInfo,
        annotation_colors = ann_colors,
        gaps_row = cumsum(table(geneCluster$geneCluster)),
        annotation_legend = T,
        show_rownames = F,
        show_colnames = F,
        cluster_rows = F,
        cluster_cols = as.hclust(dendsort(hclust(dist(t(dfheat), method=i), method = j), type=z))
      )
    }
  }
}




#######################################################################
#         intergenic novel(coding/noncoding) trans expression         #            
#######################################################################


## anova (ID)
library(DESeq2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_DESeq2_object.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")
object <- Rhesus_TransObject

## intergenic coding
coding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="coding" & TACO_gtf_trans$CPAT_potential=="coding"),]
coding_trans <- coding_trans[which(coding_trans$category_relative_detail=="intergenic"),]
dim(coding_trans) # 1701   32

## intergenic noncoding
noncoding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="noncoding" | TACO_gtf_trans$CPAT_potential=="noncoding"),]
noncoding_trans <- noncoding_trans[which(noncoding_trans$category_relative_detail=="intergenic"),]
dim(noncoding_trans) # 8002   32


####################
#  intergenic coding
####################
vsd <- varianceStabilizingTransformation(object, blind=FALSE)
vsd_data <- assay(vsd)
vsd_data_ex1 <- vsd_data[as.character(coding_trans$transcript_id), 
                         rownames(Rhesus_colData)[Rhesus_colData$SR_cortex=="cortex"]]
df <- vsd_data_ex1
dim(df) # 1701  248
df <- df[apply(df, 1, sd)!=0, ]
dim(df) # 1677  248

#### anova
colData <- Rhesus_colData[colnames(df),]
testInfo <- colData$ID
pvalue <- apply(df, 1, function(x){d <- data.frame(x,testInfo);summary(aov(x~testInfo, data=d))[[1]][1,5]})
padj <- p.adjust(pvalue, method = "fdr", n=length(pvalue))
anova_res <- data.frame(pvalue, padj)
anova_res <- anova_res[anova_res$padj<0.05,]
anova_res <- anova_res[order(anova_res$padj),]
dim(anova_res) # 1250    2
head(anova_res)
coding_anova_res <- anova_res


#### gene cluster
library(ConsensusClusterPlus)
top_df <- df[rownames(anova_res),]
ttop_df <- t(top_df)
ttop_df = scale(ttop_df)
title="revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_intergenic_coding_anova_trans_cortex"
intergenic_coding_trans_cortex = ConsensusClusterPlus(ttop_df,maxK=10,reps=100,pItem=0.8,pFeature=1,
                                                      title=title,clusterAlg="km",
                                                      innerLinkage = "ward.D2",
                                                      finalLinkage = "ward.D2",
                                                      seed = 1,
                                                      distance="euclidean",
                                                      plot="png")
# cal <- calcICL(DIfGenes_Results_SR_merge)
# cal$clusterConsensus
save(intergenic_coding_trans_cortex,
     file="revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_intergenic_coding_anova_trans_cortex.RData")

###
load("revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_intergenic_coding_anova_trans_cortex.RData")
k=2
res <- intergenic_coding_trans_cortex[[k]]
geneCluster <- data.frame(geneCluster=factor(sort(res$consensusClass)))


colInfo <- Rhesus_colData[colnames(df),c("Age_stage", "Gender", "ID")]
head(colInfo)
str(colInfo)


ann_colors <- list(ID = anno_colors$ID[unique(colInfo$ID)],
                   Age_stage = anno_colors$Age,
                   # SR = setNames(brewer.pal(5, "Set3"), unique(colInfo$SR)),
                   Gender=c(female="#D4695E", male="#C5C7A2"))
dfheat <- df[rownames(geneCluster), ]
# dfheat <- df[rownames(anova_res), ]
# dfheat <- df


##
library(dendsort)
library(pheatmap)
library(RColorBrewer)
method1 = c("euclidean", "maximum", "manhattan", "canberra", "minkowski")
method2 = c("ward.D", "ward.D2", "single",  "complete", "average")
method3 = c("average","min")
for(i in method1){
  for(j in method2){
    for(z in method3){
      i = "canberra"
      j = "single"
      z = "min"
      filename= paste("revise0530/transcripts_count_res_with_TACO_gtf/Anova_novel_intergenic_coding_trans_cortex_profile_", i, "_", j, "_",z,".pdf", sep="")
      title = paste("cortex novel&coding&intergenic trans expression\ndistance:", i, " hclust:",j, " type:",z, sep = "")
      pheatmap(
        dfheat,
        col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
        breaks = c(-3,seq(-2,2,length=20),3),
        # legend_breaks = c(-3,-1,1,3),
        filename = filename,
        scale="row",
        main = title,
        # width = 10,
        # height = 10,
        fontsize=6,
        # clustering_method = "complete",
        border_color=NA,
        # fontsize_col = 0.5,
        annotation_col = colInfo,
        annotation_colors = ann_colors,
        gaps_row = cumsum(table(geneCluster$geneCluster)),
        annotation_legend = T,
        show_rownames = F,
        show_colnames = F,
        cluster_rows = F,
        cluster_cols = as.hclust(dendsort(hclust(dist(t(dfheat), method=i), method = j), type=z))
      )
    }
  }
}



###################################
#  intergenic noncoding
###################################
vsd <- varianceStabilizingTransformation(object, blind=FALSE)
vsd_data <- assay(vsd)
vsd_data_ex1 <- vsd_data[as.character(noncoding_trans$transcript_id), 
                         rownames(Rhesus_colData)[Rhesus_colData$SR_cortex=="cortex"]]
df <- vsd_data_ex1
dim(df) #  8002  248
df <- df[apply(df, 1, sd)!=0, ]
dim(df) # 7912  248


#### anova
colData <- Rhesus_colData[colnames(df),]
testInfo <- colData$ID
pvalue <- apply(df, 1, function(x){d <- data.frame(x,testInfo);summary(aov(x~testInfo, data=d))[[1]][1,5]})
padj <- p.adjust(pvalue, method = "fdr", n=length(pvalue))
anova_res <- data.frame(pvalue, padj)
anova_res <- anova_res[anova_res$padj<0.05,]
anova_res <- anova_res[order(anova_res$padj),]
dim(anova_res) #  5737    2
head(anova_res)
noncoding_anova_res <- anova_res



#### gene cluster
library(ConsensusClusterPlus)
top_df <- df[rownames(anova_res),]
ttop_df <- t(top_df)
ttop_df = scale(ttop_df)
title="revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_intergenic_noncoding_anova_trans_cortex"
intergenic_noncoding_trans_cortex = ConsensusClusterPlus(ttop_df,maxK=10,reps=100,pItem=0.8,pFeature=1,
                                                         title=title,clusterAlg="km",
                                                         innerLinkage = "ward.D2",
                                                         finalLinkage = "ward.D2",
                                                         seed = 1,
                                                         distance="euclidean",
                                                         plot="png")
# cal <- calcICL(DIfGenes_Results_SR_merge)
# cal$clusterConsensus
save(intergenic_noncoding_trans_cortex,
     file="revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_intergenic_noncoding_anova_trans_cortex.RData")

###
load("revise0530/transcripts_count_res_with_TACO_gtf/ConsensusCluster_novel_intergenic_noncoding_anova_trans_cortex.RData")
k=7
res <- intergenic_noncoding_trans_cortex[[k]]
geneCluster <- data.frame(geneCluster=factor(sort(res$consensusClass)))


colInfo <- Rhesus_colData[colnames(df),c("Age_stage", "Gender", "ID")]
head(colInfo)
str(colInfo)


ann_colors <- list(ID = anno_colors$ID[unique(colInfo$ID)],
                   Age_stage = anno_colors$Age,
                   # SR = setNames(brewer.pal(5, "Set3"), unique(colInfo$SR)),
                   Gender=c(female="#D4695E", male="#C5C7A2"))
dfheat <- df[rownames(geneCluster), ]
# dfheat <- df[rownames(anova_res), ]
# dfheat <- df


##
library(dendsort)
library(pheatmap)
library(RColorBrewer)
method1 = c("euclidean", "maximum", "manhattan", "canberra", "minkowski")
method2 = c("ward.D", "ward.D2", "single",  "complete", "average")
method3 = c("average","min")
for(i in method1){
  for(j in method2){
    for(z in method3){
      i = "canberra"
      j = "ward.D"
      z = "min"
      filename= paste("revise0530/transcripts_count_res_with_TACO_gtf/Anova_novel_intergenic_noncoding_trans_cortex_profile_", i, "_", j, "_",z,".pdf", sep="")
      title = paste("cortex novel&noncoding&intergenic trans expression\ndistance:", i, " hclust:",j, " type:",z, sep = "")
      pheatmap(
        dfheat,
        col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
        breaks = c(-3,seq(-2,2,length=20),3),
        # legend_breaks = c(-3,-1,1,3),
        filename = filename,
        scale="row",
        main = title,
        # width = 10,
        # height = 10,
        fontsize=6,
        # clustering_method = "complete",
        border_color=NA,
        # fontsize_col = 0.5,
        annotation_col = colInfo,
        annotation_colors = ann_colors,
        # gaps_row = cumsum(table(geneCluster$geneCluster)),
        annotation_legend = T,
        show_rownames = F,
        show_colnames = F,
        # cluster_rows = F,
        cluster_cols = as.hclust(dendsort(hclust(dist(t(dfheat), method=i), method = j), type=z))
      )
    }
  }
}


#######################################################################
#                             randomForest                            #            
#######################################################################
library(randomForest)
library(mclust)
library(DESeq2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_DESeq2_object.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")
object <- Rhesus_TransObject

## intergenic coding
coding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="coding" & TACO_gtf_trans$CPAT_potential=="coding"),]
coding_trans <- coding_trans[which(coding_trans$category_relative_detail=="intergenic"),]
dim(coding_trans) # 1701   32

## intergenic noncoding
noncoding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="noncoding" | TACO_gtf_trans$CPAT_potential=="noncoding"),]
noncoding_trans <- noncoding_trans[which(noncoding_trans$category_relative_detail=="intergenic"),]
dim(noncoding_trans) # 8002   32


####################
#  intergenic coding
####################
vsd <- varianceStabilizingTransformation(object, blind=FALSE)
vsd_data <- assay(vsd)
df <- vsd_data[rownames(vsd_data) %in% rownames(coding_trans),]
dim(df) # 1589  408  
df <- df[apply(df, 1, sd)!=0, ]
dim(df) # 1589  408 



###Random Forest
library(randomForest)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

for(seedno in 1:10){
  for(mtryno in 20:50){
    set.seed(seedno)
    rf = randomForest(t(df), factor(testInfo),
                      importance=TRUE, proximity = TRUE, 
                      sampsize = nrow(t(df)),
                      mtry=mtryno,ntree = 1000)
    rf_importance = as.data.frame(importance(rf,scale=FALSE))
    rf_importance = rf_importance[order(rf_importance[,"MeanDecreaseGini"],decreasing = TRUE),]
    trans_id <- rownames(rf_importance)[1:10]
    pp <- df[trans_id,]
    pp <- melt(pp)
    colnames(pp) <- c("gene", "sample_id", "expression")
    head(pp)
    pp$SR <- colData[as.character(pp$sample_id), ]$SR_merge
    # pp$gene <- as.character(pp$gene)
    # pp$SR <- as.character(pp$SR)
    # sums <- aggregate(expression~gene, data=pp, sum)
    # sums <- sums[order(sums$expression, decreasing = T),]
    # pp$gene <- factor(pp$gene, levels = sums$gene)
    # means <- aggregate(expression~SR,data=pp, mean)
    # means <- means[order(means$expression, decreasing = T),]
    # pp$SR <- factor(pp$SR, levels = means$SR)
    
    
    ggplot(pp, aes(x=SR, y=expression, fill=SR, color=SR)) +
      geom_violin() +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face="bold"),
            axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
            # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
            text = element_text(face="bold"),
            strip.text = element_text(color = "black"),
            strip.background = element_rect(fill="white", color = "white")) +
      facet_grid(.~gene, scales="free_x") +
      coord_flip() +
      guides(fill=F, color=F) +
      xlab("") +
      ylab("expression(VST)") +
      # scale_x_discrete(limits=rev(names())) +
      scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Paired"))(16)) +
      scale_color_manual(values = colorRampPalette(brewer.pal(12,"Paired"))(16)) +
      ggtitle(paste("coding: randomForest importance top10 \n seed:", seedno, " mtry:", mtryno, sep=""))
    
    ggsave(paste("revise0530/transcripts_count_res_with_TACO_gtf/novel_coding_trans_region_rf_important_top10_expr_violin_seed", seedno, "_mtry", mtryno, ".pdf", sep=""),
           width = 10, height = 4)
  }
}



####################
#  intergenic noncoding
####################
vsd <- varianceStabilizingTransformation(object, blind=FALSE)
vsd_data <- assay(vsd)
vsd_data_ex1 <- vsd_data[as.character(noncoding_trans$transcript_id), 
                         -which(colnames(vsd_data)=="Rhesus_6_36")]
df <- vsd_data_ex1
dim(df) # 8002  407
df <- df[apply(df, 1, sd)!=0, ]
dim(df) # 7999  407

#### anova
colData <- Rhesus_colData[colnames(df),]
testInfo <- colData$SR_merge
pvalue <- apply(df, 1, function(x){d <- data.frame(x,testInfo);summary(aov(x~testInfo, data=d))[[1]][1,5]})
padj <- p.adjust(pvalue, method = "fdr", n=length(pvalue))
anova_res <- data.frame(pvalue, padj)
anova_res <- anova_res[anova_res$padj<0.05,]
anova_res <- anova_res[order(anova_res$padj),]
dim(anova_res) # 5242    2
head(anova_res)
noncoding_anova_res <- anova_res

###Random Forest
library(randomForest)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

region_order <- c("SN","HTHA","MO","PA", "AMY", "HIP", 
                  "PON",  "VTA", 
                  "GP","CC","cortex", "STR","THA", "CB","OB","PIT")
for(seedno in 2:10){
  for(mtryno in 20:50){
    seedno = 9
    mtryno = 20
    set.seed(seedno)
    rf = randomForest(t(df), factor(testInfo),
                      importance=TRUE, proximity = TRUE, 
                      sampsize = nrow(t(df)),
                      mtry=mtryno,ntree = 1000)
    rf_importance = as.data.frame(importance(rf,scale=FALSE))
    rf_importance = rf_importance[order(rf_importance[,"MeanDecreaseGini"],decreasing = TRUE),]
    trans_id <- rownames(rf_importance)[1:20]
    pp <- df[trans_id,]
    pp <- melt(pp)
    colnames(pp) <- c("gene", "sample_id", "expression")
    head(pp)
    pp$SR <- colData[as.character(pp$sample_id), ]$SR_merge
    # pp$gene <- as.character(pp$gene)
    # pp$SR <- as.character(pp$SR)
    # sums <- aggregate(expression~gene, data=pp, sum)
    # sums <- sums[order(sums$expression, decreasing = T),]
    # pp$gene <- factor(pp$gene, levels = sums$gene)
    # means <- aggregate(expression~SR,data=pp, mean)
    # means <- means[order(means$expression, decreasing = T),]
    # pp$SR <- factor(pp$SR, levels = means$SR)
    # 
    
    ggplot(pp, aes(x=SR, y=expression, fill=SR, color=SR)) +
      geom_violin() +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face="bold"),
            axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
            # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
            text = element_text(face="bold"),
            strip.text = element_text(color = "black"),
            strip.background = element_rect(fill="white", color = "white")) +
      facet_grid(.~gene, scales="free_x") +
      coord_flip() +
      guides(fill=F, color=F) +
      xlab("") +
      ylab("expression(VST)") +
      scale_x_discrete(limits=region_order) +
      scale_fill_manual(values = setNames(colorRampPalette(brewer.pal(12,"Paired"))(16), region_order)) +
      scale_color_manual(values = setNames(colorRampPalette(brewer.pal(12,"Paired"))(16), region_order)) +
      ggtitle(paste("noncoding: randomForest importance top20 \n seed:", seedno, " mtry:", mtryno, sep=""))
    
    ggsave(paste("revise0530/transcripts_count_res_with_TACO_gtf/novel_noncoding_trans_region_rf_important_top20_expr_violin_seed", seedno, "_mtry", mtryno, ".pdf", sep=""),
           width = 16, height = 4)
  }
}




#######################################################################
#                     randomForest for cortex                         #            
#######################################################################
library(randomForest)
library(mclust)
library(DESeq2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_DESeq2_object.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")
object <- Rhesus_TransObject

## intergenic coding
coding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="coding" & TACO_gtf_trans$CPAT_potential=="coding"),]
coding_trans <- coding_trans[which(coding_trans$category_relative_detail=="intergenic"),]
dim(coding_trans) # 1701   32

## intergenic noncoding
noncoding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="noncoding" | TACO_gtf_trans$CPAT_potential=="noncoding"),]
noncoding_trans <- noncoding_trans[which(noncoding_trans$category_relative_detail=="intergenic"),]
dim(noncoding_trans) # 8002   32


#############################
#  intergenic coding (cortex)
#############################
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_DESeq2_object.RData")
vsd <- varianceStabilizingTransformation(object, blind=FALSE)
vsd_data <- assay(vsd)
vsd_data_ex1 <- vsd_data[as.character(coding_trans$transcript_id), 
                         rownames(Rhesus_colData)[Rhesus_colData$SR_cortex=="cortex"]]
df <- vsd_data_ex1
dim(df) # 1701  248
df <- df[apply(df, 1, sd)!=0, ]
dim(df) # 1677  248


###Random Forest
library(randomForest)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

colData <- Rhesus_colData[colnames(df),]
testInfo <- colData$BR

for(seedno in 1:10){
  for(mtryno in 20:50){
    # seedno <- 1
    # mtryno <- 20
    set.seed(seedno)
    rf = randomForest(t(df), factor(testInfo),
                      importance=TRUE, proximity = TRUE, 
                      sampsize = nrow(t(df)),
                      mtry=mtryno,ntree = 1000)
    rf_importance = as.data.frame(importance(rf,scale=FALSE))
    rf_importance = rf_importance[order(rf_importance[,"MeanDecreaseGini"],decreasing = TRUE),]
    trans_id <- rownames(rf_importance)[1:10]
    pp <- df[trans_id,]
    pp <- melt(pp)
    colnames(pp) <- c("gene", "sample_id", "expression")
    head(pp)
    pp$BR <- colData[as.character(pp$sample_id), ]$BR
    # pp$gene <- as.character(pp$gene)
    # pp$SR <- as.character(pp$SR)
    # sums <- aggregate(expression~gene, data=pp, sum)
    # sums <- sums[order(sums$expression, decreasing = T),]
    # pp$gene <- factor(pp$gene, levels = sums$gene)
    # means <- aggregate(expression~SR,data=pp, mean)
    # means <- means[order(means$expression, decreasing = T),]
    # pp$SR <- factor(pp$SR, levels = means$SR)
    
    
    ggplot(pp, aes(x=BR, y=expression, fill=BR, color=BR)) +
      geom_violin() +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face="bold"),
            axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
            # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
            text = element_text(face="bold"),
            strip.text = element_text(color = "black"),
            strip.background = element_rect(fill="white", color = "white")) +
      facet_grid(.~gene, scales="free_x") +
      coord_flip() +
      guides(fill=F, color=F) +
      xlab("") +
      ylab("expression(VST)") +
      # scale_x_discrete(limits=rev(names())) +
      scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Paired")[(1:5)*2-1])(5)) +
      scale_color_manual(values = colorRampPalette(brewer.pal(12,"Paired")[(1:5)*2-1])(5)) +
      ggtitle(paste("cortex coding: randomForest importance top10 \n seed:", seedno, " mtry:", mtryno, sep=""))
    
    ggsave(paste("revise0530/transcripts_count_res_with_TACO_gtf/cortex_novel_coding_trans_region_rf_important_top10_expr_violin_seed", seedno, "_mtry", mtryno, ".pdf", sep=""),
           width = 10, height = 4)
  }
}


####################
#  intergenic noncoding (cortex)
####################
vsd <- varianceStabilizingTransformation(object, blind=FALSE)
vsd_data <- assay(vsd)
vsd_data_ex1 <- vsd_data[as.character(noncoding_trans$transcript_id), 
                         rownames(Rhesus_colData)[Rhesus_colData$SR_cortex=="cortex"]]
df <- vsd_data_ex1
dim(df) # 8002  248
df <- df[apply(df, 1, sd)!=0, ]
dim(df) # 7912  248


###Random Forest
library(randomForest)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

colData <- Rhesus_colData[colnames(df),]
testInfo <- colData$BR
testInfo <- colData$ID

for(seedno in 1:2){
  for(mtryno in 20:30){
    set.seed(seedno)
    rf = randomForest(t(df), factor(testInfo),
                      importance=TRUE, proximity = TRUE, 
                      sampsize = nrow(t(df)),
                      mtry=mtryno,ntree = 1000)
    rf_importance = as.data.frame(importance(rf,scale=FALSE))
    rf_importance = rf_importance[order(rf_importance[,"MeanDecreaseGini"],decreasing = TRUE),]
    trans_id <- rownames(rf_importance)[1:10]
    pp <- df[trans_id,]
    pp <- melt(pp)
    colnames(pp) <- c("gene", "sample_id", "expression")
    head(pp)
    pp$ID <- colData[as.character(pp$sample_id), ]$ID
    # pp$gene <- as.character(pp$gene)
    # pp$SR <- as.character(pp$SR)
    # sums <- aggregate(expression~gene, data=pp, sum)
    # sums <- sums[order(sums$expression, decreasing = T),]
    # pp$gene <- factor(pp$gene, levels = sums$gene)
    # means <- aggregate(expression~SR,data=pp, mean)
    # means <- means[order(means$expression, decreasing = T),]
    # pp$SR <- factor(pp$SR, levels = means$SR)
    # 
    
    ggplot(pp, aes(x=ID, y=expression, fill=ID, color=ID)) +
      geom_violin() +
      theme_bw() +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, face="bold"),
            axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),
            # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
            text = element_text(face="bold"),
            strip.text = element_text(color = "black"),
            strip.background = element_rect(fill="white", color = "white")) +
      facet_grid(.~gene, scales="free_x") +
      coord_flip() +
      guides(fill=F, color=F) +
      xlab("") +
      ylab("expression(VST)") +
      # scale_x_discrete(limits=rev(names())) +
      scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Paired")[(1:8)])(8)) +
      scale_color_manual(values = colorRampPalette(brewer.pal(12,"Paired")[(1:8)])(8)) +
      ggtitle(paste("cortex noncoding: randomForest importance top10 \n seed:", seedno, " mtry:", mtryno, sep=""))
    
    ggsave(paste("revise0530/transcripts_count_res_with_TACO_gtf/cortex_for_ID_novel_noncoding_trans_region_rf_important_top10_expr_violin_seed", seedno, "_mtry", mtryno, ".pdf", sep=""),
           width = 10, height = 4)
  }
}




#######################################################################
#                  de nove transcripts circles                        #            
#######################################################################
library(circlize)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")


fai_file = "/mnt/data1/Ref/Rhesus_macaque_8.0.1/Macaca_mulatta.Mmul_8.0.1.dna_sm.fa.fai"
fai <- read.table(fai_file, stringsAsFactors = F)
fai <- data.frame(chr=fai$V1, start=0, end=fai$V2, stringsAsFactors = F)
fai <- fai[-which(fai$chr %in% c("MT")),]
fai$chr <- as.factor(fai$chr)
head(fai)


bed <- TACO_gtf_trans[which(TACO_gtf_trans$shared_splicing=="False" | is.na(TACO_gtf_trans$shared_splicing)), c(1:3,14,20)]
sub_bed <- bed[-which(bed$seqnames=="MT"),]
sub_bed$seqnames <- as.character(sub_bed$seqnames)
sub_bed$seqnames <- factor(sub_bed$seqnames)
sub_bed$expr <- as.numeric(sub_bed$expr)
intergenic <- sub_bed[which(sub_bed$category_relative=="intergenic"), 1:3]
exonic_overlap <- sub_bed[which(sub_bed$category_relative=="exonic_overlap"), 1:3]
intragenic <- sub_bed[which(sub_bed$category_relative=="intragenic"), 1:4]
bed_list = list(intergenic, exonic_overlap, intragenic)


ppi <- 720
png("revise0530/de_novo_TACO_minExpr5.5_res/novel_trans_circles.png",
    width = 6*ppi, height = 6*ppi, res =ppi)
circos.initializeWithIdeogram(fai, 
                              plotType = c("axis", "labels"),
                              # labels.cex = 1,
                              track.height = convert_height(3, "mm")
)

# circos.genomicRainfall(bed_list, pch = 16, cex = 0.3, col = c("#F5922F", "#D1063F", "#1D72BB"))
circos.genomicDensity(intergenic, col = c("#F5922F"), track.height = 0.1)
circos.genomicDensity(exonic_overlap, col = c("#D1063F"), track.height = 0.1)
circos.genomicDensity(intragenic, col = c("#1D72BB"), track.height = 0.1)
# circos.clear()
dev.off()



##################################################
load("revise0530/de_novo_TACO_minExpr5.5_res/denovo_orth.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")

coding_trans <- TACO_gtf_trans[which(TACO_gtf_trans$CPC2_potential=="coding" & TACO_gtf_trans$CPAT_potential=="coding"),]
coding_trans <- coding_trans[which(coding_trans$shared_splicing!="True" | is.na(coding_trans$shared_splicing)),]
dim(coding_trans) # 27792    32

rownames(denovo_orth) <- denovo_orth$Rhesus_denovo_trans_id
coding_orth_trans <- denovo_orth[rownames(denovo_orth) %in% coding_trans$transcript_id,]
coding_orth_trans <- coding_orth_trans[complete.cases(coding_orth_trans$human_gene_name),]
intergenic_coding_orth_trans <- coding_orth_trans[coding_orth_trans$category_relative_detail=="intergenic",]
dim(intergenic_coding_orth_trans) # 223  38



pp_bed <- intergenic_coding_orth_trans[, c(7:9)]



circos.initializeWithIdeogram(fai, 
                              plotType = c("axis", "labels"),
                              track.height = convert_height(4, "mm"))

circos.genomicRainfall(pp_bed, pch = 16, cex = 0.3, col = c("#F5922F", "#D1063F", "#1D72BB"))





#######################################################################
#                          predicted lncRNA                           #            
#######################################################################
library(devtools)
library(rtracklayer)
library(bedr)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")

## ref
ref_gtf <- rtracklayer::import("/mnt/xdlab/ref/M.fascicularis/Mmul_8.0.1.91_ensembl/Macaca_mulatta.Mmul_8.0.1.91.chr.gtf")
ref_gtf <- as.data.frame(ref_gtf)
ref_gtf_gene <- ref_gtf[ref_gtf$type=="gene",]
ref_gtf_gene_coding <- ref_gtf_gene[ref_gtf_gene$gene_biotype == "protein_coding",]
write.table(ref_gtf_gene_coding[, 1:3], "revise0530/de_novo_TACO_minExpr5.5_res/ref_gtf_gene_coding.bed",
            row.names=F, col.name=F, quote=F, sep="\t")


## TACO
TACO_gtf <- rtracklayer::import("gtf_for_assembled/TACO_minExpr_5.5/assembly.refcomp.gtf")
TACO_gtf <- as.data.frame(TACO_gtf)
dim(TACO_gtf)


####  novel transcripts including intergenic and antisense region were reserved as the candidate lncRNAs. 
novel_trans <- TACO_gtf_trans[TACO_gtf_trans$category_relative_detail  %in% c("intergenic", "intronic_opp_strand", "opp_strand", "interleaving_opp_strand"),]
dim(novel_trans) # 14930    32
write.table(novel_trans[,c(1:3,18)], "revise0530/de_novo_TACO_minExpr5.5_res/novel_trans.bed",
            row.names=F, col.name=F, quote=F, sep="\t")

###### bedtools 
# /mnt/xdlab1/tools/bedtools2/bedtools2/bin/bedtools slop 
# -i ref_gtf_gene_coding.bed 
# -g /mnt/xdlab/ref/M.fascicularis/Mmul_8.0.1.91_ensembl/Macaca_mulatta.Mmul_8.0.1.dna_sm.fa.fai 
# -b 1000 > ref_gtf_gene_coding_slop1000.bed

# /mnt/xdlab1/tools/bedtools2/bedtools2/bin/bedtools intersect 
# -a novel_trans.bed 
# -b ref_gtf_gene_coding_slop1000.bed 
# -wao > novel_trans_intersect_ref_gtf_gene_coding_slop1000.bed

#####
intersect_res <- read.table("revise0530/de_novo_TACO_minExpr5.5_res/novel_trans_intersect_ref_gtf_gene_coding_slop1000.bed",
                            header=F, sep="\t", stringsAsFactors = F)

intersect_trans <- unique(intersect_res$V4[which(intersect_res$V8>0)])
length(intersect_trans) # 7829


###### Transcripts adjacent to known coding genes within 1000 bp were regarded as UTRs and also discarded.
novel_trans_ex_intersect <- novel_trans[-which(rownames(novel_trans) %in% intersect_trans),]
dim(novel_trans_ex_intersect) # 7101   32

###### noncoding
novel_trans_noncoding <- novel_trans_ex_intersect[-which(novel_trans_ex_intersect$CPC2_potential=="coding"&novel_trans_ex_intersect$CPAT_potential=="coding"),]
dim(novel_trans_noncoding) # 5746   32


#### with multiple exons 
TACO_gtf_multiExon <- TACO_gtf[which(as.numeric(TACO_gtf$exon_number)>0),]
multiExon_trans <- unique(TACO_gtf_multiExon$transcript_id)
novel_trans_noncoding_multiExon <- novel_trans_noncoding[which(novel_trans_noncoding$transcript_id %in% multiExon_trans), ]
dim(novel_trans_noncoding_multiExon) # 2845   32


#### and no smaller than 200 bases were reserved as lncRNAs.
novel_trans_noncoding_multiExon_detail <- TACO_gtf[TACO_gtf$transcript_id %in% rownames(novel_trans_noncoding_multiExon),]
novel_trans_noncoding_multiExon_detail <- novel_trans_noncoding_multiExon_detail[novel_trans_noncoding_multiExon_detail$type=="exon", ]
exon_width <- split(novel_trans_noncoding_multiExon_detail$width, novel_trans_noncoding_multiExon_detail$transcript_id)
trans_length <- melt(lapply(exon_width, sum))
lncRNA_trans <- trans_length$L1[trans_length$value>=200]
novel_lncRNA_trans <- novel_trans_noncoding_multiExon[lncRNA_trans,]
dim(novel_lncRNA_trans) # 2845   32
save(novel_lncRNA_trans, file="revise0530/de_novo_TACO_minExpr5.5_res/novel_lncRNA_trans.RData")

lncRNA_name <- rownames(novel_lncRNA_trans)
write.table(lncRNA_name, "revise0530/de_novo_TACO_minExpr5.5_res/lncRNA_transcripts_id_2845.txt",
            quote=F, sep="\t", row.names=F, col.names=F)


lncRNA_info <- TACO_gtf_trans[lncRNA_name,]
lncRNA_info_out <- paste("macaque_Mmul8", paste("chr", lncRNA_info$seqnames, sep=""), lncRNA_info$start, sep = "|")
lncRNA_info_out <- paste(lncRNA_info_out, lncRNA_info$end, sep="-")
lncRNA_info_out <- data.frame(transcript_id=lncRNA_name, site_info=lncRNA_info_out)
write.table(lncRNA_info_out, 
            "revise0530/de_novo_TACO_minExpr5.5_res/lncRNA_transcripts_id_info_2845.txt",
            quote=F, sep="\t", row.names=F, col.names=F)



#######################################################################
#                     predicted lncRNA expression                     #            
#######################################################################
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/novel_lncRNA_trans.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_DESeq2_object.RData")
vsd <- varianceStabilizingTransformation(Rhesus_TransObject, blind=FALSE)
vsd_data <- assay(vsd)
vsd_data_ex1 <- vsd_data[as.character(novel_lncRNA_trans$transcript_id), -which(colnames(vsd_data)=="Rhesus_6_36")]
dim(vsd_data_ex1) # 2845  407


## for region
df <- vsd_data_ex1
colData <- Rhesus_colData[colnames(vsd_data_ex1),]
testInfo <- colData$SR_merge
pvalue <- apply(df, 1, function(x){d <- data.frame(x,testInfo);summary(aov(x~testInfo, data=d))[[1]][1,5]})
padj <- p.adjust(pvalue, method = "fdr", n=length(pvalue))
anova_res <- data.frame(pvalue, padj)
anova_res <- anova_res[anova_res$padj<0.05,]
anova_res <- anova_res[order(anova_res$padj),]
dim(anova_res) # 1987    2
head(anova_res)


## for individual
df <- vsd_data_ex1[,rownames(Rhesus_colData)[Rhesus_colData$SR_cortex=="cortex"]]
dim(df) #2845  248
colData <- colData[colnames(df),]
testInfo <- colData$ID
pvalue <- apply(df, 1, function(x){d <- data.frame(x,testInfo);summary(aov(x~testInfo, data=d))[[1]][1,5]})
padj <- p.adjust(pvalue, method = "fdr", n=length(pvalue))
anova_res <- data.frame(pvalue, padj)
anova_res <- anova_res[anova_res$padj<0.05,]
anova_res <- anova_res[order(anova_res$padj),]
dim(anova_res) # 2010    2
head(anova_res)


#######################################################################
#         blastn res between lncRNA and human lncRNA database         #            
#######################################################################
library(devtools)
library(rtracklayer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/de_novo_TACO_minExpr5.5_res/novel_lncRNA_trans.RData")


#####  GENCODE
GENCODE_gtf <- rtracklayer::import("ANNOTATION/gencode.v28.long_noncoding_RNAs_ex_chr.gtf")
GENCODE_gtf <- as.data.frame(GENCODE_gtf)
GENCODE_gtf_trans <- GENCODE_gtf[GENCODE_gtf$type=="transcript", ]
rownames(GENCODE_gtf_trans) <- GENCODE_gtf_trans$transcript_id

GENCODE_res <- read.table("revise0530/de_novo_TACO_minExpr5.5_res/lncRNA_GENCODE_blastn_res.txt",
                          stringsAsFactors = F)
colnames(GENCODE_res) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                           "qstart", "qend", "sstart", "send", "evalue", "bitscore")
GENCODE_res <- GENCODE_res[order(GENCODE_res$qseqid, GENCODE_res$evalue), ]
GENCODE_res <- GENCODE_res[-which(duplicated(GENCODE_res$qseqid)),]
GENCODE_res$GENCODE_name <- GENCODE_gtf_trans[GENCODE_res$sseqid, "gene_name"]
rownames(GENCODE_res) <- GENCODE_res$qseqid
max(GENCODE_res$evalue)


#### LNCipedia
LNCipedia_res <- read.table("revise0530/de_novo_TACO_minExpr5.5_res/lncRNA_LNCipedia_blastn_res.txt",
                            stringsAsFactors = F)
colnames(LNCipedia_res) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                             "qstart", "qend", "sstart", "send", "evalue", "bitscore")
LNCipedia_res <- LNCipedia_res[order(LNCipedia_res$qseqid, LNCipedia_res$evalue), ]
LNCipedia_res <- LNCipedia_res[-which(duplicated(LNCipedia_res$qseqid)),]
LNCipedia_res$LNCipedia_name <- LNCipedia_res$sseqid
rownames(LNCipedia_res) <- LNCipedia_res$qseqid
max(LNCipedia_res$evalue)



####  NONCODE
NONCODE_res <- read.table("revise0530/de_novo_TACO_minExpr5.5_res/lncRNA_NONCODE_blastn_res.txt",
                          stringsAsFactors = F)
colnames(NONCODE_res) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                           "qstart", "qend", "sstart", "send", "evalue", "bitscore")
NONCODE_res <- NONCODE_res[order(NONCODE_res$qseqid, NONCODE_res$evalue), ]
NONCODE_res <- NONCODE_res[-which(duplicated(NONCODE_res$qseqid)),]
NONCODE_res$NONCODE_name <- NONCODE_res$sseqid
rownames(NONCODE_res) <- NONCODE_res$qseqid
max(NONCODE_res$evalue)



##############################################
dim(novel_lncRNA_trans) # 2845   32
dim(GENCODE_res) # 2146   13
dim(LNCipedia_res) # 2371   13
dim(NONCODE_res) # 2496   13


novel_lncRNA_trans$qseqid <- rownames(novel_lncRNA_trans)
merge_res <- merge(novel_lncRNA_trans, GENCODE_res[,c(1,13)], by = "qseqid", all = TRUE)
merge_res <- merge(merge_res, LNCipedia_res[,c(1,13)], by="qseqid", all = TRUE)
merge_res <- merge(merge_res, NONCODE_res[,c(1,13)], by="qseqid", all = TRUE)
dim(merge_res)


merge_res$final_name <- merge_res$NONCODE_name
merge_res$final_name[which(complete.cases(merge_res$LNCipedia_name))] <- merge_res$LNCipedia_name[which(complete.cases(merge_res$LNCipedia_name))]
merge_res$final_name[which(complete.cases(merge_res$GENCODE_name))] <- merge_res$GENCODE_name[which(complete.cases(merge_res$GENCODE_name))]
rownames(merge_res) <- merge_res$qseqid
novel_lncRNA_trans$final_name <- merge_res[rownames(novel_lncRNA_trans), "final_name"]


save(novel_lncRNA_trans, file="revise0530/de_novo_TACO_minExpr5.5_res/novel_lncRNA_trans.RData")

write.table(novel_lncRNA_trans[, c("transcript_id", "final_name", "CPC2_potential", "CPAT_potential")], "revise0530/de_novo_TACO_minExpr5.5_res/lncRNA_blast_human_lncRNA.txt", row.name=F, col.names=T, quote=F, sep="\t")


#######################################################################
#                              lncRNA-module                          #            
#######################################################################
library(WGCNA)
library(DESeq2)
library(reshape2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/novel_lncRNA_trans.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_DESeq2_object.RData")
vsd <- varianceStabilizingTransformation(Rhesus_TransObject, blind=FALSE)
vsd_data <- assay(vsd)
df <- vsd_data[as.character(novel_lncRNA_trans$transcript_id), -which(colnames(vsd_data)=="Rhesus_6_36")]
dim(df) # 2845  407

## for all region
load("revise0530/gene_count_res/WGCNA_Rhesus_net_modules_exclude_MB_ex1.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_exclude_MB_ex1.RData")

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs = net$MEs; 
dim(MEs) # 407 56
# moduleColorsWW = moduleColors
# MEs0 = moduleEigengenes(datExpr, moduleLabels)$eigengenes # MEs = net$MEs
# dim(MEs0) # [1] 407 56
MEsWW = orderMEs(MEs)
rownames(MEsWW) <- rownames(datExpr)

# cor
lncRNA_module <- cor(t(df)[rownames(MEsWW),], MEsWW)
dim(lncRNA_module) # 2845   56
lncRNA_module_pvalue = corPvalueStudent(lncRNA_module, nrow(datExpr))

pheatmap(lncRNA_module)
sub_lncRNA_module <- lncRNA_module[apply(lncRNA_module,1,function(x){any(abs(x)>0.7)}),]
dim(sub_lncRNA_module)
pheatmap(sub_lncRNA_module)


lncRNA_module_pairs <- melt(lncRNA_module)
colnames(lncRNA_module_pairs) <- c("transcript_id", "module", "cor")
rownames(lncRNA_module_pairs) <- paste(lncRNA_module_pairs$transcript_id, lncRNA_module_pairs$module, sep="_")
head(lncRNA_module_pairs)
lncRNA_module_pairs_0.7 <- lncRNA_module_pairs[abs(lncRNA_module_pairs$cor)>0.7,]
dim(lncRNA_module_pairs_0.7) # 649   3

lncRNA_module_pvalue <- melt(lncRNA_module_pvalue)
colnames(lncRNA_module_pvalue) <- c("transcript_id", "module", "pvalue")
head(lncRNA_module_pvalue)
rownames(lncRNA_module_pvalue) <- paste(lncRNA_module_pvalue$transcript_id,lncRNA_module_pvalue$module,sep="_")
lncRNA_module_pvalue_0.01 <- lncRNA_module_pvalue[lncRNA_module_pvalue$pvalue<0.01,]
dim(lncRNA_module_pvalue_0.01) # 59155     3

##
lncRNA_module_pairs_0.7 <- lncRNA_module_pairs_0.7[rownames(lncRNA_module_pairs_0.7) %in% rownames(lncRNA_module_pvalue_0.01),]
dim(lncRNA_module_pairs_0.7) # 649   3
save(lncRNA_module_pairs_0.7, file="revise0530/de_novo_TACO_minExpr5.5_res/lncRNA_module_pairs_0.7.RData")

## positive
lncRNA_module_pairs_0.7_p <- lncRNA_module_pairs_0.7[lncRNA_module_pairs_0.7$cor>0.7,]
pp <- table(lncRNA_module_pairs_0.7_p$module)[table(lncRNA_module_pairs_0.7_p$module)>0]
pp <- data.frame(pp, stringsAsFactors = F)
colnames(pp) <- c("module", "count")
pp <- pp[order(pp$count, decreasing = T),]
pp$module <- as.character(pp$module)
pp$module <- factor(pp$module, levels = pp$module)
head(pp)
str(pp)



ggplot(pp, aes(x=module, y=count)) +
  geom_bar(stat="identity",fill=brewer.pal(9,"Reds")[8]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = "black"),
        text = element_text(color="black")) +
  # scale_fill_manual(values = c(hypo=brewer.pal(9,"Blues")[8], hyper=brewer.pal(9,"Reds")[8])) +
  ggtitle("LncRNA-module positive")

ggsave("revise0530/de_novo_TACO_minExpr5.5_res/LncRNA-module_positive_cor_count_bar.pdf",
       width = 6, height = 3.5)

## negative
lncRNA_module_pairs_0.7_n <- lncRNA_module_pairs_0.7[lncRNA_module_pairs_0.7$cor< -0.7,]
pp <- table(lncRNA_module_pairs_0.7_n$module)[table(lncRNA_module_pairs_0.7_n$module)>0]
pp <- data.frame(pp, stringsAsFactors = F)
colnames(pp) <- c("module", "count")
pp <- pp[order(pp$count, decreasing = T),]
pp$module <- as.character(pp$module)
pp$module <- factor(pp$module, levels = pp$module)
head(pp)
str(pp)


ggplot(pp, aes(x=module, y=count)) +
  geom_bar(stat="identity",fill=brewer.pal(9,"Blues")[8]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, color = "black"),
        text = element_text(color="black")) +
  # scale_fill_manual(values = c(hypo=brewer.pal(9,"Blues")[8], hyper=brewer.pal(9,"Reds")[8])) +
  ggtitle("LncRNA-module negative")

ggsave("revise0530/de_novo_TACO_minExpr5.5_res/LncRNA-module_negative_cor_count_bar.pdf",
       width = 5, height = 3.5)




#######################################################################
#                            lncRNA-refRNA(coding)                    #            
#######################################################################
library(DESeq2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/novel_lncRNA_trans.RData")


### lncRNA
vsd <- varianceStabilizingTransformation(Rhesus_TransObject, blind=FALSE)
vsd_data <- assay(vsd)
lncRNA_data <- vsd_data[rownames(novel_lncRNA_trans), -which(colnames(vsd_data)=="Rhesus_6_36")]
dim(lncRNA_data) # 2845  407
lncRNA_data <- lncRNA_data[which(apply(lncRNA_data, 1, sd)>0),]
dim(lncRNA_data) # 2845  407


### ref coding RNA 
ref_coding_RNA <- rownames(gtf_ensembl_gene)[gtf_ensembl_gene$gene_biotype == "protein_coding"]
vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)
refRNA_data <- vsd_data[ref_coding_RNA, -which(colnames(vsd_data)=="Rhesus_6_36")]
dim(refRNA_data) # 20618   407
refRNA_data <- refRNA_data[which(apply(refRNA_data, 1, sd)>0),]
dim(refRNA_data) # 20257   407



### lncRNA-refcodingRNA correlation
library(reshape2)
library(WGCNA)
lncRNA_refRNA <- cor(t(lncRNA_data), t(refRNA_data))
lncRNA_refRNA_P = corPvalueStudent(lncRNA_refRNA, ncol(lncRNA_data))

lncRNA_refRNA_pairs <- melt(lncRNA_refRNA)
lncRNA_refRNA_P_pairs <- melt(lncRNA_refRNA_P)

colnames(lncRNA_refRNA_pairs) <- c("transcript_id", "gene", "correlation")
head(lncRNA_refRNA_pairs)
colnames(lncRNA_refRNA_P_pairs) <- c("transcript_id", "gene", "pvalue")
head(lncRNA_refRNA_P_pairs)


lncRNA_refRNA_pairs_0.7 <- lncRNA_refRNA_pairs[which(abs(lncRNA_refRNA_pairs$correlation)>0.7), ]
dim(lncRNA_refRNA_pairs_0.7) # 60218     3
head(lncRNA_refRNA_pairs_0.7)
lncRNA_refRNA_P_pairs_0.01 <- lncRNA_refRNA_P_pairs[which(abs(lncRNA_refRNA_P_pairs$pvalue)< 0.01), ]
dim(lncRNA_refRNA_P_pairs_0.01) # 18707942        3

a <- paste(lncRNA_refRNA_pairs_0.7$transcript_id, lncRNA_refRNA_pairs_0.7$gene, sep="_")
b <- paste(lncRNA_refRNA_P_pairs_0.01$transcript_id, lncRNA_refRNA_P_pairs_0.01$gene, sep="_")
all(a %in% b) # TRUE



lncRNA_refRNA_pairs_0.7_list <-  split(lncRNA_refRNA_pairs_0.7, factor(lncRNA_refRNA_pairs_0.7$transcript_id))
names(lncRNA_refRNA_pairs_0.7_list)

lncRNA_refRNA_pairs_0.7_list <- lapply(lncRNA_refRNA_pairs_0.7_list, 
                                       function(x){x$final_name <- novel_lncRNA_trans[as.character(x[,"transcript_id"]), "final_name"];
                                       x$gene_name <- gtf_ensembl_gene[as.character(x[,"gene"]), "gene_name"];
                                       x})

####### function
library(gProfileR)

lncRNA_coexpression_function_res <- list()
for(lncRNA in names(lncRNA_refRNA_pairs_0.7_list)[577:747]){
  genes <- as.character(lncRNA_refRNA_pairs_0.7_list[[lncRNA]][, "gene"])
  gprofiler_res <- gprofiler(
    query=genes, 
    organism = "mmulatta")
  gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),]
  gprofiler_res$gene_name <- gprofiler_res$intersection
  if(dim(gprofiler_res)[1]==0) next
  else{for(i in 1:dim(gprofiler_res)[1]){
    x <- gprofiler_res$gene_name[i]
    gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],"gene_name"], collapse=",")
  }
    lncRNA_coexpression_function_res[[lncRNA]] <- gprofiler_res
  }
}

save(lncRNA_coexpression_function_res, 
     file= "revise0530/de_novo_TACO_minExpr5.5_res/lncRNA_coexpression_function_res.RData")


#################
load("revise0530/de_novo_TACO_minExpr5.5_res/lncRNA_coexpression_function_res.RData")
unlist(lapply(lncRNA_coexpression_function_res, function(x){dim(x)[1]}))

lncRNA_coexpression_function_res2 <- list()
for(lncRNA in names(lncRNA_coexpression_function_res)){
  y <- lncRNA_coexpression_function_res[[lncRNA]]
  y <- data.frame(y, lncRNA_id=lncRNA)
  lncRNA_coexpression_function_res2[[lncRNA]] <- y
}

test <- Reduce(rbind, lncRNA_coexpression_function_res2)
test <- test[test$domain %in% c("CC", "MF", "BP", "keg"), ]
test$lncRNA_id <- as.character(test$lncRNA_id)
str(test)


test$final_name <- novel_lncRNA_trans[test$lncRNA_id, "final_name"]

lncRNA_refRNA_pairs_0.7_list[["TU11253"]]
lncRNA_refRNA_pairs_0.7_list[["TU2423"]]
lncRNA_refRNA_pairs_0.7_list[["TU9944"]]

#### igraph
library(igraph)

TU2423 <- test[test$lncRNA_id=="TU2423", ]

term_name_list <- c("synaptic signaling", "neuron differentiation", "axon development")
final_edges_list <- list()
for(term_name in term_name_list){
  genes <- strsplit(TU2423[TU2423$term.name == term_name, "gene_name"], ",")[[1]]
  if(sum(genes %in% "NA")>0){genes <- genes[-which(genes %in% "NA")]}
  a <- data.frame(term_name, genes)
  final_edges_list[[term_name]] <- a
}

r <- Reduce(rbind, final_edges_list)
genes <- unique(as.character(r$genes))
links <- matrix(data=0, nrow = length(term_name_list), ncol=length(genes))
colnames(links) <- genes
rownames(links) <- term_name_list
for(i in 1:dim(r)[1]){
  link <- r[i, ]
  links[as.character(link[, "term_name"]), as.character(link[, "genes"])] <- 1
}

net2 <- graph_from_incidence_matrix(links)
table(V(net2)$type)

V(net2)$color <- c("steel blue", "orange")[V(net2)$type+1]
V(net2)$shape <- c("square", "circle")[V(net2)$type+1]
V(net2)$label[V(net2)$type==T] <- colnames(links)
V(net2)$label[V(net2)$type==F] <- rownames(links) 

V(net2)$label.cex = .6
V(net2)$label.font = 5



plot(net2, 
     vertex.label.color="black",
     # vertex.label=t,
     vertex.label.family='Arial',
     # layout=layout_as_bipartite,
     vertex.size=(2-V(net2)$type)*10) 

t <- V(net2)$label
t <- gsub(" ", "_",t)




library(network)

routes_network <- network(r, matrix.type = "edgelist", ignore.eval = FALSE)
plot(routes_network, 
     vertex.cex = 3,
     vertex.col="blue",
     displaylabels = TRUE)



library(visNetwork)
library(networkD3)
nodes <- data.frame(id=1:length(V(net2)$label),
                    label=V(net2)$label,
                    color=rep(c("darkred", "orange"), c(3, 60)),
                    shape=V(net2)$shape)
edges <- data.frame(from=sapply(r$term_name, function(x){which(nodes$label %in% x)}),
                    to=sapply(r$genes, function(x){which(nodes$label %in% x)}))

network <- visNetwork(nodes, edges)
visSave(network, file = "/mnt/data2/Rhesus_brain/revise0530/de_novo_TACO_minExpr5.5_res/network.html", 
        background = "white")
visExport(network, type = "png",
          name="network")

nodes <- data.frame(id = 1:10, 
                    label = paste("Node", 1:10),                                 # labels
                    group = c("GrA", "GrB"),                                     # groups 
                    value = 1:10,                                                # size 
                    shape = c("square", "triangle", "box", "circle", "dot", "star",
                              "ellipse", "database", "text", "diamond"),         # shape
                    title = paste0("<p><b>", 1:10,"</b><br>Node !</p>"),         # tooltip
                    color = c("darkred", "grey", "orange", "darkblue", "purple"),# color
                    shadow = c(FALSE, TRUE, FALSE, TRUE, TRUE))                  # shadow

edges <- data.frame(from = sample(1:10,8), to = sample(1:10, 8),
                    label = paste("Edge", 1:8),                                 # labels
                    length = c(100,500),                                        # length
                    arrows = c("to", "from", "middle", "middle;to"),            # arrows
                    dashes = c(TRUE, FALSE),                                    # dashes
                    title = paste("Edge", 1:8),                                 # tooltip
                    smooth = c(FALSE, TRUE),                                    # smooth
                    shadow = c(FALSE, TRUE, FALSE, TRUE))                       # shadow

visNetwork(nodes, edges) 



#######################################################################
#                                lncRNA-ATAC                          #            
#######################################################################
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("ATAC/matrix/peak_top10000_res/Rhesus_DESeq2_object.RData") # Rhesus_PeakObject 
load("revise0530/transcripts_count_res_with_TACO_gtf/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/lncRNA_module_pairs_0.7.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_datExpr_exclude_MB_ex1.RData")
load("revise0530/gene_count_res/WGCNA_Rhesus_net_modules_exclude_MB_ex1.RData")
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")


### ATAC 
vsd <- varianceStabilizingTransformation(Rhesus_PeakObject, blind=FALSE)
vsd_data <- assay(vsd)
dim(vsd_data) # 25842    26
ATAC_data <- vsd_data


### lncRNA
lncRNA_module_pairs_0.7_p <- lncRNA_module_pairs_0.7[lncRNA_module_pairs_0.7$cor>0.7,]
dim(lncRNA_module_pairs_0.7_p) # 585   3
lncRNA_module_pairs_0.7_p_cortex <- lncRNA_module_pairs_0.7_p[lncRNA_module_pairs_0.7_p$module %in% c("ME10", "ME6"),]
dim(lncRNA_module_pairs_0.7_p_cortex) # 183   3

vsd <- varianceStabilizingTransformation(Rhesus_TransObject, blind=FALSE)
vsd_data <- assay(vsd)
lncRNA_data <- vsd_data[unique(as.character(lncRNA_module_pairs_0.7_p_cortex$transcript_id)), colnames(ATAC_data)]
dim(lncRNA_data) # 119  26

chrom_info <- TACO_gtf_trans[rownames(lncRNA_data),c("seqnames", "transcript_id"), drop=F]
chrom_info_list <- split(chrom_info, factor(as.character(chrom_info$seqnames)))
lncRNA_data_list <- lapply(chrom_info_list, function(x){lncRNA_data[x[,2],,drop=F]})


### mod RNA 
probes = names(datExpr)
moduleLabels2 <- paste("ME", net$colors, sep="")

modules <- c("ME10", "ME6")
inModule = is.finite(match(moduleLabels2, modules));
modProbes = probes[inModule];
length(modProbes) #1283

vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)
modRNA_data <- vsd_data[modProbes, colnames(ATAC_data)]
dim(modRNA_data) # 1283   26



### ref RNA 
vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)
refRNA_data <- vsd_data[, colnames(ATAC_data)]
dim(refRNA_data) # 30807    26


### lncRNA-ATAC correlation (同一染色体)
all_peaks <- rownames(ATAC_data)
all_peaks_info <- data.frame(seqnames=sapply(all_peaks, function(x){strsplit(x, ":")[[1]][1]}),
                             peak=all_peaks)
all_peaks_info_list <- split(all_peaks_info, factor(as.character(all_peaks_info$seqnames)))
ATAC_data_list <- lapply(all_peaks_info_list, function(x){ATAC_data[as.character(x[,2]),]})

lncRNA_ATAC_0.7_list <- list()
for(chrom in names(lncRNA_data_list)[1:18]){
  lncRNA_ATAC <- cor(t(lncRNA_data_list[[chrom]]), t(ATAC_data_list[[chrom]]))
  lncRNA_ATAC_pairs <- melt(lncRNA_ATAC)
  colnames(lncRNA_ATAC_pairs) <- c("transcript_id", "peak", "correlation")
  head(lncRNA_ATAC_pairs)
  lncRNA_ATAC_pairs_0.7 <- lncRNA_ATAC_pairs[which(abs(lncRNA_ATAC_pairs$correlation)>0.7), ]
  lncRNA_ATAC_0.7_list[[chrom]] <- lncRNA_ATAC_pairs_0.7
}


lncRNA_ATAC_0.7 <- Reduce(rbind, lncRNA_ATAC_0.7_list)
lncRNA_ATAC_0.7_list <- split(lncRNA_ATAC_0.7, factor(as.character(lncRNA_ATAC_0.7$transcript_id)))




### lncRNA-modRNA correlation
lncRNA_modRNA <- cor(t(lncRNA_data), t(modRNA_data))
lncRNA_modRNA_pairs <- melt(lncRNA_modRNA)
colnames(lncRNA_modRNA_pairs) <- c("transcript_id", "gene", "correlation")
head(lncRNA_modRNA_pairs)

lncRNA_modRNA_pairs_0.7 <- lncRNA_modRNA_pairs[which(abs(lncRNA_modRNA_pairs$correlation)>0.7), ]
dim(lncRNA_modRNA_pairs_0.7) # 967   3
head(lncRNA_modRNA_pairs_0.7)


### lncRNA-refRNA correlation
lncRNA_refRNA <- cor(t(lncRNA_data), t(refRNA_data))
lncRNA_refRNA_pairs <- melt(lncRNA_refRNA)
colnames(lncRNA_refRNA_pairs) <- c("transcript_id", "gene", "correlation")
head(lncRNA_refRNA_pairs)

lncRNA_refRNA_pairs_0.7 <- lncRNA_refRNA_pairs[which(abs(lncRNA_refRNA_pairs$correlation)>0.7), ]
dim(lncRNA_refRNA_pairs_0.7) # 14140     3
head(lncRNA_refRNA_pairs_0.7)



###  lncRNA-ATAC-modRNA
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
head(lncRNA_ATAC_0.7)
head(lncRNA_modRNA_pairs_0.7)

lncRNA_ATAC_modRNA <- merge(lncRNA_ATAC_0.7, lncRNA_modRNA_pairs_0.7, by="transcript_id")
dim(lncRNA_ATAC_modRNA) # 986   5

lncRNA_ATAC_modRNA$transcript_id <- as.character(lncRNA_ATAC_modRNA$transcript_id)
lncRNA_ATAC_modRNA$transcript_id <- factor(lncRNA_ATAC_modRNA$transcript_id)
lncRNA_ATAC_modRNA_list <- split(lncRNA_ATAC_modRNA, lncRNA_ATAC_modRNA$transcript_id)
unlist(lapply(lncRNA_ATAC_modRNA_list, nrow))
unlist(lapply(lncRNA_ATAC_modRNA_list, function(x){length(unique(x[,"gene"]))}))
unlist(lapply(lncRNA_ATAC_modRNA_list, function(x){length(unique(x[,"peak"]))}))
unique(lncRNA_ATAC_modRNA$gene)

lncRNA_ATAC_modRNA_list <- lapply(lncRNA_ATAC_modRNA_list, function(x){x$gene_name <- gtf_ensembl_gene[as.character(x[,"gene"]), "gene_name"];x})
gtf_ensembl_gene[as.character(unique(lncRNA_ATAC_modRNA$gene)), "gene_name"]
factor(gtf_ensembl_gene[as.character(unique(lncRNA_ATAC_modRNA$gene)), "gene_name"])


unlist(lapply(lncRNA_ATAC_modRNA_list, function(x){unique(x$gene_name)}))



###  lncRNA-ATAC-refRNA
head(lncRNA_ATAC_pairs_0.7)
head(lncRNA_refRNA_pairs_0.7)

lncRNA_ATAC_refRNA <- merge(lncRNA_ATAC_pairs_0.7, lncRNA_refRNA_pairs_0.7, by="transcript_id")
dim(lncRNA_ATAC_refRNA) # 255225      5
length(unique(lncRNA_ATAC_refRNA$transcript_id))

lncRNA_ATAC_refRNA$transcript_id <- as.character(lncRNA_ATAC_refRNA$transcript_id)
lncRNA_ATAC_refRNA$transcript_id <- factor(lncRNA_ATAC_refRNA$transcript_id)
lncRNA_ATAC_refRNA_list <- split(lncRNA_ATAC_refRNA, lncRNA_ATAC_refRNA$transcript_id)
unlist(lapply(lncRNA_ATAC_refRNA_list,nrow))




