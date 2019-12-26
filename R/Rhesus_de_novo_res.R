
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

TACO_gtf_trans_save <- data.frame(GeneID=rownames(TACO_gtf_trans), 
                            TACO_gtf_trans)
write.table(TACO_gtf_trans_save, "SourceData/SourceData_Fig.4b-d.txt",
            row.names = F, col.names = T, quote = F, sep="\t")


###
load("revise0530/de_novo_TACO_minExpr5.5_res/TACO_gtf_trans.RData")
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


####
# novel transcripts
novel <- TACO_gtf_trans[!(TACO_gtf_trans$category_relative_detail == "same_strand" & 
                            TACO_gtf_trans$shared_splicing == "True"),]
novel_id <- rownames(novel)
write.table(novel_id, "revise0530/de_novo_TACO_minExpr5.5_res/novel_transcripts_id.txt",
            row.names = F, col.names = F, quote = F, sep="\t")


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
venn(list(
  All_trans=rownames(TACO_gtf_trans),
  CPC2_coding=rownames(TACO_gtf_trans)[TACO_gtf_trans$CPC2_potential=="coding"],
  CPAT_coding=rownames(TACO_gtf_trans)[TACO_gtf_trans$CPAT_potential=="coding"]
),
# size = 30,
borders=F,
col=c("#778899","#3F60AC", "#8B2052"),
zcolor = c("#778899","#3F60AC", "#8B2052"),
cexil = 1.5,	# Character expansion for the intersection labels
cexsn = 1.5 #Character expansion for the set names,
)
dev.off()




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
intragenic <- sub_bed[which(sub_bed$category_relative=="intragenic"), 1:3]
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


