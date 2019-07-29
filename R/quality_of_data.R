####################################################
#                   RNA mapping                    #
####################################################
library(plotrix)
library(ggplot2)
library(RColorBrewer)


# .libPaths("/home/looking/R/x86_64-pc-linux-gnu-library/3.4")
setwd("/mnt/data2/Rhesus_brain")

mapInfo <- read.table("bam_stat/mapping_stat.txt", header = F)
mapInfo <- mapInfo[,c(1,2,7)]
colnames(mapInfo) <- c("sample_id", "properly_mapping_reads", "properly_mapping_rate")
mapInfo$properly_mapping_rate <- gsub("[(]|%", "", mapInfo$properly_mapping_rate)
mapInfo$properly_mapping_rate <- as.numeric(mapInfo$properly_mapping_rate)


order <- order(mapInfo$properly_mapping_reads, decreasing=T)
mapInfo <- mapInfo[order,]
mapInfo$properly_mapping_reads <- mapInfo$properly_mapping_reads/10000000
head(mapInfo)
tail(mapInfo)



### properly_mapping_reads
ggplot(mapInfo, aes(x= properly_mapping_reads)) +
  geom_histogram(color="#247BA0", fill="white") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  xlab("properly mapping reads (x10^7)") 
ggtitle("") 

ggsave("revise0530/gene_count_res/properly_mapping_reads_hist.pdf",
       width=5, height = 5)


###
ggplot(mapInfo, aes(x= properly_mapping_rate)) +
  geom_histogram(color="#F25F5C", fill="white") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  xlab("properly mapping rate (%)") 
ggtitle("") 

ggsave("revise0530/gene_count_res/properly_mapping_rate_hist.pdf",
       width=5, height = 5)

# samples <- mapInfo$sample_id
# reads <- mapInfo$properly_mapping_reads
# rate <- mapInfo$properly_mapping_rate
# 
# pdf("bam_new/mapping_stat.pdf", width=15, height=5)
# xpos <-1:416
# twoord.plot(xpos,reads,xpos,rate,xlim=c(0,416),lylim=c(0,22),rylim=c(0,100), 
#             lcol=brewer.pal(9,"Blues")[8],rcol=brewer.pal(9,"Reds")[8],xlab="sample_id",
#             ylab="properly mapping reads(x10^7)",
#             rylab="properly mapping rate(%)",type=c("bar","p"),
#             xticklab=samples,xtickpos=,
#             halfwidth=0.2,
#             mar=c(4,4,4,4))
# dev.off()




mapInfo <- mapInfo[order(mapInfo$sample_id),]
head(mapInfo)
write.table(mapInfo, "bam_new/RNA_mapping_stat.txt", row.names = F, 
            col.names = T, quote = F, sep="\t")



### histogram
rownames(mapInfo) <- mapInfo$sample_id
pdf("revise0530/gene_count_res/RNA_mapping_rate_hist.pdf")
hist(mapInfo$properly_mapping_rate)
abline(v=80.13, col="red")
abline(v=80.14, col="red")
dev.off()


pdf("revise0530/gene_count_res/RNA_mapping_reads_hist.pdf")
hist(mapInfo$properly_mapping_reads)
abline(v=5.935121, col="red")
abline(v=8.105072, col="red")
dev.off()



####################################################
#                    RNA RIN                       #
####################################################
library(ggplot2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")


Rhesus_RIN <- read.table("revise0530/RNA_RIN.txt", header = F, 
                         row.names = 1, stringsAsFactors = F)

mapInfo$RIN <- Rhesus_RIN[rownames(mapInfo),1]


pdf("revise0530/gene_count_res/correlation_between_RIN_and_mapping_rate.pdf")
plot(x=mapInfo$RIN, y=mapInfo$properly_mapping_rate)
dev.off()


###
Rhesus_FPKM <- read.table("matrix0420/Rhesus_FPKM.matrix", 
                          header=T, row.names = 1)
dim(Rhesus_FPKM) # 30807   408



n1 <- apply(Rhesus_FPKM, 2, function(x){sum(log2(x+1)>1)})
n10 <- apply(Rhesus_FPKM, 2, function(x){sum(log2(x+1)>10)})


cor.test(Rhesus_colData[colnames(Rhesus_FPKM), "RIN"], n1)
cor.test(Rhesus_colData[colnames(Rhesus_FPKM), "RIN"], n10)


pp <- data.frame(n1,n10,RIN=Rhesus_colData[colnames(Rhesus_FPKM), "RIN"])

ggplot(pp, aes(x=RIN)) +
  geom_point(aes(y=n1), color="#3385BB") +
  geom_point(aes(y=n10), color="#FE8D19") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  xlab("RIN") +
  ylab("the number of genes")



ggsave("revise0530/gene_count_res/RIN_RPKM_point.pdf")



####################################################
#                    RNA TIN                       #
####################################################
library(ggplot2)
library(reshape2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")


TIN <- read.table("revise0530/Rhesus_TIN.matrix", header=T, row.names = 1)


#### the distribution of TIN
pp <- melt(TIN)
colnames(pp) <- c("sample_id", "TIN")

ggplot(pp, aes(x= TIN, color= sample_id)) +
  geom_line(stat = "density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  # facet_wrap(~region, ncol=4) +
  # scale_color_manual(values=c("#F25F5C", "#247BA0")) +
  xlab("pvalue") +
  guides(color=F)


ggsave("revise0530/gene_count_res/TIN_density.pdf",
       height = 7, width = 11)



### medTIN
library(matrixStats)
load("revise0530/gene_count_res/Rhesus_colData.RData")


medTIN <- colMedians(as.matrix(TIN))
names(medTIN) <- colnames(TIN)

pp <- data.frame(medTIN=medTIN[rownames(Rhesus_colData)], 
                 RIN=Rhesus_colData$RIN)


ggplot(pp, aes(x= medTIN, y=RIN)) +
  geom_point(color="#8F908F") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  ggtitle("the Correlation between RIN and medTIN") +
  guides(color=F)


ggsave("revise0530/gene_count_res/the_Correlation_between_RIN_and_medTIN.pdf",
       height = 5, width = 5)



cor(pp$medTIN, pp$RIN)




####################################################
#               save RNA quality                   #
####################################################
head(Rhesus_RIN)
head(medTIN)
head(mapInfo)



quality <- Rhesus_RIN
colnames(quality) <- "RIN"
names(medTIN) <- gsub("X", "", names(medTIN))
quality$medTIN <- medTIN[rownames(quality)]
rownames(mapInfo) <- mapInfo$sample_id
quality <- data.frame(quality, mapInfo[rownames(quality),c(2:3)])


map <- setNames(c("Rhesus_1", "Rhesus_2", "Rhesus_3", "Rhesus_4", "Rhesus_5", "Rhesus_6", "Rhesus_7", "Rhesus_8"),
                c("2", "3", "4", "6", "7", "8", "10", "11"))
individual <- map[sapply(rownames(quality), function(x){strsplit(x, "_")[[1]][1]})]

quality <- data.frame(sample_id=rownames(quality), Individual=individual, quality)

write.table(quality, "revise0530/gene_count_res/RNA_qulity_stat.txt", 
            row.names = F,  col.names = T, quote = F, sep="\t")




####################################################
#                 ATAC mapping                     #
####################################################
setwd("/mnt/data2/Rhesus_brain")
mapInfo <- read.table("ATAC/bam/ATAC_mapping_info", stringsAsFactors = F)
a <- rep(c("origin", "deduplicated", "selected"), length=dim(mapInfo)[1])
mapInfo$type <- a

pp <- data.frame(sample_id=mapInfo$V1[(1:26)*3], 
                 properly_paired_reads=mapInfo$V2[(1:26)*3-2],
                 properly_paired_deduplicated_reads=mapInfo$V2[(1:26)*3-1],
                 properly_paired_selected_reads=mapInfo$V2[(1:26)*3],
                 properly_paired_mapping_rate=mapInfo$V3[(1:26)*3-2],
                 properly_paired_deduplicated_mapping_rate=mapInfo$V3[(1:26)*3-1],
                 properly_paired_selected_mapping_rate=mapInfo$V3[(1:26)*3])
pp <- pp[order(pp$sample_id),]
head(pp)
write.table(pp, "ATAC/bam/ATAC_mapping_stat.txt", row.names = F, 
            col.names = T, quote = F, sep="\t")




#####
setwd("/mnt/data2/Rhesus_brain")
mapInfo <- read.table("ATAC/bam/ATAC_mapping_stat.txt", header=T, stringsAsFactors = F)
mapInfo[,c(2,3,4)] <- mapInfo[,c(2,3,4)]/10000000



pdf("ATAC/bam/ATAC_mapping_stat.pdf", width = 8, height = 5)
par(mar = c(5, 5, 3, 4)+0.1)
plot('sample_id', 'properly mapping reads(x10^7)', xlim = c(0,6), ylim = c (0,6), type = "n", axes = F)
box()
par(new=TRUE)
sample_id<-1:26
plot(sample_id, mapInfo$properly_paired_reads, 
     pch=16, axes=FALSE, ylim=c(0,20), xlab="", ylab="", 
     type="l",col="#1B4E91",lwd=3)
par(new=TRUE)
plot(sample_id, mapInfo$properly_paired_deduplicated_reads, 
     pch=16, axes=FALSE, ylim=c(0,20), xlab="", ylab="", 
     type="o",col="#1B4E91",lwd=3)

par(new=TRUE)
plot(sample_id, mapInfo$properly_paired_selected_reads, 
     pch=16, axes=FALSE, ylim=c(0,20), xlab="", ylab="", 
     type="l",lty=3,col="#1B4E91",lwd=3)

axis(2, ylim=c(0,20),col="#1B4E91",las=1)

par(new=TRUE)
plot(sample_id, mapInfo$properly_paired_mapping_rate, pch=16, axes=FALSE, ylim=c(80,100), xlab="", ylab="", 
     type="l",col="#991B27",lwd=3)

par(new=TRUE)
plot(sample_id, 
     mapInfo$properly_paired_deduplicated_mapping_rate, 
     pch=16, axes=FALSE, ylim=c(80,100), xlab="", ylab="", 
     type="o",col="#991B27",lwd=3)
par(new=TRUE)
plot(sample_id, 
     mapInfo$properly_paired_selected_mapping_rate, 
     pch=16, axes=FALSE, ylim=c(80,100), xlab="", ylab="", 
     type="l",lty=3,col="#991B27",lwd=3)


# mtext("mapping_rate",side=4,col="#991B27",line=4)
axis(4, ylim=c(80,100), col="#991B27",col.axis="#991B27",las=1)
axis(1, sample_id, labels = mapInfo$sample_id, las=2)
# mtext("sample_id",side=1,col="#1B4E91",line=2.3)


legend("bottomleft",
       legend=c("properly_paired_reads","properly_paired_deduplicated_reads",
                "properly_paired_selected_reads","properly_paired_mapping_rate",
                "properly_paired_deduplicatecd_mapping_rate","properly_paired_selected_mapping_rate"),
       text.col=c("#1B4E91","#1B4E91","#1B4E91","#991B27","#991B27","#991B27"),
       lty=c(1,1,3,1,1,3),
       pch=c(NA,16,NA,NA,16,NA),
       col=c("#1B4E91","#1B4E91","#1B4E91","#991B27","#991B27","#991B27"),
       ncol=2)

dev.off()
