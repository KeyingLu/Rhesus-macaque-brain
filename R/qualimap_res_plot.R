library(ggplot2)
library(reshape2)
library(wesanderson) # names(wes_palettes)
options(stringAsFactors=FALSE)


setwd("/mnt/xdlab1/home/looking/brain_project")

data <- read.table("read_distribution_by_qualimap/genomic_origin.txt", 
                   header=T,sep="\t",stringsAsFactors=F)
str(data)
head(data)


p <- matrix(0, ncol=3, nrow=3)
colnames(p) <- c("Feature", "Percent", "se")
j=1
for(r in unique(data$Feature)){
    mean <- mean(data$Percent[data$Feature==r])
    se <- sd(data$Percent[data$Feature==r])
    p[j, ] <- c(r, mean, se)
    j = j+1
}

p <- data.frame(p)
p$Percent <- as.numeric(as.character(p$Percent))
p$se <- as.numeric(as.character(p$se))
head(p)
str(p)
ggplot(p, aes(x=Feature, y=Percent, fill=Feature)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text = element_text(color = "black", face="bold", size=12),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5,color = "black"),
        text = element_text(face="bold")) +
  scale_fill_manual(values = c("#CD8463", "#789EC1", "#AFC88C")) +
  geom_errorbar((aes(ymin=Percent-se, ymax=Percent+se))) +
  guides(fill=F) +
  # ggtitle("Reads Genomic origin") +
  xlab("") +
  ylab("the Percent of Reads Genomic Origin (%)") +
  scale_x_discrete(limits=c("exonic", "intronic", "intergenic"))

ggsave("read_distribution_by_qualimap/Reads_genomic_origin_percent.pdf", width = 3, height = 5)







#######################################################
#                       rRNA                          #
#######################################################
library(ggplot2)
library(reshape2)
library(wesanderson) # names(wes_palettes)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")



data <- read.table("read_rRNA_by_qualimap/reads_rRNA_fraction.txt", 
                   header=T,sep="\t",stringsAsFactors=F)
str(data)
head(data)


p <- matrix(0, ncol=3, nrow=3)
colnames(p) <- c("Feature", "Percent", "se")
j=1
for(r in unique(data$Feature)){
  mean <- mean(data$Percent[data$Feature==r])
  se <- sd(data$Percent[data$Feature==r])
  p[j, ] <- c(r, mean, se)
  j = j+1
}

p <- data.frame(p)
p$Percent <- as.numeric(as.character(p$Percent))
p$se <- as.numeric(as.character(p$se))
head(p)
str(p)



###
pp <- data[data$Feature=="exonic", ]
rownames(pp) <- pp$sample_id

pdf("revise0530/gene_count_res/RNA_mapping_rRNA_rate_hist.pdf")
hist(pp$Percent)
abline(v=0.94, col="red")
abline(v=0.81, col="red")
dev.off()



####
data <- data[data$Feature=="exonic",]
rownames(data) <- paste("X", data$sample_id, sep="")
data <- data[rownames(data) %in% rownames(Rhesus_colData),]
data <- data.frame(data, Rhesus_colData[, c("Age_Stage", "Sex")])

t.test(Percent~Age_Stage, data=data) # p-value = 0.1331
t.test(Percent~Sex, data=data) # 0.9401

ggplot(data, aes(x=Age_Stage, y=Percent, colour=Age_Stage, fill=Age_Stage)) +
  # geom_violin()+
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.7, color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(color = "black"),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white", size=1)) +
  scale_colour_manual(values = c(brewer.pal(9,"Blues")[8], brewer.pal(9,"Reds")[8])) +
  scale_fill_manual(values = c(brewer.pal(9,"Blues")[8], brewer.pal(9,"Reds")[8])) +
  xlab("") +
  guides(fill=F, colour=F)

ggsave("revise0530/gene_count_res/rRNA_rate_Age_boxplot.pdf", width = 4, height = 4)




#######################################################
#                       novel gtf                     #
#######################################################
###########  
library(ggplot2)
library(reshape2)
library(wesanderson) # names(wes_palettes)
options(stringAsFactors=FALSE)


setwd("/mnt/data2/Rhesus_brain")

data <- read.table("read_novel_gtf_by_qualimap/reads_novel_gtf_fraction.txt", 
                   header=T,sep="\t",stringsAsFactors=F)
str(data)
head(data)


p <- matrix(0, ncol=3, nrow=3)
colnames(p) <- c("Feature", "Percent", "se")
j=1
for(r in unique(data$Feature)){
  mean <- mean(data$Percent[data$Feature==r])
  se <- sd(data$Percent[data$Feature==r])
  p[j, ] <- c(r, mean, se)
  j = j+1
}

p <- data.frame(p)
p$Percent <- as.numeric(as.character(p$Percent))
p$se <- as.numeric(as.character(p$se))
head(p)
str(p)

label <- paste(round(p$Percent,2), "±", round(p$se,2), sep="")
ggplot(p, aes(x=Feature, y=Percent, fill=Feature)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text = element_text(color = "black", face="bold", size=12),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5,color = "black"),
        text = element_text(face="bold")) +
  scale_fill_manual(values = c("#CD8463", "#789EC1", "#AFC88C")) +
  geom_errorbar((aes(ymin=Percent-se, ymax=Percent+se))) +
  guides(fill=F) +
  # ggtitle("Reads Genomic origin") +
  xlab("") +
  geom_text(aes(label=label), vjust=-1) + 
  ylab("the Percent of Reads Genomic Origin (%)") +
  scale_x_discrete(limits=c("exonic", "intronic", "intergenic"))

ggsave("read_novel_gtf_by_qualimap/Reads_genomic_origin_percent_for_novel_gtf.pdf",
       width = 5, height = 5)











