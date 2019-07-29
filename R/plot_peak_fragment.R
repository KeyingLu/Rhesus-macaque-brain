library(ggplot2)
library(reshape2)

setwd("/mnt/data2/Rhesus_brain")

sample_list <- read.table("ATAC/fragment/sample_list.txt",header=F, stringsAsFactors = F)
sample_list <- sample_list$V1

fragment_list <- list()
for(sample_id in sample_list){
  filename = paste("ATAC/fragment/", sample_id, ".fragment", sep="")
  system(paste("awk '$4=$3-$2 {print $4}' ", filename, " > ATAC/fragment/tmp", sep=""))
  fragment <- read.table("ATAC/fragment/tmp", header = F, stringsAsFactors = F)
  fragment_list[[sample_id]] <- fragment
}

pp <- melt(fragment_list)
pp <- pp[,2:3]
colnames(pp) <- c("length", "sample_id")
str(pp)

ggplot(pp, aes(x=length, colour=sample_id))+
  geom_line(stat="density") + 
  theme_bw() +
  scale_x_continuous(limits=c(0, 1000))


ggsave("ATAC/fragment/all_sample_fragment_density.pdf",
       width = 5, height = 4.5)




