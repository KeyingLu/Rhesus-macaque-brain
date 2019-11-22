# library("lme4")
library(DESeq2)
# library(MASS)
library(mgcv)
library(BiocParallel)
options(stringAsFactors=FALSE)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Age_res.RData")
load("revise0530/gene_count_res/Sex_res.RData")
str(Rhesus_colData)



#######################################################
#                   gam results                       #
#######################################################
colData <- Rhesus_colData
gam_func <- function(r){
  library(mgcv)
  load("revise0530/gene_count_res/RhesusGeneCount.RData")
  load("revise0530/gene_count_res/Rhesus_colData.RData")
  colData <- Rhesus_colData
  sub_colData <- colData[colData$SR_merge == r,]
  df <- RhesusGeneCount[, rownames(sub_colData)]
  df <- df[apply(df, 1, function(x){sum(x>1)>3}), ]
  ####
  Sex_pvalue <- c()
  Age_pvalue <- c()
  Sex_coef <- c()
  Age_coef <- c()
  for(j in 1:dim(df)[1]){
    gene_id <- rownames(df)[j]
    GeneExpr = t(df[j,])
    colnames(GeneExpr) <- c("GeneExpr")           
    dd <- data.frame(GeneExpr, sub_colData)
    dd$ID <- factor(dd$ID)
    ####
    coef=tryCatch({
      if(length(unique(sub_colData$SR))>1){
        m.nb <- gam(GeneExpr ~ Sex + Age_Stage + SR + s(ID, bs = "re"),
                    data=dd,
                    family = "nb")
      }else{ m.nb <- gam(GeneExpr ~ Sex + Age_Stage ,
                         data=dd,  
                         family = "nb" )
      }
      coef <- c(summary(m.nb)$p.pv, summary(m.nb)$p.coeff)
    },error=function(e){
      cat("ERROR :",conditionMessage(e), "\n")
      cat(j,unique(colData$SR_merge))
      coef <- c(NA, NA, NA, NA, NA, NA)
      return(coef)})
    ####
    Sex_pvalue <- c(Sex_pvalue, coef[2])
    Age_pvalue <- c(Age_pvalue, coef[3])
    Sex_coef <- c(Sex_coef, coef[5])
    Age_coef <- c(Age_coef, coef[6])
  }
  value = data.frame(Sex=Sex_pvalue, Age=Age_pvalue, 
                     Sex_coef=Sex_coef, Age_coef=Age_coef)
  rownames(value) <- rownames(df)
  return(value)
}


param <- SnowParam(workers = 16, type = "SOCK")
gam_res  <- bplapply(unique(colData$SR_merge), gam_func, BPPARAM = param)
names(gam_res) <- unique(colData$SR_merge)


save(gam_res, 
     file = "revise0530/gene_count_res/gam_res.RData")



pp <- melt(gam_res)
ggplot(pp, aes(x= value, fill= variable, color=variable)) +
  geom_histogram(position="identity", alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  facet_wrap(~L1, ncol=4, scales = "free") +
  scale_color_manual(values=c("#1B9E77", "#D95F02")) +
  scale_fill_manual(values=c("#1B9E77",  "#D95F02")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue (gam)") 


ggsave("revise0530/gene_count_res/gam_pvalue_hist.pdf",
       height = 7, width = 11)






##########
load("revise0530/gene_count_res/Age_res.RData")
load("revise0530/gene_count_res/Sex_res.RData")

gam_genes <- list()
for(r in names(gam_res)){
  pvalue_df <- gam_res[[r]]
  for(condition in colnames(pvalue_df)){
    if(condition == "Sex"){
      res <- results(Sex_res[[r]],
                     contrast = c("Sex", c("male", "female")))
      res <- data.frame(pvalue=pvalue_df[, "Sex"], 
                        res[rownames(pvalue_df), "log2FoldChange",drop=F])
    }
    if(condition == "Age"){
      res <- results(Age_res[[r]],
                     contrast = c("Age_Stage", c("Mid", "Young")))
      res <- data.frame(pvalue=pvalue_df[, "Age"], 
                        res[rownames(pvalue_df), "log2FoldChange",drop=F])
    }
    #####
    res <- res[complete.cases(res),]
    hyper <- res[res$pvalue<0.05 & res$log2FoldChange>1, ]
    hypo <- res[res$pvalue<0.05 & res$log2FoldChange < -1, ]
    if(condition == "Sex"){
      gam_genes[[r]][[condition]] <- list(male=hyper, female=hypo)
    }
    if(condition == "Age"){
      gam_genes[[r]][[condition]] <- list(Mid=hyper, Young=hypo)
    }
}
}


####
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")

Sex_genes_InEachSR <- lapply(gam_genes, function(x){y <- x[["Sex"]]; list(male=rownames(y[["male"]]), female=rownames(y[["female"]]))})
Sex_genes_count <- melt(Sex_genes_InEachSR)
colnames(Sex_genes_count) <- c("gene", "Sex", "SR")
Sex_genes_count1 <- table(Sex_genes_count$SR, Sex_genes_count$Sex)
Sex_genes_count1 <- melt(Sex_genes_count1)
colnames(Sex_genes_count1) <- c("SR", "Sex", "counts")
head(Sex_genes_count1)

Sex_genes_count2 <- table(Sex_genes_count$SR)
Sex_genes_count2 <- Sex_genes_count2[order(Sex_genes_count2, decreasing=T)]
head(Sex_genes_count2)
Sex_genes_count1$SR <- factor(Sex_genes_count1$SR, levels=(names(Sex_genes_count2)))


### plot
pp <- Sex_genes_count1
ggplot(pp,aes(x=SR, y=counts, fill=Sex)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  scale_fill_manual(values = anno_colors$Sex) +
  ggtitle("the counts of DEG for Sex")

ggsave("revise0530/gene_count_res/Sex_DEG_counts_bar_InSR_gam.pdf",
       width=6, height = 5)


Age_genes_InEachSR <- lapply(gam_genes, function(x){y <- x[["Age"]]; list(Mid=rownames(y[["Mid"]]), Young=rownames(y[["Young"]]))})
Age_genes_count <- melt(Age_genes_InEachSR)
colnames(Age_genes_count) <- c("gene", "Age", "SR")
Age_genes_count1 <- table(Age_genes_count$SR, Age_genes_count$Age)
Age_genes_count1 <- melt(Age_genes_count1)
colnames(Age_genes_count1) <- c("SR", "Age", "counts")
head(Age_genes_count1)

Age_genes_count2 <- table(Age_genes_count$SR)
Age_genes_count2 <- Age_genes_count2[order(Age_genes_count2, decreasing=T)]
head(Age_genes_count2)
Age_genes_count1$SR <- factor(Age_genes_count1$SR, levels=(names(Age_genes_count2)))


### plot
pp <- Age_genes_count1
ggplot(pp,aes(x=SR, y=counts, fill=Age)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = anno_colors$Age_Stage) +
  ggtitle("the counts of DEG for age")


ggsave("revise0530/gene_count_res/Age_DEG_counts_bar_InSR_gam.pdf",
       width=6, height = 5)


############
load("revise0530/gene_count_res/gam_res.RData")
pvalue <- gam_res[["MO"]]
pvalue <- data.frame(gene_id=rownames(pvalue), pvalue)

write.table(pvalue, "revise0530/gene_count_res/MO_pvalue.txt",
            col.names = T, row.names = F, quote = F, sep = "\t")




#######################################################
#                   glmmTMB results                   #
#######################################################
rm(list=ls())
library(snow)
library(DESeq2)
options(stringAsFactors=FALSE)
set.seed(100)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
str(Rhesus_colData)
clus <- makeCluster(16)


glmmTMB_res_func <- function(r){
  library(glmmTMB)
  load("revise0530/gene_count_res/Rhesus_colData.RData")
  load("revise0530/gene_count_res/RhesusGeneCount.RData")
  colData <- Rhesus_colData
  sub_colData <- colData[colData$SR_merge == r,]
  df <- RhesusGeneCount[, rownames(sub_colData)]
  df <- df[apply(df, 1, function(x){sum(x>1)>3}), ]
  ####
  Sex_pvalue <- c()
  Age_pvalue <- c()
  ###
  for(i in 1:dim(df)[1]){
    gene_id <- rownames(df)[i]
    GeneExpr = t(df[i,])
    colnames(GeneExpr) <- c("GeneExpr")           
    dd <- data.frame(GeneExpr, sub_colData, stringsAsFactors = F)
    dd$ID <- factor(dd$ID)
    coef=tryCatch({
      if(length(unique(dd$SR))>1){
        m.nb <- glmmTMB(GeneExpr ~ Sex + Age_Stage + SR + (1|ID),  family = nbinom2,  data=dd)
      }else{
        m.nb <- glmmTMB(GeneExpr ~ Sex + Age_Stage ,  family = nbinom2,  data=dd)
      }
      coef <- summary(m.nb)$coefficients$cond[,4]
    },error=function(e){
      cat("ERROR :",conditionMessage(e), "\n")
      coef <- c("NA", "NA", "NA")
      return(coef)})
    Sex_pvalue <- c(Sex_pvalue, coef[2])
    Age_pvalue <- c(Age_pvalue, coef[3])
  }
  ##
  pvalue = data.frame(Sex=Sex_pvalue, Age=Age_pvalue)
  rownames(pvalue) <- rownames(df)
  return(pvalue) 
}



clusterExport(clus, "glmmTMB_res_func")
region_list <- list()
for(r in unique(Rhesus_colData$SR_merge)){region_list[[r]] <- r}

glmmTMB_res <- parLapply(clus,
                         region_list,
                         function(x) glmmTMB_res_func(x))



save(glmmTMB_res, 
     file = "revise0530/gene_count_res/glmmTMB_res.RData")




################################################################
#                    locfdr + fdrTbL                           #
################################################################
rm(list = ls())
library(BiocParallel)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/gam_res.RData")


########## adjusted
adjust_func <- function(r){
  library(locfdr)
  library(fdrci)
  load("revise0530/gene_count_res/Permutaion18_gam_res.RData")
  load("revise0530/gene_count_res/gam_res.RData")
  permute_list <- Permutaion18_gam_res[[r]]
  #### Sex
  Sex_null <- list()
  for(i in names(permute_list)[1:18]){
    p <- permute_list[[i]][,"Sex"]
    p[p > 0.9999999999] <- 0.9999999999
    p[p < 1e-10] <- 1e-10
    z <- qnorm(p)
    # z <- permute_list[[i]][,"Sex_zvalue"]
    fdr <- locfdr(z, bre = 1000)
    a <- fdr$mat
    b <- sample(a[,1], 100000, prob = a[,6], replace = T)
    Sex_null[[i]] <- data.frame(pvalue=pnorm(b))
  }
  observe_pvalue <- gam_res[[r]][, "Sex"]
  Sex_tbl <- fdrTbl(observe_pvalue, Sex_null,
                    "pvalue", length(observe_pvalue), 1, 50)
  #### Age
  Age_null <- list()
  for(i in names(permute_list)[19:36]){
    p <- permute_list[[i]][,"Age"]
    p[p > 0.9999999999] <- 0.9999999999
    p[p < 1e-10] <- 1e-10
    z <- qnorm(p)
    # z <- permute_list[[i]][,"Age_zvalue"]
    fdr <- locfdr(z, bre = 1000)
    a <- fdr$mat
    b <- sample(a[,1], 100000, prob = a[,6], replace = T)
    Age_null[[i]] <- data.frame(pvalue=pnorm(b))
  }
  observe_pvalue <- gam_res[[r]][, "Age"]
  Age_tbl <- fdrTbl(observe_pvalue, Age_null,
                    "pvalue", length(observe_pvalue), 1, 50)
  #####
  return(list(Sex=Sex_tbl, Age=Age_tbl))
}


param <- SnowParam(workers = 16, type = "SOCK")
locfdr_fdrTbL_res  <- bplapply(names(gam_res), adjust_func, BPPARAM = param)
names(locfdr_fdrTbL_res) <- names(gam_res)


save(locfdr_fdrTbL_res, 
     file="revise0530/gene_count_res/locfdr_fdrTbL_res.RData")




##########    gene filter
library(DESeq2)
load("revise0530/gene_count_res/Age_res.RData")
load("revise0530/gene_count_res/Sex_res.RData")
load("revise0530/gene_count_res/gam_res.RData")


Age_locfdr_fdrTbL_genes <- list()
Sex_locfdr_fdrTbL_genes <- list()
for(r in names(locfdr_fdrTbL_res)){
  fdrres <- locfdr_fdrTbL_res[[r]]
  for(condition in names(fdrres)){
    sub_fdrres <- fdrres[[condition]]
    if(condition == "Sex"){
      res <- results(Sex_res[[r]],
                     contrast = c("Sex", c("male", "female")))
      res <- data.frame(pvalue=gam_res[[r]][rownames(res), "Sex"], 
                        res[, "log2FoldChange",drop=F])
    }
    if(condition == "Age"){
      res <- results(Age_res[[r]],
                     contrast = c("Age_Stage", c("Mid", "Young")))
      res <- data.frame(pvalue=gam_res[[r]][rownames(res), "Age"], 
                        res[, "log2FoldChange",drop=F])
    }
    #####
    res <- res[complete.cases(res),]
    res$logP <- -log10(res$pvalue)
    res <- res[order(res$logP, decreasing = T),]
    #####
    sub_fdrres <- sub_fdrres[complete.cases(sub_fdrres),]
    sub_fdrres <- sub_fdrres[sub_fdrres$fdr < 0.1,]
    if(dim(sub_fdrres)[1]==0){print(c(r, condition)); next}
    else{
      logP <- sub_fdrres[1, "threshold"]
      ##
      hyper <- res[res$logP>logP & res$log2FoldChange>1, ]
      hypo <- res[res$logP>logP & res$log2FoldChange < -1, ]
      if(condition == "Sex"){
        Sex_locfdr_fdrTbL_genes[[r]] <- list(male=hyper, female=hypo)
      }
      if(condition == "Age"){
        Age_locfdr_fdrTbL_genes[[r]] <- list(Mid=hyper, Young=hypo)
      }
    }
  }
}

save(Sex_locfdr_fdrTbL_genes, Age_locfdr_fdrTbL_genes, 
     file="revise0530/gene_count_res/locfdr_fdrTbL_genes.RData")


### write
load("revise0530/gene_count_res/locfdr_fdrTbL_genes.RData")
load("revise0530/gene_count_res/gam_res.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")


### Age
Age_table <- list()
for(r in names(Age_locfdr_fdrTbL_genes)){
  res <- gam_res[[r]]
  Mid <- Age_locfdr_fdrTbL_genes[[r]][["Mid"]][, 1:2]
  Young <- Age_locfdr_fdrTbL_genes[[r]][["Young"]][, 1:2]
  coef <- res[, c("Sex_coef", "Age_coef")]
  if(dim(Mid)[1]>0){Mid <- data.frame(region=r, hyper="Mid", Mid, coef[rownames(Mid),, drop=F])}
  if(dim(Young)[1]>0){Young <- data.frame(region=r, hyper="Young", Young, coef[rownames(Young),, drop=F])}
  dif <- rbind(Mid, Young)
  dif <- data.frame(gene=rownames(dif), gene_name=gtf_ensembl_gene[rownames(dif), "gene_name"], dif)
  Age_table[[r]] <- dif
}


Age_table <- Reduce(rbind, Age_table)
write.table(Age_table, "revise0530/gene_count_res/Age_table.txt",
            row.names = F, col.names = T, sep="\t", quote = F)


### Sex
##
Sex_table <- list()
for(r in names(Sex_locfdr_fdrTbL_genes)){
  res <- gam_res[[r]]
  male <- Sex_locfdr_fdrTbL_genes[[r]][["male"]][, 1:2]
  female <- Sex_locfdr_fdrTbL_genes[[r]][["female"]][, 1:2]
  coef <- res[, c("Sex_coef", "Age_coef")]
  if(dim(male)[1]>0){male <- data.frame(region=r, hyper="male", male, coef[rownames(male),, drop=F])}
  if(dim(female)[1]>0){female <- data.frame(region=r, hyper="female", female, coef[rownames(female),, drop=F])}
  dif <- rbind(male, female)
  dif <- data.frame(gene=rownames(dif), gene_name=gtf_ensembl_gene[rownames(dif), "gene_name"], dif)
  Sex_table[[r]] <- dif
}


Sex_table <- Reduce(rbind, Sex_table)
write.table(Sex_table, "revise0530/gene_count_res/Sex_table.txt",
            row.names = F, col.names = T, sep="\t", quote = F)








####
library(ggplot2)
library(reshape2)
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")

## Sex
Sex_genes_InEachSR <- lapply(Sex_locfdr_fdrTbL_genes, function(x){list(male=rownames(x[["male"]]), female=rownames(x[["female"]]))})
Sex_genes_count <- melt(Sex_genes_InEachSR)
colnames(Sex_genes_count) <- c("gene", "Sex", "SR")
Sex_genes_count1 <- table(Sex_genes_count$SR, Sex_genes_count$Sex)
Sex_genes_count1 <- melt(Sex_genes_count1)
colnames(Sex_genes_count1) <- c("SR", "Sex", "counts")
head(Sex_genes_count1)

Sex_genes_count2 <- table(Sex_genes_count$SR)
Sex_genes_count2 <- Sex_genes_count2[order(Sex_genes_count2, decreasing=T)]
head(Sex_genes_count2)
Sex_genes_count1$SR <- factor(Sex_genes_count1$SR, levels=(names(Sex_genes_count2)))


### plot
pp <- Sex_genes_count1
ggplot(pp,aes(x=SR, y=counts, fill=Sex)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  scale_fill_manual(values = anno_colors$Sex) +
  ggtitle("the counts of DEG for Sex")

ggsave("revise0530/gene_count_res/Sex_DEG_counts_bar_InSR_locfdr_fdrTbL.pdf",
       width=6, height = 5)

#### Age
Age_genes_InEachSR <- lapply(Age_locfdr_fdrTbL_genes, function(x){list(Mid=rownames(x[["Mid"]]), Young=rownames(x[["Young"]]))})
Age_genes_count <- melt(Age_genes_InEachSR)
colnames(Age_genes_count) <- c("gene", "Age", "SR")
Age_genes_count1 <- table(Age_genes_count$SR, Age_genes_count$Age)
Age_genes_count1 <- melt(Age_genes_count1)
colnames(Age_genes_count1) <- c("SR", "Age", "counts")
head(Age_genes_count1)

Age_genes_count2 <- table(Age_genes_count$SR)
Age_genes_count2 <- Age_genes_count2[order(Age_genes_count2, decreasing=T)]
head(Age_genes_count2)
Age_genes_count1$SR <- factor(Age_genes_count1$SR, levels=(names(Age_genes_count2)))


### plot
pp <- Age_genes_count1
ggplot(pp,aes(x=SR, y=counts, fill=Age)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = anno_colors$Age_Stage) +
  ggtitle("the counts of DEG for age")


ggsave("revise0530/gene_count_res/Age_DEG_counts_bar_InSR_locfdr_fdrTbL.pdf",
       width=6, height = 5)



#############################################################
#         Age(Sex) genes function (after permutations)      #
#############################################################
library(gProfileR)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive2/r1750_e91_eg38/web")
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/locfdr_fdrTbL_genes.RData")
# Sex_locfdr_fdrTbL_genes, Age_locfdr_fdrTbL_genes


locfdr_fdrTbL_genes <- Age_locfdr_fdrTbL_genes
gprofiler_res_list <- list()
for(region in c("MO", "PIT")){
  for(age in c("Mid", "Young")){
    genes <- rownames(locfdr_fdrTbL_genes[[region]][[age]])
    gprofiler_res <- gprofiler(query=genes, 
                               organism = "mmulatta",
                               correction_method="fdr"
    )
    gprofiler_res$gene_name <- gprofiler_res$intersection
    if(dim(gprofiler_res)[1]==0) next
    else{for(i in 1:dim(gprofiler_res)[1]){
      x <- gprofiler_res$gene_name[i]
      gprofiler_res$gene_name[i] <- paste(gtf_ensembl_gene[strsplit(x,",")[[1]],"gene_name"], collapse=",")
    }
      gprofiler_res <- gprofiler_res[order(gprofiler_res$p.value),]
      filename=paste("revise0530/gene_count_res/", region,"_", age, "_hyper_function_permutation.txt", sep="")
      # write.table(gprofiler_res, filename,
      #             quote = F, sep="\t", row.names = F, col.names = T)
      gprofiler_res_list[[region]][[age]] <- gprofiler_res
    }
  }
}


#####################################################################
#         Age specific genes pheatmap (after permutations)          #
#####################################################################
library(DESeq2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/gene_count_res/locfdr_fdrTbL_genes.RData")
load("revise0530/gene_count_res/Age_res.RData")


#################
## scale function
#################
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


#####

for(region in names(Age_locfdr_fdrTbL_genes)){
  Object <- Age_res[[region]]
  vsd <- varianceStabilizingTransformation(Object, blind=FALSE)
  vsd_data <- assay(vsd)
  ####
  rowInfo <-  melt(list(Mid=rownames(Age_locfdr_fdrTbL_genes[[region]][["Mid"]]), 
                        Young=rownames(Age_locfdr_fdrTbL_genes[[region]][["Young"]])))
  if(dim(rowInfo)[1]<=1) next
  colnames(rowInfo) <- c("gene", "Age")
  rownames(rowInfo) <- as.character(rowInfo$gene)
  rowInfo$gene_name <- gtf_ensembl_gene[as.character(rowInfo$gene), "gene_name"]
  rowInfo$gene_name[which(is.na(rowInfo$gene_name))] <- as.character(rowInfo$gene)[which(is.na(rowInfo$gene_name))]
  head(rowInfo)
  ####
  colData <- colData(Object)
  colInfo <- data.frame(colData[, c("Sex", "Age_Stage")])
  head(colInfo)
  str(colInfo)
  
  ann_colors <- list(Age_Stage = anno_colors$Age_Stage,
                     Sex =anno_colors$Sex)
  
  dfheat <- vsd_data[rownames(rowInfo), rownames(colInfo)]
  rownames(dfheat) <- rowInfo$gene_name
  scale_df <- scale_mat(dfheat, "row")
  scale_df[scale_df>4] <- 3
  filename = paste("revise0530/gene_count_res/Age_DEG_profile_in_", region, "_after_permutation.pdf", sep="")
  pheatmap(
    dfheat,
    col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
    breaks = c(-3,seq(-2,2,length=20),3),
    # legend_breaks = c(-2,-1,1,2),
    filename = filename,
    main = region,
    scale="row",
    # main = main,
    width = 10,
    height = 10,
    # clustering_method = "single",
    border_color=NA,
    # fontsize = 0.7,
    annotation_col = colInfo,
    annotation_colors = ann_colors,
    gaps_row = cumsum(table(rowInfo$Age)[c("Mid", "Young")])[complete.cases(cumsum(table(rowInfo$Age)[c("Mid", "Young")]))],
    annotation_legend = T,
    # show_rownames = F,
    show_colnames = T,
    cluster_rows = F 
  )
}


########################
#   Age selected genes   
########################
library(ggplot2)
library(DESeq2)
library(reshape2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_DESeq2_object.RData")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/gtf_ensembl_gene.RData")
load("revise0530/gene_count_res/Age_res.RData")
load("revise0530/gene_count_res/locfdr_fdrTbL_genes.RData")


vsd <- varianceStabilizingTransformation(Rhesus_GeneObject, blind=FALSE)
vsd_data <- assay(vsd)

genes <- c("CHRNA1", "EBF3")
gene_ids <- gtf_ensembl_gene[gtf_ensembl_gene$gene_name %in% genes, ] 
data <- vsd_data[rownames(gene_ids),]
rownames(data) <- gene_ids$gene_name
data <- data[genes, ]
pp <- data.frame(t(data), Rhesus_colData[colnames(data),c("Age_Stage", "SR_merge")])
pp <- pp[pp$SR_merge %in% c("CB", "CC", "HTHA", "MO", "PA", "PON", "SN", "THA"),]
p <- melt(pp)
colnames(p)[3:4] <- c("gene", "expression")

ggplot(p, aes(x=SR_merge, y=expression, colour=Age_Stage, fill=Age_Stage)) +
  # geom_violin()+
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=0.7, color = "black"),
        axis.text.y = element_text(color = "black"),
        text = element_text(color = "black"),
        strip.text = element_text(size=rel(1)),
        strip.background = element_rect(fill="white", size=1)) +
  scale_colour_manual(values = c(Mid=brewer.pal(9,"Blues")[8], Young=brewer.pal(9,"Reds")[8])) +
  scale_fill_manual(values = c(Mid=brewer.pal(9,"Blues")[8], Young=brewer.pal(9,"Reds")[8])) +
  facet_wrap(~gene, ncol=3, scales="free_y") +
  scale_x_discrete(limits=c("MO", "SN", "CC", "PON", "HTHA", "THA", "PA", "CB"))+
  xlab("") 
  # guides(fill=F, colour=F)

ggsave("revise0530/gene_count_res/Age_selected_genes_expression_boxplot.pdf", width = 7, height = 3)




#####################################################################
#                  Sex specific genes pheatmap                      #
#####################################################################
library(DESeq2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")
load("revise0530/gene_count_res/locfdr_fdrTbL_genes.RData")
load("revise0530/gene_count_res/Sex_res.RData")


#################
## scale function
#################
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


#####
for(region in names(Sex_locfdr_fdrTbL_genes)){
  Object <- Sex_res[[region]]
  vsd <- varianceStabilizingTransformation(Object, blind=FALSE)
  vsd_data <- assay(vsd)
  ####
  # rowInfo <-  melt(list(male=rownames(Sex_permutation_genes[[region]][["male"]]), 
  #                       female=rownames(Sex_permutation_genes[[region]][["female"]])))
  rowInfo <-  melt(list(male=rownames(Sex_locfdr_fdrTbL_genes[[region]][["male"]]),
                        female=rownames(Sex_locfdr_fdrTbL_genes[[region]][["female"]])))
  if(dim(rowInfo)[1]<=1) next
  colnames(rowInfo) <- c("gene", "Sex")
  rownames(rowInfo) <- as.character(rowInfo$gene)
  rowInfo$gene_name <- gtf_ensembl_gene[as.character(rowInfo$gene), "gene_name"]
  rowInfo$gene_name[which(is.na(rowInfo$gene_name))] <- as.character(rowInfo$gene)[which(is.na(rowInfo$gene_name))]
  head(rowInfo)
  ####
  colData <- colData(Object)
  colInfo <- data.frame(colData[, c("Sex", "Age_Stage")])
  head(colInfo)
  str(colInfo)
  
  ann_colors <- list(Age_Stage = anno_colors$Age_Stage,
                     Sex =anno_colors$Sex)
  
  dfheat <- vsd_data[rownames(rowInfo), rownames(colInfo)]
  rownames(dfheat) <- rowInfo$gene_name
  scale_df <- scale_mat(dfheat, "row")
  scale_df[scale_df>4] <- 3
  filename = paste("revise0530/gene_count_res/Sex_DEG_profile_in_", region, "_after_permutation.pdf", sep="")
  pheatmap(
    dfheat,
    col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(21),
    breaks = c(-3,seq(-2,2,length=20),3),
    # legend_breaks = c(-2,-1,1,2),
    filename = filename,
    main = region,
    scale="row",
    # main = main,
    width = 10,
    height = 10,
    # clustering_method = "single",
    border_color=NA,
    # fontsize = 0.7,
    annotation_col = colInfo,
    annotation_colors = ann_colors,
    gaps_row = cumsum(table(rowInfo$Sex)[c("male", "female")])[complete.cases(cumsum(table(rowInfo$Sex)[c("male", "female")]))],
    annotation_legend = T,
    show_rownames = T,
    show_colnames = F,
    cluster_rows = F 
  )
}





#######################################################
#                   gam Permutation                   #
#######################################################
rm(list=ls())
set.seed(100)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")



gam_function <- function(r){
  library(BiocParallel)
  setwd("/mnt/data2/Rhesus_brain")
  load("revise0530/gene_count_res/Rhesus_colData.RData")
  load("revise0530/gene_count_res/RhesusGeneCount.RData")
  # clus <- makeCluster(50)
  ####
  colData <- Rhesus_colData[Rhesus_colData$SR_merge == r, ]
  df <- RhesusGeneCount[, rownames(colData)]
  df <- df[apply(df, 1, function(x){sum(x>1)>3}), ]
  ##########
  ### for ID
  ##########
  testData <- Rhesus_colData[, c("ID", "Sex", "Age_Stage")]
  testData <- testData[!duplicated(testData$ID),]
  test_tmp <- testData
  i <- 1
  random_list <- list()
  while(i <= 1000){
    test_tmp$ID_test <- sample(test_tmp$ID)
    rownames(test_tmp) <- test_tmp$ID_test
    test_tmp <- test_tmp[order(test_tmp$Sex,test_tmp$Age_Stage,test_tmp$ID_test),]
    if(length(random_list)==0){
      random_list[[i]] <- test_tmp
      # test_tmp$IDnum <- as.numeric(sapply(test_tmp$ID_test, function(x){strsplit(x, "_")[[1]][2]}))
      # uniquePermID = as.numeric(paste(test_tmp$IDnum,collapse = ""))
      # print(uniquePermID)
      i <- i + 1
    }else{
      a <- lapply(random_list, function(x){all(x$ID_test == test_tmp$ID_test)})
      if(all(a==F)){
        random_list[[i]] <- test_tmp
        i <- i + 1
      } else{i <- i}
    }
  }
  # random_list <- lapply(random_list, function(x) list(x, df))
  #################
  random_func <- function(random_data){
    library(mgcv)
    setwd("/mnt/data2/Rhesus_brain")
    load("revise0530/gene_count_res/Rhesus_colData.RData")
    test_tmp <- random_data
    # test_tmp <- random_data[[1]]
    # df <- random_data[[2]]
    colData <- Rhesus_colData[Rhesus_colData$SR_merge == r, ]
    colData$Age_test <- test_tmp[colData$ID, "Age_Stage"]
    colData$Sex_test <- test_tmp[colData$ID, "Sex"]
    Sex_pvalue <- c()
    Age_pvalue <- c()
    for(j in 1:dim(df)[1]){
      gene_id <- rownames(df)[j]
      GeneExpr = t(df[j,])
      colnames(GeneExpr) <- c("GeneExpr")           
      dd <- data.frame(GeneExpr, colData)
      dd$ID <- factor(dd$ID)
      ####
      coef=tryCatch({
        if(length(unique(colData$SR))>1){
          m.nb <- gam(GeneExpr ~ Sex_test + Age_test + SR + s(ID, bs = "re"),
                      data=dd,
                      family = "nb")
        }else{ m.nb <- gam(GeneExpr ~ Sex_test + Age_test ,
                           data=dd,  
                           family = "nb" )
        }
        coef <- summary(m.nb)$p.pv
      },error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
        cat(j,unique(colData$SR_merge))
        coef <- c(NA, NA, NA)
        return(coef)})
      ####
      Sex_pvalue <- c(Sex_pvalue, coef[2])
      Age_pvalue <- c(Age_pvalue, coef[3])
    }
    pvalue = data.frame(Sex=Sex_pvalue, Age=Age_pvalue)
    rownames(pvalue) <- rownames(df)
    return(pvalue)
  }
  ##
  param <- SnowParam(workers = 50, type = "SOCK")
  random_res <- bplapply(random_list, random_func, BPPARAM = param)
  return(random_res)
}
  
###
# region_list <- list()
# for(r in unique(Rhesus_colData$SR_merge)){region_list[[r]] <- r}
# Permutaion_gam_res <- lapply(region_list, 
#                              function(x) gam_function(x))

Permutaion_gam_res <- list()
for(r in unique(Rhesus_colData$SR_merge)[11:16]){
  Permutaion_gam_res[[r]] <- gam_function(r)
}


Permutaion1000_gam_res_ID <- Permutaion_gam_res


save(Permutaion1000_gam_res_ID,
     file="revise0530/gene_count_res/Permutaion1000_gam_res_ID.RData")




################
###   for SN
################
rm(list=ls())
library(snow)
set.seed(100)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
clus <- makeCluster(100)


random_gamfunction <- function(test_tmp){
  r <- "SN"
  library(mgcv)
  setwd("/mnt/data2/Rhesus_brain")
  load("revise0530/gene_count_res/Rhesus_colData.RData")
  load("revise0530/gene_count_res/RhesusGeneCount.RData")
  colData <- Rhesus_colData[Rhesus_colData$SR_merge == r, ]
  df <- RhesusGeneCount[, rownames(colData)]
  df <- df[apply(df, 1, function(x){sum(x>1)>3}), ]
  colData$Age_test <- test_tmp[colData$ID, "Age_Stage"]
  colData$Sex_test <- test_tmp[colData$ID, "Sex"]
  Sex_pvalue <- c()
  Age_pvalue <- c()
  for(j in 1:dim(df)[1]){
    gene_id <- rownames(df)[j]
    GeneExpr = t(df[j,])
    colnames(GeneExpr) <- c("GeneExpr")           
    dd <- data.frame(GeneExpr, colData)
    dd$ID <- factor(dd$ID)
    ####
    coef=tryCatch({
      if(length(unique(colData$SR))>1){
        m.nb <- gam(GeneExpr ~ Sex_test + Age_test + SR + s(ID, bs = "re"),
                    data=dd,
                    family = "nb")
      }else{ m.nb <- gam(GeneExpr ~ Sex_test + Age_test ,
                         data=dd,  
                         family = "nb" )
      }
      coef <- summary(m.nb)$p.pv
    },error=function(e){
      cat("ERROR :",conditionMessage(e), "\n")
      cat(j,unique(colData$SR_merge))
      coef <- c(NA, NA, NA)
      return(coef)})
    ####
    Sex_pvalue <- c(Sex_pvalue, coef[2])
    Age_pvalue <- c(Age_pvalue, coef[3])
  }
  pvalue = data.frame(Sex=Sex_pvalue, Age=Age_pvalue)
  rownames(pvalue) <- rownames(df)
  return(pvalue)
}



load("revise0530/gene_count_res/Rhesus_colData.RData")
testData <- Rhesus_colData[, c("ID", "Sex", "Age_Stage")]
testData <- testData[!duplicated(testData$ID),]
##########
### random ID
##########
test_tmp <- testData
i <- 1
random_list <- list()
while(i <= 100){
  test_tmp$ID_test <- sample(test_tmp$ID)
  rownames(test_tmp) <- test_tmp$ID_test
  test_tmp <- test_tmp[order(test_tmp$Sex,test_tmp$Age_Stage,test_tmp$ID_test),]
  if(length(random_list)==0){
    random_list[[i]] <- test_tmp
    i <- i + 1
  }else{
    a <- lapply(random_list, function(x){all(x$ID_test == test_tmp$ID_test)})
    if(all(a==F)){
      random_list[[i]] <- test_tmp
      i <- i + 1
    } else{i <- i}
  }
}
##
clusterExport(clus, "random_gamfunction")
region_gam_res <- parLapply(clus, random_list, function(x)
  random_gamfunction(x))

stopCluster(clus)

Permutation100_SN_gam_res <- region_gam_res 

save(Permutation100_SN_gam_res, 
     file="revise0530/gene_count_res/Permutation100_SN_gam_res.RData")


#### histgram
library(ggplot2)
library(reshape2)
load("revise0530/gene_count_res/Permutation20_HIP_gam_res.RData")
load("revise0530/gene_count_res/gam_res.RData")

permutation_pvalue <- melt(Permutation100_SN_gam_res)

### Sex
pvalue <- gam_res[["HIP"]][,"Sex"]
pp <- data.frame(pvalue)
ggplot(pp, aes(x= pvalue)) +
  geom_histogram(position="identity", alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  # facet_wrap(~region, ncol=4, scales = "free") +
  # scale_color_manual(values=c("#1B9E77", "#D95F02")) +
  # scale_fill_manual(values=c("#1B9E77", "#D95F02")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue for Sex") 
ggsave("revise0530/gene_count_res/HIP_gam_res_pvalue.pdf")


pp <- permutation_pvalue[permutation_pvalue$variable=="Sex",]
ggplot(pp, aes(x= value)) +
  geom_histogram(position="identity", alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  # facet_wrap(~L1, ncol=10, scales = "free") +
  # scale_color_manual(values=c("#1B9E77", "#D95F02")) +
  # scale_fill_manual(values=c("#1B9E77", "#D95F02")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue for Sex (permutation)") 
ggsave("revise0530/gene_count_res/SN_gam_pvalue_permutation_Sex.pdf",
       # width=20, height = 20
       )



### Age
pvalue <- HIP_glmmTMB_res[["HIP"]][,"Age"]
pp <- data.frame(pvalue)
ggplot(pp, aes(x= pvalue)) +
  geom_histogram(position="identity", alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  # facet_wrap(~region, ncol=4, scales = "free") +
  # scale_color_manual(values=c("#1B9E77", "#D95F02")) +
  # scale_fill_manual(values=c("#1B9E77", "#D95F02")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue for Age") 
ggsave("revise0530/gene_count_res/HIP_gam_pvalue_Age.pdf")


pp <- permutation_pvalue[permutation_pvalue$variable=="Age",]
ggplot(pp, aes(x= value)) +
  geom_histogram(position="identity", alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  facet_wrap(~L1, ncol=5, scales = "free") +
  # scale_color_manual(values=c("#1B9E77", "#D95F02")) +
  # scale_fill_manual(values=c("#1B9E77", "#D95F02")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue for Age (permutation)") 
ggsave("revise0530/gene_count_res/HIP_gam_pvalue_permutation_Age.pdf",
       width=10, height = 8
       )



############ SN
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
load("revise0530/gene_count_res/RhesusGeneCount.RData")
colData <- Rhesus_colData[Rhesus_colData$SR_merge == "SN", ]
df <- RhesusGeneCount[, rownames(colData)]
df <- df[apply(df, 1, function(x){sum(x>1)>3}), ]

colData <- data.frame(sample_id=rownames(colData), colData)
write.table(colData, "revise0530/gene_count_res/SN_colData.txt",
            row.names = F, col.names = T, sep="\t", quote = F)
df <- data.frame(gene_id=rownames(df), df)
write.table(df, "revise0530/gene_count_res/SN_expr.txt",
            row.names = F, col.names = T, sep="\t", quote = F)




#######################################################
#          gam Permutation (固定Sex或Age)             #
#######################################################
rm(list=ls())
set.seed(100)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")



gam_function <- function(r){
  library(BiocParallel)
  setwd("/mnt/data2/Rhesus_brain")
  load("revise0530/gene_count_res/Rhesus_colData.RData")
  load("revise0530/gene_count_res/RhesusGeneCount.RData")
  ####
  colData <- Rhesus_colData[Rhesus_colData$SR_merge == r, ]
  df <- RhesusGeneCount[, rownames(colData)]
  df <- df[apply(df, 1, function(x){sum(x>1)>3}), ]
  ##########
  ### for Sex
  ##########
  library(gtools)
  testData <- Rhesus_colData[, c("ID", "Sex", "Age_Stage")]
  testData <- testData[!duplicated(testData$ID),]
  testData <- testData[order(testData$Age_Stage, testData$Sex, testData$ID),]
  comb1 <- combinations(4,2,testData$ID[1:4])
  comb2 <- combinations(4,2,testData$ID[5:8])
  Sex_ID_test <- data.frame(testData$ID)
  for(i in 1:6){
    for(j in 1:6){
      ID_order <- c(comb1[i,], sort(setdiff(testData$ID[1:4], comb1[i,])),
                    comb2[j,], sort(setdiff(testData$ID[5:8], comb2[j,])))
      ID_order_tmp <- c(ID_order[3:4], ID_order[1:2], ID_order[7:8], ID_order[5:6])
      TF <- apply(Sex_ID_test,2,function(x){all(x==ID_order_tmp)})
      if(any(TF==T)) next
      Sex_ID_test <- data.frame(Sex_ID_test, ID_order)
    }
  }
  Sex_ID_test <- Sex_ID_test[,-1]
  Sex_random_list <- apply(Sex_ID_test, 2, function(x){test_tmp <- testData
                                    test_tmp$ID_test <- x;
                                    rownames(test_tmp) <- test_tmp$ID_test;
                                    test_tmp})
  ##########
  ### for Age
  ##########
  library(gtools)
  testData <- Rhesus_colData[, c("ID", "Sex", "Age_Stage")]
  testData <- testData[!duplicated(testData$ID),]
  testData <- testData[order(testData$Sex, testData$Age_Stage, testData$ID),]
  comb1 <- combinations(4,2,testData$ID[1:4])
  comb2 <- combinations(4,2,testData$ID[5:8])
  Age_ID_test <- data.frame(testData$ID)
  for(i in 1:6){
    for(j in 1:6){
      ID_order <- c(comb1[i,], sort(setdiff(testData$ID[1:4], comb1[i,])),
                    comb2[j,], sort(setdiff(testData$ID[5:8], comb2[j,])))
      ID_order_tmp <- c(ID_order[3:4], ID_order[1:2], ID_order[7:8], ID_order[5:6])
      TF <- apply(Age_ID_test,2,function(x){all(x==ID_order_tmp)})
      if(any(TF==T)) next
      Age_ID_test <- data.frame(Age_ID_test, ID_order)
    }
  }
  Age_ID_test <- Age_ID_test[,-1]
  Age_random_list <- apply(Age_ID_test, 2, function(x){test_tmp <- testData
                                           test_tmp$ID_test <- x;
                                           rownames(test_tmp) <- test_tmp$ID_test;
                                           test_tmp})
  #####################
  random_func <- function(random_data){
    library(mgcv)
    setwd("/mnt/data2/Rhesus_brain")
    load("revise0530/gene_count_res/Rhesus_colData.RData")
    test_tmp <- random_data
    # test_tmp <- random_data[[1]]
    # df <- random_data[[2]]
    colData <- Rhesus_colData[Rhesus_colData$SR_merge == r, ]
    colData$Age_test <- test_tmp[colData$ID, "Age_Stage"]
    colData$Sex_test <- test_tmp[colData$ID, "Sex"]
    Sex_pvalue <- c()
    Age_pvalue <- c()
    Sex_zvalue <- c()
    Age_zvalue <- c()
    for(j in 1:dim(df)[1]){
      gene_id <- rownames(df)[j]
      GeneExpr = t(df[j,])
      colnames(GeneExpr) <- c("GeneExpr")           
      dd <- data.frame(GeneExpr, colData)
      dd$ID <- factor(dd$ID)
      ####
      coef=tryCatch({
        if(length(unique(colData$SR))>1){
          m.nb <- gam(GeneExpr ~ Sex_test + Age_test + SR + s(ID, bs = "re"),
                      data=dd,
                      family = "nb")
        }else{ m.nb <- gam(GeneExpr ~ Sex_test + Age_test ,
                           data=dd,  
                           family = "nb" )
        }
        coef <- summary(m.nb)$p.pv
        zvalue <- summary(m.nb)$p.t
        coef <- c(coef, zvalue)
      },error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
        cat(j,unique(colData$SR_merge))
        coef <- c(NA, NA, NA, NA, NA, NA)
        return(coef)})
      ####
      Sex_pvalue <- c(Sex_pvalue, coef[2])
      Age_pvalue <- c(Age_pvalue, coef[3])
      Sex_zvalue <- c(Sex_zvalue, coef[5])
      Age_zvalue <- c(Age_zvalue, coef[6])
    }
    value = data.frame(Sex=Sex_pvalue, Age=Age_pvalue, 
                       Sex_zvalue=Sex_zvalue, Age_zvalue=Age_zvalue)
    rownames(value) <- rownames(df)
    return(value)
  }
  ##
  random_list <- c(Sex_random_list, Age_random_list)
  names(random_list) <- c(paste("Sex", 1:18, sep=""), paste("Age", 1:18, sep=""))
  param <- SnowParam(workers = 36, type = "SOCK")
  random_res <- bplapply(random_list, random_func, BPPARAM = param)
  return(random_res)
}


Permutaion18_gam_res <- list()
for(r in unique(Rhesus_colData$SR_merge)){
  Permutaion18_gam_res[[r]] <- gam_function(r)
}



save(Permutaion18_gam_res,
     file="revise0530/gene_count_res/Permutaion18_gam_res.RData")




#### /mnt/data1/Tools/R-3.6.0/R
# saveRDS(Permutaion_gam_res, file="revise0530/gene_count_res/test.Rds")
# test <- readRDS("revise0530/gene_count_res/test.Rds")
# source("revise0530/gene_count_res/test.RData")
# 
# for(r in names(Permutaion_gam_res)){
#   Permutaion1000_gam_res_ID[[r]] <- Permutaion_gam_res[[r]]
# }

# 
# save(Permutaion1000_gam_res_ID,
#      file="revise0530/gene_count_res/Permutaion1000_gam_res_ID.RData")
# 


######## Permutaion18_gam_res pvalue histogram
library(reshape2)
library(ggplot2)
load("revise0530/gene_count_res/Permutaion18_gam_res.RData")
load("revise0530/gene_count_res/gam_res.RData")


### 
test_list <- list()
for(r in names(Permutaion18_gam_res)){
  for(i in names(Permutaion18_gam_res[[r]])[1:18]){
    test_list[[r]][[i]] <- Permutaion18_gam_res[[r]][[i]][, "Sex"]
  }
  for(i in names(Permutaion18_gam_res[[r]])[19:36]){
    test_list[[r]][[i]] <- Permutaion18_gam_res[[r]][[i]][, "Age"]
  }
}

pp <- melt(test_list)
pp$class <- gsub("[0-9]+", "" , pp$L2)

## Age pvalue distribution
p <- pp[pp$class=="Age", ]
p <- data.frame(p[, c("value", "L1")], state="Permutation")
observe <- melt(lapply(gam_res, function(x){x[, "Age"]}))
observe <- data.frame(observe, state="observe")
p <- p[sample(1:dim(p)[1], dim(observe)[1]), ]
p <- rbind(observe, p)
ggplot(p, aes(x= value, fill= state, color=state)) +
  geom_histogram(position="identity", alpha=0.3) +
  # geom_line(stat = "density") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        # legend.position = c(0.8,0.8),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  facet_wrap(~L1, ncol=4, scales = "free") +
  scale_color_manual(values=c("#F25F5C", "#247BA0")) +
  scale_fill_manual(values=c("#F25F5C", "#247BA0")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue for Age") 


ggsave("revise0530/gene_count_res/Age_permutation18_pvalue_hist.pdf",
       height = 7, width = 11)


## Sex pvalue distribution
p <- pp[pp$class=="Sex", ]
p <- data.frame(p[, c("value", "L1")], state="Permutation")
observe <- melt(lapply(gam_res, function(x){x[, "Sex"]}))
observe <- data.frame(observe, state="observe")
p <- p[sample(1:dim(p)[1], dim(observe)[1]), ]
p <- rbind(observe, p)
ggplot(p, aes(x= value, fill= state, color=state)) +
  geom_histogram(position="identity", alpha=0.3) +
  # geom_line(stat = "density") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        # legend.position = c(0.8,0.8),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  facet_wrap(~L1, ncol=4, scales = "free") +
  scale_color_manual(values=c("#1B9E77", "#D95F02")) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue for Sex") 


ggsave("revise0530/gene_count_res/Sex_permutation18_pvalue_hist.pdf",
       height = 7, width = 11)












### pp
pdf("revise0530/gene_count_res/permutation18_pvalue_hist.pdf",
    width = 10, height = 10)
for(r in names(Permutaion18_gam_res)){
  p <- pp[pp$L1==r,]
  plot=ggplot(p, aes(x= value)) +
    geom_histogram(position="identity", alpha=0.3) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                      colour = "black", size = 11),
          axis.text.y  = element_text(colour = "black"),
          strip.background = element_rect(fill="white", colour = "white", size=1)) +
    facet_wrap(~L2, ncol=6, scales = "free") +
    # scale_color_manual(values=c("#1B9E77", "#D95F02")) +
    # scale_fill_manual(values=c("#1B9E77", "#D95F02")) +
    xlab("pvalue") +
    ggtitle(r)
  plot(plot)
}

dev.off()






#######################################################
#                 glmmTMB Permutation                 #
#######################################################
rm(list=ls())
### permutation Age and Sex
library(snow)
set.seed(100)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Rhesus_colData.RData")
clus <- makeCluster(16)

glmmTMB_function <- function(r){
  library(glmmTMB)
  load("revise0530/gene_count_res/Rhesus_colData.RData")
  load("revise0530/gene_count_res/RhesusGeneCount.RData")
  testData <- Rhesus_colData[, c("ID", "Sex", "Age_Stage")]
  testData <- testData[!duplicated(testData$ID),]
  colData <- Rhesus_colData[Rhesus_colData$SR_merge == r, ]
  df <- RhesusGeneCount[, rownames(colData)]
  df <- df[apply(df, 1, function(x){sum(x>1)>3}), ]
  ###
  # for Age and Sex
  ###
  test_tmp <- testData
  i <- 1
  random_list <- list()
  while(i <= 20){
    test_tmp$Age_test <- sample(test_tmp$Age_Stage)
    test_tmp$Sex_test <- sample(test_tmp$Sex)
    rownames(test_tmp) <- test_tmp$ID
    pair <- paste(test_tmp$Age_test, test_tmp$Sex_test)
    len <- length(unique(pair))
    if(len!=2){
      if(length(random_list)==0){
        random_list[[i]] <- test_tmp
        i <- i + 1
      }else{
        a <- lapply(random_list, function(x){all(x == test_tmp)})
        if(all(a==F)){
          random_list[[i]] <- test_tmp
          i <- i + 1
        } else{i <- i}
      }
    }else{i <- i}
  }
  ##########
  ### for ID
  ##########
  test_tmp <- testData
  i <- 1
  random_list <- list()
  while(i <= 20){
    test_tmp$ID_test <- sample(test_tmp$ID)
    rownames(test_tmp) <- test_tmp$ID_test
    test_tmp <- test_tmp[order(test_tmp$Sex,test_tmp$Age_Stage,test_tmp$ID_test),]
    if(length(random_list)==0){
      random_list[[i]] <- test_tmp
            i <- i + 1
    }else{
      a <- lapply(random_list, function(x){all(x$ID_test == test_tmp$ID_test)})
      if(all(a==F)){
        random_list[[i]] <- test_tmp
                i <- i + 1
      } else{i <- i}
    }
  }
  ##
  ##
  region_glmmTMB_res <- list()
  for(i in 1:20){
    test_tmp <- random_list[[i]]
    colData$Age_test <- test_tmp[colData$ID, "Age_test"]
    colData$Sex_test <- test_tmp[colData$ID, "Sex_test"]
    Sex_pvalue <- c()
    Age_pvalue <- c()
    ##
    for(j in 1:dim(df)[1]){
      GeneExpr = t(df[j,])
      colnames(GeneExpr) <- c("GeneExpr")           
      dd <- data.frame(GeneExpr, colData, stringsAsFactors = F)
      dd$ID <- factor(dd$ID)
      coef=tryCatch({
        if(length(unique(dd$SR))>1){
          m.nb <- glmmTMB(GeneExpr ~ Sex_test + Age_test + SR + (1|ID),  family = nbinom2,  data=dd)
        }else{
          m.nb <- glmmTMB(GeneExpr ~ Sex_test + Age_test ,  family = nbinom2,  data=dd)
        }
        coef <- summary(m.nb)$coefficients$cond[,4]
      },error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
        coef <- c(NA, NA, NA)
        return(coef)})
      Sex_pvalue <- c(Sex_pvalue, coef[2])
      Age_pvalue <- c(Age_pvalue, coef[3])
    }
    ##
    pvalue = data.frame(Sex=Sex_pvalue, Age=Age_pvalue)
    rownames(pvalue) <- rownames(df)
    region_glmmTMB_res[[i]] <- pvalue
  }
  return(region_glmmTMB_res)  
}



clusterExport(clus,"glmmTMB_function")
region_list <- list()
for(r in unique(Rhesus_colData$SR_merge)){region_list[[r]] <- r}

Permutaion_glmmTMB_res <- parLapply(clus, 
                                region_list, 
                                function(x) glmmTMB_function(x))


save(Permutaion_glmmTMB_res,
     file="revise0530/gene_count_res/Permutaion20_glmmTMB_res.RData")




################################################################
#     pvalue adjusted using permutations pvalue (per gene)     #
################################################################
rm(list=ls())
library(BiocParallel)
library(reshape2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/gam_res.RData")


adjust_func <- function(r){
  library(reshape2)
  load("revise0530/gene_count_res/gam_res.RData")
  load("revise0530/gene_count_res/Permutaion1000_gam_res_ID.RData")
  region_gam_res <- gam_res[[r]]
  Permutation1000_gam_res <- Permutaion1000_gam_res_ID[[r]]
  pvalue <- list()
  for(i in 1:dim(region_gam_res)[1]){
    ## Sex
    Sex_per <- melt(lapply(Permutation1000_gam_res, function(x){x[i,"Sex"]}))
    per <- sum(Sex_per$value < region_gam_res[i, "Sex"])/1000
    pvalue[["Sex"]] <- c(pvalue[["Sex"]], per)
    ## Age
    Sex_per <- melt(lapply(Permutation1000_gam_res, function(x){x[i,"Age"]}))
    per <- sum(Sex_per$value < region_gam_res[i, "Age"])/1000
    pvalue[["Age"]] <- c(pvalue[["Age"]], per)
  }
  ###
  qvalue <- lapply(pvalue, function(x){padj <- p.adjust(x, method = "fdr", n=length(x))})
  qvalue <- data.frame(Sex=qvalue[["Sex"]], Age=qvalue[["Age"]])
  ###
  pvalue <- data.frame(Sex=pvalue[["Sex"]], Age=pvalue[["Age"]])
  rownames(pvalue) <- rownames(region_gam_res)
  rownames(qvalue) <- rownames(region_gam_res)
  return(list(pvalue, qvalue))
}


param <- SnowParam(workers = 16, type = "SOCK")
permutation_pqvalue  <- bplapply(names(gam_res), adjust_func, BPPARAM = param)
names(permutation_pqvalue) <- names(gam_res)
permutation_pvalue <- lapply(permutation_pqvalue, function(x) x[[1]])

permutation1000_pqvalue <- permutation_pqvalue
save(permutation1000_pqvalue, 
     file="revise0530/gene_count_res/permutation1000_pqvalue.RData")



########### permutation_pvalue distirbution
rm(list=ls())
library(reshape2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/permutation_pqvalue.RData")
load("revise0530/gene_count_res/gam_res.RData")

permutation_pvalue <- lapply(permutation_pqvalue, function(x) x[[1]])

permute <- melt(permutation_pvalue)
observe <- melt(gam_res)
pp <- rbind(data.frame(observe, type="observed"),
            data.frame(permute, type="adjusted")
            )

ggplot(pp, aes(x= value, fill= type, color=type)) +
  geom_histogram(position="identity", alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  facet_wrap(variable~L1, ncol=8, scales = "free") +
  scale_color_manual(values=c("#F25F5C", "#247BA0")) +
  scale_fill_manual(values=c("#F25F5C", "#247BA0")) +
  xlab("pvalue") +
  ggtitle("the histogram of adjusted pvalue ") 


ggsave("revise0530/gene_count_res/Permutation_adjusted_pvalue_hist.pdf",
       height = 7, width = 20)







##########
# //mnt/data1/Tools/R-3.6.0/bin/R
library("DESeq2")
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Age_res.RData")
load("revise0530/gene_count_res/Sex_res.RData")
load("revise0530/gene_count_res/gam_res.RData")
load("revise0530/gene_count_res/permutation1000_pqvalue.RData")


permutation_qvalue <- lapply(permutation1000_pqvalue, function(x) x[[2]])


Permutation_genes <- list()
for(r in names(permutation_qvalue)){
  qvalue <- permutation_qvalue[[r]]
  for(condition in names(qvalue)){
    if(condition == "Sex"){
      res <- results(Sex_res[[r]],
                     contrast = c("Sex", c("male", "female")))
      res <- data.frame(qvalue=qvalue[rownames(res), "Sex"], 
                        res[, "log2FoldChange",drop=F])
    }
    if(condition == "Age"){
      res <- results(Age_res[[r]],
                     contrast = c("Age_Stage", c("Mid", "Young")))
      res <- data.frame(qvalue=qvalue[, "Age"], 
                        res[rownames(gam_res[[r]]), "log2FoldChange",drop=F])
    }
    #####
    res <- res[complete.cases(res),]
    res <- res[order(res$qvalue, decreasing = T),]
    hyper <- res[res$qvalue<0.1 & res$log2FoldChange>1, ]
    hypo <- res[res$qvalue<0.1 & res$log2FoldChange < -1, ]
    if(condition == "Sex"){
      Permutation_genes[[r]][[condition]] <- list(male=hyper, female=hypo)
    }
    if(condition == "Age"){
      Permutation_genes[[r]][[condition]] <- list(Mid=hyper, Young=hypo)
    }
  }
}

save(Permutation_genes, 
     file="revise0530/gene_count_res/Permutation_genes.RData")




#############################################
# plot the count of differential genes for Sex
library(reshape2)
library(ggplot2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Permutation_genes.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")


Sex_genes_InEachSR <- lapply(Permutation_genes, function(x){y <- x[["Sex"]]; list(male=rownames(y[["male"]]), female=rownames(y[["female"]]))})
Sex_genes_count <- melt(Sex_genes_InEachSR)
colnames(Sex_genes_count) <- c("gene", "Sex", "SR")
Sex_genes_count1 <- table(Sex_genes_count$SR, Sex_genes_count$Sex)
Sex_genes_count1 <- melt(Sex_genes_count1)
colnames(Sex_genes_count1) <- c("SR", "Sex", "counts")
head(Sex_genes_count1)

Sex_genes_count2 <- table(Sex_genes_count$SR)
Sex_genes_count2 <- Sex_genes_count2[order(Sex_genes_count2, decreasing=T)]
head(Sex_genes_count2)
Sex_genes_count1$SR <- factor(Sex_genes_count1$SR, levels=(names(Sex_genes_count2)))


### plot
pp <- Sex_genes_count1
ggplot(pp,aes(x=SR, y=counts, fill=Sex)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  scale_fill_manual(values = anno_colors$Sex) +
  ggtitle("the counts of DEG for Sex (after permutation)")

ggsave("revise0530/gene_count_res/Sex_DEG_counts_bar_InSR_permutation_adjust.pdf",
       width=6, height = 5)


Age_genes_InEachSR <- lapply(Permutation_genes, function(x){y <- x[["Age"]]; list(Mid=rownames(y[["Mid"]]), Young=rownames(y[["Young"]]))})
Age_genes_count <- melt(Age_genes_InEachSR)
colnames(Age_genes_count) <- c("gene", "Age", "SR")
Age_genes_count1 <- table(Age_genes_count$SR, Age_genes_count$Age)
Age_genes_count1 <- melt(Age_genes_count1)
colnames(Age_genes_count1) <- c("SR", "Age", "counts")
head(Age_genes_count1)

Age_genes_count2 <- table(Age_genes_count$SR)
Age_genes_count2 <- Age_genes_count2[order(Age_genes_count2, decreasing=T)]
head(Age_genes_count2)
Age_genes_count1$SR <- factor(Age_genes_count1$SR, levels=(names(Age_genes_count2)))


### plot
pp <- Age_genes_count1
ggplot(pp,aes(x=SR, y=counts, fill=Age)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = c(Mid=brewer.pal(9,"Blues")[8], Young=brewer.pal(9,"Reds")[8])) +
  ggtitle("the counts of DEG for age (after permutation)")


ggsave("revise0530/gene_count_res/Age_DEG_counts_bar_InSR_permutation_adjust.pdf",
       width=6, height = 5)





################################################################
#            pvalue adjusted using permutations pvalue         #
################################################################
rm(list=ls())
library(fdrci)
library(DESeq2)
library(reshape2)
library(ggplot2)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/gam_res.RData")
load("revise0530/gene_count_res/Permutaion100_gam_res_ID.RData")


Permutaion_gam_res <- Permutaion_gam_res_ID

fdr_od_res <- list()
fdrTbl_res <- list()
pvalue_res <- list()
for(r in names(Permutaion_gam_res)){
  perm_list <- Permutaion_gam_res[[r]]
  for(condition in c("Age", "Sex")){
    observe <- gam_res[[r]]
    observe_pvalue <- observe[,condition, drop=F]
    sub_perm_list <- lapply(perm_list, function(x){a <- x[rownames(observe_pvalue),condition, drop=F];
                                                   colnames(a) <- c("pvalue");
                                                   a})
    pp <- data.frame(observe_pvalue, sub_perm_list)
    pp <- pp[complete.cases(pp),]
    pp <- melt(pp)
    pp$variable <- as.character(pp$variable)
    pp$variable[1:dim(observe_pvalue)[1]] <- "observe_pvalue"
    pp$variable[pp$variable != "observe_pvalue"] <- "random_pvalue"
    pvalue_res[[r]][[condition]] <- data.frame(pp, region=r)
    ######
    test <- fdrTbl(observe_pvalue, sub_perm_list, "pvalue", length(observe_pvalue), 1, 50)
    test2 <- fdr_od(observe_pvalue, sub_perm_list, "pvalue", length(observe_pvalue), 
                    thres=0.05, cl=.95)
    fdr_od_res[[r]][[condition]] <- test2
    fdrTbl_res[[r]][[condition]] <- test
  }
}

save(pvalue_res, fdr_od_res, fdrTbl_res,
     file="revise0530/gene_count_res/fdrci_gam_res.RData")



#####  Age pvalue distribution
sample_pvalue_res <- list()
for(r in names(pvalue_res)){
  x <- pvalue_res[[r]][["Age"]]
  a <- x[c(1:sum(x$variable=="observe_pvalue"), sample(which(x$variable=="random_pvalue"), sum(x$variable=="observe_pvalue"))), ]
  sample_pvalue_res[[r]] <- a
  }

pp <- Reduce(rbind, sample_pvalue_res)

ggplot(pp, aes(x= value, fill= variable, color=variable)) +
  geom_histogram(position="identity", alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  facet_wrap(~region, ncol=4, scales = "free") +
  scale_color_manual(values=c("#F25F5C", "#247BA0")) +
  scale_fill_manual(values=c("#F25F5C", "#247BA0")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue for Age") 


ggsave("revise0530/gene_count_res/Age_permutation_pvalue_hist.pdf",
       height = 7, width = 11)


#####  Sex pvalue distribution
sample_pvalue_res <- list()
for(r in names(pvalue_res)){
  x <- pvalue_res[[r]][["Sex"]]
  a <- x[c(1:sum(x$variable=="observe_pvalue"), sample(which(x$variable=="random_pvalue"), sum(x$variable=="observe_pvalue"))), ]
  sample_pvalue_res[[r]] <- a
}
pp <- Reduce(rbind, sample_pvalue_res)

ggplot(pp, aes(x= value, fill= variable, color=variable)) +
  geom_histogram(position="identity", alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  facet_wrap(~region, ncol=4, scales = "free") +
  scale_color_manual(values=c("#1B9E77", "#D95F02")) +
  scale_fill_manual(values=c("#1B9E77", "#D95F02")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue for Sex") 


ggsave("revise0530/gene_count_res/Sex_permutation_pvalue_hist.pdf",
       height = 7, width = 11)



##########
load("revise0530/gene_count_res/Age_res.RData")
load("revise0530/gene_count_res/Sex_res.RData")

Permutation_genes <- list()
for(r in names(fdrTbl_res)){
  fdrres <- fdrTbl_res[[r]]
  for(condition in names(fdrres)){
    sub_fdrres <- fdrres[[condition]]
    if(condition == "Sex"){
      res <- results(Sex_res[[r]],
                      contrast = c("Sex", c("male", "female")))
      res <- data.frame(pvalue=gam_res[[r]][rownames(res), "Sex"], 
                        res[, "log2FoldChange",drop=F])
    }
    if(condition == "Age"){
      res <- results(Age_res[[r]],
                     contrast = c("Age_Stage", c("Mid", "Young")))
      res <- data.frame(pvalue=gam_res[[r]][, "Age"], 
                        res[rownames(gam_res[[r]]), "log2FoldChange",drop=F])
    }
    #####
    res <- res[complete.cases(res),]
    res$logP <- -log10(res$pvalue)
    res <- res[order(res$logP, decreasing = T),]
    #####
    sub_fdrres <- sub_fdrres[complete.cases(sub_fdrres),]
    sub_fdrres <- sub_fdrres[sub_fdrres$fdr < 0.1,]
    if(dim(sub_fdrres)[1]==0) next
    else{
      logP <- sub_fdrres[1, "threshold"]
      ##
      hyper <- res[res$logP>logP & res$log2FoldChange>1, ]
      hypo <- res[res$logP>logP & res$log2FoldChange < -1, ]
      if(condition == "Sex"){
        Permutation_genes[[r]][[condition]] <- list(male=hyper, female=hypo)
      }
      if(condition == "Age"){
        Permutation_genes[[r]][[condition]] <- list(Mid=hyper, Young=hypo)
      }
    }
  }
 }

save(Permutation_genes, 
     file="revise0530/gene_count_res/Permutation_genes.RData")




#############################################
########### plot the count of differential genes for Sex
library(reshape2)
library(ggplot2)
library(RColorBrewer)
setwd("/mnt/data2/Rhesus_brain")
load("revise0530/gene_count_res/Permutation_genes.RData")
load("revise0530/gene_count_res/Rhesus_anno_colors.RData")


Sex_genes_InEachSR <- lapply(Permutation_genes, function(x){y <- x[["Sex"]]; list(male=rownames(y[["male"]]), female=rownames(y[["female"]]))})
Sex_genes_count <- melt(Sex_genes_InEachSR)
colnames(Sex_genes_count) <- c("gene", "Sex", "SR")
Sex_genes_count1 <- table(Sex_genes_count$SR, Sex_genes_count$Sex)
Sex_genes_count1 <- melt(Sex_genes_count1)
colnames(Sex_genes_count1) <- c("SR", "Sex", "counts")
head(Sex_genes_count1)

Sex_genes_count2 <- table(Sex_genes_count$SR)
Sex_genes_count2 <- Sex_genes_count2[order(Sex_genes_count2, decreasing=T)]
head(Sex_genes_count2)
Sex_genes_count1$SR <- factor(Sex_genes_count1$SR, levels=(names(Sex_genes_count2)))


### plot
pp <- Sex_genes_count1
ggplot(pp,aes(x=SR, y=counts, fill=Sex)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  scale_fill_manual(values = anno_colors$Sex) +
  ggtitle("the counts of DEG for Sex (after permutation)")

ggsave("revise0530/gene_count_res/Sex_DEG_counts_bar_InSR_permutation_adjust.pdf",
       width=6, height = 5)


Age_genes_InEachSR <- lapply(Permutation_genes, function(x){y <- x[["Age"]]; list(Mid=rownames(y[["Mid"]]), Young=rownames(y[["Young"]]))})
Age_genes_count <- melt(Age_genes_InEachSR)
colnames(Age_genes_count) <- c("gene", "Age", "SR")
Age_genes_count1 <- table(Age_genes_count$SR, Age_genes_count$Age)
Age_genes_count1 <- melt(Age_genes_count1)
colnames(Age_genes_count1) <- c("SR", "Age", "counts")
head(Age_genes_count1)

Age_genes_count2 <- table(Age_genes_count$SR)
Age_genes_count2 <- Age_genes_count2[order(Age_genes_count2, decreasing=T)]
head(Age_genes_count2)
Age_genes_count1$SR <- factor(Age_genes_count1$SR, levels=(names(Age_genes_count2)))


### plot
pp <- Age_genes_count1
ggplot(pp,aes(x=SR, y=counts, fill=Age)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = c(Mid=brewer.pal(9,"Blues")[8], Young=brewer.pal(9,"Reds")[8])) +
  ggtitle("the counts of DEG for age (after permutation)")


ggsave("revise0530/gene_count_res/Age_DEG_counts_bar_InSR_permutation_adjust.pdf",
       width=6, height = 5)







#######################################################################
expr <- c(132,138,153,235,255,265,124,122,103,231,239,217,222,231,216,323,308,283,190,203,198,281,264,285)
sex1 <- rep("male",3); sex2 <- rep("female", 3)
sexind <- c(sex1, sex1, sex1, sex1, sex2, sex2, sex2, sex2)
age1 <- rep("young",3); age2 <- rep("midage",3)
ageind <- c(age1, age2, age1, age2, age1, age2, age1, age2)
ind = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8)


dd <- data.frame(expr, sexind, ageind, ind)


lrtest(gam(expr ~ sexind + s(ind, bs = "re"), data=dd, family = "nb"), 
       gam(expr ~ sexind + ageind + s(ind, bs = "re"), data=dd, family = "nb"))
lrtest(gam(expr ~ ageind + s(ind, bs = "re"), data=dd, family = "nb"), 
       gam(expr ~ sexind + ageind + s(ind, bs = "re"), data=dd, family = "nb"))





#################################################################
#                           simulation
#################################################################
rm(list=ls())
library(mgcv)
setwd("/mnt/data2/Rhesus_brain")
annotation <- read.csv("revise0530/gene_count_res/annotation_ind60rep30.csv")
annotation$ind <- paste("R", annotation$ind, sep="")
colnames(annotation) <- c("ID", "Sex", "Age_Stage")

gene50 <- read.csv("revise0530/gene_count_res/ind60rep30(1).csv")
##########
### for ID
##########
testData <- annotation[!duplicated(annotation$ID),]
test_tmp <- testData
i <- 1
random_list <- list()
while(i <= 20){
  test_tmp$ID_test <- sample(test_tmp$ID)
  rownames(test_tmp) <- test_tmp$ID_test
  test_tmp <- test_tmp[order(test_tmp$Sex,test_tmp$Age_Stage,test_tmp$ID_test),]
  if(length(random_list)==0){
    random_list[[i]] <- test_tmp
     i <- i + 1
  }else{
    a <- lapply(random_list, function(x){all(x$ID_test == test_tmp$ID_test)})
    if(all(a==F)){
      random_list[[i]] <- test_tmp
           i <- i + 1
    } else{i <- i}
  }
}
##
test_gam_res <- list()
for(i in 1:20){
  test_tmp <- random_list[[i]]
  annotation$Age_test <- test_tmp[annotation$ID, "Age_Stage"]
  annotation$Sex_test <- test_tmp[annotation$ID, "Sex"]
  Sex_pvalue <- c()
  Age_pvalue <- c()
  for(j in 1:dim(gene50)[1]){
    GeneExpr = t(gene50[j,])
    colnames(GeneExpr) <- c("GeneExpr")           
    dd <- data.frame(GeneExpr, annotation)
    dd$ID <- factor(dd$ID)
    ####
    coef=tryCatch({
      m.nb <- gam(GeneExpr ~ Sex_test + Age_test + s(ID, bs = "re"),
                  data=dd,
                  family = "nb")
      coef <- summary(m.nb)$p.pv
    },error=function(e){
      cat("ERROR :",conditionMessage(e), "\n")
      cat(j,unique(colData$SR_merge))
      coef <- c("NA", "NA", "NA")
      return(coef)})
    ####
    Sex_pvalue <- c(Sex_pvalue, coef[2])
    Age_pvalue <- c(Age_pvalue, coef[3])
  }
  pvalue = data.frame(Sex=Sex_pvalue, Age=Age_pvalue)
  rownames(pvalue) <- rownames(df)
  test_gam_res[[i]] <- pvalue
}



####  observe
Sex_pvalue <- c()
Age_pvalue <- c()
for(j in 1:dim(gene50)[1]){
  GeneExpr = t(gene50[j,])
  colnames(GeneExpr) <- c("GeneExpr")           
  dd <- data.frame(GeneExpr, annotation)
  dd$ID <- factor(dd$ID)
  ####
  coef=tryCatch({
    m.nb <- gam(GeneExpr ~ Sex + Age_Stage + s(ID, bs = "re"),
                data=dd,
                family = "nb")
    coef <- summary(m.nb)$p.pv
  },error=function(e){
    cat("ERROR :",conditionMessage(e), "\n")
    cat(j,unique(colData$SR_merge))
    coef <- c(NA, NA, NA)
    return(coef)})
  ####
  Sex_pvalue <- c(Sex_pvalue, coef[2])
  Age_pvalue <- c(Age_pvalue, coef[3])
}
pvalue = data.frame(Sex=Sex_pvalue, Age=Age_pvalue)
rownames(pvalue) <- rownames(df)


pdf("revise0530/gene_count_res/test_Sex_pvalue_hist.pdf")
hist(Sex_pvalue)
dev.off()




#################################################
#### histgram
library(ggplot2)
library(reshape2)

permutation_pvalue <- melt(test_gam_res)

### Sex
pp <- permutation_pvalue[permutation_pvalue$variable=="Sex",]
ggplot(pp, aes(x= value)) +
  geom_histogram(position="identity", alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  # facet_wrap(~region, ncol=4, scales = "free") +
  # scale_color_manual(values=c("#1B9E77", "#D95F02")) +
  # scale_fill_manual(values=c("#1B9E77", "#D95F02")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue for Sex (permutation)") 
ggsave("revise0530/gene_count_res/gene50_gam_pvalue_permutation_Sex.pdf")



### Age
pp <- permutation_pvalue[permutation_pvalue$variable=="Age",]
ggplot(pp, aes(x= value)) +
  geom_histogram(position="identity", alpha=0.3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x  = element_text(angle = 0, hjust=0.5, vjust=0.5, 
                                    colour = "black", size = 11),
        axis.text.y  = element_text(colour = "black"),
        strip.background = element_rect(fill="white", colour = "white", size=1)) +
  # facet_wrap(~region, ncol=4, scales = "free") +
  # scale_color_manual(values=c("#1B9E77", "#D95F02")) +
  # scale_fill_manual(values=c("#1B9E77", "#D95F02")) +
  xlab("pvalue") +
  ggtitle("the histogram of pvalue for Age (permutation)") 
ggsave("revise0530/gene_count_res/gene50_gam_pvalue_permutation_Age.pdf")






#################
rm(list=ls())
library(mgcv)
library(lmtest)
setwd("/mnt/data2/Rhesus_brain")

annotation <- read.csv("revise0530/gene_count_res/annotation_ind60rep30.csv")
annotation$ind <- paste("R", annotation$ind, sep="")
colnames(annotation) <- c("ID", "Sex", "Age_Stage")

dat <- read.csv("revise0530/gene_count_res/gene1.csv")
dat$ind <- factor(dat$ind)

m.full <- gam(y ~ sex + age + s(ind, bs = "re"), family = "nb", data = dat)
m.reduce <- gam(y ~ age + s(ind, bs = "re"), family = "nb", data = dat)
summary.gam(m.full)
summary(m.full)
anova(m.full, m.reduce, test = "LRT")
lrtest(m.reduce, m.full)


#########
library(MapGAM)


m.full <- modgam(y ~ sex + age + s(ind, bs = "re"), 
                 family = "nb", 
                 data = dat)




####################
### for Age and Sex
###################
# test_tmp <- testData
# i <- 1
# random_list <- list()
# while(i <= 20){
#   test_tmp$Age_test <- sample(test_tmp$Age_Stage)
#   test_tmp$Sex_test <- sample(test_tmp$Sex)
#   rownames(test_tmp) <- test_tmp$ID
#   pair <- paste(test_tmp$Age_test, test_tmp$Sex_test)
#   len <- length(unique(pair))
#   if(len!=2){
#     if(length(random_list)==0){
#       random_list[[i]] <- test_tmp
#       i <- i + 1
#     }else{
#       a <- lapply(random_list, function(x){all(x == test_tmp)})
#       if(all(a==F)){
#         random_list[[i]] <- test_tmp
#         i <- i + 1
#         } else{i <- i}
#       }
#   }else{i <- i}
# }


