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







#######################################################
#          gam Permutation (fix Sex or Age)           #
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

colnames(p) <- c("Pvalue", "Region", "state")
write.table(p, "SourceData/Fig.S3a.txt",
            col.names = T, row.names = F, quote=F, sep="\t")


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
  geom_boxplot(alpha=0.5, outlier.colour = NA)+
  geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.1), color="black", size=0.2) +
  # geom_dotplot(aes(fill=Age_Stage), binaxis='y', stackdir='center', stackratio=1.5, dotsize=1.2) +
  # geom_jitter(, position = position_jitter(width = .5)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5, color = "black"),
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

ggsave("revise0530/gene_count_res/Age_selected_genes_expression_boxplot.pdf", 
       width = 7, height = 2)


colnames(p)[2] <- "Region"
write.table(p, "SourceData/Fig.S3c.txt",
            col.names = T, row.names = F, quote=F, sep="\t")




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







