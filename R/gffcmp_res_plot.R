library(ggplot2)
library(reshape2)
library(wesanderson) # names(wes_palettes)
options(stringAsFactors=FALSE)


setwd("/mnt/xdlab1/home/looking/brain_project")


stat1 <- read.table("gtf_for_assembled/gffcmp.stats_for_plot_1.txt", sep="\t", header = T)
stat2 <- read.table("gtf_for_assembled/gffcmp.stats_for_plot_2.txt", sep="\t", header = T)

pp <- melt(stat1)
colnames(pp) <- c("Type", "Accuracy", "Percent")
ggplot(pp, aes(x=Type, y=Percent, fill=Type)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = wes_palette("Royal2")) +
  facet_wrap(~Accuracy,ncol=2) +
  xlab("") +
  ylab("Percent(%)") +
  guides(fill=F) +
  geom_text(aes(y=Percent, label=Percent), vjust=-0.2, size=3) 

ggsave("gtf_for_assembled/gffcmp_sensitivity_precision_bar.pdf", 
       width=8, height = 5)    
## 


##
pp <- stat2[,c(1,4)]
ggplot(pp, aes(x=Type, y=Percent, fill=Type)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = colorRampPalette(wes_palette("Royal2"))(6)) +
  # facet_wrap(~Accuracy,ncol=2) +
  xlab("") +
  ylab("Percent(%)") +
  guides(fill=F) +
  geom_text(aes(y=Percent, label=Percent), vjust=-0.2, size=3) 

ggsave("gtf_for_assembled/gffcmp_Missed_Novel_bar.pdf", 
       width=5.5, height = 5)  



########################################################################
#                           TACO minExpr                               #
########################################################################

### sensitivity_precision
stat1 <- read.table("gtf_for_assembled/TACO_minExpr_gffcmp.stats_sensitivity_precision.txt", 
                    sep="\t", header = T)

pp <- stat1[which(stat1$Type=="Exon_level" | stat1$Type=="Intron_level"),]
F1_exon <-  pp[pp$Sensitivity_Precision == "F1_score" &pp$Type== "Exon_level",]
F1_exon <- F1_exon[order(F1_exon$Percent, decreasing = T), ]
F1_Intron <-  pp[pp$Sensitivity_Precision == "F1_score" &pp$Type== "Intron_level",]
F1_Intron <- F1_Intron[order(F1_Intron$Percent, decreasing = T), ]
head(F1_exon)
head(F1_Intron)

ggplot(pp, aes(x=minExpr, y=Percent, color=Sensitivity_Precision)) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 0.5, vjust=0.5, colour = "black"),
        strip.background = element_rect(fill = "white"),
        text = element_text(face="bold")) +
  scale_color_manual(values = wes_palette("GrandBudapest")) +
  facet_wrap(~Type,ncol=1, scales = "free_x") +
  xlab("TACO parameter minExpr") +
  ylab("Percent(%)") +
  geom_vline(xintercept=4.5, colour="red", linetype="dashed") +
  geom_vline(xintercept=6.5, colour="red", linetype="dashed") +
  scale_x_continuous(limits = c(0.5,9.5)) +
  guides(color=guide_legend(title = NULL))
  

  ggsave("gtf_for_assembled/TACO_minExpr_gffcmp_gffcmp_sensitivity_precision_F1.pdf",
       width=4, height = 5)
## 


### Novel_Missed
minExpr = 5.5
stat2 <- read.table("gtf_for_assembled/TACO_minExpr_gffcmp.stats_Novel_Missed.txt", sep="\t", header = T)
head(stat2)
stat2 <- stat2[stat2$minExpr==minExpr,]

pp <- stat2
ggplot(pp, aes(x=Type, y=Percent, fill=Type)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = colorRampPalette(wes_palette("Royal2"))(6)) +
  # facet_wrap(~Accuracy,ncol=2) +
  xlab("") +
  ylab("Percent(%)") +
  guides(fill=F) +
  geom_text(aes(y=Percent, label=Percent), vjust=-0.2, size=3) 

ggsave("gtf_for_assembled/TACO_minExpr5.5_gffcmp_Missed_Novel_bar.pdf", 
       width=5.5, height = 5)  



## transcripts count
library(ggplot2)
setwd("/mnt/xdlab1/home/looking/brain_project")
stat3 <-  read.table("gtf_for_assembled/TACO_minExpr_gffcmp.stats_transcripts_count.txt", 
                             sep="\t", header = T)

ggplot(stat3, aes(x=minExpr, y=transcripts_count)) +
  geom_line(color=wes_palette("GrandBudapest")[2]) +
  geom_point(color=wes_palette("GrandBudapest")[3]) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black")) +
  scale_x_continuous(limits = c(0.5,9.5)) +
  geom_vline(xintercept=5.5, colour="red", linetype="dashed") +
  geom_text(aes(y=transcripts_count, label=transcripts_count), vjust=-0.2, size=2)+
  xlab("TACO parameter minExpr") +
  ylab("the Number of Transcripts")

ggsave("gtf_for_assembled/TACO_minExpr5.5_gffcmp_transcripts_line.pdf", 
       width=5.2, height = 5)  




 ########################################################################
#                             TACO FRAC                                #
########################################################################

### sensitivity_precision
stat1 <- read.table("gtf_for_assembled/TACO_FRAC_gffcmp.stats_sensitivity_precision.txt", 
                    sep="\t", header = T)

pp <- stat1[which(stat1$Type=="Exon_level" | stat1$Type=="Intron_level"),]
F1_exon <-  pp[pp$Sensitivity_Precision == "F1_score" &pp$Type== "Exon_level",]
F1_exon <- F1_exon[order(F1_exon$Percent, decreasing = T), ]
F1_Intron <-  pp[pp$Sensitivity_Precision == "F1_score" &pp$Type== "Intron_level",]
F1_Intron <- F1_Intron[order(F1_Intron$Percent, decreasing = T), ]
head(F1_exon)
head(F1_Intron)

ggplot(pp, aes(x=FRAC, y=Percent, color=Sensitivity_Precision)) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        text = element_text(face="bold")) +
  scale_color_manual(values = wes_palette("GrandBudapest")) +
  facet_wrap(~Type,ncol=1, scales = "free_x") +
  xlab("TACO parameter FRAC") +
  ylab("Percent(%)") +
  # scale_x_continuous(breaks = seq(1:30)/2, limits = c(1,8)) +
  guides(color=guide_legend(title = NULL))


ggsave("gtf_for_assembled/TACO_FRAC_gffcmp_sensitivity_precision_F1.pdf",
       width=6, height = 5)
#


### Novel_Missed
stat2 <- read.table("gtf_for_assembled/TACO_FRAC_gffcmp.stats_Novel_Missed.txt", sep="\t", header = T)
head(stat2)

pp <- stat2[,c(1,2, 5)]
ggplot(pp, aes(x=minExpr, y=Percent, color=Type)) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = colorRampPalette(wes_palette("Royal2"))(6)) +
  # facet_wrap(~Accuracy,ncol=2) +
  xlab("TACO parameter minExpr") +
  ylab("Percent(%)") +
  scale_x_continuous(breaks = seq(1:30)/2, limits = c(1,8)) +
  # geom_text(aes(y=Percent, label=Percent), vjust=-0.2, size=3) +
  guides(fill=F) 


# ggsave("gtf_for_assembled/gffcmp_Missed_Novel_bar.pdf", 
#        width=5.5, height = 5)  


########################################################################
#                           stringTie FRAC                             #
########################################################################

### sensitivity_precision
stat1 <- read.table("gtf_for_assembled/stringTie_FRAC_gffcmp.stats_sensitivity_precision.txt", 
                    sep="\t", header = T)

pp <- stat1[which(stat1$Type=="Exon_level" | stat1$Type=="Intron_level"),]
F1_exon <-  pp[pp$Sensitivity_Precision == "F1_score" &pp$Type== "Exon_level",]
F1_exon <- F1_exon[order(F1_exon$Percent, decreasing = T), ]
F1_Intron <-  pp[pp$Sensitivity_Precision == "F1_score" &pp$Type== "Intron_level",]
F1_Intron <- F1_Intron[order(F1_Intron$Percent, decreasing = T), ]
head(F1_exon)
head(F1_Intron)

ggplot(pp, aes(x=FRAC, y=Percent, color=Sensitivity_Precision)) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        text = element_text(face="bold")) +
  scale_color_manual(values = wes_palette("GrandBudapest")) +
  facet_wrap(~Type,ncol=1, scales = "free_x") +
  xlab("stringTie parameter FRAC") +
  ylab("Percent(%)") +
  # scale_x_continuous(breaks = seq(1:30)/2, limits = c(1,8)) +
  guides(color=guide_legend(title = NULL))


ggsave("gtf_for_assembled/stringTie_FRAC_gffcmp_sensitivity_precision_F1.pdf",
       width=6, height = 5)
#


### Novel_Missed
stat2 <- read.table("gtf_for_assembled/TACO_FRAC_gffcmp.stats_Novel_Missed.txt", sep="\t", header = T)
head(stat2)

pp <- stat2[,c(1,2, 5)]
ggplot(pp, aes(x=minExpr, y=Percent, color=Type)) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = colorRampPalette(wes_palette("Royal2"))(6)) +
  # facet_wrap(~Accuracy,ncol=2) +
  xlab("TACO parameter minExpr") +
  ylab("Percent(%)") +
  scale_x_continuous(breaks = seq(1:30)/2, limits = c(1,8)) +
  # geom_text(aes(y=Percent, label=Percent), vjust=-0.2, size=3) +
  guides(fill=F) 


# ggsave("gtf_for_assembled/gffcmp_Missed_Novel_bar.pdf", 
#        width=5.5, height = 5)  



########################################################################
#                           stringTie cov                              #
########################################################################

### sensitivity_precision
stat1 <- read.table("training_stringTie_cov/stringTie_cov_gffcmp.stats_sensitivity_precision.txt", 
                    sep="\t", header = T)

pp <- stat1[which(stat1$Type=="Exon_level" | stat1$Type=="Intron_level"),]
F1_exon <-  pp[pp$Sensitivity_Precision == "F1_score" &pp$Type== "Exon_level",]
F1_exon <- F1_exon[order(F1_exon$Percent, decreasing = T), ]
F1_Intron <-  pp[pp$Sensitivity_Precision == "F1_score" &pp$Type== "Intron_level",]
F1_Intron <- F1_Intron[order(F1_Intron$Percent, decreasing = T), ]
head(F1_exon)
head(F1_Intron)

ggplot(pp, aes(x=cov, y=Percent, color=Sensitivity_Precision)) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        text = element_text(face="bold")) +
  scale_color_manual(values = wes_palette("GrandBudapest")) +
  facet_wrap(~Type,ncol=1, scales = "free_x") +
  xlab("stringTie parameter cov") +
  ylab("Percent(%)") +
  # scale_x_continuous(breaks = seq(1:30)/2, limits = c(1,8)) +
  guides(color=guide_legend(title = NULL))


ggsave("gtf_for_assembled/stringTie_cov_gffcmp_sensitivity_precision_F1.pdf",
       width=6, height = 5)
#


### Novel_Missed
stat2 <- read.table("gtf_for_assembled/TACO_FRAC_gffcmp.stats_Novel_Missed.txt", sep="\t", header = T)
head(stat2)

pp <- stat2[,c(1,2, 5)]
ggplot(pp, aes(x=minExpr, y=Percent, color=Type)) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
        text = element_text(face="bold")) +
  scale_fill_manual(values = colorRampPalette(wes_palette("Royal2"))(6)) +
  # facet_wrap(~Accuracy,ncol=2) +
  xlab("TACO parameter minExpr") +
  ylab("Percent(%)") +
  scale_x_continuous(breaks = seq(1:30)/2, limits = c(1,8)) +
  # geom_text(aes(y=Percent, label=Percent), vjust=-0.2, size=3) +
  guides(fill=F) 


# ggsave("gtf_for_assembled/gffcmp_Missed_Novel_bar.pdf", 
#        width=5.5, height = 5)  

