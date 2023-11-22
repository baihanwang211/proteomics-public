rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(tidyverse)
library(ggplot2)
library(ggpattern)
library(ggpubr)
# library(ppcor)
library(gtsummary)
library(Hmisc)
library(ggthemes)
library(datawizard)
library(scales)
library(ckbplotr)


######################## load raw, transformed and cleaned somascan data

overlap <- read.csv("overlap.csv")
olink <- read.csv("olink.csv")
somascan_normalised_log <- read.csv("somascan_normalised_log.csv")
somascan_non_normalised_log <- read.csv("somascan_non_normalised_log.csv")

# get ids for subcohort participants

ascertainment <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/data_baseline_ascertainments.csv")
subcohort_id <- ascertainment$csid[ascertainment$olinkexp1536_chd_b1_subcohort==1]

# only keep subcohort

olink <- olink[which(olink$csid %in% subcohort_id),]
somascan_normalised_log <- somascan_normalised_log[which(somascan_normalised_log$csid %in% subcohort_id),]
somascan_non_normalised_log <- somascan_non_normalised_log[which(somascan_non_normalised_log$csid %in% subcohort_id),]

########## check correlation

overlap_cor <- overlap

for (i in 1:nrow(overlap_cor)) {
  
  df <- merge(olink[,c("csid",overlap_cor$olink_id[i])],somascan_normalised_log[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$rho_olink_soma_normal[i] <- cor.test(df[[2]], df[[3]], method = 'spearman')$estimate
  
  df <- merge(olink[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised_log[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$rho_olink_soma_non_normal[i] <- cor.test(df[[2]], df[[3]], method = 'spearman')$estimate
  
}

write.csv(overlap_cor,"overlap_cor_subcohort.csv", quote=F, row.names=F)

# keep one-to-one

protein_dup <- unique(overlap_cor$uniprot_id[duplicated(overlap_cor$uniprot_id)]) # get proteins targeted by multiple aptamers

aptamer_dup <- unique(overlap_cor$somascan_id[duplicated(overlap_cor$somascan_id)]) # get aptamers targeted by multiple proteins

overlap_1_to_1_cor <- overlap_cor[-c(which(overlap_cor$uniprot_id %in% protein_dup),which(overlap_cor$somascan_id %in% aptamer_dup)), ]

write.csv(overlap_1_to_1_cor,"overlap_1_to_1_cor_subcohort.csv", quote=F, row.names=F)

# overlap_1_to_1_cor <- read.csv("overlap_1_to_1_cor_subcohort.csv")

##################################### get median

names(overlap_1_to_1_cor)

coeff_sum <- data.frame(coeff = names(overlap_1_to_1_cor)[9:10],
                        median = numeric(2), max = numeric(2), min = numeric(2))

coeff_sum$median[1:2] <- lapply(overlap_1_to_1_cor[,c(9:10)],median)
coeff_sum$max[1:2] <- lapply(overlap_1_to_1_cor[,c(9:10)],max)
coeff_sum$min[1:2] <- lapply(overlap_1_to_1_cor[,c(9:10)],min)

print(coeff_sum, row.names=F)

##################################### plot

# spearman olink vs soma

plot_rho_olink_soma_normal <- ggplot(overlap_1_to_1_cor, aes(x=rho_olink_soma_normal)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  annotate("text", x = 0.4, y = 200, label = "Median rho = 0.19", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("OLINK vs SomaScan-ANML", subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12))  

plot_rho_olink_soma_normal

plot_rho_olink_soma_non_normal <- ggplot(overlap_1_to_1_cor, aes(x=rho_olink_soma_non_normal)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  annotate("text", x = 0.5, y = 200, label = "Median rho = 0.24", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("OLINK vs SomaScan (non-ANML)", subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12))  

plot_rho_olink_soma_non_normal

plot_rho_olink_soma <- ggarrange(plot_rho_olink_soma_normal,plot_rho_olink_soma_non_normal,ncol=2,nrow=1)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_rho_olink_soma_1_to_1_subcohort.png",plot_rho_olink_soma,width=14,height=6)
