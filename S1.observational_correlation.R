##### this script estimates the observational correlations between Olink and SomaScan protein measurements

rm(list = ls())

setwd("")

library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpattern)
library(ggpubr)
library(gtsummary)
library(Hmisc)
library(ggthemes)
library(datawizard)
library(ckbplotr)
library(scales)
library(reshape2)

##### read overlapping proteins and find those that are one-to-one matched

overlap <- read_excel("") # read file that lists overlapping proteins

protein_dup <- unique(overlap$uniprot_id[duplicated(overlap$uniprot_id)]) # get proteins targeted by multiple aptamers

length(protein_dup) # how many proteins targeted by multiple aptamers

overlap_protein_dup <- overlap[which(overlap$uniprot_id %in% protein_dup), ] # list of proteins 

length(unique(overlap_protein_dup$somascan_id))

write.csv(overlap_protein_dup,"overlap_protein_dup.csv", quote=F, row.names=F)

aptamer_dup <- unique(overlap$somascan_id[duplicated(overlap$somascan_id)]) # get aptamers targeted by multiple proteins

length(aptamer_dup) # how many aptamers targeted by multiple proteins

overlap_aptamer_dup <- overlap[which(overlap$somascan_id %in% aptamer_dup), ] # list of proteins 

length(unique(overlap_aptamer_dup$uniprot_id))

write.csv(overlap_aptamer_dup,"overlap_aptamer_dup.csv", quote=F, row.names=F)

# remove duplicates

overlap_1_to_1 <- overlap[-c(which(overlap$uniprot_id %in% protein_dup),which(overlap$somascan_id %in% aptamer_dup)), ]

write.csv(overlap_1_to_1,"overlap_1_to_1.csv", quote=F, row.names=F)

##### load raw olink and somascan data

olink <- read.csv("olink_overlap.csv")
somascan_non_ANML <- read.csv("somascan_non_ANML_overlap.csv")

# log transform somascan data
somascan_non_ANML_log <- somascan_non_ANML
somascan_non_ANML_log[2:length(somascan_non_ANML_log)] <- lapply(somascan_non_ANML_log[2:length(somascan_non_ANML_log)],log)

##### check correlation

overlap_cor <- overlap

# spearman

for (i in 1:nrow(overlap_cor)) {
  
  df <- merge(olink[,c("csid",overlap_cor$olink_id[i])],somascan_non_ANML[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$rho_olink_soma_non_ANML[i] <- cor.test(df[[2]], df[[3]], method = 'spearman')$estimate
  
  print(i)
  
}

# pearson

for (i in 1:nrow(overlap_cor)) {

  df <- merge(olink[,c("csid",overlap_cor$olink_id[i])],somascan_non_ANML_log[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_non_ANML_log[i] <- cor.test(df[[2]], df[[3]], method = 'pearson')$estimate
  
  print(i)
  
}


# save

write.csv(overlap_cor,"overlap_cor.csv", quote=F, row.names=F)

##### plot distribution

# calculate median rho

median_rho <- format(round(median(overlap_cor$rho_olink_soma_non_ANML),2),nsmall=2)

# spearman olink vs soma

plot_rho_olink_soma_non_ANML <- ggplot(overlap_cor, aes(x=rho_olink_soma_non_ANML)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,450) +
  annotate("text", x = 0.55, y = 200, label = paste("Median rho =", median_rho), size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_non_ANML)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18),text = element_text(size = 14))  

plot_rho_olink_soma_non_ANML

# pearson olink vs soma

median_r <- format(round(median(overlap_cor$r_olink_soma_non_ANML),2),nsmall=2)

plot_r_olink_soma_non_ANML <- ggplot(overlap_cor, aes(x=r_olink_soma_non_ANML_log)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,500) +
  annotate("text", x = 0.5, y = 200, label = paste("Median r =", median_r), size = 12/.pt) +
  geom_vline(aes(xintercept=median(r_olink_soma_non_ANML_log)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") + xlab("Pearson's r") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18),text = element_text(size = 14))  

plot_r_olink_soma_non_ANML


##### only keep 1_to_1 matched pairs

overlap_1_to_1_cor <- overlap_cor[-c(which(overlap_cor$uniprot_id %in% protein_dup),which(overlap_cor$somascan_id %in% aptamer_dup)), ]

write.csv(overlap_1_to_1_cor,"overlap_1_to_1_cor.csv", quote=F, row.names=F)

##### plot distribution

# calculate median rho

median(overlap_1_to_1_cor$rho_olink_soma_non_ANML)

# spearman olink vs soma

median_rho <- format(round(median(overlap_1_to_1_cor$rho_olink_soma_non_ANML),2),nsmall=2)

plot_rho_olink_soma_non_ANML <- ggplot(overlap_1_to_1_cor, aes(x=rho_olink_soma_non_ANML)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,450) +
  annotate("text", x = 0.55, y = 200, label = paste("Median rho =", median_rho), size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_non_ANML)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18),text = element_text(size = 14))  

plot_rho_olink_soma_non_ANML

# pearson olink vs soma

median_r <- format(round(median(overlap_1_to_1_cor$r_olink_soma_non_ANML),2),nsmall=2)

plot_r_olink_soma_non_ANML <- ggplot(overlap_1_to_1_cor, aes(x=r_olink_soma_non_ANML_log)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,500) +
  annotate("text", x = 0.5, y = 200, label = paste("Median r =", median_r), size = 12/.pt) +
  geom_vline(aes(xintercept=median(r_olink_soma_non_ANML_log)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") + xlab("Pearson's r") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18),text = element_text(size = 14))  

plot_r_olink_soma_non_ANML
