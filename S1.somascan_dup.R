rm(list = ls())

library(tidyverse)
library(dplyr)
library(ggplot2)
library(scales)
library(ckbplotr)
library(ggpubr)

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

somascan_normalised_log <- read.csv("somascan_normalised_log.csv")
somascan_non_normalised_log <- read.csv("somascan_non_normalised_log.csv")

somascan_protein <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/somalogic_meta.csv")[,c("aptname","uniprot","organism")]

# human aptamers n = 7335
somascan_protein <- somascan_protein[somascan_protein$organism=="Human",]
somascan_protein <- somascan_protein[,-3]

# aptamers with uniprot id n=7301
somascan_protein <- somascan_protein[somascan_protein$uniprot!="",]

# split different uniprot ids into different rows
somascan_protein <- separate_rows(somascan_protein, "uniprot", sep = "\\|", convert = T)

# count how many aptamers per protein
n_aptamer <- data.frame(table(somascan_protein$uniprot))
names(n_aptamer) <- c("uniprot","n_aptamer")

# merge
somascan_protein <- merge(somascan_protein,n_aptamer,by="uniprot")

# count
table(somascan_protein$n_aptamer)

#1    2    3    4    5    6    7    8    9 
#5562 1450  252   64   10   18    7    8   18

# calculate correlations between each pair for the same protein
# pairwise correlations

protein_dup <- unique(somascan_protein$uniprot[somascan_protein$n_aptamer>1])

pairwise_cor_normal <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("aptamer_1", "aptamer_2", "rho", "r", "uniprot"))))

for (i in 1:length(protein_dup)) {
  aptamer <- somascan_protein$aptname[somascan_protein$uniprot==protein_dup[i]]
  df <- select(somascan_normalised_log,aptamer)
  rho <- cor(df,method = "spearman")
  r <- cor(df,method = "pearson")
  cor_results <- data.frame(aptamer_1=rownames(rho)[row(rho)[upper.tri(rho)]], 
                            aptamer_2=colnames(rho)[col(rho)[upper.tri(rho)]], 
                            rho=rho[upper.tri(rho)],
                            r=r[upper.tri(r)])
  cor_results$uniprot <- protein_dup[i]
  pairwise_cor_normal <- rbind(pairwise_cor_normal,cor_results)
}

write.csv(pairwise_cor_normal,"pairwise_cor_soma_normal_all.csv",quote=F, row.names=F)

pairwise_cor_non_normal <- data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("aptamer_1", "aptamer_2", "rho", "r", "uniprot"))))

for (i in 1:length(protein_dup)) {
  aptamer <- somascan_protein$aptname[somascan_protein$uniprot==protein_dup[i]]
  df <- select(somascan_non_normalised_log,aptamer)
  rho <- cor(df,method = "spearman")
  r <- cor(df,method = "pearson")
  cor_results <- data.frame(aptamer_1=rownames(rho)[row(rho)[upper.tri(rho)]], 
                            aptamer_2=colnames(rho)[col(rho)[upper.tri(rho)]], 
                            rho=rho[upper.tri(rho)],
                            r=r[upper.tri(r)])
  cor_results$uniprot <- protein_dup[i]
  pairwise_cor_non_normal <- rbind(pairwise_cor_non_normal,cor_results)
}

write.csv(pairwise_cor_non_normal,"pairwise_cor_soma_non_normal_all.csv",quote=F, row.names=F)

# pairwise_cor_normal <- read.csv("pairwise_cor_normal.csv") 
# pairwise_cor_non_normal <- read.csv("pairwise_cor_non_normal.csv") 

# plot

min(pairwise_cor_normal$rho)
max(pairwise_cor_normal$rho)
median(pairwise_cor_normal$rho)

plot_protein_dup_normal <- ggplot(pairwise_cor_normal, aes(x=rho)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,150) +
  annotate("text", x = 0.4, y = 150, label = "Median rho = 0.15", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("SomaScan (ANML)") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12))  

plot_protein_dup_normal

min(pairwise_cor_non_normal$rho)
max(pairwise_cor_non_normal$rho)
median(pairwise_cor_non_normal$rho)

plot_protein_dup_non_normal <- ggplot(pairwise_cor_non_normal, aes(x=rho)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,150) +
  annotate("text", x = 0.7, y = 150, label = "Median rho = 0.44", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("SomaScan (non-ANML)") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12))  

plot_protein_dup_non_normal

plot_protein_dup <- ggarrange(plot_protein_dup_normal,plot_protein_dup_non_normal,ncol=2,nrow=1)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_soma_protein_dup.png",plot_protein_dup,width=14,height=6)
