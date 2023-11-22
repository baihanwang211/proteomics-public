rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(scales)
library(dplyr)
library(ckbplotr)

# read

overlap_cor <- read.csv("overlap_cor.csv")
olink <- read.csv("olink.csv")
somascan_normalised_log <- read.csv("somascan_normalised_log.csv")
somascan_non_normalised_log <- read.csv("somascan_non_normalised_log.csv")

################### proteins targeted by multiple aptamers

protein_dup <- unique(overlap_cor$uniprot_id[duplicated(overlap_cor$uniprot_id)]) # get proteins targeted by multiple aptamers

length(protein_dup) # how many proteins targeted by multiple aptamers

overlap_cor_protein_dup <- overlap_cor[which(overlap_cor$uniprot_id %in% protein_dup), ] # list of proteins 

length(unique(overlap_cor_protein_dup$somascan_id))

write.csv(overlap_cor_protein_dup,"overlap_cor_protein_dup.csv", quote=F, row.names=F)

# overlap_cor_protein_dup <- read.csv("overlap_cor_protein_dup.csv")

# count how many proteins there are

n_aptamer <- data.frame(table(overlap_cor_protein_dup$uniprot_id))
names(n_aptamer) <- c("uniprot_id","n_aptamer")
table(n_aptamer$n_aptamer)

# 2   3   4   5   6   8   9 
# 399  56  11   1   2   1   2

# merge with original dataset

overlap_cor_protein_dup <- merge(overlap_cor_protein_dup,n_aptamer,by="uniprot_id")

## now plot olink vs SomaScan correlation

min(overlap_cor_protein_dup$rho_olink_soma_normal)
max(overlap_cor_protein_dup$rho_olink_soma_normal)
median(overlap_cor_protein_dup$rho_olink_soma_normal)

plot_rho_olink_soma_protein_dup_normal <- ggplot(overlap_cor_protein_dup, aes(x=rho_olink_soma_normal)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,150) +
  annotate("text", x = 0.6, y = 100, label = "Median rho = 0.31", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("OLINK vs SomaScan-ANML",subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12))  

plot_rho_olink_soma_protein_dup_normal

min(overlap_cor_protein_dup$rho_olink_soma_non_normal)
max(overlap_cor_protein_dup$rho_olink_soma_non_normal)
median(overlap_cor_protein_dup$rho_olink_soma_non_normal)

plot_rho_olink_soma_protein_dup_non_normal <- ggplot(overlap_cor_protein_dup, aes(x=rho_olink_soma_non_normal)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,150) +
  annotate("text", x = 0.6, y = 100, label = "Median rho = 0.37", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("OLINK vs SomaScan-non-ANML",subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12))  

plot_rho_olink_soma_protein_dup_non_normal

plot_rho_olink_soma_protein_dup <- ggarrange(plot_rho_olink_soma_protein_dup_normal,plot_rho_olink_soma_protein_dup_non_normal,ncol=2,nrow=1)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_rho_olink_soma_protein_dup.png",plot_rho_olink_soma_protein_dup,width=14,height=6)


## plot correlations between aptamers

pairwise_cor_soma_normal <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("aptamer_1", "aptamer_2", "rho", "r", "uniprot_id", "olink_id"))))

for (i in 1:length(protein_dup)) {
  df <- select(somascan_normalised_log,overlap_cor$somascan_id[overlap_cor$uniprot_id==protein_dup[i]])
  rho <- cor(df,method = "spearman")
  r <- cor(df,method = "pearson")
  cor_results_soma <- data.frame(aptamer_1=rownames(rho)[row(rho)[upper.tri(rho)]], 
                            aptamer_2=colnames(rho)[col(rho)[upper.tri(rho)]], 
                            rho=rho[upper.tri(rho)],
                            r=r[upper.tri(r)])
  cor_results_soma$uniprot_id <- protein_dup[i]
  cor_results_soma$olink_id <- unique(overlap_cor$olink_id[overlap_cor$uniprot_id==protein_dup[i]])
  pairwise_cor_soma_normal <- rbind(pairwise_cor_soma_normal,cor_results_soma)
}

write.csv(pairwise_cor_soma_normal,"pairwise_cor_soma_normal.csv", quote=F, row.names=F)
# pairwise_cor_soma_normal <- read.csv("pairwise_cor_soma_normal.csv")

pairwise_cor_soma_non_normal <- data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("aptamer_1", "aptamer_2", "rho", "r", "uniprot_id", "olink_id"))))

for (i in 1:length(protein_dup)) {
  df <- select(somascan_non_normalised_log,overlap_cor$somascan_id[overlap_cor$uniprot_id==protein_dup[i]])
  rho <- cor(df,method = "spearman")
  r <- cor(df,method = "pearson")
  cor_results_soma <- data.frame(aptamer_1=rownames(rho)[row(rho)[upper.tri(rho)]], 
                                 aptamer_2=colnames(rho)[col(rho)[upper.tri(rho)]], 
                                 rho=rho[upper.tri(rho)],
                                 r=r[upper.tri(r)])
  cor_results_soma$uniprot_id <- protein_dup[i]
  cor_results_soma$olink_id <- unique(overlap_cor$olink_id[overlap_cor$uniprot_id==protein_dup[i]])
  pairwise_cor_soma_non_normal <- rbind(pairwise_cor_soma_non_normal,cor_results_soma)
}

write.csv(pairwise_cor_soma_non_normal,"pairwise_cor_soma_non_normal.csv", quote=F, row.names=F)
# pairwise_cor_soma_non_normal <- read.csv("pairwise_cor_soma_non_normal.csv")

# plot

min(pairwise_cor_soma_normal$rho)
max(pairwise_cor_soma_normal$rho)
median(pairwise_cor_soma_normal$rho)

plot_pairwise_cor_soma_normal <- ggplot(pairwise_cor_soma_normal, aes(x=rho)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,100) +
  annotate("text", x = 0.4, y = 70, label = "Median rho = 0.18", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("SomaScan-ANML",subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12))  

plot_pairwise_cor_soma_normal

min(pairwise_cor_soma_non_normal$rho)
max(pairwise_cor_soma_non_normal$rho)
median(pairwise_cor_soma_non_normal$rho)

plot_pairwise_cor_soma_non_normal <- ggplot(pairwise_cor_soma_non_normal, aes(x=rho)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,100) +
  annotate("text", x = 0.7, y = 70, label = "Median rho = 0.46", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("SomaScan-non-ANML",subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12))  

plot_pairwise_cor_soma_non_normal

plot_pairwise_cor_soma <- ggarrange(plot_pairwise_cor_soma_normal,plot_pairwise_cor_soma_non_normal,ncol=2,nrow=1)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_pairwise_cor_soma.png",plot_pairwise_cor_soma,width=14,height=6)



################### aptamers matched to more than one olink

aptamer_dup <- unique(overlap_cor$somascan_id[duplicated(overlap_cor$somascan_id)]) # get aptamers targeted by multiple proteins

length(aptamer_dup) # how many aptamers targeted by multiple proteins

overlap_cor_aptamer_dup <- overlap_cor[which(overlap_cor$somascan_id %in% aptamer_dup), ] # list of proteins 

length(unique(overlap_cor_aptamer_dup$uniprot_id))

write.csv(overlap_cor_aptamer_dup,"overlap_cor_aptamer_dup.csv", quote=F, row.names=F)

# overlap_cor_aptamer_dup <- read.csv("overlap_cor_aptamer_dup.csv")

## plot olink vs soma

min(overlap_cor_aptamer_dup$rho_olink_soma_normal)
max(overlap_cor_aptamer_dup$rho_olink_soma_normal)
median(overlap_cor_aptamer_dup$rho_olink_soma_normal)

plot_rho_olink_soma_aptamer_dup_normal <- ggplot(overlap_cor_aptamer_dup, aes(x=rho_olink_soma_normal)) + 
  geom_histogram(binwidth = 0.2,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.2, 1, 0.2), lim = c(-0.25, 1.05), labels = label_number(accuracy = 0.1)) +
  ylim(0,20) +
  annotate("text", x = 0.3, y = 15, label = "Median rho = 0.07", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("OLINK vs SomaScan-ANML",subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12))  

plot_rho_olink_soma_aptamer_dup_normal

min(overlap_cor_aptamer_dup$rho_olink_soma_non_normal)
max(overlap_cor_aptamer_dup$rho_olink_soma_non_normal)
median(overlap_cor_aptamer_dup$rho_olink_soma_non_normal)

plot_rho_olink_soma_aptamer_dup_non_normal <- ggplot(overlap_cor_aptamer_dup, aes(x=rho_olink_soma_non_normal)) + 
  geom_histogram(binwidth = 0.2,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.2, 1, 0.2), lim = c(-0.25, 1.05), labels = label_number(accuracy = 0.1)) +
  ylim(0,20) +
  annotate("text", x = 0.4, y = 15, label = "Median rho = 0.16", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("OLINK vs SomaScan-non-ANML",subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12))  

plot_rho_olink_soma_aptamer_dup_non_normal

# plot everything together

plot_rho_olink_soma_aptamer_dup <- ggarrange(plot_rho_olink_soma_aptamer_dup_normal,plot_rho_olink_soma_aptamer_dup_non_normal,ncol=2,nrow=1)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_rho_olink_soma_aptamer_dup.png",plot_rho_olink_soma_aptamer_dup,width=14,height=6)


## plot correlation between olink reagents

pairwise_cor_olink <- data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c("olink_1", "olink_2", "rho", "r", "uniprot_id_1", "uniprot_id_2", "somascan_id"))))

for (i in 1:length(aptamer_dup)) {
  df <- select(olink,overlap_cor$olink_id[overlap_cor$somascan_id==aptamer_dup[i]])
  rho <- cor(df,method = "spearman")
  r <- cor(df,method = "pearson")
  cor_results_olink <- data.frame(olink_1=rownames(rho)[row(rho)[upper.tri(rho)]], 
                                 olink_2=colnames(rho)[col(rho)[upper.tri(rho)]], 
                                 rho=rho[upper.tri(rho)],
                                 r=r[upper.tri(r)])
  # this is only correct for 1 amptamer matched to 2 olink
  cor_results_olink$uniprot_1 <- unique(overlap_cor$uniprot_id[overlap_cor$olink_id==cor_results_olink$olink_1])
  cor_results_olink$uniprot_2 <- unique(overlap_cor$uniprot_id[overlap_cor$olink_id==cor_results_olink$olink_2])
  cor_results_olink$somascan_id <- aptamer_dup[i]
  pairwise_cor_olink <- rbind(pairwise_cor_olink,cor_results_olink)
}

write.csv(pairwise_cor_olink,"pairwise_cor_olink.csv", quote=F, row.names=F)
# pairwise_cor_olink <- read.csv("pairwise_cor_olink.csv")

# plot

min(pairwise_cor_olink$rho)
max(pairwise_cor_olink$rho)
median(pairwise_cor_olink$rho)

plot_pairwise_cor_olink <- ggplot(pairwise_cor_olink, aes(x=rho)) + 
  geom_histogram(binwidth = 0.1,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2), lim = c(0, 0.8), labels = label_number(accuracy = 0.1)) +
  ylim(0,6) +
  annotate("text", x = 0.45, y = 4, label = "Median r = 0.31", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("OLINK",subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12))  

plot_pairwise_cor_olink

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_pairwise_cor_olink.png",plot_pairwise_cor_olink,width=7,height=6)
