

rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(scales)

# read

overlap_cor <- read.csv("overlap_cor_subcohort.csv")

################### proteins targeted by multiple aptamers

protein_dup <- unique(overlap_cor$uniprot_id[duplicated(overlap_cor$uniprot_id)]) # get proteins targeted by multiple aptamers

length(protein_dup) # how many proteins targeted by multiple aptamers

overlap_cor_protein_dup <- overlap_cor[which(overlap_cor$uniprot_id %in% protein_dup), ] # list of proteins 

length(unique(overlap_cor_protein_dup$somascan_id))

write.csv(overlap_cor_protein_dup,"overlap_cor_protein_dup_subcohort.csv", quote=F, row.names=F)

# count how many proteins there are

n_aptamer <- data.frame(table(overlap_cor_protein_dup$uniprot_id))
names(n_aptamer) <- c("uniprot_id","n_aptamer")
table(n_aptamer$n_aptamer)

# 2   3   4   5   6   8   9 
# 399  56  11   1   2   1   2

# merge with original dataset

overlap_cor_protein_dup <- merge(overlap_cor_protein_dup,n_aptamer,by="uniprot_id")

# now plot 

min(overlap_cor_protein_dup$rho_olink_soma_normal)
max(overlap_cor_protein_dup$rho_olink_soma_normal)
median(overlap_cor_protein_dup$rho_olink_soma_normal)

plot_rho_olink_soma_protein_dup_normal <- ggplot(overlap_cor_protein_dup, aes(x=rho_olink_soma_normal)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,150) +
  annotate("text", x = 0.6, y = 100, label = "Median rho = 0.29", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("OLINK vs SomaScan (ANML)") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_few() +
  theme(text = element_text(size = 12))  

plot_rho_olink_soma_protein_dup_normal

min(overlap_cor_protein_dup$rho_olink_soma_non_normal)
max(overlap_cor_protein_dup$rho_olink_soma_non_normal)
median(overlap_cor_protein_dup$rho_olink_soma_non_normal)

plot_rho_olink_soma_protein_dup_non_normal <- ggplot(overlap_cor_protein_dup, aes(x=rho_olink_soma_non_normal)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,150) +
  annotate("text", x = 0.6, y = 100, label = "Median rho = 0.35", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("OLINK vs SomaScan (non-ANML)") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_few() +
  theme(text = element_text(size = 12))  

plot_rho_olink_soma_protein_dup_non_normal

# plot everything together

plot_rho_olink_soma_protein_dup <- ggarrange(plot_rho_olink_soma_protein_dup_normal,plot_rho_olink_soma_protein_dup_non_normal,ncol=2,nrow=1)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_rho_olink_soma_protein_dup_subcohort.png",plot_rho_olink_soma_protein_dup,width=14,height=6)


################### aptamers matched to more than one olink

aptamer_dup <- unique(overlap_cor$somascan_id[duplicated(overlap_cor$somascan_id)]) # get aptamers targeted by multiple proteins

length(aptamer_dup) # how many aptamers targeted by multiple proteins

overlap_cor_aptamer_dup <- overlap_cor[which(overlap_cor$somascan_id %in% aptamer_dup), ] # list of proteins 

length(unique(overlap_cor_aptamer_dup$uniprot_id))

write.csv(overlap_cor_aptamer_dup,"overlap_cor_aptamer_dup_subcohort.csv", quote=F, row.names=F)

# plot

min(overlap_cor_aptamer_dup$rho_olink_soma_normal)
max(overlap_cor_aptamer_dup$rho_olink_soma_normal)
median(overlap_cor_aptamer_dup$rho_olink_soma_normal)

plot_rho_olink_soma_aptamer_dup_normal <- ggplot(overlap_cor_aptamer_dup, aes(x=rho_olink_soma_normal)) + 
  geom_histogram(binwidth = 0.2,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.2, 1, 0.2), lim = c(-0.25, 1.05), labels = label_number(accuracy = 0.1)) +
  ylim(0,20) +
  annotate("text", x = 0.6, y = 15, label = "Median rho = 0.09", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("OLINK vs SomaScan (ANML)") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_few() +
  theme(text = element_text(size = 12))  

plot_rho_olink_soma_aptamer_dup_normal

min(overlap_cor_aptamer_dup$rho_olink_soma_non_normal)
max(overlap_cor_aptamer_dup$rho_olink_soma_non_normal)
median(overlap_cor_aptamer_dup$rho_olink_soma_non_normal)

plot_rho_olink_soma_aptamer_dup_non_normal <- ggplot(overlap_cor_aptamer_dup, aes(x=rho_olink_soma_non_normal)) + 
  geom_histogram(binwidth = 0.2,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.2, 1, 0.2), lim = c(-0.25, 1.05), labels = label_number(accuracy = 0.1)) +
  ylim(0,20) +
  annotate("text", x = 0.6, y = 15, label = "Median rho = 0.14", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("OLINK vs SomaScan (non-ANML)") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_few() +
  theme(text = element_text(size = 12))  

plot_rho_olink_soma_aptamer_dup_non_normal

# plot everything together

plot_rho_olink_soma_aptamer_dup <- ggarrange(plot_rho_olink_soma_aptamer_dup_normal,plot_rho_olink_soma_aptamer_dup_non_normal,ncol=2,nrow=1)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_rho_olink_soma_aptamer_dup_subcohort.png",plot_rho_olink_soma_aptamer_dup,width=14,height=6)

