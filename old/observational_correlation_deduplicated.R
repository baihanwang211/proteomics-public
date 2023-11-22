rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpattern)
library(ppcor)
library(RColorBrewer)
library(gtsummary)
library(ggpubr)
library(Hmisc)
library(VennDiagram)
library(ggthemes)

########################## shortcut

overlap <- read.csv("overlap.csv")

overlap_deduplicated_normal <- overlap[order(overlap$rho_olink_somascan_normal,decreasing=TRUE),]

overlap_deduplicated_normal <- overlap_deduplicated_normal[!duplicated(overlap_deduplicated_normal$uniprot_id),]

overlap_deduplicated_non_normal <- overlap[order(overlap$rho_olink_somascan_non_normal,decreasing=TRUE),]

overlap_deduplicated_non_normal <- overlap_deduplicated_non_normal[!duplicated(overlap_deduplicated_non_normal$uniprot_id),]

median(overlap_deduplicated_normal$rho_olink_somascan_normal)

median(overlap_deduplicated_non_normal$rho_olink_somascan_non_normal)

# plot

olink_somascan_normal <- ggplot(overlap_deduplicated_normal, aes(x=rho_olink_somascan_normal)) + 
  geom_histogram(color="black", fill="white") +
  xlim(-0.5,1) +
  ylim(0,450) +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan (normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  theme_few()

olink_somascan_non_normal <- ggplot(overlap_deduplicated_non_normal, aes(x=rho_olink_somascan_non_normal)) + 
  geom_histogram(color="black", fill="white") +
  xlim(-0.5,1) +
  ylim(0,450) +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan (non-normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  theme_few()

olink_somascan <- ggarrange(olink_somascan_normal,olink_somascan_non_normal,labels=c("A","B"),ncol=1,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_deduplicated.png",olink_somascan,width=9,height=9)

## plot separate

batch <- sort(unique(overlap$batch))
dilution <- sort(unique(overlap$dilution))
dilution_label <- c("0.005%","0.5%","20%")

# normalised

olink_somascan_separate_normal_list <- list()

i <- 1

for (x in 1:2) {
  for (y in 1:3){
    olink_somascan_separate_normal_list[[i]] <- ggplot(overlap_deduplicated_normal[overlap_deduplicated_normal$batch==batch[x] & overlap_deduplicated_normal$dilution==dilution[y],], aes(x=rho_olink_somascan_normal)) + 
      geom_histogram(color="black", fill="white") +
      geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
      xlim(-0.5,1) +
      ylim(0,250) +
      ggtitle(paste0("Olink (batch ",batch[x],") vs SomaScan (",dilution_label[y], ", normalised)"))  + 
      xlab("Spearman correlation coefficient") + 
      ylab("Frequency") +
      labs(fill="SomaScan dilution") +
      theme_few()
    i <- i+1
  }
}

olink_somascan_separate_normal <- ggarrange(olink_somascan_separate_normal_list[[1]],olink_somascan_separate_normal_list[[2]],olink_somascan_separate_normal_list[[3]],olink_somascan_separate_normal_list[[4]],olink_somascan_separate_normal_list[[5]],olink_somascan_separate_normal_list[[6]],labels=c("A","B","C","D","E","F"),ncol=3,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_separate_normal_deduplicated.png",olink_somascan_separate_normal,width=16,height=9)

# non-normalised

olink_somascan_separate_non_normal_list <- list()

i <- 1

for (x in 1:2) {
  for (y in 1:3){
    olink_somascan_separate_non_normal_list[[i]] <- ggplot(overlap_deduplicated_non_normal[overlap_deduplicated_non_normal$batch==batch[x] & overlap_deduplicated_non_normal$dilution==dilution[y],], aes(x=rho_olink_somascan_non_normal)) + 
      geom_histogram(color="black", fill="white") +
      geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
      xlim(-0.5,1) +
      ylim(0,250) +
      ggtitle(paste0("Olink (batch ",batch[x],") vs SomaScan (",dilution_label[y], ", non-normalised)"))  + 
      xlab("Spearman correlation coefficient") + 
      ylab("Frequency") +
      labs(fill="SomaScan dilution") +
      theme_few()
    i <- i+1
  }
}

olink_somascan_separate_non_normal <- ggarrange(olink_somascan_separate_non_normal_list[[1]],olink_somascan_separate_non_normal_list[[2]],olink_somascan_separate_non_normal_list[[3]],olink_somascan_separate_non_normal_list[[4]],olink_somascan_separate_non_normal_list[[5]],olink_somascan_separate_non_normal_list[[6]],labels=c("A","B","C","D","E","F"),ncol=3,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_separate_non_normal_deduplicated.png",olink_somascan_separate_non_normal,width=16,height=9)

## plot (separate by batch)

# cols <- brewer.pal(n = 3, name = "Set2")

olink_1_soma_normal <- ggplot(overlap_deduplicated_normal[overlap_deduplicated_normal$batch==1,], aes(x=rho_olink_somascan_normal,
                                                              fill=factor(dilution,levels=c("0.2","0.005","5e-05"),labels=c("20%","0.5%","0.005%")))) + 
  geom_histogram(colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,300) +
  ggtitle("Olink (batch 1) vs SomaScan (normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="SomaScan dilution") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_2_soma_normal <- ggplot(overlap_deduplicated_normal[overlap_deduplicated_normal$batch==2,], aes(x=rho_olink_somascan_normal,
                                                              fill=factor(dilution,levels=c("0.2","0.005","5e-05"),labels=c("20%","0.5%","0.005%")))) + 
  geom_histogram(colour="black") +
  xlim(-0.5,1) +
  ylim(0,300) +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink (batch 2) vs SomaScan (normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="SomaScan dilution") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_1_soma_non_normal <- ggplot(overlap_deduplicated_non_normal[overlap_deduplicated_non_normal$batch==1,], aes(x=rho_olink_somascan_non_normal,
                                                                  fill=factor(dilution,levels=c("0.2","0.005","5e-05"),labels=c("20%","0.5%","0.005%")))) + 
  geom_histogram(colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,300) +
  ggtitle("Olink (batch 1) vs SomaScan (non-normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="SomaScan dilution") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_2_soma_non_normal <- ggplot(overlap_deduplicated_non_normal[overlap_deduplicated_non_normal$batch==2,], aes(x=rho_olink_somascan_non_normal,
                                                                  fill=factor(dilution,levels=c("0.2","0.005","5e-05"),labels=c("20%","0.5%","0.005%")))) + 
  geom_histogram(colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,300) +
  ggtitle("Olink (batch 2) vs SomaScan (non-normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="SomaScan dilution") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_somascan_separate_by_batch <- ggarrange(olink_1_soma_normal,olink_2_soma_normal,olink_1_soma_non_normal,olink_2_soma_non_normal,labels=c("A","B","C","D"),ncol=2,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_separate_by_batch_deduplicated.png",olink_somascan_separate_by_batch,width=19.2,height=10.8)

## plot (separate by dilution)

# cols <- brewer.pal(n = 3, name = "Pastel2")[c(1,2)]

olink_soma_5e05_normal <- ggplot(overlap_deduplicated_normal[overlap_deduplicated_normal$dilution==5e-05,], aes(x=rho_olink_somascan_normal,
                                                                        fill=factor(batch,levels=c("2","1")))) + 
  geom_histogram(bins=20,colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,500) +
  ggtitle("Olink vs SomaScan (0.005%, normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Olink Batch") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_soma_0.005_normal <- ggplot(overlap_deduplicated_normal[overlap_deduplicated_normal$dilution==0.005,], aes(x=rho_olink_somascan_normal,
                                                                         fill=factor(batch,levels=c("2","1")))) + 
  geom_histogram(bins=20,colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,500) +
  ggtitle("Olink vs SomaScan (0.5%, normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Olink Batch") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_soma_0.2_normal <- ggplot(overlap_deduplicated_normal[overlap_deduplicated_normal$dilution==0.2,], aes(x=rho_olink_somascan_normal,
                                                                     fill=factor(batch,levels=c("2","1")))) + 
  geom_histogram(bins=20,colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,500) +
  ggtitle("Olink vs SomaScan (20%, normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Olink Batch") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_soma_5e05_non_normal <- ggplot(overlap_deduplicated_non_normal[overlap_deduplicated_non_normal$dilution==5e-05,], aes(x=rho_olink_somascan_non_normal,
                                                                            fill=factor(batch,levels=c("2","1")))) + 
  geom_histogram(bins=20,colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,500) +
  ggtitle("Olink vs SomaScan (0.005%, non-normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Olink Batch") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_soma_0.005_non_normal <- ggplot(overlap_deduplicated_non_normal[overlap_deduplicated_non_normal$dilution==0.005,], aes(x=rho_olink_somascan_non_normal,
                                                                             fill=factor(batch,levels=c("2","1")))) + 
  geom_histogram(bins=20,colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,500) +
  ggtitle("Olink vs SomaScan (0.5%, non-normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Olink Batch") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_soma_0.2_non_normal <- ggplot(overlap_deduplicated_non_normal[overlap_deduplicated_non_normal$dilution==0.2,], aes(x=rho_olink_somascan_non_normal,
                                                                         fill=factor(batch,levels=c("2","1")))) + 
  geom_histogram(bins=20,colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,500) +
  ggtitle("Olink vs SomaScan (20%, non-normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Olink Batch") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_somascan_separate_by_dilution <- ggarrange(olink_soma_5e05_normal,olink_soma_0.005_normal,olink_soma_0.2_normal,olink_soma_5e05_non_normal,olink_soma_0.005_non_normal,olink_soma_0.2_non_normal,labels=c("A","B","C","D","E","F"),ncol=3,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_separate_by_dilution_deduplicated.png",olink_somascan_separate_by_dilution,width=19.2,height=10.8)

## plot (combined by batch)

# get colour palette

# cols <- brewer.pal(n = 3, name = "Pastel1")

# normalised

olink_soma_normal_batch <- ggplot(overlap_deduplicated_normal) + 
  geom_histogram_pattern(aes(x=rho_olink_somascan_normal,
                             pattern = factor(batch,levels=c("2","1")),
                             fill=factor(dilution,levels=c("0.2","0.005","5e-05"),labels=c("20%","0.5%","0.005%"))), 
                         colour="black", 
                         pattern_angle = 45, 
                         pattern_fill = "black",
                         pattern_density=0.001,
                         pattern_spacing=0.01) +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  # scale_fill_manual(values=cols) +
  scale_pattern_manual(values=c('stripe',"none")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = c("none", "none", "none")))) +
  theme_few() +
  ggtitle("Olink vs SomaScan (normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") + labs(pattern ='Olink batch',fill="SomaScan dilution")

# non-normalised

olink_soma_non_normal_batch <- ggplot(overlap_deduplicated_non_normal) + 
  geom_histogram_pattern(aes(x=rho_olink_somascan_non_normal,
                             pattern = factor(batch,levels=c("2","1")),
                             fill=factor(dilution,levels=c("0.2","0.005","5e-05"),labels=c("20%","0.5%","0.005%"))), 
                         colour="black", 
                         pattern_angle = 45, 
                         pattern_fill = "black",
                         pattern_density=0.001,
                         pattern_spacing=0.01) +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  # scale_fill_manual(values=cols) +
  scale_pattern_manual(values=c('stripe',"none")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = c("none", "none", "none")))) +
  theme_few() +
  ggtitle("Olink vs SomaScan (non-normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") + labs(pattern ='Olink batch',fill="SomaScan dilution")

olink_somascan_combined_batch <- ggarrange(olink_soma_normal_batch,olink_soma_non_normal_batch,labels=c("A","B"),ncol=1,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_combined_batch_deduplicated.png",olink_somascan_combined_batch,width=9.6,height=10.8)

## plot (combined by dilution)

# cols <- brewer.pal(n = 3, name = "Pastel1")

# normalised

olink_soma_normal_dilution <- ggplot(overlap_deduplicated_normal) + 
  geom_histogram_pattern(aes(x=rho_olink_somascan_normal,
                             fill=factor(dilution,levels=c("0.2","0.005","5e-05"),labels=c("20%","0.5%","0.005%")),
                             pattern = factor(batch,levels=c("2","1"))), 
                         colour="black", 
                         pattern_angle = 45, 
                         pattern_fill = "black",
                         pattern_density=0.001,
                         pattern_spacing=0.01) +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,450) +
  # scale_fill_manual(values=cols) +
  scale_pattern_manual(values=c('stripe',"none")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = c("none", "none", "none")))) +
  theme_few() +
  ggtitle("Olink vs SomaScan (normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") + labs(pattern ='Olink batch',fill="SomaScan dilution")

# non-normalised

olink_soma_non_normal_dilution <- ggplot(overlap_deduplicated_non_normal) + 
  geom_histogram_pattern(aes(x=rho_olink_somascan_non_normal,
                             fill=factor(dilution,levels=c("0.2","0.005","5e-05"),labels=c("20%","0.5%","0.005%")),
                             pattern = factor(batch,levels=c("2","1"))), 
                         colour="black", 
                         pattern_angle = 45, 
                         pattern_fill = "black",
                         pattern_density=0.001,
                         pattern_spacing=0.01) +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,450) +
  # scale_fill_manual(values=cols) +
  scale_pattern_manual(values=c('stripe',"none")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = c("none", "none", "none")))) +
  theme_few() +
  ggtitle("Olink vs SomaScan (non-normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") + labs(pattern ='Olink batch',fill="SomaScan dilution")

olink_somascan_combined_dilution <- ggarrange(olink_soma_normal_dilution,olink_soma_non_normal_dilution,labels=c("A","B"),ncol=1,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_combined_dilution_deduplicated.png",olink_somascan_combined_dilution,width=9.6,height=10.8)
