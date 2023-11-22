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
library(datawizard)

## load olink data

# firstly load linkage file for uniprot id

olink_protein <- read_excel("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00190-V1/olink_uniprot.xlsx")[c(1,3,4)]

names(olink_protein) <- c("uniprot_id","olink_id","panel")

olink_protein <- olink_protein[!duplicated(olink_protein$olink_id),]

olink_protein$olink_id <- tolower(olink_protein$olink_id)

olink_protein$olink_id <- paste0("ol_", olink_protein$olink_id)

olink_protein$olink_id <- gsub("-","_",olink_protein$olink_id)

olink_protein$batch <- 1

olink_protein[grep("II",olink_protein$panel),]$batch <-2

olink_protein$panel <- gsub("_II","",olink_protein$panel)

## load olink ckb data

olink_cardiometabolic <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00190-V1/data_baseline_olink_cardiometabolic.csv")

olink_inflammation <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00190-V1/data_baseline_olink_inflammation.csv")

olink_neurology <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00190-V1/data_baseline_olink_neurology.csv")

olink_oncology <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00190-V1/data_baseline_olink_oncology.csv")

olink_all <- olink_cardiometabolic %>% inner_join(olink_inflammation, by='csid', suffix=c("","_dup")) %>% inner_join(olink_neurology, by='csid', suffix=c("","_dup")) %>% inner_join(olink_oncology, by='csid', suffix=c("","_dup"))

# find difference between proteins in the assay and ckb data

setdiff(olink_protein$olink_id,names(olink_all)[-1])

setdiff(names(olink_all)[-1],olink_protein$olink_id)

# only keep those proteins that are also included in assay list

olink_unique <- olink_all[,-grep("dup",names(olink_all))]

olink <- olink_all[, c(1, which(names(olink_all) %in% olink_protein$olink_id))]

# find duplicated proteins

olink_protein_dup <- names(olink_all)[grep("_dup",names(olink_all))]

olink_protein_dup <- gsub("_dup", "", olink_protein_dup[1:6])

# check where they are from

olink_protein_dup_df <- as.data.frame(olink_protein_dup)

for (i in 1:length(olink_protein_dup)){
  olink_protein_dup_df$cardiometabolic[i] <- olink_protein_dup[i] %in% names(olink_cardiometabolic)
  olink_protein_dup_df$inflammation[i] <- olink_protein_dup[i] %in% names(olink_inflammation)
  olink_protein_dup_df$neurology[i] <- olink_protein_dup[i] %in% names(olink_neurology)
  olink_protein_dup_df$oncology[i] <- olink_protein_dup[i] %in% names(olink_oncology)
}

# get data for each duplicated protein

olink_dup_list <- list()

for (i in 1:length(olink_protein_dup)){
  olink_dup_list[[i]] <- as.data.frame(olink_all$csid)
  names(olink_dup_list[[i]]) <- "csid"
  olink_dup_list[[i]]$cardiometabolic <- olink_cardiometabolic[,grep(olink_protein_dup[i],names(olink_cardiometabolic))]
  olink_dup_list[[i]]$inflammation <- olink_inflammation[,grep(olink_protein_dup[i],names(olink_inflammation))]
  olink_dup_list[[i]]$neurology <- olink_neurology[,grep(olink_protein_dup[i],names(olink_neurology))]
  olink_dup_list[[i]]$oncology <- olink_oncology[,grep(olink_protein_dup[i],names(olink_oncology))]
  olink_dup_list[[i]]$protein <- olink_protein_dup[i]
}

# they are actually all the same

# update the linkage file

olink_protein <- olink_protein[which(olink_protein$olink_id %in% names(olink)),]

## load somascan data

# load linkage file first

somascan_protein <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00190-V1/somalogic_meta.csv")[,c("aptname","uniprot","dilution2")]

names(somascan_protein) <- c("somascan_id","uniprot_id","dilution")

# split uniprot column into 2 because some  aptamers have more than 1 uniprot id

somascan_protein <- separate(data = somascan_protein, col = uniprot_id, into = c("uniprot_id_1", "uniprot_id_2","uniprot_id_3"), sep = "\\|")

# load normalised somascan data and only keep overlapping proteins

somascan_normalised <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00190-V1/data_baseline_somalogic.csv")

# somascan_normalised <- somascan_normalised[which(somascan_normalised$variable %in% somascan_protein$somascan_id),]

# change to wide format

somascan_normalised <- pivot_wider(somascan_normalised, names_from = variable, values_from = value)

somascan_normalised <- somascan_normalised[somascan_normalised$qc==0,-c(2,3)]

## load non-normalised somascan data and only keep overlapping proteins

somascan_non_normalised <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00190-V1/data_baseline_somalogic_qc.csv")

# somascan_non_normalised <- somascan_non_normalised[which(somascan_non_normalised$variable %in% somascan_protein$somascan_id),]

# change to wide format

somascan_non_normalised <- pivot_wider(somascan_non_normalised, names_from = variable, values_from = value)

somascan_non_normalised <- somascan_non_normalised[somascan_non_normalised$qc==0,-c(2,3)]


## merge 

overlap_1 <- merge(olink_protein, somascan_protein, by.x="uniprot_id", by.y="uniprot_id_1")
names(overlap_1)[c(6,7)] <- c("somascan_uniprot_2","somascan_uniprot_3")

overlap_2 <- merge(olink_protein, somascan_protein, by.x="uniprot_id", by.y="uniprot_id_2")
names(overlap_2)[c(6,7)] <- c("somascan_uniprot_2","somascan_uniprot_3")

overlap_3 <- merge(olink_protein, somascan_protein, by.x="uniprot_id", by.y="uniprot_id_3")

overlap <- rbind(overlap_1,overlap_2)

olink_soma_normal <- merge(olink, somascan_normalised, by="csid")

olink_soma_normal<- olink_soma_normal[,c(1,which(names(olink_soma_normal) %in% c(overlap$olink_id,overlap$somascan_id)))]

olink_soma_non_normal <- merge(olink, somascan_non_normalised, by="csid")

olink_soma_non_normal<- olink_soma_non_normal[,c(1,which(names(olink_soma_non_normal) %in% c(overlap$olink_id,overlap$somascan_id)))]

# load demographic variables and merge

baseline <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00190-V1/data_baseline_questionnaires.csv")

baseline_var <- baseline[,c("csid","age_at_study_date_x100","is_female", "region_code", "region_is_urban")]

baseline_var$age_at_study_date_x100 <- baseline_var$age_at_study_date_x100/100

names(baseline_var)[2] <- "age"

olink_soma_normal <- merge(baseline_var, olink_soma_normal, by="csid")

olink_soma_non_normal <- merge(baseline_var, olink_soma_non_normal, by="csid") 

# load plate id and merge

olink_plate <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00190-V1/data_baseline_olink_plates.csv")

olink_plate$batch <- 1

olink_plate[grep("II",olink_plate$panel_full),]$batch <- 2

olink_plate_1 <- olink_plate[olink_plate$batch==1,]
                             
olink_plate_2 <- olink_plate[olink_plate$batch==2,]

olink_plate_1 <- pivot_wider(olink_plate_1, names_from = panel_full, values_from = plateid)

olink_plate_2 <- pivot_wider(olink_plate_2, names_from = panel_full, values_from = plateid)

olink_plate <- merge(olink_plate_1,olink_plate_2,by=c("csid","index"))

olink_plate <- olink_plate[,-which(names(olink_plate) %in% c("index","batch.x","batch.y"))]

colnames(olink_plate)[2:9] <- paste("olink_plt", colnames(olink_plate)[2:9], sep = "_")

somascan_plate <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00190-V1/data_baseline_somalogic_samples.csv")[c("csid","plateid")]

names(somascan_plate)[2] <- "somascan_plt_id"

olink_soma_normal <- merge(olink_plate, olink_soma_normal, by="csid")

olink_soma_normal <- merge(somascan_plate, olink_soma_normal, by="csid")

olink_soma_non_normal <- merge(olink_plate, olink_soma_non_normal, by="csid")

olink_soma_non_normal <- merge(somascan_plate, olink_soma_non_normal, by="csid")

## plot a Venn diagram

display.brewer.all()

cols <- brewer.pal(n = 3, name = "Set1")[c(1,2)]

length(unique(overlap$uniprot_id))
length(unique(c(somascan_protein$uniprot_id_1,somascan_protein$uniprot_id_2,somascan_protein$uniprot_id_3)))

venn.diagram(
  list(Olink = 1:2921, SomaScan = 755:7170),
  filename = 'K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/venn_diagram.png',
  fontfamily = "sans",	cat.fontfamily = "sans",
  fill=cols)

## run correlation between olink and somascan

# unadjusted

for (i in 1:nrow(overlap)) {
  overlap$rho_olink_somascan_normal[i] <- cor.test(olink_soma_normal[[overlap$olink_id[i]]], olink_soma_normal[[overlap$somascan_id[i]]], method = 'spearman')$estimate
  overlap$rho_olink_somascan_non_normal[i] <- cor.test(olink_soma_non_normal[[overlap$olink_id[i]]], olink_soma_non_normal[[overlap$somascan_id[i]]], method = 'spearman')$estimate
}

# adjusted

olink_soma_normal_adj <- data_adjust(olink_soma_normal,effect="somascan_plt_id",select="seq")
olink_soma_normal_adj <- data_adjust(olink_soma_normal,effect=names(olink_soma_normal)[3:10],select="ol_")

olink_soma_non_normal_adj <- data_adjust(olink_soma_non_normal,effect="somascan_plt_id",select="seq")
olink_soma_non_normal_adj <- data_adjust(olink_soma_non_normal,effect=names(olink_soma_non_normal)[3:10],select="ol_")

for (i in 1:nrow(overlap)) {
  overlap$rho_olink_somascan_normal_adj[i] <- cor.test(olink_soma_normal_adj[[overlap$olink_id[i]]], olink_soma_normal_adj[[overlap$somascan_id[i]]], method = 'spearman')$estimate
  overlap$rho_olink_somascan_non_normal_adj[i] <- cor.test(olink_soma_non_normal_adj[[overlap$olink_id[i]]], olink_soma_non_normal_adj[[overlap$somascan_id[i]]], method = 'spearman')$estimate
}

# normalised vs non_normalised

for (i in 1:nrow(overlap)) {
  overlap$rho_somascan_normal_non_normal[i] <- cor.test(olink_soma_normal[[overlap$somascan_id[i]]], olink_soma_non_normal[[overlap$somascan_id[i]]], method = 'spearman')$estimate
}

write.csv(overlap,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/overlap.csv", quote=F, row.names=F)

# save

write.csv(olink_soma_normal,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/olink_soma_normal.csv", quote=F, row.names=F)

write.csv(olink_soma_non_normal,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/olink_soma_non_normal.csv", quote=F, row.names=F)

write.csv(overlap,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/overlap.csv", quote=F, row.names=F)

################################################################################ shortcut

olink_soma_normal <- read.csv("olink_soma_normal.csv")

olink_soma_non_normal <- read.csv("olink_soma_non_normal.csv")

overlap <- read.csv("overlap.csv")

## get summary and save plots
  
summary(overlap$rho_olink_somascan_normal)

summary(overlap$rho_olink_somascan_non_normal)

summary(overlap$rho_somascan_normal_non_normal)

length(unique(overlap$uniprot_id))

length(unique(overlap$olink_id))

length(unique(overlap$somascan_id))

## summary table

table_1 <- overlap %>% 
  select(rho_olink_somascan_normal,rho_olink_somascan_non_normal,batch) %>%
  tbl_summary(by = batch, 
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ 
                c("{median} ({p25}, {p75})","{min}, {max}")) 
# %>%
#   modify_header(label = "**Variable**",stat_1 = "**Batch 1, N = 1,619**",stat_2 = "**Batch 2, N = 1,064**")

table_1

table_2 <- overlap %>% 
  select(rho_olink_somascan_normal,rho_olink_somascan_non_normal,dilution) %>%
  tbl_summary(by = dilution, 
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ 
                c("{median} ({p25}, {p75})","{min}, {max}")) 
# %>%
#   modify_header(label = "**Variable**",stat_1 = "**0.005%, N = 128**",stat_2 = "**0.5%, N = 555**",stat_3 = "**20%, N = 2,000**")

table_2

table_3 <- overlap %>% 
  select(rho_olink_somascan_normal,rho_olink_somascan_non_normal,panel) %>%
  tbl_summary(by = panel, 
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ 
                c("{median} ({p25}, {p75})","{min}, {max}")) 
table_3

# check number of proteins

length(unique(overlap$uniprot_id))

protein_count <- as.data.frame(table(overlap$uniprot_id))

# check the number of proteins targeted by multiple aptamers

table(protein_count$Freq)

overlap_dup <- overlap[duplicated(overlap$uniprot_id) | duplicated(overlap$uniprot_id, fromLast=TRUE),]

length(unique(overlap_dup$uniprot_id))

length(unique(overlap_dup$somascan_id))

# check correlation between aptamers

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# normalised data

somascan_cor_aptamer_normal <- as.data.frame(matrix(ncol = 5, nrow = 0))

names(somascan_cor_aptamer_normal) <-  c("row","column","cor","p","uniprot_id")

for (i in 1:length(unique(overlap_dup$uniprot_id))) {
  protein_multi_amptamer <- unique(overlap_dup$uniprot_id)[i]
  aptamer_multi_amptamer <- overlap_dup$somascan_id[overlap_dup$uniprot_id == unique(overlap_dup$uniprot_id)[i]]
  somascan_multi_aptamer <- olink_soma_normal[,aptamer_multi_amptamer]
  cor_aptamer_matrix <- rcorr(as.matrix(somascan_multi_aptamer),type="spearman")
  cor_aptamer <- flattenCorrMatrix(cor_aptamer_matrix$r,cor_aptamer_matrix$P)
  cor_aptamer$uniprot_id <- unique(overlap_dup$uniprot_id)[i]
  somascan_cor_aptamer_normal <- rbind(somascan_cor_aptamer_normal,cor_aptamer)
}

summary(somascan_cor_aptamer_normal$cor)

somascan_aptamer_normal <- ggplot(somascan_cor_aptamer_normal, aes(x=cor)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=median(cor)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,100) +
  ggtitle("SomaScan aptamers targeting the same protein (normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  theme_few()

# non-normalised

somascan_cor_aptamer_non_normal <- as.data.frame(matrix(ncol = 5, nrow = 0))

names(somascan_cor_aptamer_non_normal) <-  c("row","column","cor","p","uniprot_id")

for (i in 1:length(unique(overlap_dup$uniprot_id))) {
  protein_multi_amptamer <- unique(overlap_dup$uniprot_id)[i]
  aptamer_multi_amptamer <- overlap_dup$somascan_id[overlap_dup$uniprot_id == unique(overlap_dup$uniprot_id)[i]]
  somascan_multi_aptamer <- olink_soma_non_normal[,aptamer_multi_amptamer]
  cor_aptamer_matrix <- rcorr(as.matrix(somascan_multi_aptamer),type="spearman")
  cor_aptamer <- flattenCorrMatrix(cor_aptamer_matrix$r,cor_aptamer_matrix$P)
  cor_aptamer$uniprot_id <- unique(overlap_dup$uniprot_id)[i]
  somascan_cor_aptamer_non_normal <- rbind(somascan_cor_aptamer_non_normal,cor_aptamer)
}

summary(somascan_cor_aptamer_non_normal$cor)

somascan_aptamer_non_normal <- ggplot(somascan_cor_aptamer_non_normal, aes(x=cor)) + 
  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=median(cor)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,100) +
  ggtitle("SomaScan aptamers targeting the same protein (non-normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  theme_few()
somascan_aptamer <- ggarrange(somascan_aptamer_normal,somascan_aptamer_non_normal,labels=c("A","B"),ncol=1,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/somascan_aptamer.png",somascan_aptamer,width=9,height=9)

# plot

olink_somascan_normal <- ggplot(overlap, aes(x=rho_olink_somascan_normal)) + 
  geom_histogram(color="black", fill="white") +
  xlim(-0.5,1) +
  ylim(0,550) +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan (normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  theme_few()

olink_somascan_non_normal <- ggplot(overlap, aes(x=rho_olink_somascan_non_normal)) + 
  geom_histogram(color="black", fill="white") +
  xlim(-0.5,1) +
  ylim(0,550) +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan (non-normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  theme_few()

olink_somascan <- ggarrange(olink_somascan_normal,olink_somascan_non_normal,labels=c("A","B"),ncol=1,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan.png",olink_somascan,width=9,height=9)

somascan_normal_non_normal <- ggplot(overlap, aes(x=rho_somascan_normal_non_normal)) + 
  geom_histogram(color="black", fill="white") +
  xlim(0,1) +
  geom_vline(aes(xintercept=median(rho_somascan_normal_non_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("SomaScan normalised vs non-normalised") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  theme_few()

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/somascan_normal_non_normal.png",somascan_normal_non_normal,width=9,height=6)

# # adjusted
# 
summary(overlap$rho_olink_somascan_normal_adj)

olink_somascan_normal_adj <- ggplot(overlap, aes(x=rho_olink_somascan_normal_adj)) + 
  geom_histogram(color="black", fill="white") +
  xlim(-0.5,1) +
  ylim(0,550) +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal_adj)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan (normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  theme_few()

summary(overlap$rho_olink_somascan_non_normal_adj)

olink_somascan_non_normal_adj <- ggplot(overlap, aes(x=rho_olink_somascan_non_normal_adj)) + 
  geom_histogram(color="black", fill="white") +
  xlim(-0.5,1) +
  ylim(0,550) +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal_adj)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan (non-normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  theme_few()

olink_somascan_adj <- ggarrange(olink_somascan_normal_adj,olink_somascan_non_normal_adj,labels=c("A","B"),ncol=1,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_adj.png",olink_somascan_adj,width=9.6,height=10.8)


## plot separate

batch <- sort(unique(overlap$batch))
dilution <- sort(unique(overlap$dilution))
dilution_label <- c("0.005%","0.5%","20%")

# normalised

olink_somascan_separate_normal_list <- list()

i <- 1

for (x in 1:2) {
  for (y in 1:3){
    olink_somascan_separate_normal_list[[i]] <- ggplot(overlap[overlap$batch==batch[x] & overlap$dilution==dilution[y],], aes(x=rho_olink_somascan_normal)) + 
      geom_histogram(bins=20,color="black", fill="white") +
      geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
      xlim(-0.5,1) +
      ylim(0,350) +
      ggtitle(paste0("Olink (batch ",batch[x],") vs SomaScan (",dilution_label[y], ", normalised)"))  + 
      xlab("Spearman correlation coefficient") + 
      ylab("Frequency") +
      labs(fill="SomaScan dilution") +
      theme_few()
    i <- i+1
  }
}

olink_somascan_separate_normal <- ggarrange(olink_somascan_separate_normal_list[[1]],olink_somascan_separate_normal_list[[2]],olink_somascan_separate_normal_list[[3]],olink_somascan_separate_normal_list[[4]],olink_somascan_separate_normal_list[[5]],olink_somascan_separate_normal_list[[6]],labels=c("A","B","C","D","E","F"),ncol=3,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_separate_normal.png",olink_somascan_separate_normal,width=16,height=9)

# non-normalised

olink_somascan_separate_non_normal_list <- list()

i <- 1

for (x in 1:2) {
  for (y in 1:3){
    olink_somascan_separate_non_normal_list[[i]] <- ggplot(overlap[overlap$batch==batch[x] & overlap$dilution==dilution[y],], aes(x=rho_olink_somascan_non_normal)) + 
      geom_histogram(bins=20,color="black", fill="white") +
      geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
      xlim(-0.5,1) +
      ylim(0,350) +
      ggtitle(paste0("Olink (batch ",batch[x],") vs SomaScan (",dilution_label[y], ", non-normalised)"))  + 
      xlab("Spearman correlation coefficient") + 
      ylab("Frequency") +
      labs(fill="SomaScan dilution") +
      theme_few()
    i <- i+1
  }
}

olink_somascan_separate_non_normal <- ggarrange(olink_somascan_separate_non_normal_list[[1]],olink_somascan_separate_non_normal_list[[2]],olink_somascan_separate_non_normal_list[[3]],olink_somascan_separate_non_normal_list[[4]],olink_somascan_separate_non_normal_list[[5]],olink_somascan_separate_non_normal_list[[6]],labels=c("A","B","C","D","E","F"),ncol=3,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_separate_non_normal.png",olink_somascan_separate_non_normal,width=16,height=9)

## plot (separate by batch)

# cols <- brewer.pal(n = 3, name = "Set2")

olink_1_soma_normal <- ggplot(overlap[overlap$batch==1,], aes(x=rho_olink_somascan_normal,
                                                              fill=factor(dilution,levels=c("0.2","0.005","5e-05"),labels=c("20%","0.5%","0.005%")))) + 
  geom_histogram(colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,320) +
  ggtitle("Olink (batch 1) vs SomaScan (normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="SomaScan dilution") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_2_soma_normal <- ggplot(overlap[overlap$batch==2,], aes(x=rho_olink_somascan_normal,
                                                              fill=factor(dilution,levels=c("0.2","0.005","5e-05"),labels=c("20%","0.5%","0.005%")))) + 
  geom_histogram(colour="black") +
  xlim(-0.5,1) +
  ylim(0,320) +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink (batch 2) vs SomaScan (normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="SomaScan dilution") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_1_soma_non_normal <- ggplot(overlap[overlap$batch==1,], aes(x=rho_olink_somascan_non_normal,
                                                                  fill=factor(dilution,levels=c("0.2","0.005","5e-05"),labels=c("20%","0.5%","0.005%")))) + 
  geom_histogram(colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,320) +
  ggtitle("Olink (batch 1) vs SomaScan (non-normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="SomaScan dilution") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_2_soma_non_normal <- ggplot(overlap[overlap$batch==2,], aes(x=rho_olink_somascan_non_normal,
                                                                  fill=factor(dilution,levels=c("0.2","0.005","5e-05"),labels=c("20%","0.5%","0.005%")))) + 
  geom_histogram(colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,320) +
  ggtitle("Olink (batch 2) vs SomaScan (non-normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="SomaScan dilution") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_somascan_separate_by_batch <- ggarrange(olink_1_soma_normal,olink_2_soma_normal,olink_1_soma_non_normal,olink_2_soma_non_normal,labels=c("A","B","C","D"),ncol=2,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_separate_by_batch.png",olink_somascan_separate_by_batch,width=19.2,height=10.8)

## plot (separate by dilution)

# cols <- brewer.pal(n = 3, name = "Pastel2")[c(1,2)]

olink_soma_5e05_normal <- ggplot(overlap[overlap$dilution==5e-05,], aes(x=rho_olink_somascan_normal,
                                                                        fill=factor(batch,levels=c("2","1")))) + 
  geom_histogram(bins=20,colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,620) +
  ggtitle("Olink vs SomaScan (0.005%, normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Olink Batch") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_soma_0.005_normal <- ggplot(overlap[overlap$dilution==0.005,], aes(x=rho_olink_somascan_normal,
                                                                         fill=factor(batch,levels=c("2","1")))) + 
  geom_histogram(bins=20,colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,620) +
  ggtitle("Olink vs SomaScan (0.5%, normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Olink Batch") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_soma_0.2_normal <- ggplot(overlap[overlap$dilution==0.2,], aes(x=rho_olink_somascan_normal,
                                                                     fill=factor(batch,levels=c("2","1")))) + 
  geom_histogram(bins=20,colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,620) +
  ggtitle("Olink vs SomaScan (20%, normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Olink Batch") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_soma_5e05_non_normal <- ggplot(overlap[overlap$dilution==5e-05,], aes(x=rho_olink_somascan_non_normal,
                                                                            fill=factor(batch,levels=c("2","1")))) + 
  geom_histogram(bins=20,colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,620) +
  ggtitle("Olink vs SomaScan (0.005%, non-normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Olink Batch") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_soma_0.005_non_normal <- ggplot(overlap[overlap$dilution==0.005,], aes(x=rho_olink_somascan_non_normal,
                                                                             fill=factor(batch,levels=c("2","1")))) + 
  geom_histogram(bins=20,colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,620) +
  ggtitle("Olink vs SomaScan (0.5%, non-normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Olink Batch") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_soma_0.2_non_normal <- ggplot(overlap[overlap$dilution==0.2,], aes(x=rho_olink_somascan_non_normal,
                                                                         fill=factor(batch,levels=c("2","1")))) + 
  geom_histogram(bins=20,colour="black") +
  geom_vline(aes(xintercept=median(rho_olink_somascan_non_normal)),color="black", linetype="dashed", linewidth=1) +
  xlim(-0.5,1) +
  ylim(0,620) +
  ggtitle("Olink vs SomaScan (20%, non-normalised)") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Olink Batch") +
  # scale_fill_manual(values=cols) +
  theme_few()

olink_somascan_separate_by_dilution <- ggarrange(olink_soma_5e05_normal,olink_soma_0.005_normal,olink_soma_0.2_normal,olink_soma_5e05_non_normal,olink_soma_0.005_non_normal,olink_soma_0.2_non_normal,labels=c("A","B","C","D","E","F"),ncol=3,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_separate_by_dilution.png",olink_somascan_separate_by_dilution,width=19.2,height=10.8)

## plot (combined by batch)

# get colour palette

# cols <- brewer.pal(n = 3, name = "Pastel1")

# normalised

olink_soma_normal_batch <- ggplot(overlap) + 
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
      ylim(0,520) +
      # scale_fill_manual(values=cols) +
      scale_pattern_manual(values=c('stripe',"none")) +
      guides(pattern = guide_legend(override.aes = list(fill = "white")),
             fill = guide_legend(override.aes = list(pattern = c("none", "none", "none")))) +
      theme_few() +
      ggtitle("Olink vs SomaScan (normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") + labs(pattern ='Olink batch',fill="SomaScan dilution")

# non-normalised

olink_soma_non_normal_batch <- ggplot(overlap) + 
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
  ylim(0,520) +
  # scale_fill_manual(values=cols) +
  scale_pattern_manual(values=c('stripe',"none")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = c("none", "none", "none")))) +
  theme_few() +
  ggtitle("Olink vs SomaScan (non-normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") + labs(pattern ='Olink batch',fill="SomaScan dilution")

olink_somascan_combined_batch <- ggarrange(olink_soma_normal_batch,olink_soma_non_normal_batch,labels=c("A","B"),ncol=1,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_combined_batch.png",olink_somascan_combined_batch,width=9.6,height=10.8)

## plot (combined by dilution)

# cols <- brewer.pal(n = 3, name = "Pastel1")

# normalised

olink_soma_normal_dilution <- ggplot(overlap) + 
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
  ylim(0,550) +
  # scale_fill_manual(values=cols) +
  scale_pattern_manual(values=c('stripe',"none")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = c("none", "none", "none")))) +
  theme_few() +
  ggtitle("Olink vs SomaScan (normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") + labs(pattern ='Olink batch',fill="SomaScan dilution")

# non-normalised

olink_soma_non_normal_dilution <- ggplot(overlap) + 
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
  ylim(0,550) +
  # scale_fill_manual(values=cols) +
  scale_pattern_manual(values=c('stripe',"none")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = c("none", "none", "none")))) +
  theme_few() +
  ggtitle("Olink vs SomaScan (non-normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") + labs(pattern ='Olink batch',fill="SomaScan dilution")

olink_somascan_combined_dilution <- ggarrange(olink_soma_normal_dilution,olink_soma_non_normal_dilution,labels=c("A","B"),ncol=1,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_combined_dilution.png",olink_somascan_combined_dilution,width=9.6,height=10.8)
