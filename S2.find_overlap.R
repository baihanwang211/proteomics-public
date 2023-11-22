rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(readxl)
library(tidyverse)
library(VennDiagram)
library(RColorBrewer)

############################################################## load olink

# firstly load linkage file for uniprot id

olink_protein <- read_excel("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/olink_uniprot.xlsx")[c(1,3,4)]

names(olink_protein) <- c("uniprot_id","olink_id","panel")

olink_protein <- olink_protein[!duplicated(olink_protein$olink_id),]

olink_protein$olink_id <- gsub("-",".",olink_protein$olink_id)

# olink_protein$olink_id <- tolower(olink_protein$olink_id)
# 
# olink_protein$olink_id <- paste0("ol_", olink_protein$olink_id)
# 
# olink_protein$olink_id <- gsub("-","_",olink_protein$olink_id)

olink_protein$batch <- 1

olink_protein[grep("II",olink_protein$panel),]$batch <-2

## load olink ckb data

# long format to get lod

# olink_all <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/data_baseline_olink_explore.csv")
# 
# olink_lod <- olink_all[!duplicated(olink_all$csid),c("csid","lod")]
# 
# names(olink_lod)[2] <- "olink_lod"

# wide format

olink <- read.csv("olink.csv")

# update the linkage file

setdiff(olink_protein$olink_id,names(olink))
setdiff(names(olink),olink_protein$olink_id)

olink_protein <- olink_protein[which(olink_protein$olink_id %in% names(olink)),]

########################################################################################### load somascan data

# load linkage file first

somascan_protein <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/somalogic_meta.csv")[,c("aptname","uniprot","dilution2","organism")]

names(somascan_protein)[1:3] <- c("somascan_id","uniprot_id","dilution")

# split uniprot column into 3 because some aptamers have more than 1 uniprot id

somascan_protein <- separate(data = somascan_protein, col = uniprot_id, into = c("uniprot_id_1", "uniprot_id_2","uniprot_id_3"), sep = "\\|")

# restrict to human proteins (and exclude aptamers that don't target proteins)

somascan_protein_human <- somascan_protein[somascan_protein$organism=="Human" & somascan_protein$uniprot_id_1!="",]

somascan_protein_human_long <- pivot_longer(somascan_protein_human, cols=2:4, values_to = "uniprot_id")

length(unique(somascan_protein_human_long$uniprot_id))

n_aptamer <- data.frame(table(somascan_protein_human_long$uniprot_id))

n_aptamer_sum <- data.frame(table(n_aptamer$Freq))

# somascan_protein_1_to_1 <- somascan_protein_1_to_1[,-2]
# 
# somascan_protein_1_to_1 <- somascan_protein_1_to_1[!is.na(somascan_protein_1_to_1$uniprot_id),]
# 
# somascan_protein_1_to_1 <- somascan_protein_1_to_1[!(duplicated(somascan_protein_1_to_1$somascan_id) | duplicated(somascan_protein_1_to_1$somascan_id, fromLast = TRUE) |
#                                                      duplicated(somascan_protein_1_to_1$uniprot_id) | duplicated(somascan_protein_1_to_1$uniprot_id, fromLast = TRUE)), ]

## merge 

overlap_1 <- merge(olink_protein, somascan_protein_human, by.x="uniprot_id", by.y="uniprot_id_1")
names(overlap_1)[c(6,7)] <- c("somascan_uniprot_2","somascan_uniprot_3")

overlap_2 <- merge(olink_protein, somascan_protein_human, by.x="uniprot_id", by.y="uniprot_id_2")
names(overlap_2)[c(6,7)] <- c("somascan_uniprot_2","somascan_uniprot_3")

overlap_3 <- merge(olink_protein, somascan_protein_human, by.x="uniprot_id", by.y="uniprot_id_3")

overlap <- rbind(overlap_1,overlap_2)

## check numbers; there were 2923 proteins profiled by Olink in CKB

length(unique(overlap$uniprot_id))

# ncol(olink_all[,-grep("dup",names(olink_all))])-1

length(unique(c(somascan_protein_human$uniprot_id_1,somascan_protein_human$uniprot_id_2,somascan_protein_human$uniprot_id_3)))

## plot a Venn diagram

display.brewer.all()

cols <- brewer.pal(n = 3, name = "Set1")[c(1,2)]

length(unique(overlap$uniprot_id))
length(unique(c(somascan_protein$uniprot_id_1,somascan_protein$uniprot_id_2,somascan_protein$uniprot_id_3)))

venn.diagram(
  list(Olink = 1:2923, SomaScan = 756:7152),
  filename = 'K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/venn_diagram.png',
  fontfamily = "sans",	cat.fontfamily = "sans",
  category.names = c("Olink Explore 3072" , "SomaScan Assay v4.1"),
  height = 2500, 
  width = 2500, 
  cat.pos = c(-25, 30),
  cat.dist = c(0.04, 0.04),
  fill=cols)

venn.diagram(
  list(Olink = 1:2923, SomaScan = 756:7152),
  filename = 'K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/venn_diagram_big_font.png',
  fontfamily = "sans",	cat.fontfamily = "sans",
  category.names = c("" , ""),
  height = 1250, 
  width = 1250, 
  cat.pos = c(-25, 30),
  cat.dist = c(0.04, 0.04),
  fill=cols)


# save

table(overlap$organism) # check if they are human proteins (they all are)

overlap <- overlap[,-which(names(overlap) %in% "organism")]

write.csv(overlap,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/overlap.csv", quote=F, row.names=F)


# load 

# overlap <- read.csv("overlap.csv")

## find one to one proteins

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

# test

test1 <- overlap[-which(overlap$uniprot_id %in% protein_dup), ] # remove

test2 <- overlap[-which(overlap$somascan_id %in% aptamer_dup), ] # remove

# remove

overlap_1_to_1 <- overlap[-c(which(overlap$uniprot_id %in% protein_dup),which(overlap$somascan_id %in% aptamer_dup)), ]

# further checks

sum(duplicated(overlap_1_to_1$uniprot_id))

sum(duplicated(overlap_1_to_1$somascan_id))

write.csv(overlap_1_to_1,"overlap_1_to_1.csv", quote=F, row.names=F)

# overlap_1_to_1 <- read.csv("./overlap_1_to_1.csv")

## compare old and new data

# overlap_old <- read.csv("./old/overlap.csv")
# 
# setdiff(overlap$uniprot_id,overlap_old$uniprot_id)
# 
# names(olink)[grep("HLA",names(olink))]