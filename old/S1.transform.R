rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(readxl)
library(tidyverse)
library(datawizard)
library(RNOmni)

# firstly load linkage file for uniprot id

olink_protein <- read_excel("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/olink_uniprot.xlsx")[c(1,3,4)]

names(olink_protein) <- c("uniprot_id","olink_id","panel")

olink_protein <- olink_protein[!duplicated(olink_protein$olink_id),]

olink_protein$olink_id <- tolower(olink_protein$olink_id)

olink_protein$olink_id <- paste0("ol_", olink_protein$olink_id)

olink_protein$olink_id <- gsub("-","_",olink_protein$olink_id)

olink_protein$batch <- 1

olink_protein[grep("II",olink_protein$panel),]$batch <-2

olink_protein$panel <- gsub("_II","",olink_protein$panel)

## load olink ckb data

olink_all <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/data_baseline_olink_explore.csv")

# check duplicates

olink_all %>%
  dplyr::group_by(csid, assay) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 

olink_unique <- olink_all[!duplicated(olink_all[c("csid","assay")]),]

olink <- pivot_wider(olink_unique[,c(1,3,5)], id_cols = csid, names_from = assay, values_from = npx)

names(olink)

write.csv(olink,"olink.csv", quote=F, row.names=F)


olink_cardiometabolic <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/data_baseline_olink_cardiometabolic.csv")

olink_inflammation <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/data_baseline_olink_inflammation.csv")

olink_neurology <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/data_baseline_olink_neurology.csv")

olink_oncology <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/data_baseline_olink_oncology.csv")

olink_all <- olink_cardiometabolic %>% inner_join(olink_inflammation, by='csid', suffix=c("","_dup")) %>% inner_join(olink_neurology, by='csid', suffix=c("","_dup")) %>% inner_join(olink_oncology, by='csid', suffix=c("","_dup"))

# find difference between proteins in the assay and ckb data

setdiff(olink_protein$olink_id,names(olink_all)[-1])

setdiff(names(olink_all)[-1],olink_protein$olink_id)

# find unique proteins measured in ckb

olink_unique <- olink_all[,-grep("dup",names(olink_all))]

# only keep those proteins that are also included in assay list

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
  olink_dup_list[[i]]$cardiometabolic <- olink_cardiometabolic[,which(names(olink_cardiometabolic)==olink_protein_dup[i])]
  olink_dup_list[[i]]$inflammation <- olink_inflammation[,which(names(olink_inflammation)==olink_protein_dup[i])]
  olink_dup_list[[i]]$neurology <- olink_neurology[,which(names(olink_neurology)==olink_protein_dup[i])]
  olink_dup_list[[i]]$oncology <- olink_oncology[,which(names(olink_oncology)==olink_protein_dup[i])]
  olink_dup_list[[i]]$protein <- olink_protein_dup[i]
}

cor.test(olink_dup_list[[4]]$cardiometabolic,olink_dup_list[[4]]$inflammation)

# they are actually all the same

# save

write.csv(olink,"olink.csv", quote=F, row.names=F)


## load somascan data

# load normalised somascan data

somascan_normalised <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/data_baseline_somalogic.csv")

# change to wide format

somascan_normalised <- pivot_wider(somascan_normalised, names_from = variable, values_from = value)

somascan_normalised <- somascan_normalised[somascan_normalised$qc==0,-c(2,3)]

## load non-normalised somascan data 

somascan_non_normalised <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/data_baseline_somalogic_qc.csv")

# change to wide format

somascan_non_normalised <- pivot_wider(somascan_non_normalised, names_from = variable, values_from = value)

somascan_non_normalised <- somascan_non_normalised[somascan_non_normalised$qc==0,-c(2,3)]

## log

somascan_normalised_log <- somascan_normalised
somascan_normalised_log[2:length(somascan_normalised_log)] <- lapply(somascan_normalised_log[2:length(somascan_normalised_log)],log)

somascan_non_normalised_log <- somascan_non_normalised
somascan_non_normalised_log[2:length(somascan_non_normalised_log)] <- lapply(somascan_non_normalised_log[2:length(somascan_non_normalised_log)],log)

# save

write.csv(somascan_normalised,"somascan_normalised.csv", quote=F, row.names=F)
write.csv(somascan_non_normalised,"somascan_non_normalised.csv", quote=F, row.names=F)

write.csv(somascan_normalised_log,"somascan_normalised_log.csv", quote=F, row.names=F)
write.csv(somascan_non_normalised_log,"somascan_non_normalised_log.csv", quote=F, row.names=F)

