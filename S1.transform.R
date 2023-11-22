rm(list = ls())

setwd("")

library(readxl)
library(tidyverse)
library(datawizard)
library(RNOmni)

## load olink ckb data

olink_all <- read.csv("")

# olink_test <- pivot_wider(olink_all[,c(1,3,5)], id_cols = csid, names_from = assay, values_from = npx)

# check duplicates

olink_protein_dup <- unique(olink_all$assay[duplicated(olink_all[c("csid","assay")])])

olink_dup <- olink_all[which(olink_all$assay %in% olink_protein_dup),]

olink_dup <- pivot_wider(olink_dup[,c(1,3,5)], id_cols = csid, names_from = assay, values_from = npx)

# olink_all %>%
#   dplyr::group_by(csid, assay) %>%
#   dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
#   dplyr::filter(n > 1L) 

# unique olink proteins

olink_unique <- olink_all[!duplicated(olink_all[c("csid","assay")]),]

olink <- pivot_wider(olink_unique[,c(1,3,5)], id_cols = csid, names_from = assay, values_from = npx)

names(olink)

write.csv(olink,"olink.csv", quote=F, row.names=F)

## load somascan data

# load normalised somascan data

somascan_normalised <- read.csv("")

# change to wide format

somascan_normalised <- pivot_wider(somascan_normalised, names_from = variable, values_from = value)

somascan_normalised <- somascan_normalised[somascan_normalised$qc==0,-c(2,3)]

## load non-normalised somascan data 

somascan_non_normalised <- read.csv("")

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

