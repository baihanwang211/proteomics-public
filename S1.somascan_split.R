rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(readxl)
library(tidyverse)
library(datawizard)
library(RNOmni)

## load somascan data

# load normalised somascan data

somascan_normalised <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/data_baseline_somalogic.csv")

# change to wide format

somascan_normalised <- pivot_wider(somascan_normalised, names_from = variable, values_from = value)

dup <- somascan_normalised$csid[somascan_normalised$qc==1]

normalised_dup_1 <- somascan_normalised[somascan_normalised$csid %in% dup & somascan_normalised$qc==0, ]

normalised_dup_2 <- somascan_normalised[somascan_normalised$csid %in% dup & somascan_normalised$qc==1, ]

rho_normal <- data.frame(somascan_id=names(normalised_dup_1)[4:ncol(normalised_dup_1)],
                         rho=NA)

for (i in 4:ncol(normalised_dup_1)) {
  rho_normal$rho[i-3] <- cor.test(normalised_dup_1[[i]],normalised_dup_2[[i]],method="spearman")$estimate
}

median(rho_normal$rho)

## load non-normalised somascan data 

somascan_non_normalised <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/data_baseline_somalogic_qc.csv")

# change to wide format

somascan_non_normalised <- pivot_wider(somascan_non_normalised, names_from = variable, values_from = value)

dup <- somascan_non_normalised$csid[somascan_non_normalised$qc==1]

non_normalised_dup_1 <- somascan_non_normalised[somascan_non_normalised$csid %in% dup & somascan_non_normalised$qc==0, ]

non_normalised_dup_2 <- somascan_non_normalised[somascan_non_normalised$csid %in% dup & somascan_non_normalised$qc==1, ]

rho_non_normal <- data.frame(somascan_id=names(non_normalised_dup_1)[4:ncol(non_normalised_dup_1)],
                         rho=NA)

for (i in 4:ncol(non_normalised_dup_1)) {
  rho_non_normal$rho[i-3] <- cor.test(non_normalised_dup_1[[i]],non_normalised_dup_2[[i]],method="spearman")$estimate
}

median(rho_non_normal$rho)