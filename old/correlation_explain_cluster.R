rm(list = ls())

setwd("/gpfs3/well/ckb/users/dma206/proteomics/data")

library(Boruta)

####### run Random-forest

overlap_annot <- readRDS("overlap_annot.RDS")

overlap_annot_normal <- overlap_annot[,c(9,3,4,8,14,15,17:length(overlap_annot))]

overlap_annot_non_normal <- overlap_annot[,c(10,3,4,8,14,16,17:length(overlap_annot))]

## run machine learning

# normalised

set.seed(47)

boruta_normal <- Boruta(rho_olink_somascan_normal~., data = overlap_annot_normal, doTrace = 2, maxRuns = 5000)

saveRDS(boruta_normal,"boruta_normal.RDS")

# non-normalised

set.seed(47)

boruta_non_normal <- Boruta(rho_olink_somascan_non_normal~., data = overlap_annot_non_normal, doTrace = 2, maxRuns = 5000)

saveRDS(boruta_non_normal,"boruta_non_normal.RDS")
