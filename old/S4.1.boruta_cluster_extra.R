rm(list = ls())

setwd("/gpfs3/well/ckb/users/dma206/proteomics/data")

library(Boruta)

# overlap_1_to_1_annot_normal <- readRDS("overlap_1_to_1_annot_normal.RDS")
overlap_1_to_1_annot_non_normal <- readRDS("overlap_1_to_1_annot_non_normal.RDS")

# define levels of batch

# overlap_1_to_1_annot_normal$batch <- factor(overlap_1_to_1_annot_normal$batch, levels = c("1","2"))
overlap_1_to_1_annot_non_normal$batch <- factor(overlap_1_to_1_annot_non_normal$batch, levels = c("1","2"))

# normalised

# set.seed(47)
# 
# boruta_normal <- Boruta(rho_olink_soma_normal~., data = overlap_1_to_1_annot_normal, doTrace = 2, maxRuns = 5000000)
# 
# saveRDS(boruta_normal,"boruta_normal.RDS")

# non_normalised

set.seed(47)

boruta_non_normal <- Boruta(rho_olink_soma_non_normal~., data = overlap_1_to_1_annot_non_normal, doTrace = 2, maxRuns = 5000000)

saveRDS(boruta_non_normal,"boruta_non_normal.RDS")
