rm(list = ls())

setwd("/gpfs3/well/ckb/users/dma206/proteomics/data")

library(Boruta)

####### run Random-forest

overlap_annot <- readRDS("overlap_annot.RDS")

names(overlap_annot)

overlap_annot_deduplicated_normal <-  overlap_annot[order(overlap_annot$rho_olink_somascan_normal,decreasing=TRUE),]

overlap_annot_deduplicated_normal <- overlap_annot_deduplicated_normal[!duplicated(overlap_annot_deduplicated_normal$uniprot_id),]

overlap_annot_deduplicated_normal <- overlap_annot_deduplicated_normal[,c(9,3,4,8,14,15,17:length(overlap_annot_deduplicated_normal))]

overlap_annot_deduplicated_non_normal <-  overlap_annot[order(overlap_annot$rho_olink_somascan_non_normal,decreasing=TRUE),]

overlap_annot_deduplicated_non_normal <- overlap_annot_deduplicated_non_normal[!duplicated(overlap_annot_deduplicated_non_normal$uniprot_id),]

overlap_annot_deduplicated_non_normal <- overlap_annot_deduplicated_non_normal[,c(10,3,4,8,14,16,17:length(overlap_annot_deduplicated_non_normal))]

## run machine learning

# normalised

set.seed(47)

boruta_normal <- Boruta(rho_olink_somascan_normal~., data = overlap_annot_deduplicated_normal, doTrace = 2, maxRuns = 5000)

saveRDS(boruta_normal,"boruta_normal_deduplicated.RDS")

# non-normalised

set.seed(47)

boruta_non_normal <- Boruta(rho_olink_somascan_non_normal~., data = overlap_annot_deduplicated_non_normal, doTrace = 2, maxRuns = 5000)

saveRDS(boruta_non_normal,"boruta_non_normal_deduplicated.RDS")
