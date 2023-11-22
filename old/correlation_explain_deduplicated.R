rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(UniProt.ws)
library(ckbplotr)
library(hausekeep)
library(stringr)
library(Boruta)
library(ggplot2)
library(reshape)
library(ggpubr)
library(RColorBrewer)

####### run Random-forest

overlap_annot <- readRDS("overlap_annot.RDS")

names(overlap_annot)

overlap_annot_deduplicated_normal <-  overlap_annot[order(overlap_annot$rho_olink_somascan_normal,decreasing=TRUE),]

overlap_annot_deduplicated_normal <- overlap_annot_deduplicated_normal[!duplicated(overlap_annot_deduplicated_normal$uniprot_id),]

overlap_annot_deduplicated_normal <- overlap_annot_deduplicated_normal[,c(8,3,7,13,14,16:length(overlap_annot_deduplicated_normal))]

overlap_annot_deduplicated_non_normal <-  overlap_annot[order(overlap_annot$rho_olink_somascan_non_normal,decreasing=TRUE),]

overlap_annot_deduplicated_non_normal <- overlap_annot_deduplicated_non_normal[!duplicated(overlap_annot_deduplicated_non_normal$uniprot_id),]

overlap_annot_deduplicated_non_normal <- overlap_annot_deduplicated_non_normal[,c(9,3,7,13,15,16:length(overlap_annot))]

## run machine learning

# normalised

set.seed(47)

boruta_normal <- Boruta(rho_olink_somascan_normal~., data = overlap_annot_deduplicated_normal, doTrace = 2, maxRuns = 1000)

saveRDS(boruta_normal,"boruta_normal_deduplicated.RDS")

print(boruta_normal)

boruta_normal_df <- attStats(boruta_normal)

print(boruta_normal_df)

# reshape data

boruta_normal_melt <- melt(as.data.frame(boruta_normal$ImpHistory))

boruta_normal_melt <- boruta_normal_melt[is.finite(boruta_normal_melt$value),]

# remove shadow variables

boruta_normal_melt <- boruta_normal_melt[-grep("shadow",boruta_normal_melt$variable),]

names(boruta_normal_melt) <- c("variable","importance")

# color by decision

boruta_normal_decision <- data.frame(variable=row.names(boruta_normal_df),decision=as.vector(boruta_normal_df$decision))

boruta_normal_melt <- merge(boruta_normal_melt,boruta_normal_decision,by="variable")

boruta_normal_melt$decision <- factor(boruta_normal_melt$decision,levels=c("Confirmed","Tentative","Rejected"))

# plot

plot(boruta_normal)

boruta_normal_plot <- ggplot(boruta_normal_melt, aes(x=reorder(variable, importance, FUN = median), y=importance)) + 
  geom_boxplot(aes(fill=decision)) + 
  coord_flip() +
  xlab("Variable") +
  ylab("Importance") +
  labs(fill="Decision") +
  theme_bw() +
  theme(text = element_text(size = 17))

boruta_normal_plot

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/boruta_normal_plot_deduplicated.png",boruta_normal_plot,width=19.2,height=10.8)

# non-normalised

set.seed(47)

boruta_non_normal <- Boruta(rho_olink_somascan_non_normal~., data = overlap_annot_deduplicated_non_normal, doTrace = 2, maxRuns = 1000)

saveRDS(boruta_non_normal,"boruta_non_normal_deduplicated.RDS")

print(boruta_non_normal)

boruta_non_normal_df <- attStats(boruta_non_normal)

print(boruta_non_normal_df)

# reshape data

boruta_non_normal_melt <- melt(as.data.frame(boruta_non_normal$ImpHistory))

boruta_non_normal_melt <- boruta_non_normal_melt[is.finite(boruta_non_normal_melt$value),]

# remove shadow variables

boruta_non_normal_melt <- boruta_non_normal_melt[-grep("shadow",boruta_non_normal_melt$variable),]

names(boruta_non_normal_melt) <- c("variable","importance")

# color by decision

boruta_non_normal_decision <- data.frame(variable=row.names(boruta_non_normal_df),decision=as.vector(boruta_non_normal_df$decision))

boruta_non_normal_melt <- merge(boruta_non_normal_melt,boruta_non_normal_decision,by="variable")

boruta_non_normal_melt$decision <- factor(boruta_non_normal_melt$decision,levels=c("Confirmed","Tentative","Rejected"))

# plot

plot(boruta_non_normal)

boruta_non_normal_plot <- ggplot(boruta_non_normal_melt, aes(x=reorder(variable, importance, FUN = median), y=importance)) + 
  geom_boxplot(aes(fill=decision)) + 
  coord_flip() +
  xlab("Variable") +
  ylab("Importance") +
  labs(fill="Decision") +
  theme_bw() +
  theme(text = element_text(size = 17))

boruta_non_normal_plot

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/boruta_non_normal_plot_deduplicated.png",boruta_non_normal_plot,width=19.2,height=10.8)

boruta_plot <- ggarrange(boruta_normal_plot,boruta_non_normal_plot,labels=c("A","B"),ncol=2,nrow=1)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/boruta_plot_deduplicated.png",boruta_plot,width=19.2*1.5,height=10.8*1.5)
