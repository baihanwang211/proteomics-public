rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(hausekeep)
library(stringr)
library(Boruta)
library(ggplot2)
library(reshape)
library(ggpubr)
library(RColorBrewer)

# normalised

boruta_normal <- readRDS("boruta_normal_deduplicated.RDS")

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

# change cases and special characters

boruta_normal_melt$variable <- gsub("`", "", boruta_normal_melt$variable)

boruta_normal_melt$variable <- gsub("\\.", " ", boruta_normal_melt$variable)

boruta_normal_melt$variable <- gsub("outlier_olink", "%outlier Olink", boruta_normal_melt$variable)

boruta_normal_melt$variable <- gsub("outlier_soma_normal", "%outlier SomaScan", boruta_normal_melt$variable)

boruta_normal_melt$variable <- gsub("n_isoform", "Number of isoforms", boruta_normal_melt$variable)

boruta_normal_melt$variable <- gsub("dilution", "SomaScan dilution", boruta_normal_melt$variable)

boruta_normal_melt$variable <- gsub("batch", "Olink batch", boruta_normal_melt$variable)

boruta_normal_melt$variable <- gsub("panel", "Olink panel", boruta_normal_melt$variable)

boruta_normal_melt$variable <- paste(toupper(substr(boruta_normal_melt$variable, 1, 1)), substr(boruta_normal_melt$variable, 2, nchar(boruta_normal_melt$variable)), sep="")

unique(boruta_normal_melt$variable)

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

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/boruta_normal_plot_deduplicated.png",boruta_normal_plot,width=15,height=10)

# non_normalised

boruta_non_normal <- readRDS("boruta_non_normal_deduplicated.RDS")

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

# change cases and special characters

boruta_non_normal_melt$variable <- gsub("`", "", boruta_non_normal_melt$variable)

boruta_non_normal_melt$variable <- gsub("\\.", " ", boruta_non_normal_melt$variable)

boruta_non_normal_melt$variable <- gsub("outlier_olink", "%outlier Olink", boruta_non_normal_melt$variable)

boruta_non_normal_melt$variable <- gsub("outlier_soma_non_normal", "%outlier SomaScan", boruta_non_normal_melt$variable)

boruta_non_normal_melt$variable <- gsub("n_isoform", "Number of isoforms", boruta_non_normal_melt$variable)

boruta_non_normal_melt$variable <- gsub("dilution", "SomaScan dilution", boruta_non_normal_melt$variable)

boruta_non_normal_melt$variable <- gsub("batch", "Olink batch", boruta_non_normal_melt$variable)

boruta_non_normal_melt$variable <- gsub("panel", "Olink panel", boruta_non_normal_melt$variable)

boruta_non_normal_melt$variable <- paste(toupper(substr(boruta_non_normal_melt$variable, 1, 1)), substr(boruta_non_normal_melt$variable, 2, nchar(boruta_non_normal_melt$variable)), sep="")

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

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/boruta_non_normal_plot_deduplicated.png",boruta_non_normal_plot,width=15,height=10)
