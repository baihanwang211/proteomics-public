rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(tidyverse)
library(ggplot2)
library(ggpattern)
# library(ppcor)
library(gtsummary)
library(ggpubr)
library(Hmisc)
library(ggthemes)
library(datawizard)

######################## load raw, transformed and cleaned somascan data


overlap <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/overlap.csv")

olink_clean <- read.csv("olink_clean.csv")
olink_clean_adj_rint <- read.csv("olink_clean_adj_rint.csv")

somascan_normalised <- read.csv("somascan_normalised.csv")
somascan_normalised_log_clean <- read.csv("somascan_normalised_log_clean.csv")
somascan_normalised_1x_clean <- read.csv("somascan_normalised_1x_clean.csv")
somascan_normalised_log_clean_adj_rint <- read.csv("somascan_normalised_log_clean_adj_rint.csv")
somascan_normalised_1x_clean_adj_rint <- read.csv("somascan_normalised_1x_clean_adj_rint.csv")

somascan_non_normalised <- read.csv("somascan_non_normalised.csv")
somascan_non_normalised_log_clean <- read.csv("somascan_non_normalised_log_clean.csv")
somascan_non_normalised_1x_clean <- read.csv("somascan_non_normalised_1x_clean.csv")
somascan_non_normalised_log_clean_adj_rint <- read.csv("somascan_non_normalised_log_clean_adj_rint.csv")
somascan_non_normalised_1x_clean_adj_rint <- read.csv("somascan_non_normalised_1x_clean_adj_rint.csv")

################## some quick random checks on distribution

hist(olink_clean[[47]])
hist(olink_clean_adj_rint[[47]])

hist(somascan_normalised[[47]])
hist(somascan_normalised_log_clean[[47]])
hist(somascan_normalised_1x_clean[[47]])
hist(somascan_normalised_log_clean_adj_rint[[47]])
hist(somascan_normalised_1x_clean_adj_rint[[47]])

hist(somascan_non_normalised[[47]])
hist(somascan_non_normalised_log_clean[[47]])
hist(somascan_non_normalised_1x_clean[[47]])
hist(somascan_non_normalised_log_clean_adj_rint[[47]])
hist(somascan_non_normalised_1x_clean_adj_rint[[47]])

########## check correlation

overlap_cor <- overlap


## spearman

for (i in 1:nrow(overlap_cor)) {
  df_normal <- merge(olink_clean[,c("csid",overlap_cor$olink_id[i])],somascan_normalised[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$rho_olink_soma_normal[i] <- cor.test(df_normal[[2]], df_normal[[3]], method = 'spearman')$estimate
  
  df_non_normal <- merge(olink_clean[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$rho_olink_soma_non_normal[i] <- cor.test(df_non_normal[[2]], df_non_normal[[3]], method = 'spearman')$estimate
}



## normalised

# olink raw

for (i in 1:nrow(overlap_cor)) {
  df_normal <- merge(olink_clean[,c("csid",overlap_cor$olink_id[i])],somascan_normalised[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_normal[i] <- cor.test(df_normal[[2]], df_normal[[3]], method = 'pearson')$estimate
  
  df_normal <- merge(olink_clean[,c("csid",overlap_cor$olink_id[i])],somascan_normalised_log_clean[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_normal_log[i] <- cor.test(df_normal[[2]], df_normal[[3]], method = 'pearson')$estimate
  
  df_normal <- merge(olink_clean[,c("csid",overlap_cor$olink_id[i])],somascan_normalised_1x_clean[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_normal_1x_op[i] <- -cor.test(df_normal[[2]], df_normal[[3]], method = 'pearson')$estimate
  
  df_normal <- merge(olink_clean[,c("csid",overlap_cor$olink_id[i])],somascan_normalised_log_clean_adj_rint[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_normal_log_rint[i] <- cor.test(df_normal[[2]], df_normal[[3]], method = 'pearson')$estimate
  
  df_normal <- merge(olink_clean[,c("csid",overlap_cor$olink_id[i])],somascan_normalised_1x_clean_adj_rint[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_normal_1x_rint_op[i] <- -cor.test(df_normal[[2]], df_normal[[3]], method = 'pearson')$estimate
}

# olink rint

for (i in 1:nrow(overlap_cor)) {
  df_normal <- merge(olink_clean_adj_rint[,c("csid",overlap_cor$olink_id[i])],somascan_normalised[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_rint_soma_normal[i] <- cor.test(df_normal[[2]], df_normal[[3]], method = 'pearson')$estimate
  
  df_normal <- merge(olink_clean_adj_rint[,c("csid",overlap_cor$olink_id[i])],somascan_normalised_log_clean[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_rint_soma_normal_log[i] <- cor.test(df_normal[[2]], df_normal[[3]], method = 'pearson')$estimate
  
  df_normal <- merge(olink_clean_adj_rint[,c("csid",overlap_cor$olink_id[i])],somascan_normalised_1x_clean[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_rint_soma_normal_1x_op[i] <- -cor.test(df_normal[[2]], df_normal[[3]], method = 'pearson')$estimate
  
  df_normal <- merge(olink_clean_adj_rint[,c("csid",overlap_cor$olink_id[i])],somascan_normalised_log_clean_adj_rint[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_rint_soma_normal_log_rint[i] <- cor.test(df_normal[[2]], df_normal[[3]], method = 'pearson')$estimate
  
  df_normal <- merge(olink_clean_adj_rint[,c("csid",overlap_cor$olink_id[i])],somascan_normalised_1x_clean_adj_rint[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_rint_soma_normal_1x_rint_op[i] <- -cor.test(df_normal[[2]], df_normal[[3]], method = 'pearson')$estimate
}


## non-normalised

# olink raw

for (i in 1:nrow(overlap_cor)) {
  df_non_normal <- merge(olink_clean[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_non_normal[i] <- cor.test(df_non_normal[[2]], df_non_normal[[3]], method = 'pearson')$estimate
  
  df_non_normal <- merge(olink_clean[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised_log_clean[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_non_normal_log[i] <- cor.test(df_non_normal[[2]], df_non_normal[[3]], method = 'pearson')$estimate
  
  df_non_normal <- merge(olink_clean[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised_1x_clean[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_non_normal_1x_op[i] <- -cor.test(df_non_normal[[2]], df_non_normal[[3]], method = 'pearson')$estimate
  
  df_non_normal <- merge(olink_clean[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised_log_clean_adj_rint[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_non_normal_log_rint[i] <- cor.test(df_non_normal[[2]], df_non_normal[[3]], method = 'pearson')$estimate
  
  df_non_normal <- merge(olink_clean[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised_1x_clean_adj_rint[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_non_normal_1x_rint_op[i] <- -cor.test(df_non_normal[[2]], df_non_normal[[3]], method = 'pearson')$estimate
}


# olink rint

for (i in 1:nrow(overlap_cor)) {
  df_non_normal <- merge(olink_clean_adj_rint[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_rint_soma_non_normal[i] <- cor.test(df_non_normal[[2]], df_non_normal[[3]], method = 'pearson')$estimate
  
  df_non_normal <- merge(olink_clean_adj_rint[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised_log_clean[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_rint_soma_non_normal_log[i] <- cor.test(df_non_normal[[2]], df_non_normal[[3]], method = 'pearson')$estimate
  
  df_non_normal <- merge(olink_clean_adj_rint[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised_1x_clean[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_rint_soma_non_normal_1x_op[i] <- -cor.test(df_non_normal[[2]], df_non_normal[[3]], method = 'pearson')$estimate
  
  df_non_normal <- merge(olink_clean_adj_rint[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised_log_clean_adj_rint[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_rint_soma_non_normal_log_rint[i] <- cor.test(df_non_normal[[2]], df_non_normal[[3]], method = 'pearson')$estimate
  
  df_non_normal <- merge(olink_clean_adj_rint[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised_1x_clean_adj_rint[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_rint_soma_non_normal_1x_rint_op[i] <- -cor.test(df_non_normal[[2]], df_non_normal[[3]], method = 'pearson')$estimate
}

# save

write.csv(overlap_cor,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/overlap_cor.csv", quote=F, row.names=F)


##################################### get median

coeff_sum <- data.frame(coeff = names(overlap_cor)[9:30],
                        median = numeric(22))

coeff_sum$median <- lapply(overlap_cor[,c(9:30)],median)

print(coeff_sum, row.names=F)


##################################### plot

# spearman olink vs soma

plot_rho_olink_soma_normal <- ggplot(overlap_cor, aes(x=rho_olink_soma_normal)) + 
  geom_histogram(color="black", fill="white") +
  xlim(-0.5,1) +
  ylim(0,1100) +
  geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
  # ggtitle("Olink vs SomaScan (normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  theme_few()

plot_rho_olink_soma_non_normal <- ggplot(overlap_cor, aes(x=rho_olink_soma_non_normal)) + 
  geom_histogram(color="black", fill="white") +
  xlim(-0.5,1) +
  ylim(0,1100) +
  geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
  # ggtitle("Olink vs SomaScan (non-normalised)") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  theme_few()

plot_rho_olink_soma <- ggarrange(plot_rho_olink_soma_normal,plot_rho_olink_soma_non_normal,labels=c("A","B"),ncol=1,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_rho_olink_soma.png",plot_rho_olink_soma,width=9,height=9)

# pearson olink vs soma

plot_list <- list()

coeff_name <- names(overlap_cor)[11:30]

plot_list <- list()

for (i in 1:20) {
  plot_list[[i]] <- ggplot(overlap_cor, aes(x=.data[[coeff_name[i]]])) + 
    geom_histogram(color="black", fill="white") +
    xlim(-0.55,1) +
    ylim(0,1100) +
    geom_vline(aes(xintercept=median(.data[[coeff_name[i]]])),color="black", linetype="dashed", linewidth=1) +
    xlab(coeff_name[i]) +
    theme_few()
}

plot_r_olink_soma <- ggarrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
                               plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],
                               plot_list[[9]],plot_list[[10]],plot_list[[11]],plot_list[[12]],
                               plot_list[[13]],plot_list[[14]],plot_list[[15]],plot_list[[16]],
                               plot_list[[17]],plot_list[[18]],plot_list[[19]],plot_list[[20]],
                               ncol=5,nrow=4)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_r_olink_soma.png",plot_r_olink_soma,width=19.2,height=10.8)

