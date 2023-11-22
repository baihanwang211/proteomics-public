rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(egg)
library(forcats)
library(ggpattern)
library(ckbplotr)
library(scales)
library(RColorBrewer)

overlap_1_to_1_assoc <- read.csv("overlap_1_to_1_assoc.csv")

## number of hits

# define variables 

variable <- c("region_mean_temp","hours_since_last_ate","age","is_female","region_is_urban",
              "smoking_ever_regular","alcohol_regular_vs_occasion",
              "married", "school", "bmi_calc", "dbp_mean", "sbp_mean", "heart_rate_mean",
              "random_glucose", "poor_health", "diabetes_diag", "kidney_dis_diag", "cancer_diag")

for (i in 1:length(variable)) {
  overlap_1_to_1_assoc$olink_sig_raw <-F
  overlap_1_to_1_assoc$olink_sig_raw[overlap_1_to_1_assoc[10+i*6]<0.05] <- T

  overlap_1_to_1_assoc$soma_normal_sig_raw <- F
  overlap_1_to_1_assoc$soma_normal_sig_raw[overlap_1_to_1_assoc[12+i*6]<0.05] <- T

  overlap_1_to_1_assoc$soma_non_normal_sig_raw <- F
  overlap_1_to_1_assoc$soma_non_normal_sig_raw[overlap_1_to_1_assoc[14+i*6]<0.05] <- T
  
  names(overlap_1_to_1_assoc)[(ncol(overlap_1_to_1_assoc)-2):ncol(overlap_1_to_1_assoc)] <-
    paste(names(overlap_1_to_1_assoc)[(ncol(overlap_1_to_1_assoc)-2):ncol(overlap_1_to_1_assoc)],variable[i],sep="_")
}

# fdr

for (i in 1:length(variable)) {
  overlap_1_to_1_assoc$olink_p_fdr <- F
  overlap_1_to_1_assoc$olink_p_fdr <- p.adjust(overlap_1_to_1_assoc[[10+i*6]],method = "fdr")
  overlap_1_to_1_assoc$olink_sig_fdr <-F
  overlap_1_to_1_assoc$olink_sig_fdr[overlap_1_to_1_assoc$olink_p_fdr<0.05] <- T
  
  overlap_1_to_1_assoc$soma_normal_p_fdr <- F
  overlap_1_to_1_assoc$soma_normal_p_fdr <- p.adjust(overlap_1_to_1_assoc[[12+i*6]],method = "fdr")
  overlap_1_to_1_assoc$soma_normal_sig_fdr <-F
  overlap_1_to_1_assoc$soma_normal_sig_fdr[overlap_1_to_1_assoc$soma_normal_p_fdr<0.05] <- T
  
  overlap_1_to_1_assoc$soma_non_normal_p_fdr <- F
  overlap_1_to_1_assoc$soma_non_normal_p_fdr <- p.adjust(overlap_1_to_1_assoc[[14+i*6]],method = "fdr")
  overlap_1_to_1_assoc$soma_non_normal_sig_fdr <-F
  overlap_1_to_1_assoc$soma_non_normal_sig_fdr[overlap_1_to_1_assoc$soma_non_normal_p_fdr<0.05] <- T
  
  names(overlap_1_to_1_assoc)[(ncol(overlap_1_to_1_assoc)-5):ncol(overlap_1_to_1_assoc)] <-
    paste(names(overlap_1_to_1_assoc)[(ncol(overlap_1_to_1_assoc)-5):ncol(overlap_1_to_1_assoc)],variable[i],sep="_")
}

## count uncorrected

n_hit <- data.frame(colSums(overlap_1_to_1_assoc[grep("sig_raw",names(overlap_1_to_1_assoc))]))
n_hit$variable <- row.names(n_hit)
names(n_hit)[1] <- "n"

# put them into a table

label <- c("Ambient temperature", "Hours since last ate", "Age", "Female", "Urban area",
           "Smoking", "Alcohol drinking", "Married", "Went to high school", "BMI", "DBP", "SBP", "Heart rate",
           "Random glucose", "Poor self-rated health", "Diabetes", "Kidney disease", "Cancer")

level <- c("Ambient temperature", "Hours since last ate", 
           "Age", "Female", "Urban area", "Married", "Went to high school", 
           "BMI", "DBP", "SBP", "Heart rate", "Random glucose", 
           "Poor self-rated health", "Diabetes", "Kidney disease", "Cancer",
           "Smoking", "Alcohol drinking")

var_sample <- c("Ambient temperature","Hours since last ate")
var_demo <- c("Age","Female","Urban area","Married", "Went to high school")
var_phy <- c("BMI", "DBP", "SBP", "Heart rate", "Random glucose")
var_med <- c("Poor self-rated health", "Diabetes", "Kidney disease", "Cancer")
var_life <- c("Smoking","Alcohol drinking")

for (i in 1:length(variable)) {
  n_hit$label[grep(variable[i],n_hit$variable)] <- label[i]
}

n_hit$label <- factor(n_hit$label, levels = level)

n_hit$platform <- NA
n_hit$platform[grep("olink",n_hit$variable)] <- "OLINK"
n_hit$platform[grep("soma_normal",n_hit$variable)] <- "SomaScan (ANML)"
n_hit$platform[grep("soma_non_normal",n_hit$variable)] <- "SomaScan (non-ANML)"

n_hit$platform <- factor(n_hit$platform, c("OLINK","SomaScan (ANML)","SomaScan (non-ANML)"))

n_hit$cat <- NA
n_hit$cat[n_hit$label %in% var_sample] <- "Sample-\nrelated"
n_hit$cat[n_hit$label %in% var_demo] <- "Socio-demographics"
n_hit$cat[n_hit$label %in% var_phy] <- "Physical\nmeasurements"
n_hit$cat[n_hit$label %in% var_med] <- "Medical\nhistory"
n_hit$cat[n_hit$label %in% var_life ] <- "Lifestyle"

n_hit$cat <- factor(n_hit$cat, levels=c("Sample-\nrelated","Socio-demographics","Physical\nmeasurements","Medical\nhistory","Lifestyle"))

hex <- rev(hue_pal()(3))

plot_n_hit <- ggplot(n_hit, aes(y=n, x=fct_rev(label), fill=fct_rev(platform))) + 
  geom_bar(stat="identity",position="dodge",colour="black") +
  ylim (0, 1694) +
  xlab("") +
  ylab("Number of significant associations") +
  labs(fill = "Platform") +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(strip.placement = "outside") +
  facet_grid(cat~.,scales = "free_y",space = "free",switch = "y") +
  coord_flip() +
  scale_fill_manual(values=hex)

plot_n_hit

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_n_hit.png",plot_n_hit,width=9,height=6)

## count fdr

names(overlap_1_to_1_assoc)
n_hit_fdr <- data.frame(colSums(overlap_1_to_1_assoc[grep("sig_fdr",names(overlap_1_to_1_assoc))]))
n_hit_fdr$variable <- row.names(n_hit_fdr)
names(n_hit_fdr)[1] <- "n"

for (i in 1:length(variable)) {
  n_hit_fdr$label[grep(variable[i],n_hit_fdr$variable)] <- label[i]
}

n_hit_fdr$label <- factor(n_hit_fdr$label, levels = level)

n_hit_fdr$platform <- NA
n_hit_fdr$platform[grep("olink",n_hit_fdr$variable)] <- "OLINK"
n_hit_fdr$platform[grep("soma_normal",n_hit_fdr$variable)] <- "SomaScan (ANML)"
n_hit_fdr$platform[grep("soma_non_normal",n_hit_fdr$variable)] <- "SomaScan (non-ANML)"

n_hit_fdr$platform <- factor(n_hit_fdr$platform, c("OLINK","SomaScan (ANML)","SomaScan (non-ANML)"))

n_hit_fdr$cat <- NA
n_hit_fdr$cat[n_hit_fdr$label %in% var_sample] <- "Sample-\nrelated"
n_hit_fdr$cat[n_hit_fdr$label %in% var_demo] <- "Socio-demographics"
n_hit_fdr$cat[n_hit_fdr$label %in% var_phy] <- "Physical\nmeasurements"
n_hit_fdr$cat[n_hit_fdr$label %in% var_med] <- "Medical\nhistory"
n_hit_fdr$cat[n_hit_fdr$label %in% var_life ] <- "Lifestyle"

n_hit_fdr$cat <- factor(n_hit_fdr$cat, levels=c("Sample-\nrelated","Socio-demographics","Physical\nmeasurements","Medical\nhistory","Lifestyle"))

plot_n_hit_fdr <- ggplot(n_hit_fdr, aes(y=n, x=fct_rev(label), fill=fct_rev(platform))) + 
  geom_bar(stat="identity",position="dodge",colour="black") +
  ylim (0, 1694) +
  xlab("") +
  ylab("Number of significant associations") +
  labs(fill = "Platform") +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(strip.placement = "outside") +
  facet_grid(cat~.,scales = "free_y",space = "free",switch = "y") +
  scale_fill_manual(values=hex)

plot_n_hit_fdr

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_n_hit_fdr.png",plot_n_hit_fdr,width=9,height=6)









## whether hits are shared for uncorrected

trans_normal <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(variable)) {
  df <- overlap_1_to_1_assoc[grep(variable[i],names(overlap_1_to_1_assoc))]
  
  df$concord[df[2]<0.05 & df[4]<0.05 & sign(df[1])==sign(df[3])] <- "Shared"
  df$concord[df[2]<0.05 & df[4]>=0.05] <- "Olink-specific"
  df$concord[df[2]>=0.05 & df[4]<0.05] <- "SomaScan-specific"
  df$concord[df[2]<0.05 & df[4]<0.05 & sign(df[1])!=sign(df[3])] <- "Opposite direction"
  df$concord[df[2]>=0.05 & df[4]>=0.05] <- "Neither"
  
  df$concord <- factor(df$concord, levels=c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
  # table(df$concord)
  
  df_sum <- data.frame(table(df$concord))
  names(df_sum) <- c("concord","freq")
  df_sum$variable <- variable[i]
  df_sum$label <- label[i]
  df_sum$soma <- "OLINK vs SomaScan (ANML)"
  
  trans_normal <- rbind(trans_normal,df_sum)
}

trans_non_normal <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(variable)) {
  df <- overlap_1_to_1_assoc[grep(variable[i],names(overlap_1_to_1_assoc))]
  
  df$concord[df[2]<0.05 & df[6]<0.05 & sign(df[1])==sign(df[5])] <- "Shared"
  df$concord[df[2]<0.05 & df[6]>=0.05] <- "Olink-specific"
  df$concord[df[2]>=0.05 & df[6]<0.05] <- "SomaScan-specific"
  df$concord[df[2]<0.05 & df[6]<0.05 & sign(df[1])!=sign(df[5])] <- "Opposite direction"
  df$concord[df[2]>=0.05 & df[6]>=0.05] <- "Neither"
  
  df$concord <- factor(df$concord, levels=c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
  # table(df$concord)
  
  df_sum <- data.frame(table(df$concord))
  names(df_sum) <- c("concord","freq")
  df_sum$variable <- variable[i]
  df_sum$label <- label[i]
  df_sum$soma <- "OLINK vs SomaScan (non-ANML)"
  
  trans_non_normal <- rbind(trans_non_normal,df_sum)
}

trans <- rbind(trans_normal,trans_non_normal)

# plot

hex <- rev(hue_pal()(4))

trans$label <- factor(trans$label, levels = level)

trans$concord <- factor(trans$concord, levels=c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))

trans$percentage <- trans$freq/1694*100

trans$cat <- NA
trans$cat[trans$label %in% var_sample] <- "Sample-\nrelated"
trans$cat[trans$label %in% var_demo] <- "Socio-demographics"
trans$cat[trans$label %in% var_phy] <- "Physical\nmeasurements"
trans$cat[trans$label %in% var_med] <- "Medical\nhistory"
trans$cat[trans$label %in% var_life ] <- "Lifestyle"

trans$cat <- factor(trans$cat, levels=c("Sample-\nrelated","Socio-demographics","Physical\nmeasurements","Medical\nhistory","Lifestyle"))

plot_hit_concordance <- ggplot(trans, aes(y=percentage, x=fct_rev(label), fill=fct_rev(concord))) + 
  geom_bar(stat="identity",position="stack",colour="black") +
  xlab("") +
  ylab("Percentage (%)") +
  labs(fill = "Concordance") +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(strip.placement = "outside") +
  facet_grid(cat~soma,scales = "free_y",space = "free",switch = "y") +
  scale_fill_manual(values=c("gray80",hex))

plot_hit_concordance

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_hit_concordance.png",plot_hit_concordance,width=9,height=6)


## whether hits are shared for fdr

trans_normal_fdr <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(variable)) {
  df <- overlap_1_to_1_assoc[grep(variable[i],names(overlap_1_to_1_assoc))]
  
  df$concord[df[10]<0.05 & df[12]<0.05 & sign(df[1])==sign(df[3])] <- "Shared"
  df$concord[df[10]<0.05 & df[12]>=0.05] <- "Olink-specific"
  df$concord[df[10]>=0.05 & df[12]<0.05] <- "SomaScan-specific"
  df$concord[df[10]<0.05 & df[12]<0.05 & sign(df[1])!=sign(df[3])] <- "Opposite direction"
  df$concord[df[10]>=0.05 & df[12]>=0.05] <- "Neither"
  
  df$concord <- factor(df$concord, levels=c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
  # table(df$concord)
  
  df_sum <- data.frame(table(df$concord))
  names(df_sum) <- c("concord","freq")
  df_sum$variable <- variable[i]
  df_sum$label <- label[i]
  df_sum$soma <- "OLINK vs SomaScan (ANML)"
  
  trans_normal_fdr <- rbind(trans_normal_fdr,df_sum)
}

trans_non_normal_fdr <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(variable)) {
  df <- overlap_1_to_1_assoc[grep(variable[i],names(overlap_1_to_1_assoc))]
  
  df$concord[df[10]<0.05 & df[14]<0.05 & sign(df[1])==sign(df[5])] <- "Shared"
  df$concord[df[10]<0.05 & df[14]>=0.05] <- "Olink-specific"
  df$concord[df[10]>=0.05 & df[14]<0.05] <- "SomaScan-specific"
  df$concord[df[10]<0.05 & df[14]<0.05 & sign(df[1])!=sign(df[5])] <- "Opposite direction"
  df$concord[df[10]>=0.05 & df[14]>=0.05] <- "Neither"
  
  df$concord <- factor(df$concord, levels=c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
  # table(df$concord)
  
  df_sum <- data.frame(table(df$concord))
  names(df_sum) <- c("concord","freq")
  df_sum$variable <- variable[i]
  df_sum$label <- label[i]
  df_sum$soma <- "OLINK vs SomaScan (non-ANML)"
  
  trans_non_normal_fdr <- rbind(trans_non_normal_fdr,df_sum)
}

trans_fdr <- rbind(trans_normal_fdr,trans_non_normal_fdr)

# plot

trans_fdr$label <- factor(trans_fdr$label, levels = level)

trans_fdr$concord <- factor(trans_fdr$concord, levels=c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))

trans_fdr$percentage <- trans_fdr$freq/1694*100

trans_fdr$cat <- NA
trans_fdr$cat[trans_fdr$label %in% var_sample] <- "Sample-\nrelated"
trans_fdr$cat[trans_fdr$label %in% var_demo] <- "Socio-demographics"
trans_fdr$cat[trans_fdr$label %in% var_phy] <- "Physical\nmeasurements"
trans_fdr$cat[trans_fdr$label %in% var_med] <- "Medical\nhistory"
trans_fdr$cat[trans_fdr$label %in% var_life ] <- "Lifestyle"

trans_fdr$cat <- factor(trans_fdr$cat, levels=c("Sample-\nrelated","Socio-demographics","Physical\nmeasurements","Medical\nhistory","Lifestyle"))

plot_hit_concordance_fdr <- ggplot(trans_fdr, aes(y=percentage, x=fct_rev(label), fill=fct_rev(concord))) + 
  geom_bar(stat="identity",position="stack",colour="black") +
  facet_grid(~ soma) +
  xlab("") +
  ylab("Percentage (%)") +
  labs(fill = "Concordance") +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(text = element_text(size = 10)) +
  theme(strip.placement = "outside") +
  facet_grid(cat~soma,scales = "free_y",space = "free",switch = "y") +
  scale_fill_manual(values=c("gray80",hex))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_hit_concordance_fdr.png",plot_hit_concordance_fdr,width=9,height=6)


## whether hits are shared between ANML and non-ANML (uncorrected)

trans_soma <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(variable)) {
  df <- overlap_1_to_1_assoc[grep(variable[i],names(overlap_1_to_1_assoc))]
  
  df$concord[df[4]<0.05 & df[6]<0.05 & sign(df[3])==sign(df[5])] <- "Shared"
  df$concord[df[4]<0.05 & df[6]>=0.05] <- "ANML-specific"
  df$concord[df[4]>=0.05 & df[6]<0.05] <- "Non-ANML-specific"
  df$concord[df[4]<0.05 & df[6]<0.05 & sign(df[3])!=sign(df[5])] <- "Opposite direction"
  df$concord[df[4]>=0.05 & df[6]>=0.05] <- "Neither"
  
  df$concord <- factor(df$concord, levels=c("Shared","Opposite direction","ANML-specific","Non-ANML-specific","Neither"))
  # table(df$concord)
  
  df_sum <- data.frame(table(df$concord))
  names(df_sum) <- c("concord","freq")
  df_sum$variable <- variable[i]
  df_sum$label <- label[i]

  trans_soma <- rbind(trans_soma,df_sum)
}

# plot

trans_soma$label <- factor(trans_soma$label, levels = level)

trans_soma$concord <- factor(trans_soma$concord, levels=c("Shared","Opposite direction","ANML-specific","Non-ANML-specific","Neither"))

trans_soma$percentage <- trans_soma$freq/1694*100

trans_soma$cat <- NA
trans_soma$cat[trans_soma$label %in% var_sample] <- "Sample-\nrelated"
trans_soma$cat[trans_soma$label %in% var_demo] <- "Socio-demographics"
trans_soma$cat[trans_soma$label %in% var_phy] <- "Physical\nmeasurements"
trans_soma$cat[trans_soma$label %in% var_med] <- "Medical\nhistory"
trans_soma$cat[trans_soma$label %in% var_life ] <- "Lifestyle"

trans_soma$cat <- factor(trans_soma$cat, levels=c("Sample-\nrelated","Socio-demographics","Physical\nmeasurements","Medical\nhistory","Lifestyle"))

plot_hit_concordance_soma <- ggplot(trans_soma, aes(y=percentage, x=fct_rev(label), fill=fct_rev(concord))) + 
  geom_bar(stat="identity",position="stack",colour="black") +
  xlab("") +
  ylab("Percentage (%)") +
  labs(fill = "Concordance") +
  ggtitle("SomaScan ANML vs non-ANML") +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(strip.placement = "outside") + 
  facet_grid(cat~.,scales = "free_y",space = "free",switch = "y") +
  scale_fill_manual(values=c("gray80",hex))

plot_hit_concordance_soma

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_hit_concordance_soma.png",plot_hit_concordance_soma,width=6,height=6)

## whether hits are shared between ANML and non-ANML (fdr)

trans_soma_fdr <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(variable)) {
  df <- overlap_1_to_1_assoc[grep(variable[i],names(overlap_1_to_1_assoc))]
  
  df$concord[df[12]<0.05 & df[14]<0.05 & sign(df[3])==sign(df[5])] <- "Shared"
  df$concord[df[12]<0.05 & df[14]>=0.05] <- "ANML-specific"
  df$concord[df[12]>=0.05 & df[14]<0.05] <- "Non-ANML-specific"
  df$concord[df[12]<0.05 & df[14]<0.05 & sign(df[3])!=sign(df[5])] <- "Opposite direction"
  df$concord[df[12]>=0.05 & df[14]>=0.05] <- "Neither"
  
  df$concord <- factor(df$concord, levels=c("Shared","Opposite direction","ANML-specific","Non-ANML-specific","Neither"))
  # table(df$concord)
  
  df_sum <- data.frame(table(df$concord))
  names(df_sum) <- c("concord","freq")
  df_sum$variable <- variable[i]
  df_sum$label <- label[i]
  
  trans_soma_fdr <- rbind(trans_soma_fdr,df_sum)
}

# plot

trans_soma_fdr$label <- factor(trans_soma_fdr$label, levels = level)

trans_soma_fdr$concord <- factor(trans_soma_fdr$concord, levels=c("Shared","Opposite direction","ANML-specific","Non-ANML-specific","Neither"))

trans_soma_fdr$percentage <- trans_soma_fdr$freq/1694*100

trans_soma_fdr$cat <- NA
trans_soma_fdr$cat[trans_soma_fdr$label %in% var_sample] <- "Sample-\nrelated"
trans_soma_fdr$cat[trans_soma_fdr$label %in% var_demo] <- "Socio-demographics"
trans_soma_fdr$cat[trans_soma_fdr$label %in% var_phy] <- "Physical\nmeasurements"
trans_soma_fdr$cat[trans_soma_fdr$label %in% var_med] <- "Medical\nhistory"
trans_soma_fdr$cat[trans_soma_fdr$label %in% var_life ] <- "Lifestyle"

trans_soma_fdr$cat <- factor(trans_soma_fdr$cat, levels=c("Sample-\nrelated","Socio-demographics","Physical\nmeasurements","Medical\nhistory","Lifestyle"))

plot_hit_concordance_soma_fdr <- ggplot(trans_soma_fdr, aes(y=percentage, x=fct_rev(label), fill=fct_rev(concord))) + 
  geom_bar(stat="identity",position="stack",colour="black") +
  xlab("Variable") +
  ylab("Percentage (%)") +
  ggtitle("SomaScan ANML vs non-ANML") +
  labs(fill = "Concordance") +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(strip.placement = "outside") +
  facet_grid(cat~.,scales = "free_y",space = "free",switch = "y") +
  scale_fill_manual(values=c("gray80",hex))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_hit_concordance_soma_fdr.png",plot_hit_concordance_soma_fdr,width=6,height=6)




# # percentage shared out of pairs with at least one sig
# 
# trans$percent_shared <- 0
# 
# for (i in 1:(length(variable)*2)) {
#   trans$percent_shared[((5*i-4)):(5*i)] <- trans$freq[5*i-4]/(1694-trans$freq[5*i])
# }




## correlation coefficient between beta

cor_es_normal <- data.frame(matrix(ncol = 5, nrow = 0))
cor_es_non_normal <- data.frame(matrix(ncol = 5, nrow = 0))
cor_es_soma <- data.frame(matrix(ncol = 5, nrow = 0))

for (i in 1:length(variable)) {
  
  df <- overlap_1_to_1_assoc[grep(paste0("es_",variable[i]),names(overlap_1_to_1_assoc))]
  
  df_1 <- data.frame(variable=variable[i],label=label[i],
                     r=cor(df[[1]],df[[2]]),
                     lci=cor.test(df[[1]],df[[2]])$conf.int[1],
                     uci=cor.test(df[[1]],df[[2]])$conf.int[2])
  
  cor_es_normal <- rbind(cor_es_normal,df_1)
  
  df_2 <- data.frame(variable=variable[i],label=label[i],
                     r=cor(df[[1]],df[[3]]),
                     lci=cor.test(df[[1]],df[[3]])$conf.int[1],
                     uci=cor.test(df[[1]],df[[3]])$conf.int[2])
  
  cor_es_non_normal <- rbind(cor_es_non_normal,df_2)
  
  df_3 <- data.frame(variable=variable[i],label=label[i],
                     r=cor(df[[2]],df[[3]]),
                     lci=cor.test(df[[2]],df[[3]])$conf.int[1],
                     uci=cor.test(df[[2]],df[[3]])$conf.int[2])
  
  cor_es_soma <- rbind(cor_es_soma,df_3)
}

cor_es_normal$label <- factor(cor_es_normal$label, levels = level)
shared_normal <- trans_normal_fdr[trans_normal_fdr$concord=="Shared",c("label","concord","freq")]
cor_es_normal <- merge(cor_es_normal,shared_normal,by="label")

cor_es_non_normal$label <- factor(cor_es_non_normal$label, levels = level)
shared_non_normal <- trans_non_normal_fdr[trans_non_normal_fdr$concord=="Shared",c("label","concord","freq")]
cor_es_non_normal <- merge(cor_es_non_normal,shared_non_normal,by="label")

cor_es_soma$label <- factor(cor_es_soma$label, levels = level)
shared_soma <- trans_soma_fdr[trans_soma_fdr$concord=="Shared",c("label","concord","freq")]
cor_es_soma <- merge(cor_es_soma,shared_soma,by="label")

cor_es <- list(cor_es_normal,cor_es_non_normal)

row_labels <- data.frame(
  subgroup = level,
  group    = c("Sample-related","Sample-related",
               "Socio-demographics","Socio-demographics","Socio-demographics","Socio-demographics","Socio-demographics",
               "Physical measurements","Physical measurements","Physical measurements","Physical measurements","Physical measurements",
               "Medical history","Medical history","Medical history","Medical history",
               "Lifestyle","Lifestyle"),
  label = level)

plot_cor_es <- forest_plot(panels = cor_es,
            exponentiate = F,
            panel.headings = c("OLINK vs SomaScan (ANML)","OLINK vs SomaScan (non-ANML)"),
            row.labels = row_labels,
            row.labels.levels = c("group","subgroup"),
            rows = unique(row_labels$group),
            col.key = "label",
            col.estimate = "r",
            col.lci = "lci",
            col.uci = "uci",
            col.left = "freq",
            col.left.heading = "No. shared hits",
            xlab = "Pearson's r (95% CI)",
            xlim = c(0.2,0.7),
            # estcolumn = F,
            pointsize = 2,
            base_size = 12,
            col.right.heading = "Pearson's r (95% CI)"
)

plot_cor_es_bg <- plot_cor_es$plot + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_cor_es.png",plot_cor_es_bg,width=14,height=6)

plot_cor_es_soma <- forest_plot(panels = list(cor_es_soma),
                           exponentiate = F,
                           panel.headings = c("SomaScan ANML vs non-ANML"),
                           row.labels = row_labels,
                           row.labels.levels = c("group","subgroup"),
                           rows = unique(row_labels$group),
                           col.key = "label",
                           col.estimate = "r",
                           col.lci = "lci",
                           col.uci = "uci",
                           col.left = "freq",
                           col.left.heading = "No. shared hits",
                           xlab = "Pearson's r (95% CI)",
                           # xlim = c(0.2,0.7),
                           # estcolumn = F,
                           pointsize = 2,
                           base_size = 12,
                           col.right.heading = "Pearson's r (95% CI)"
)

plot_cor_es_soma_bg <- plot_cor_es_soma$plot + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_cor_es_soma.png",plot_cor_es_soma_bg,width=7,height=6)



## correlation coefficient between beta, only for shared associations

cor_es_normal_sig <- data.frame(matrix(ncol = 5, nrow = 0))
cor_es_non_normal_sig <- data.frame(matrix(ncol = 5, nrow = 0))
cor_es_soma_sig <- data.frame(matrix(ncol = 5, nrow = 0))

# remove variables with two few shared hits

variable_sig <- variable[! variable %in% c("married","school","kidney_dis_diag","cancer_diag")]
label_sig <- label[! label %in% c("Married","Went to high school","Kidney disease","Cancer")]
level_sig <- level[! level %in% c("Married","Went to high school","Kidney disease","Cancer")]

for (i in 1:length(variable_sig)) {
  
  df <- overlap_1_to_1_assoc[c(grep(paste0("es_",variable_sig[i]),names(overlap_1_to_1_assoc)),
                               grep(paste0("p_fdr_",variable_sig[i]),names(overlap_1_to_1_assoc)))]
  
  df_1 <- df[df[4]<0.05 & df[5]<0.05 & sign(df[1])==sign(df[2]),]

  df_1_result <- data.frame(variable=variable_sig[i],label=label_sig[i],
                     r=cor(df_1[[1]],df_1[[2]]),
                     lci=cor.test(df_1[[1]],df_1[[2]])$conf.int[1],
                     uci=cor.test(df_1[[1]],df_1[[2]])$conf.int[2])
  
  cor_es_normal_sig <- rbind(cor_es_normal_sig,df_1_result)
  
  df_2 <- df[df[4]<0.05 & df[6]<0.05 & sign(df[1])==sign(df[3]),]
  
  df_2_result <- data.frame(variable=variable_sig[i],label=label_sig[i],
                     r=cor(df_2[[1]],df_2[[3]]),
                     lci=cor.test(df_2[[1]],df_2[[3]])$conf.int[1],
                     uci=cor.test(df_2[[1]],df_2[[3]])$conf.int[2])
  
  cor_es_non_normal_sig <- rbind(cor_es_non_normal_sig,df_2_result)
  
  df_3 <- df[df[5]<0.05 & df[6]<0.05 & sign(df[2])==sign(df[3]),]
  
  df_3_result <- data.frame(variable=variable_sig[i],label=label_sig[i],
                     r=cor(df_3[[2]],df_3[[3]]),
                     lci=cor.test(df_3[[2]],df_3[[3]])$conf.int[1],
                     uci=cor.test(df_3[[2]],df_3[[3]])$conf.int[2])
  
  cor_es_soma_sig <- rbind(cor_es_soma_sig,df_3_result)
}

cor_es_normal_sig$label <- factor(cor_es_normal_sig$label, levels = level_sig)
shared_normal_sig <- trans_normal_fdr[trans_normal_fdr$concord=="Shared",c("label","concord","freq")]
cor_es_normal_sig <- merge(cor_es_normal_sig,shared_normal,by="label")

cor_es_non_normal_sig$label <- factor(cor_es_non_normal_sig$label, levels = level_sig)
shared_non_normal_sig <- trans_non_normal_fdr[trans_non_normal_fdr$concord=="Shared",c("label","concord","freq")]
cor_es_non_normal_sig <- merge(cor_es_non_normal_sig,shared_non_normal,by="label")

cor_es_soma_sig$label <- factor(cor_es_soma_sig$label, levels = level_sig)
shared_soma_sig <- trans_soma_fdr[trans_soma_fdr$concord=="Shared",c("label","concord","freq")]
cor_es_soma_sig <- merge(cor_es_soma_sig,shared_soma,by="label")

cor_es_sig <- list(cor_es_normal_sig,cor_es_non_normal_sig)

row_labels_sig <- data.frame(
  subgroup = level_sig,
  group    = c("Sample-related","Sample-related",
               "Socio-demographics","Socio-demographics","Socio-demographics",
               "Physical measurements","Physical measurements","Physical measurements","Physical measurements","Physical measurements",
               "Medical history","Medical history",
               "Lifestyle","Lifestyle"),
  label = level_sig)

plot_cor_es_sig <- forest_plot(panels = cor_es_sig,
                           exponentiate = F,
                           panel.headings = c("OLINK vs SomaScan (ANML)","OLINK vs SomaScan (non-ANML)"),
                           row.labels = row_labels_sig,
                           row.labels.levels = c("group","subgroup"),
                           rows = unique(row_labels$group),
                           col.key = "label",
                           col.estimate = "r",
                           col.lci = "lci",
                           col.uci = "uci",
                           col.left = "freq",
                           col.left.heading = "No. shared hits",
                           xlab = "Pearson's r (95% CI)",
                           # xlim = c(0.2,0.7),
                           # estcolumn = F,
                           pointsize = 2,
                           base_size = 12,
                           col.right.heading = "Pearson's r (95% CI)"
)

plot_cor_es_sig_bg <- plot_cor_es_sig$plot + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_cor_es_sig.png",plot_cor_es_sig_bg,width=14,height=6)

plot_cor_es_soma_sig <- forest_plot(panels = list(cor_es_soma_sig),
                                exponentiate = F,
                                panel.headings = c("SomaScan ANML vs non-ANML"),
                                row.labels = row_labels_sig,
                                row.labels.levels = c("group","subgroup"),
                                rows = unique(row_labels$group),
                                col.key = "label",
                                col.estimate = "r",
                                col.lci = "lci",
                                col.uci = "uci",
                                col.left = "freq",
                                col.left.heading = "No. shared hits",
                                xlab = "Pearson's r (95% CI)",
                                # xlim = c(0.2,0.7),
                                # estcolumn = F,
                                pointsize = 2,
                                base_size = 12,
                                col.right.heading = "Pearson's r (95% CI)"
)

plot_cor_es_soma_sig_bg <- plot_cor_es_soma_sig$plot + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_cor_es_soma_sig.png",plot_cor_es_soma_sig_bg,width=7,height=6)



write.csv(overlap_1_to_1_assoc,"overlap_1_to_1_assoc_sig.csv", quote=F, row.names=F)



# ### temperature
# 
# ## scatter
# 
# # beta
# 
# scatter_temp_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_temp, y=soma_normal_es_temp,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.07,0.07) +
#   ylim(-0.07,0.07) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_temp_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_temp,overlap_1_to_1_assoc$soma_normal_es_temp)
# 
# # temp non-normal
# 
# scatter_temp_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_temp, y=soma_non_normal_es_temp,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.07,0.07) +
#   ylim(-0.07,0.07) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("Non-ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   geom_abline(linetype = "dashed") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_temp_non_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_temp,overlap_1_to_1_assoc$soma_non_normal_es_temp)
# 
# scatter_temp <- ggarrange(scatter_temp_normal,scatter_temp_non_normal,ncol=2,nrow=1) 
# 
# scatter_temp <- annotate_figure(scatter_temp, top = text_grob("Ambient temperature", face = "bold", size = 28))
# 
# scatter_temp
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/scatter_temp.png",scatter_temp,width=16,height=8, bg = "white")
# 
# # histogram temp
# 
# overlap_1_to_1_assoc$sig_temp_normal <- "n.s."
# overlap_1_to_1_assoc$sig_temp_normal[overlap_1_to_1_assoc$olink_sig_temp==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_temp_normal[overlap_1_to_1_assoc$soma_normal_sig_temp==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_temp_normal[overlap_1_to_1_assoc$olink_sig_temp==T & overlap_1_to_1_assoc$soma_normal_sig_temp==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_temp_normal <- factor(overlap_1_to_1_assoc$sig_temp_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_temp_normal)
# 
# hist_temp_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_normal,fill=sig_temp_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_temp_normal
# 
# overlap_1_to_1_assoc$sig_temp_non_normal <- "n.s."
# overlap_1_to_1_assoc$sig_temp_non_normal[overlap_1_to_1_assoc$olink_sig_temp==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_temp_non_normal[overlap_1_to_1_assoc$soma_non_normal_sig_temp==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_temp_non_normal[overlap_1_to_1_assoc$olink_sig_temp==T & overlap_1_to_1_assoc$soma_non_normal_sig_temp==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_temp_non_normal <- factor(overlap_1_to_1_assoc$sig_temp_non_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_temp_non_normal)
# 
# hist_temp_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_non_normal,fill=sig_temp_non_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("Non-ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_temp_non_normal
# 
# hist_temp <- ggarrange(hist_temp_normal,hist_temp_non_normal,ncol=2,nrow=1)
# 
# hist_temp <- annotate_figure(hist_temp, top = text_grob("Ambient temperature", face = "bold", size = 28))
# 
# hist_temp
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_temp.png",hist_temp,width=16,height=6,bg="white")
# 
# 
# 
# # hours since laste ate
# 
# scatter_ate_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_ate, y=soma_normal_es_ate,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.12,0.12) +
#   ylim(-0.12,0.12) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_ate_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_ate,overlap_1_to_1_assoc$soma_normal_es_ate)
# 
# # ate non-normal
# 
# scatter_ate_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_ate, y=soma_non_normal_es_ate,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.12,0.12) +
#   ylim(-0.12,0.12) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("Non-ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   geom_abline(linetype = "dashed") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_ate_non_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_ate,overlap_1_to_1_assoc$soma_non_normal_es_ate)
# 
# scatter_ate <- ggarrange(scatter_ate_normal,scatter_ate_non_normal,ncol=2,nrow=1) 
# 
# scatter_ate <- annotate_figure(scatter_ate, top = text_grob("Hours since last ate", face = "bold", size = 28))
# 
# scatter_ate
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/scatter_ate.png",scatter_ate,width=16,height=8, bg = "white")
# 
# # histogram ate
# 
# overlap_1_to_1_assoc$sig_ate_normal <- "n.s."
# overlap_1_to_1_assoc$sig_ate_normal[overlap_1_to_1_assoc$olink_sig_ate==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_ate_normal[overlap_1_to_1_assoc$soma_normal_sig_ate==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_ate_normal[overlap_1_to_1_assoc$olink_sig_ate==T & overlap_1_to_1_assoc$soma_normal_sig_ate==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_ate_normal <- factor(overlap_1_to_1_assoc$sig_ate_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_ate_normal)
# 
# hist_ate_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_normal,fill=sig_ate_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_ate_normal
# 
# overlap_1_to_1_assoc$sig_ate_non_normal <- "n.s."
# overlap_1_to_1_assoc$sig_ate_non_normal[overlap_1_to_1_assoc$olink_sig_ate==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_ate_non_normal[overlap_1_to_1_assoc$soma_non_normal_sig_ate==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_ate_non_normal[overlap_1_to_1_assoc$olink_sig_ate==T & overlap_1_to_1_assoc$soma_non_normal_sig_ate==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_ate_non_normal <- factor(overlap_1_to_1_assoc$sig_ate_non_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_ate_non_normal)
# 
# hist_ate_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_non_normal,fill=sig_ate_non_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("Non-ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_ate_non_normal
# 
# hist_ate <- ggarrange(hist_ate_normal,hist_ate_non_normal,ncol=2,nrow=1)
# 
# hist_ate <- annotate_figure(hist_ate, top = text_grob("Hours since last ate", face = "bold", size = 28))
# 
# hist_ate
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_ate.png",hist_ate,width=16,height=6,bg="white")
# 
# 
# 
# 
# # age
# 
# scatter_age_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_age, y=soma_normal_es_age,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.07,0.07) +
#   ylim(-0.07,0.07) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_age_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_age,overlap_1_to_1_assoc$soma_normal_es_age)
# 
# # age non-normal
# 
# scatter_age_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_age, y=soma_non_normal_es_age,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.07,0.07) +
#   ylim(-0.07,0.07) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("Non-ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   geom_abline(linetype = "dashed") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_age_non_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_age,overlap_1_to_1_assoc$soma_non_normal_es_age)
# 
# scatter_age <- ggarrange(scatter_age_normal,scatter_age_non_normal,ncol=2,nrow=1) 
# 
# scatter_age <- annotate_figure(scatter_age, top = text_grob("Age", face = "bold", size = 28))
# 
# scatter_age
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/scatter_age.png",scatter_age,width=16,height=8, bg = "white")
# 
# # histogram age
# 
# overlap_1_to_1_assoc$sig_age_normal <- "n.s."
# overlap_1_to_1_assoc$sig_age_normal[overlap_1_to_1_assoc$olink_sig_age==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_age_normal[overlap_1_to_1_assoc$soma_normal_sig_age==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_age_normal[overlap_1_to_1_assoc$olink_sig_age==T & overlap_1_to_1_assoc$soma_normal_sig_age==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_age_normal <- factor(overlap_1_to_1_assoc$sig_age_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_age_normal)
# 
# hist_age_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_normal,fill=sig_age_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_age_normal
# 
# overlap_1_to_1_assoc$sig_age_non_normal <- "n.s."
# overlap_1_to_1_assoc$sig_age_non_normal[overlap_1_to_1_assoc$olink_sig_age==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_age_non_normal[overlap_1_to_1_assoc$soma_non_normal_sig_age==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_age_non_normal[overlap_1_to_1_assoc$olink_sig_age==T & overlap_1_to_1_assoc$soma_non_normal_sig_age==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_age_non_normal <- factor(overlap_1_to_1_assoc$sig_age_non_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_age_non_normal)
# 
# hist_age_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_non_normal,fill=sig_age_non_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("Non-ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_age_non_normal
# 
# hist_age <- ggarrange(hist_age_normal,hist_age_non_normal,ncol=2,nrow=1)
# 
# hist_age <- annotate_figure(hist_age, top = text_grob("Age", face = "bold", size = 28))
# 
# hist_age
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_age.png",hist_age,width=16,height=6,bg="white")
# 
# 
# 
# 
# # sex
# 
# scatter_sex_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_sex, y=soma_normal_es_sex,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-2,2) +
#   ylim(-2,2) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_sex_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_sex,overlap_1_to_1_assoc$soma_normal_es_sex)
# 
# # sex non-normal
# 
# scatter_sex_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_sex, y=soma_non_normal_es_sex,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-2,2) +
#   ylim(-2,2) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("Non-ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   geom_abline(linetype = "dashed") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_sex_non_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_sex,overlap_1_to_1_assoc$soma_non_normal_es_sex)
# 
# scatter_sex <- ggarrange(scatter_sex_normal,scatter_sex_non_normal,ncol=2,nrow=1) 
# 
# scatter_sex <- annotate_figure(scatter_sex, top = text_grob("Sex", face = "bold", size = 28))
# 
# scatter_sex
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/scatter_sex.png",scatter_sex,width=16,height=8, bg = "white")
# 
# # histogram sex
# 
# overlap_1_to_1_assoc$sig_sex_normal <- "n.s."
# overlap_1_to_1_assoc$sig_sex_normal[overlap_1_to_1_assoc$olink_sig_sex==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_sex_normal[overlap_1_to_1_assoc$soma_normal_sig_sex==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_sex_normal[overlap_1_to_1_assoc$olink_sig_sex==T & overlap_1_to_1_assoc$soma_normal_sig_sex==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_sex_normal <- factor(overlap_1_to_1_assoc$sig_sex_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_sex_normal)
# 
# hist_sex_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_normal,fill=sig_sex_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_sex_normal
# 
# overlap_1_to_1_assoc$sig_sex_non_normal <- "n.s."
# overlap_1_to_1_assoc$sig_sex_non_normal[overlap_1_to_1_assoc$olink_sig_sex==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_sex_non_normal[overlap_1_to_1_assoc$soma_non_normal_sig_sex==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_sex_non_normal[overlap_1_to_1_assoc$olink_sig_sex==T & overlap_1_to_1_assoc$soma_non_normal_sig_sex==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_sex_non_normal <- factor(overlap_1_to_1_assoc$sig_sex_non_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_sex_non_normal)
# 
# hist_sex_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_non_normal,fill=sig_sex_non_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("Non-ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_sex_non_normal
# 
# hist_sex <- ggarrange(hist_sex_normal,hist_sex_non_normal,ncol=2,nrow=1)
# 
# hist_sex <- annotate_figure(hist_sex, top = text_grob("Sex", face = "bold", size = 28))
# 
# hist_sex
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_sex.png",hist_sex,width=16,height=6,bg="white")
# 
# 
# 
# 
# 
# 
# 
# # bmi
# 
# scatter_bmi_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_bmi, y=soma_normal_es_bmi,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.16,0.16) +
#   ylim(-0.16,0.16) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_bmi_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_bmi,overlap_1_to_1_assoc$soma_normal_es_bmi)
# 
# # bmi non-normal
# 
# scatter_bmi_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_bmi, y=soma_non_normal_es_bmi,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.16,0.16) +
#   ylim(-0.16,0.16) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("Non-ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   geom_abline(linetype = "dashed") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_bmi_non_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_bmi,overlap_1_to_1_assoc$soma_non_normal_es_bmi)
# 
# scatter_bmi <- ggarrange(scatter_bmi_normal,scatter_bmi_non_normal,ncol=2,nrow=1) 
# 
# scatter_bmi <- annotate_figure(scatter_bmi, top = text_grob("BMI", face = "bold", size = 28))
# 
# scatter_bmi
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/scatter_bmi.png",scatter_bmi,width=16,height=8, bg = "white")
# 
# # histogram bmi
# 
# overlap_1_to_1_assoc$sig_bmi_normal <- "n.s."
# overlap_1_to_1_assoc$sig_bmi_normal[overlap_1_to_1_assoc$olink_sig_bmi==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_bmi_normal[overlap_1_to_1_assoc$soma_normal_sig_bmi==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_bmi_normal[overlap_1_to_1_assoc$olink_sig_bmi==T & overlap_1_to_1_assoc$soma_normal_sig_bmi==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_bmi_normal <- factor(overlap_1_to_1_assoc$sig_bmi_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_bmi_normal)
# 
# hist_bmi_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_normal,fill=sig_bmi_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_bmi_normal
# 
# overlap_1_to_1_assoc$sig_bmi_non_normal <- "n.s."
# overlap_1_to_1_assoc$sig_bmi_non_normal[overlap_1_to_1_assoc$olink_sig_bmi==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_bmi_non_normal[overlap_1_to_1_assoc$soma_non_normal_sig_bmi==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_bmi_non_normal[overlap_1_to_1_assoc$olink_sig_bmi==T & overlap_1_to_1_assoc$soma_non_normal_sig_bmi==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_bmi_non_normal <- factor(overlap_1_to_1_assoc$sig_bmi_non_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_bmi_non_normal)
# 
# hist_bmi_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_non_normal,fill=sig_bmi_non_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("Non-ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_bmi_non_normal
# 
# hist_bmi <- ggarrange(hist_bmi_normal,hist_bmi_non_normal,ncol=2,nrow=1)
# 
# hist_bmi <- annotate_figure(hist_bmi, top = text_grob("BMI", face = "bold", size = 28))
# 
# hist_bmi
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_bmi.png",hist_bmi,width=16,height=6,bg="white")
# 
# 
# 
# 
# 
# 
# 
# 
# # sbp
# 
# scatter_sbp_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_sbp, y=soma_normal_es_sbp,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.015,0.015) +
#   ylim(-0.015,0.015) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_sbp_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_sbp,overlap_1_to_1_assoc$soma_normal_es_sbp)
# 
# # sbp non-normal
# 
# scatter_sbp_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_sbp, y=soma_non_normal_es_sbp,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.015,0.015) +
#   ylim(-0.015,0.015) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("Non-ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   geom_abline(linetype = "dashed") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_sbp_non_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_sbp,overlap_1_to_1_assoc$soma_non_normal_es_sbp)
# 
# scatter_sbp <- ggarrange(scatter_sbp_normal,scatter_sbp_non_normal,ncol=2,nrow=1) 
# 
# scatter_sbp <- annotate_figure(scatter_sbp, top = text_grob("Systolic blood pressure", face = "bold", size = 28))
# 
# scatter_sbp
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/scatter_sbp.png",scatter_sbp,width=16,height=8, bg = "white")
# 
# # histogram sbp
# 
# overlap_1_to_1_assoc$sig_sbp_normal <- "n.s."
# overlap_1_to_1_assoc$sig_sbp_normal[overlap_1_to_1_assoc$olink_sig_sbp==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_sbp_normal[overlap_1_to_1_assoc$soma_normal_sig_sbp==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_sbp_normal[overlap_1_to_1_assoc$olink_sig_sbp==T & overlap_1_to_1_assoc$soma_normal_sig_sbp==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_sbp_normal <- factor(overlap_1_to_1_assoc$sig_sbp_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_sbp_normal)
# 
# hist_sbp_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_normal,fill=sig_sbp_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_sbp_normal
# 
# overlap_1_to_1_assoc$sig_sbp_non_normal <- "n.s."
# overlap_1_to_1_assoc$sig_sbp_non_normal[overlap_1_to_1_assoc$olink_sig_sbp==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_sbp_non_normal[overlap_1_to_1_assoc$soma_non_normal_sig_sbp==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_sbp_non_normal[overlap_1_to_1_assoc$olink_sig_sbp==T & overlap_1_to_1_assoc$soma_non_normal_sig_sbp==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_sbp_non_normal <- factor(overlap_1_to_1_assoc$sig_sbp_non_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_sbp_non_normal)
# 
# hist_sbp_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_non_normal,fill=sig_sbp_non_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("Non-ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_sbp_non_normal
# 
# hist_sbp <- ggarrange(hist_sbp_normal,hist_sbp_non_normal,ncol=2,nrow=1)
# 
# hist_sbp <- annotate_figure(hist_sbp, top = text_grob("Systolic blood pressure", face = "bold", size = 28))
# 
# hist_sbp
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_sbp.png",hist_sbp,width=16,height=6,bg="white")
# 
# 
# 
# 
# 
# 
# # heart rate
# 
# scatter_hr_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_hr, y=soma_normal_es_hr,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.02,0.02) +
#   ylim(-0.02,0.02) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_hr_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_hr,overlap_1_to_1_assoc$soma_normal_es_hr)
# 
# # hr non-normal
# 
# scatter_hr_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_hr, y=soma_non_normal_es_hr,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.02,0.02) +
#   ylim(-0.02,0.02) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("Non-ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   geom_abline(linetype = "dashed") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_hr_non_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_hr,overlap_1_to_1_assoc$soma_non_normal_es_hr)
# 
# scatter_hr <- ggarrange(scatter_hr_normal,scatter_hr_non_normal,ncol=2,nrow=1) 
# 
# scatter_hr <- annotate_figure(scatter_hr, top = text_grob("Heart rate", face = "bold", size = 28))
# 
# scatter_hr
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/scatter_hr.png",scatter_hr,width=16,height=8, bg = "white")
# 
# # histogram hr
# 
# overlap_1_to_1_assoc$sig_hr_normal <- "n.s."
# overlap_1_to_1_assoc$sig_hr_normal[overlap_1_to_1_assoc$olink_sig_hr==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_hr_normal[overlap_1_to_1_assoc$soma_normal_sig_hr==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_hr_normal[overlap_1_to_1_assoc$olink_sig_hr==T & overlap_1_to_1_assoc$soma_normal_sig_hr==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_hr_normal <- factor(overlap_1_to_1_assoc$sig_hr_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_hr_normal)
# 
# hist_hr_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_normal,fill=sig_hr_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_hr_normal
# 
# overlap_1_to_1_assoc$sig_hr_non_normal <- "n.s."
# overlap_1_to_1_assoc$sig_hr_non_normal[overlap_1_to_1_assoc$olink_sig_hr==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_hr_non_normal[overlap_1_to_1_assoc$soma_non_normal_sig_hr==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_hr_non_normal[overlap_1_to_1_assoc$olink_sig_hr==T & overlap_1_to_1_assoc$soma_non_normal_sig_hr==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_hr_non_normal <- factor(overlap_1_to_1_assoc$sig_hr_non_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_hr_non_normal)
# 
# hist_hr_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_non_normal,fill=sig_hr_non_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("Non-ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_hr_non_normal
# 
# hist_hr <- ggarrange(hist_hr_normal,hist_hr_non_normal,ncol=2,nrow=1)
# 
# hist_hr <- annotate_figure(hist_hr, top = text_grob("Heart rate", face = "bold", size = 28))
# 
# hist_hr
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_hr.png",hist_hr,width=16,height=6,bg="white")
# 
# 
# 
# 
# 
# 
# # random glucose
# 
# scatter_glu_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_glu, y=soma_normal_es_glu,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.2,0.2) +
#   ylim(-0.2,0.2) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_glu_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_glu,overlap_1_to_1_assoc$soma_normal_es_glu)
# 
# # glu non-normal
# 
# scatter_glu_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_glu, y=soma_non_normal_es_glu,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-0.2,0.2) +
#   ylim(-0.2,0.2) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("Non-ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   geom_abline(linetype = "dashed") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_glu_non_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_glu,overlap_1_to_1_assoc$soma_non_normal_es_glu)
# 
# scatter_glu <- ggarrange(scatter_glu_normal,scatter_glu_non_normal,ncol=2,nrow=1) 
# 
# scatter_glu <- annotate_figure(scatter_glu, top = text_grob("Random glucose", face = "bold", size = 28))
# 
# scatter_glu
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/scatter_glu.png",scatter_glu,width=16,height=8, bg = "white")
# 
# # histogram glu
# 
# overlap_1_to_1_assoc$sig_glu_normal <- "n.s."
# overlap_1_to_1_assoc$sig_glu_normal[overlap_1_to_1_assoc$olink_sig_glu==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_glu_normal[overlap_1_to_1_assoc$soma_normal_sig_glu==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_glu_normal[overlap_1_to_1_assoc$olink_sig_glu==T & overlap_1_to_1_assoc$soma_normal_sig_glu==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_glu_normal <- factor(overlap_1_to_1_assoc$sig_glu_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_glu_normal)
# 
# hist_glu_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_normal,fill=sig_glu_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_glu_normal
# 
# overlap_1_to_1_assoc$sig_glu_non_normal <- "n.s."
# overlap_1_to_1_assoc$sig_glu_non_normal[overlap_1_to_1_assoc$olink_sig_glu==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_glu_non_normal[overlap_1_to_1_assoc$soma_non_normal_sig_glu==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_glu_non_normal[overlap_1_to_1_assoc$olink_sig_glu==T & overlap_1_to_1_assoc$soma_non_normal_sig_glu==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_glu_non_normal <- factor(overlap_1_to_1_assoc$sig_glu_non_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_glu_non_normal)
# 
# hist_glu_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_non_normal,fill=sig_glu_non_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("Non-ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_glu_non_normal
# 
# hist_glu <- ggarrange(hist_glu_normal,hist_glu_non_normal,ncol=2,nrow=1)
# 
# hist_glu <- annotate_figure(hist_glu, top = text_grob("Random glucose", face = "bold", size = 28))
# 
# hist_glu
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_glu.png",hist_glu,width=16,height=6,bg="white")
# 
# 
# 
# 
# 
# 
# 
# 
# # smoking
# 
# scatter_smoking_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_smoking, y=soma_normal_es_smoking,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-1,1) +
#   ylim(-1,1) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_smoking_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_smoking,overlap_1_to_1_assoc$soma_normal_es_smoking)
# 
# # smoking non-normal
# 
# scatter_smoking_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_smoking, y=soma_non_normal_es_smoking,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-1,1) +
#   ylim(-1,1) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("Non-ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   geom_abline(linetype = "dashed") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_smoking_non_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_smoking,overlap_1_to_1_assoc$soma_non_normal_es_smoking)
# 
# scatter_smoking <- ggarrange(scatter_smoking_normal,scatter_smoking_non_normal,ncol=2,nrow=1) 
# 
# scatter_smoking <- annotate_figure(scatter_smoking, top = text_grob("Ever regular smoking", face = "bold", size = 28))
# 
# scatter_smoking
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/scatter_smoking.png",scatter_smoking,width=16,height=8, bg = "white")
# 
# # histogram smoking
# 
# overlap_1_to_1_assoc$sig_smoking_normal <- "n.s."
# overlap_1_to_1_assoc$sig_smoking_normal[overlap_1_to_1_assoc$olink_sig_smoking==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_smoking_normal[overlap_1_to_1_assoc$soma_normal_sig_smoking==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_smoking_normal[overlap_1_to_1_assoc$olink_sig_smoking==T & overlap_1_to_1_assoc$soma_normal_sig_smoking==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_smoking_normal <- factor(overlap_1_to_1_assoc$sig_smoking_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_smoking_normal)
# 
# hist_smoking_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_normal,fill=sig_smoking_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_smoking_normal
# 
# overlap_1_to_1_assoc$sig_smoking_non_normal <- "n.s."
# overlap_1_to_1_assoc$sig_smoking_non_normal[overlap_1_to_1_assoc$olink_sig_smoking==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_smoking_non_normal[overlap_1_to_1_assoc$soma_non_normal_sig_smoking==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_smoking_non_normal[overlap_1_to_1_assoc$olink_sig_smoking==T & overlap_1_to_1_assoc$soma_non_normal_sig_smoking==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_smoking_non_normal <- factor(overlap_1_to_1_assoc$sig_smoking_non_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_smoking_non_normal)
# 
# hist_smoking_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_non_normal,fill=sig_smoking_non_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("Non-ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_smoking_non_normal
# 
# hist_smoking <- ggarrange(hist_smoking_normal,hist_smoking_non_normal,ncol=2,nrow=1)
# 
# hist_smoking <- annotate_figure(hist_smoking, top = text_grob("Ever regular smoking", face = "bold", size = 28))
# 
# hist_smoking
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_smoking.png",hist_smoking,width=16,height=6,bg="white")
# 
# 
# 
# 
# 
# 
# 
# # alcohol
# 
# scatter_alcohol_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_alcohol, y=soma_normal_es_alcohol,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-1,1) +
#   ylim(-1,1) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_alcohol_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_alcohol,overlap_1_to_1_assoc$soma_normal_es_alcohol)
# 
# # alcohol non-normal
# 
# scatter_alcohol_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_alcohol, y=soma_non_normal_es_alcohol,stroke=NA)) + 
#   geom_point(size=3,alpha = 0.3) +
#   xlim(-1,1) +
#   ylim(-1,1) +
#   geom_abline(linetype = "dashed") +
#   ggtitle("Non-ANML") +
#   xlab("Olink beta") +
#   ylab("SomaScan beta") +
#   geom_abline(linetype = "dashed") +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(), 
#         axis.line = element_line(colour = "black"),
#         text = element_text(size = 20))
# 
# scatter_alcohol_non_normal
# 
# cor.test(overlap_1_to_1_assoc$olink_es_alcohol,overlap_1_to_1_assoc$soma_non_normal_es_alcohol)
# 
# scatter_alcohol <- ggarrange(scatter_alcohol_normal,scatter_alcohol_non_normal,ncol=2,nrow=1) 
# 
# scatter_alcohol <- annotate_figure(scatter_alcohol, top = text_grob("Current regular vs occasional alcohol drinking", face = "bold", size = 28))
# 
# scatter_alcohol
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/scatter_alcohol.png",scatter_alcohol,width=16,height=8, bg = "white")
# 
# # histogram alcohol
# 
# overlap_1_to_1_assoc$sig_alcohol_normal <- "n.s."
# overlap_1_to_1_assoc$sig_alcohol_normal[overlap_1_to_1_assoc$olink_sig_alcohol==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_alcohol_normal[overlap_1_to_1_assoc$soma_normal_sig_alcohol==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_alcohol_normal[overlap_1_to_1_assoc$olink_sig_alcohol==T & overlap_1_to_1_assoc$soma_normal_sig_alcohol==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_alcohol_normal <- factor(overlap_1_to_1_assoc$sig_alcohol_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_alcohol_normal)
# 
# hist_alcohol_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_normal,fill=sig_alcohol_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_alcohol_normal
# 
# overlap_1_to_1_assoc$sig_alcohol_non_normal <- "n.s."
# overlap_1_to_1_assoc$sig_alcohol_non_normal[overlap_1_to_1_assoc$olink_sig_alcohol==T] <- "sig. Olink"
# overlap_1_to_1_assoc$sig_alcohol_non_normal[overlap_1_to_1_assoc$soma_non_normal_sig_alcohol==T] <- "sig. SomaScan"
# overlap_1_to_1_assoc$sig_alcohol_non_normal[overlap_1_to_1_assoc$olink_sig_alcohol==T & overlap_1_to_1_assoc$soma_non_normal_sig_alcohol==T] <- "sig. both"
# overlap_1_to_1_assoc$sig_alcohol_non_normal <- factor(overlap_1_to_1_assoc$sig_alcohol_non_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
# table(overlap_1_to_1_assoc$sig_alcohol_non_normal)
# 
# hist_alcohol_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_non_normal,fill=sig_alcohol_non_normal)) + 
#   geom_histogram(color="black") +
#   xlim(-0.4,1) +
#   ylim(0,400) +
#   ggtitle("Non-ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
#   scale_fill_discrete(name = "Significance") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA))
# 
# hist_alcohol_non_normal
# 
# hist_alcohol <- ggarrange(hist_alcohol_normal,hist_alcohol_non_normal,ncol=2,nrow=1)
# 
# hist_alcohol <- annotate_figure(hist_alcohol, top = text_grob("Current regular vs occasional alcohol drinking", face = "bold", size = 28))
# 
# hist_alcohol
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_alcohol.png",hist_alcohol,width=16,height=6,bg="white")
# 
# 
