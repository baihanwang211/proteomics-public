rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(ggplot2)
library(ggpubr)
library(forcats)
library(egg)
library(ckbplotr)
library(scales)
library(ggpubr)
library(egg)
library(ggpattern)
library(VennDiagram)
library(RColorBrewer)

# number of pQTLs

pqtl <- read.delim("overlap_pqtls.txt")
pqtl <- pqtl[,-c(2,3)]
names(pqtl) <- c("somascan_id","olink_cis","olink_trans","soma_normal_cis","soma_normal_trans","soma_non_normal_cis","soma_non_normal_trans")

# load somamer matched to multiple olink

overlap_multiple_olink <- read.csv("overlap_cor_one_soma_to_multiple_olink.csv")
overlap_multiple_olink_pqtl <- merge(overlap_multiple_olink,pqtl,by=c("somascan_id"))
length(unique(overlap_multiple_olink_pqtl$somascan_id))
length(unique(overlap_multiple_olink_pqtl$olink_id))

# numbers

multiple_olink_pqtl_tf_count <- data.frame(Platform = c("OLINK","SomaScan-ANML","SomaScan-non-ANML","OLINK","SomaScan-ANML","SomaScan-non-ANML"),
  pQTL = c("Cis-pQTL","Cis-pQTL","Cis-pQTL","Trans-pQTL","Trans-pQTL","Trans-pQTL"),
  Frequency = c(length(unique(overlap_multiple_olink_pqtl$olink_id[overlap_multiple_olink_pqtl$olink_cis>0])),
                length(unique(overlap_multiple_olink_pqtl$somascan_id[overlap_multiple_olink_pqtl$soma_normal_cis>0])),
                length(unique(overlap_multiple_olink_pqtl$somascan_id[overlap_multiple_olink_pqtl$soma_non_normal_cis>0])),
                length(unique(overlap_multiple_olink_pqtl$olink_id[overlap_multiple_olink_pqtl$olink_trans>0])),
                length(unique(overlap_multiple_olink_pqtl$somascan_id[overlap_multiple_olink_pqtl$soma_normal_trans>0])),
                length(unique(overlap_multiple_olink_pqtl$somascan_id[overlap_multiple_olink_pqtl$soma_non_normal_trans>0]))))


multiple_olink_pqtl_count <- data.frame(Platform = c("OLINK","SomaScan-ANML","SomaScan-non-ANML","OLINK","SomaScan-ANML","SomaScan-non-ANML"),
                                        pQTL = c("Cis-pQTL","Cis-pQTL","Cis-pQTL","Trans-pQTL","Trans-pQTL","Trans-pQTL"),
                                        Frequency = c(sum(overlap_multiple_olink_pqtl$olink_cis[!duplicated(overlap_multiple_olink_pqtl$olink_id)]),
                                                      sum(overlap_multiple_olink_pqtl$soma_normal_cis[!duplicated(overlap_multiple_olink_pqtl$somascan_id)]),
                                                      sum(overlap_multiple_olink_pqtl$soma_non_normal_cis[!duplicated(overlap_multiple_olink_pqtl$somascan_id)]),
                                                      sum(overlap_multiple_olink_pqtl$olink_trans[!duplicated(overlap_multiple_olink_pqtl$olink_id)]),
                                                      sum(overlap_multiple_olink_pqtl$soma_normal_trans[!duplicated(overlap_multiple_olink_pqtl$somascan_id)]),
                                                      sum(overlap_multiple_olink_pqtl$soma_non_normal_trans[!duplicated(overlap_multiple_olink_pqtl$somascan_id)])))


multiple_olink_pqtl_tf_count_add <- multiple_olink_pqtl_tf_count
multiple_olink_pqtl_tf_count_add$Frequency <- multiple_olink_pqtl_count$Frequency - multiple_olink_pqtl_tf_count$Frequency

multiple_olink_pqtl_tf_count$cat<- "original"
multiple_olink_pqtl_tf_count_add$cat<- "additional"


multiple_olink_pqtl_count_all <- rbind(multiple_olink_pqtl_tf_count,multiple_olink_pqtl_tf_count_add)

multiple_olink_pqtl_count_all$Frequency_all <- NA
multiple_olink_pqtl_count_all$Frequency_all[multiple_olink_pqtl_count_all$cat=="additional"] <- multiple_olink_pqtl_count$Frequency

multiple_olink_pqtl_count_all$Frequency_shade <- NA
multiple_olink_pqtl_count_all$Frequency_shade[multiple_olink_pqtl_count_all$cat=="original"] <- multiple_olink_pqtl_count_all$Frequency[multiple_olink_pqtl_count_all$cat=="original"]

# plot stacked bar together

hex <- hue_pal()(3)

plot_multiple_olink_pqtl_hit_all <- ggplot(multiple_olink_pqtl_count_all,aes(y=Frequency, x=fct_rev(Platform), fill=fct_rev(Platform), pattern=cat)) +
  geom_bar_pattern(stat = "identity",
                   position = "stack",
                   colour = 'black',
                   pattern_fill = "black",
                   pattern_spacing = 0.03,
                   pattern_frequency = 5,
                   pattern_angle = 45,
                   pattern_density=0.01) +
  ylim(0,25) +
  geom_text(aes(label=Frequency_all),position="stack",hjust=-0.1) +
  geom_label(aes(label=Frequency_shade),position=position_stack(vjust = 0.5),label.size = NA) +
  scale_fill_manual(values = rev(hex)) +
  scale_pattern_manual(values=c("none","stripe")) +
  xlab("") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(legend.position = "none") +
  coord_flip() +
  facet_grid(pQTL~.)

plot_multiple_olink_pqtl_hit_all

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_multiple_olink_pqtl_hit_all.png",plot_multiple_olink_pqtl_hit_all,width=9,height=6)


# load olink matched to multiple somamers

overlap_multiple_soma <- read.csv("overlap_cor_one_olink_to_multiple_soma.csv")
overlap_multiple_soma_pqtl <- merge(overlap_multiple_soma,pqtl,by=c("somascan_id"))
length(unique(overlap_multiple_soma_pqtl$somascan_id))
length(unique(overlap_multiple_soma_pqtl$olink_id))

# numbers

multiple_soma_pqtl_tf_count <- data.frame(Platform = c("OLINK","SomaScan-ANML","SomaScan-non-ANML","OLINK","SomaScan-ANML","SomaScan-non-ANML"),
                                           pQTL = c("Cis-pQTL","Cis-pQTL","Cis-pQTL","Trans-pQTL","Trans-pQTL","Trans-pQTL"),
                                           Frequency = c(length(unique(overlap_multiple_soma_pqtl$olink_id[overlap_multiple_soma_pqtl$olink_cis>0])),
                                                         length(unique(overlap_multiple_soma_pqtl$somascan_id[overlap_multiple_soma_pqtl$soma_normal_cis>0])),
                                                         length(unique(overlap_multiple_soma_pqtl$somascan_id[overlap_multiple_soma_pqtl$soma_non_normal_cis>0])),
                                                         length(unique(overlap_multiple_soma_pqtl$olink_id[overlap_multiple_soma_pqtl$olink_trans>0])),
                                                         length(unique(overlap_multiple_soma_pqtl$somascan_id[overlap_multiple_soma_pqtl$soma_normal_trans>0])),
                                                         length(unique(overlap_multiple_soma_pqtl$somascan_id[overlap_multiple_soma_pqtl$soma_non_normal_trans>0]))))



multiple_soma_pqtl_count <- data.frame(Platform = c("OLINK","SomaScan-ANML","SomaScan-non-ANML","OLINK","SomaScan-ANML","SomaScan-non-ANML"),
                                        pQTL = c("Cis-pQTL","Cis-pQTL","Cis-pQTL","Trans-pQTL","Trans-pQTL","Trans-pQTL"),
                                        Frequency = c(sum(overlap_multiple_soma_pqtl$olink_cis[!duplicated(overlap_multiple_soma_pqtl$olink_id)]),
                                                      sum(overlap_multiple_soma_pqtl$soma_normal_cis[!duplicated(overlap_multiple_soma_pqtl$somascan_id)]),
                                                      sum(overlap_multiple_soma_pqtl$soma_non_normal_cis[!duplicated(overlap_multiple_soma_pqtl$somascan_id)]),
                                                      sum(overlap_multiple_soma_pqtl$olink_trans[!duplicated(overlap_multiple_soma_pqtl$olink_id)]),
                                                      sum(overlap_multiple_soma_pqtl$soma_normal_trans[!duplicated(overlap_multiple_soma_pqtl$somascan_id)]),
                                                      sum(overlap_multiple_soma_pqtl$soma_non_normal_trans[!duplicated(overlap_multiple_soma_pqtl$somascan_id)])))


multiple_soma_pqtl_tf_count_add <- multiple_soma_pqtl_tf_count
multiple_soma_pqtl_tf_count_add$Frequency <- multiple_soma_pqtl_count$Frequency - multiple_soma_pqtl_tf_count$Frequency

multiple_soma_pqtl_tf_count$cat<- "original"
multiple_soma_pqtl_tf_count_add$cat<- "additional"


multiple_soma_pqtl_count_all <- rbind(multiple_soma_pqtl_tf_count,multiple_soma_pqtl_tf_count_add)

multiple_soma_pqtl_count_all$Frequency_all <- NA
multiple_soma_pqtl_count_all$Frequency_all[multiple_soma_pqtl_count_all$cat=="additional"] <- multiple_soma_pqtl_count$Frequency

multiple_soma_pqtl_count_all$Frequency_shade <- NA
multiple_soma_pqtl_count_all$Frequency_shade[multiple_soma_pqtl_count_all$cat=="original"] <- multiple_soma_pqtl_count_all$Frequency[multiple_soma_pqtl_count_all$cat=="original"]

# plot stacked bar together

hex <- hue_pal()(3)

plot_multiple_soma_pqtl_hit_all <- ggplot(multiple_soma_pqtl_count_all,aes(y=Frequency, x=fct_rev(Platform), fill=fct_rev(Platform), pattern=cat)) +
  geom_bar_pattern(stat = "identity",
                   position = "stack",
                   colour = 'black',
                   pattern_fill = "black",
                   pattern_spacing = 0.03,
                   pattern_frequency = 5,
                   pattern_angle = 45,
                   pattern_density=0.01) +
  ylim(0,1200) +
  geom_text(aes(label=Frequency_all),position="stack",hjust=-0.1) +
  geom_label(aes(label=Frequency_shade),position=position_stack(vjust = 0.5),label.size = NA) +
  scale_fill_manual(values = rev(hex)) +
  scale_pattern_manual(values=c("none","stripe")) +
  xlab("") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(legend.position = "none") +
  coord_flip() +
  facet_grid(pQTL~.)

plot_multiple_soma_pqtl_hit_all

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_multiple_soma_pqtl_hit_all.png",plot_multiple_soma_pqtl_hit_all,width=9,height=6)




###### load coloc results and merge with somamers matched to multiple olink

coloc_normal <- read.delim("coloc_anml_non_one_to_one.txt")
names(coloc_normal)
coloc_normal <- coloc_normal[,c("olink","soma","shared_cis","shared_trans")]
names(coloc_normal) <- c("olink_id","somascan_id","coloc_normal_cis","coloc_normal_trans")
coloc_normal <- coloc_normal[!duplicated(coloc_normal[,c("olink_id","somascan_id")]),]
overlap_multiple_olink_coloc <- merge(overlap_multiple_olink_pqtl,coloc_normal,by=c("somascan_id","olink_id"),all.x=T)

coloc_non_normal <- read.delim("coloc_nanml_non_one_to_one.txt")
names(coloc_non_normal)
coloc_non_normal <- coloc_non_normal[,c("olink","soma","shared_cis","shared_trans")]
names(coloc_non_normal) <- c("olink_id","somascan_id","coloc_non_normal_cis","coloc_non_normal_trans")
coloc_non_normal <- coloc_non_normal[!duplicated(coloc_non_normal[,c("olink_id","somascan_id")]),]
overlap_multiple_olink_coloc <- merge(overlap_multiple_olink_coloc,coloc_non_normal,by=c("somascan_id","olink_id"),all.x=T)

write.csv(overlap_multiple_olink_coloc,"overlap_multiple_olink_coloc.csv", quote=F, row.names=F)

# replace NA with 0

overlap_multiple_olink_coloc[is.na(overlap_multiple_olink_coloc)] <- 0

# pairwise comparison (% coloc)
# normal

overlap_multiple_olink_coloc_percent_normal <- data.frame(somascan_id=unique(overlap_multiple_olink_coloc$somascan_id))

for (i in 1:nrow(overlap_multiple_olink_coloc_percent_normal)) {
  df <- overlap_multiple_olink_coloc[overlap_multiple_olink_coloc$somascan_id==overlap_multiple_olink_coloc_percent_normal$somascan_id[i],]
  overlap_multiple_olink_coloc_percent_normal$olink[i] <- sum(df$olink_cis!=0)
  overlap_multiple_olink_coloc_percent_normal$soma[i] <- df$soma_normal_cis!=0
  overlap_multiple_olink_coloc_percent_normal$coloc[i] <- sum(df$coloc_normal_cis!=0)
  overlap_multiple_olink_coloc_percent_normal$total[i] <- nrow(df)
  # overlap_multiple_olink_coloc_percent_normal$percent[i] <- sum(df$coloc_normal_cis!=0)/nrow(df)
}

tb2 <- table(overlap_multiple_olink_coloc_percent_normal$total,overlap_multiple_olink_coloc_percent_normal$olink)
tb1 <- table(overlap_multiple_olink_coloc_percent_normal$total,overlap_multiple_olink_coloc_percent_normal$soma)
n_pqtl_multiple_olink_normal <- cbind(tb1,tb2)
n_pqtl_multiple_olink_normal

table(overlap_multiple_olink_coloc_percent_normal$total,overlap_multiple_olink_coloc_percent_normal$coloc)

# non-normal

overlap_multiple_olink_coloc_percent_non_normal <- data.frame(somascan_id=unique(overlap_multiple_olink_coloc$somascan_id))

for (i in 1:nrow(overlap_multiple_olink_coloc_percent_non_normal)) {
  df <- overlap_multiple_olink_coloc[overlap_multiple_olink_coloc$somascan_id==overlap_multiple_olink_coloc_percent_non_normal$somascan_id[i],]
  overlap_multiple_olink_coloc_percent_non_normal$olink[i] <- sum(df$olink_cis!=0)
  overlap_multiple_olink_coloc_percent_non_normal$soma[i] <- df$soma_normal_cis!=0
  overlap_multiple_olink_coloc_percent_non_normal$coloc[i] <- sum(df$coloc_non_normal_cis!=0)
  overlap_multiple_olink_coloc_percent_non_normal$total[i] <- nrow(df)
  # overlap_multiple_olink_coloc_percent_non_normal$percent[i] <- sum(df$coloc_non_normal_cis!=0)/nrow(df)
}

tb2 <- table(overlap_multiple_olink_coloc_percent_non_normal$total,overlap_multiple_olink_coloc_percent_non_normal$olink)
tb1 <- table(overlap_multiple_olink_coloc_percent_non_normal$total,overlap_multiple_olink_coloc_percent_non_normal$soma)
n_pqtl_multiple_olink_non_normal <- cbind(tb1,tb2)
n_pqtl_multiple_olink_non_normal

table(overlap_multiple_olink_coloc_percent_non_normal$total,overlap_multiple_olink_coloc_percent_non_normal$coloc)

###### merge with olink reagents matched to multiple somamers

overlap_multiple_soma_coloc <- merge(overlap_multiple_soma_pqtl,coloc_normal,by=c("somascan_id","olink_id"),all.x=T)
overlap_multiple_soma_coloc <- merge(overlap_multiple_soma_coloc,coloc_non_normal,by=c("somascan_id","olink_id"),all.x=T)

length(unique(overlap_multiple_soma_coloc$olink_id[overlap_multiple_soma_coloc$coloc_normal_cis!=0]))

length(unique(overlap_multiple_soma_coloc$olink_id[overlap_multiple_soma_coloc$coloc_non_normal_cis!=0]))

write.csv(overlap_multiple_soma_coloc,"overlap_multiple_soma_coloc.csv", quote=F, row.names=F)

# replace NA with 0

overlap_multiple_soma_coloc[is.na(overlap_multiple_soma_coloc)] <- 0

# pairwise comparison (% coloc)
# normal

overlap_multiple_soma_coloc_percent_normal <- data.frame(olink_id=unique(overlap_multiple_soma_coloc$olink_id))

for (i in 1:nrow(overlap_multiple_soma_coloc_percent_normal)) {
  df <- overlap_multiple_soma_coloc[overlap_multiple_soma_coloc$olink_id==overlap_multiple_soma_coloc_percent_normal$olink_id[i],]
  overlap_multiple_soma_coloc_percent_normal$olink[i] <- df$olink_cis!=0
  overlap_multiple_soma_coloc_percent_normal$soma[i] <- sum(df$soma_normal_cis!=0)
  overlap_multiple_soma_coloc_percent_normal$coloc[i] <- sum(df$coloc_normal_cis!=0)
  overlap_multiple_soma_coloc_percent_normal$total[i] <- nrow(df)
  # overlap_multiple_soma_coloc_percent_normal$percent[i] <- sum(df$coloc_normal_cis!=0)/nrow(df)
}

tb1 <- table(overlap_multiple_soma_coloc_percent_normal$total,overlap_multiple_soma_coloc_percent_normal$olink)
tb2 <- table(overlap_multiple_soma_coloc_percent_normal$total,overlap_multiple_soma_coloc_percent_normal$soma)
n_pqtl_multiple_soma_normal <- cbind(tb1,tb2)
n_pqtl_multiple_soma_normal

table(overlap_multiple_soma_coloc_percent_normal$total,overlap_multiple_soma_coloc_percent_normal$coloc)

# non-normal

overlap_multiple_soma_coloc_percent_non_normal <- data.frame(olink_id=unique(overlap_multiple_soma_coloc$olink_id))

for (i in 1:nrow(overlap_multiple_soma_coloc_percent_non_normal)) {
  df <- overlap_multiple_soma_coloc[overlap_multiple_soma_coloc$olink_id==overlap_multiple_soma_coloc_percent_non_normal$olink_id[i],]
  overlap_multiple_soma_coloc_percent_non_normal$olink[i] <- df$olink_cis!=0
  overlap_multiple_soma_coloc_percent_non_normal$soma[i] <- sum(df$soma_non_normal_cis!=0)
  overlap_multiple_soma_coloc_percent_non_normal$coloc[i] <- sum(df$coloc_non_normal_cis!=0)
  overlap_multiple_soma_coloc_percent_non_normal$total[i] <- nrow(df)
  # overlap_multiple_soma_coloc_percent_non_normal$percent[i] <- sum(df$coloc_non_normal_cis!=0)/nrow(df)
}

tb1 <- table(overlap_multiple_soma_coloc_percent_non_normal$total,overlap_multiple_soma_coloc_percent_non_normal$olink)
tb2 <- table(overlap_multiple_soma_coloc_percent_non_normal$total,overlap_multiple_soma_coloc_percent_non_normal$soma)
n_pqtl_multiple_soma_non_normal <- cbind(tb1,tb2)
n_pqtl_multiple_soma_non_normal

table(overlap_multiple_soma_coloc_percent_non_normal$total,overlap_multiple_soma_coloc_percent_non_normal$coloc)


# 
# 
# # histogram
# 
# ## plot shared only using shade
# 
# overlap_multiple_olink_coloc$coloc_normal_cis_tf <- F
# overlap_multiple_olink_coloc$coloc_normal_cis_tf[overlap_multiple_olink_coloc$coloc_normal_cis!=0] <- T
# table(overlap_multiple_olink_coloc$coloc_normal_cis_tf)
# 
# median_rho_all_normal <- format(round(median(overlap_multiple_olink_coloc$rho_olink_soma_normal),2),nsmall=2)
# median_rho_coloc_normal <- format(round(median(overlap_multiple_olink_coloc$rho_olink_soma_normal[overlap_multiple_olink_coloc$coloc_normal_cis_tf==T]),2),nsmall=2)
# 
# hist_pqtl_multiple_olink_coloc_normal <- ggplot(overlap_multiple_olink_coloc, aes(x=rho_olink_soma_normal,pattern=coloc_normal_cis_tf)) + 
#   geom_histogram_pattern(binwidth = 0.05,
#                          boundary=0,
#                          pattern_fill = "black", 
#                          fill    = 'white',
#                          colour  = 'black',
#                          pattern_spacing = 0.015,
#                          pattern_frequency = 5, 
#                          pattern_angle = 45,
#                          pattern_density=0.01) +
#   scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
#   ylim(0,300) +
#   ggtitle("OLINK vs SomaScan-ANML", subtitle ="") +
#   xlab("Spearman's rho") + 
#   ylab("Frequency") +
#   labs(fill="Shared") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
#   theme(text = element_text(size = 12)) +
#   scale_pattern_manual(values=c("none","stripe")) +
#   theme(legend.position = "none") +
#   annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_normal, "(all)"), size = 12/.pt, hjust = 0) +
#   annotate("text", x = -0.4, y = 285, label = paste("Median rho =", median_rho_coloc_normal, "(coloc cis)"), size = 12/.pt, hjust = 0)
# 
# hist_pqtl_multiple_olink_coloc_normal
# 
# # non-normal
# 
# overlap_multiple_olink_coloc$coloc_non_normal_cis_tf <- F
# overlap_multiple_olink_coloc$coloc_non_normal_cis_tf[overlap_multiple_olink_coloc$coloc_non_normal_cis!=0] <- T
# table(overlap_multiple_olink_coloc$coloc_non_normal_cis_tf)
# 
# median_rho_all_non_normal <- format(round(median(overlap_multiple_olink_coloc$rho_olink_soma_non_normal),2),nsmall=2)
# median_rho_coloc_non_normal <- format(round(median(overlap_multiple_olink_coloc$rho_olink_soma_non_normal[overlap_multiple_olink_coloc$coloc_non_normal_cis_tf==T]),2),nsmall=2)
# 
# hist_pqtl_multiple_olink_coloc_non_normal <- ggplot(overlap_multiple_olink_coloc, aes(x=rho_olink_soma_non_normal,pattern=coloc_non_normal_cis_tf)) + 
#   geom_histogram_pattern(binwidth = 0.05,
#                          boundary=0,
#                          pattern_fill = "black", 
#                          fill    = 'white',
#                          colour  = 'black',
#                          pattern_spacing = 0.015,
#                          pattern_frequency = 5, 
#                          pattern_angle = 45,
#                          pattern_density=0.01) +
#   scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
#   ylim(0,300) +
#   ggtitle("OLINK vs SomaScan-non-ANML", subtitle ="") +
#   xlab("Spearman's rho") + 
#   ylab("Frequency") +
#   labs(fill="Shared") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
#   theme(text = element_text(size = 12)) +
#   scale_pattern_manual(values=c("none","stripe")) +
#   theme(legend.position = "none") +
#   annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_non_normal, "(all)"), size = 12/.pt, hjust = 0) +
#   annotate("text", x = -0.4, y = 285, label = paste("Median rho =", median_rho_coloc_non_normal, "(coloc cis)"), size = 12/.pt, hjust = 0)
# 
# hist_pqtl_multiple_olink_coloc_non_normal
# 
# hist_pqtl_multiple_olink_coloc <- ggarrange(hist_pqtl_multiple_olink_coloc_normal,hist_pqtl_multiple_olink_coloc_non_normal,ncol=2,nrow=1)
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_pqtl_multiple_olink_coloc.png",hist_pqtl_multiple_olink_coloc,width=14,height=6)
# 
