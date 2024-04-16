### this script loads pQTL and colocalisation results from GWAS and plots results for non-one-to-one matched proteins

rm(list = ls())

setwd("")

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

overlap_coloc <- read.csv("overlap_coloc.csv")

# load somamer matched to multiple olink

overlap_multiple_olink <- read.csv("overlap_cor_one_soma_to_multiple_olink.csv")
overlap_multiple_olink_pqtl <- merge(overlap_multiple_olink,overlap_coloc,by=c("olink_id","uniprot_id","somascan_id"))
length(unique(overlap_multiple_olink_pqtl$somascan_id))
length(unique(overlap_multiple_olink_pqtl$olink_id))

# numbers

multiple_olink_pqtl_tf_count <- data.frame(Platform = c("Olink","SomaScan-ANML","SomaScan-non-ANML","Olink","SomaScan-ANML","SomaScan-non-ANML"),
  pQTL = c("Cis-pQTL","Cis-pQTL","Cis-pQTL","Trans-pQTL","Trans-pQTL","Trans-pQTL"),
  Frequency = c(length(unique(overlap_multiple_olink_pqtl$olink_id[overlap_multiple_olink_pqtl$olink_cis>0])),
                length(unique(overlap_multiple_olink_pqtl$somascan_id[overlap_multiple_olink_pqtl$soma_normal_cis>0])),
                length(unique(overlap_multiple_olink_pqtl$somascan_id[overlap_multiple_olink_pqtl$soma_non_normal_cis>0])),
                length(unique(overlap_multiple_olink_pqtl$olink_id[overlap_multiple_olink_pqtl$olink_trans>0])),
                length(unique(overlap_multiple_olink_pqtl$somascan_id[overlap_multiple_olink_pqtl$soma_normal_trans>0])),
                length(unique(overlap_multiple_olink_pqtl$somascan_id[overlap_multiple_olink_pqtl$soma_non_normal_trans>0]))))

multiple_olink_pqtl_count <- data.frame(Platform = c("Olink","SomaScan-ANML","SomaScan-non-ANML","Olink","SomaScan-ANML","SomaScan-non-ANML"),
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

ggsave("",plot_multiple_olink_pqtl_hit_all,width=9,height=6)


# load olink matched to multiple somamers

overlap_multiple_soma <- read.csv("overlap_cor_one_olink_to_multiple_soma.csv")
overlap_multiple_soma_pqtl <- merge(overlap_multiple_soma,overlap_coloc,by=c("olink_id","uniprot_id","somascan_id"))
length(unique(overlap_multiple_soma_pqtl$somascan_id))
length(unique(overlap_multiple_soma_pqtl$olink_id))

# numbers

multiple_soma_pqtl_tf_count <- data.frame(Platform = c("Olink","SomaScan-ANML","SomaScan-non-ANML","Olink","SomaScan-ANML","SomaScan-non-ANML"),
                                           pQTL = c("Cis-pQTL","Cis-pQTL","Cis-pQTL","Trans-pQTL","Trans-pQTL","Trans-pQTL"),
                                           Frequency = c(length(unique(overlap_multiple_soma_pqtl$olink_id[overlap_multiple_soma_pqtl$olink_cis>0])),
                                                         length(unique(overlap_multiple_soma_pqtl$somascan_id[overlap_multiple_soma_pqtl$soma_normal_cis>0])),
                                                         length(unique(overlap_multiple_soma_pqtl$somascan_id[overlap_multiple_soma_pqtl$soma_non_normal_cis>0])),
                                                         length(unique(overlap_multiple_soma_pqtl$olink_id[overlap_multiple_soma_pqtl$olink_trans>0])),
                                                         length(unique(overlap_multiple_soma_pqtl$somascan_id[overlap_multiple_soma_pqtl$soma_normal_trans>0])),
                                                         length(unique(overlap_multiple_soma_pqtl$somascan_id[overlap_multiple_soma_pqtl$soma_non_normal_trans>0]))))



multiple_soma_pqtl_count <- data.frame(Platform = c("Olink","SomaScan-ANML","SomaScan-non-ANML","Olink","SomaScan-ANML","SomaScan-non-ANML"),
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

ggsave("",plot_multiple_soma_pqtl_hit_all,width=9,height=6)




###### load coloc results and merge with somamers matched to multiple olink


# pairwise comparison (% coloc)
# normal

overlap_multiple_olink_pqtl_percent_normal <- data.frame(somascan_id=unique(overlap_multiple_olink_pqtl$somascan_id))

for (i in 1:nrow(overlap_multiple_olink_pqtl_percent_normal)) {
  df <- overlap_multiple_olink_pqtl[overlap_multiple_olink_pqtl$somascan_id==overlap_multiple_olink_pqtl_percent_normal$somascan_id[i],]
  overlap_multiple_olink_pqtl_percent_normal$olink_cis[i] <- sum(df$olink_cis!=0)
  overlap_multiple_olink_pqtl_percent_normal$soma_cis[i] <- df$soma_normal_cis!=0
  overlap_multiple_olink_pqtl_percent_normal$olink_trans[i] <- sum(df$olink_trans!=0)
  overlap_multiple_olink_pqtl_percent_normal$soma_trans[i] <- df$soma_normal_trans!=0
  overlap_multiple_olink_pqtl_percent_normal$coloc_cis[i] <- sum(df$coloc_normal_cis!=0)
  overlap_multiple_olink_pqtl_percent_normal$coloc_trans[i] <- sum(df$coloc_normal_trans!=0)
  overlap_multiple_olink_pqtl_percent_normal$total[i] <- nrow(df)
}

tb2 <- table(overlap_multiple_olink_pqtl_percent_normal$total,overlap_multiple_olink_pqtl_percent_normal$olink_cis)
tb1 <- table(overlap_multiple_olink_pqtl_percent_normal$total,overlap_multiple_olink_pqtl_percent_normal$soma_cis)
n_pqtl_multiple_olink_normal <- cbind(tb1,tb2)
n_pqtl_multiple_olink_normal

table(overlap_multiple_olink_pqtl_percent_normal$total,overlap_multiple_olink_pqtl_percent_normal$coloc_cis)
table(overlap_multiple_olink_pqtl_percent_normal$total,overlap_multiple_olink_pqtl_percent_normal$coloc_trans)

# non-normal

overlap_multiple_olink_pqtl_percent_non_normal <- data.frame(somascan_id=unique(overlap_multiple_olink_pqtl$somascan_id))

for (i in 1:nrow(overlap_multiple_olink_pqtl_percent_non_normal)) {
  df <- overlap_multiple_olink_pqtl[overlap_multiple_olink_pqtl$somascan_id==overlap_multiple_olink_pqtl_percent_non_normal$somascan_id[i],]
  overlap_multiple_olink_pqtl_percent_non_normal$olink_cis[i] <- sum(df$olink_cis!=0)
  overlap_multiple_olink_pqtl_percent_non_normal$soma_cis[i] <- df$soma_normal_cis!=0
  overlap_multiple_olink_pqtl_percent_non_normal$olink_trans[i] <- sum(df$olink_trans!=0)
  overlap_multiple_olink_pqtl_percent_non_normal$soma_trans[i] <- df$soma_normal_trans!=0
  overlap_multiple_olink_pqtl_percent_non_normal$coloc_cis[i] <- sum(df$coloc_non_normal_cis!=0)
  overlap_multiple_olink_pqtl_percent_non_normal$coloc_trans[i] <- sum(df$coloc_non_normal_trans!=0)
  overlap_multiple_olink_pqtl_percent_non_normal$total[i] <- nrow(df)
  # overlap_multiple_olink_pqtl_percent_non_normal$percent[i] <- sum(df$coloc_non_normal_cis!=0)/nrow(df)
}

tb2 <- table(overlap_multiple_olink_pqtl_percent_non_normal$total,overlap_multiple_olink_pqtl_percent_non_normal$olink_cis)
tb1 <- table(overlap_multiple_olink_pqtl_percent_non_normal$total,overlap_multiple_olink_pqtl_percent_non_normal$soma_cis)
n_pqtl_multiple_olink_non_normal <- cbind(tb1,tb2)
n_pqtl_multiple_olink_non_normal

table(overlap_multiple_olink_pqtl_percent_non_normal$total,overlap_multiple_olink_pqtl_percent_non_normal$coloc_cis)
table(overlap_multiple_olink_pqtl_percent_non_normal$total,overlap_multiple_olink_pqtl_percent_non_normal$coloc_trans)

###### merge with olink reagents matched to multiple somamers

# pairwise comparison (% coloc)
# normal

overlap_multiple_soma_pqtl_percent_normal <- data.frame(olink_id=unique(overlap_multiple_soma_pqtl$olink_id))

for (i in 1:nrow(overlap_multiple_soma_pqtl_percent_normal)) {
  df <- overlap_multiple_soma_pqtl[overlap_multiple_soma_pqtl$olink_id==overlap_multiple_soma_pqtl_percent_normal$olink_id[i],]
  # how many have olink cis
  overlap_multiple_soma_pqtl_percent_normal$olink_cis[i] <- df$olink_cis[1]!=0
  # how many have soma cis
  overlap_multiple_soma_pqtl_percent_normal$soma_cis[i] <- sum(df$soma_normal_cis!=0)
  # how many colocalise
  overlap_multiple_soma_pqtl_percent_normal$coloc_cis[i] <- sum(df$coloc_normal_cis!=0)
  # how many have olink trans
  overlap_multiple_soma_pqtl_percent_normal$olink_trans[i] <- df$olink_trans[1]!=0
  # how many have soma trans
  overlap_multiple_soma_pqtl_percent_normal$soma_trans[i] <- sum(df$soma_normal_trans!=0)
  # how many colocalise
  overlap_multiple_soma_pqtl_percent_normal$coloc_trans[i] <- sum(df$coloc_normal_trans!=0)
  # how many pairs in total/how many SOMAmers are matched to one olink reagent
  overlap_multiple_soma_pqtl_percent_normal$total[i] <- nrow(df)
  # overlap_multiple_soma_pqtl_percent_normal$percent[i] <- sum(df$coloc_normal_cis!=0)/nrow(df)
}

# how many have olink cis, grouped by number of matched SOMAmers
table(overlap_multiple_soma_pqtl_percent_normal$total,overlap_multiple_soma_pqtl_percent_normal$olink_cis)

# how many have soma cis, grouped by number of matched SOMAmers
table(overlap_multiple_soma_pqtl_percent_normal$total,overlap_multiple_soma_pqtl_percent_normal$soma_cis)

# how many have coloc cis-pqtls, grouped by matched SOMAmers
table(overlap_multiple_soma_pqtl_percent_normal$total,overlap_multiple_soma_pqtl_percent_normal$coloc_cis)
table(overlap_multiple_soma_pqtl_percent_normal$total,overlap_multiple_soma_pqtl_percent_normal$coloc_trans)

sum(overlap_multiple_soma_pqtl_percent_normal$coloc_cis>0)
sum(overlap_multiple_soma_pqtl_percent_normal$coloc_trans>0)

# how many have matched SOMAmers
table(overlap_multiple_soma_pqtl_percent_normal$total)

# non-normal

overlap_multiple_soma_pqtl_percent_non_normal <- data.frame(olink_id=unique(overlap_multiple_soma_pqtl$olink_id))

for (i in 1:nrow(overlap_multiple_soma_pqtl_percent_non_normal)) {
  df <- overlap_multiple_soma_pqtl[overlap_multiple_soma_pqtl$olink_id==overlap_multiple_soma_pqtl_percent_non_normal$olink_id[i],]
  # how many have olink cis
  overlap_multiple_soma_pqtl_percent_non_normal$olink_cis[i] <- df$olink_cis[1]!=0
  # how many have soma cis
  overlap_multiple_soma_pqtl_percent_non_normal$soma_cis[i] <- sum(df$soma_non_normal_cis!=0)
  # how many colocalise
  overlap_multiple_soma_pqtl_percent_non_normal$coloc_cis[i] <- sum(df$coloc_non_normal_cis!=0)
  # how many have olink trans
  overlap_multiple_soma_pqtl_percent_non_normal$olink_trans[i] <- df$olink_trans[1]!=0
  # how many have soma trans
  overlap_multiple_soma_pqtl_percent_non_normal$soma_trans[i] <- sum(df$soma_non_normal_trans!=0)
  # how many colocalise
  overlap_multiple_soma_pqtl_percent_non_normal$coloc_trans[i] <- sum(df$coloc_non_normal_trans!=0)
  # how many pairs in total/how many SOMAmers are matched to one olink reagent
  overlap_multiple_soma_pqtl_percent_non_normal$total[i] <- nrow(df)
  # overlap_multiple_soma_pqtl_percent_non_normal$percent[i] <- sum(df$coloc_non_normal_cis!=0)/nrow(df)
}

# how many have olink cis, grouped by number of matched SOMAmers
table(overlap_multiple_soma_pqtl_percent_non_normal$total,overlap_multiple_soma_pqtl_percent_non_normal$olink_cis)

# how many have soma cis, grouped by number of matched SOMAmers
table(overlap_multiple_soma_pqtl_percent_non_normal$total,overlap_multiple_soma_pqtl_percent_non_normal$soma_cis)

# how many have coloc cis-pqtls, grouped by matched SOMAmers
table(overlap_multiple_soma_pqtl_percent_non_normal$total,overlap_multiple_soma_pqtl_percent_non_normal$coloc_cis)
table(overlap_multiple_soma_pqtl_percent_non_normal$total,overlap_multiple_soma_pqtl_percent_non_normal$coloc_trans)

sum(overlap_multiple_soma_pqtl_percent_non_normal$coloc_cis>0)
sum(overlap_multiple_soma_pqtl_percent_non_normal$coloc_cis)
sum(overlap_multiple_soma_pqtl_percent_non_normal$coloc_trans>0)
sum(overlap_multiple_soma_pqtl_percent_non_normal$coloc_trans)

# how many have matched SOMAmers
table(overlap_multiple_soma_pqtl_percent_non_normal$total)
