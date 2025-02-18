### this script estimates the observational correlations between Olink and SomaScan protein measurements

rm(list = ls())

setwd("")

library(tidyverse)
library(ggplot2)
library(ggpattern)
library(ggpubr)
# library(ppcor)
library(gtsummary)
library(Hmisc)
library(ggthemes)
library(datawizard)
library(ckbplotr)
library(scales)
library(reshape2)

######################## load raw, transformed and cleaned somascan data

overlap <- read.csv("overlap.csv")
overlap_1_to_1 <- read.csv("overlap_1_to_1.csv")
olink <- read.csv("olink.csv")
somascan_normalised <- read.csv("somascan_normalised.csv")
somascan_non_normalised <- read.csv("somascan_non_normalised.csv")
somascan_normalised_log <- read.csv("somascan_normalised_log.csv")
somascan_non_normalised_log <- read.csv("somascan_non_normalised_log.csv")


########## check correlation

overlap_cor <- overlap

## spearman

for (i in 1:nrow(overlap_cor)) {
  
  df <- merge(olink[,c("csid",overlap_cor$olink_id[i])],somascan_normalised[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$rho_olink_soma_normal[i] <- cor.test(df[[2]], df[[3]], method = 'spearman')$estimate
  
  df <- merge(olink[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$rho_olink_soma_non_normal[i] <- cor.test(df[[2]], df[[3]], method = 'spearman')$estimate
  
}


## pearson

for (i in 1:nrow(overlap_cor)) {
  
  df <- merge(olink[,c("csid",overlap_cor$olink_id[i])],somascan_normalised[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_normal[i] <- cor.test(df[[2]], df[[3]], method = 'pearson')$estimate
  
  df <- merge(olink[,c("csid",overlap_cor$olink_id[i])],somascan_normalised_log[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_normal_log[i] <- cor.test(df[[2]], df[[3]], method = 'pearson')$estimate
  
  df <- merge(olink[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_non_normal[i] <- cor.test(df[[2]], df[[3]], method = 'pearson')$estimate
  
  df <- merge(olink[,c("csid",overlap_cor$olink_id[i])],somascan_non_normalised_log[,c("csid",overlap_cor$somascan_id[i])],by="csid")
  overlap_cor$r_olink_soma_non_normal_log[i] <- cor.test(df[[2]], df[[3]], method = 'pearson')$estimate
  
}


# save

write.csv(overlap_cor,"overlap_cor.csv", quote=F, row.names=F)

# check correlation with CI

overlap_cor_ci <- overlap

## spearman

for (i in 1:nrow(overlap_cor_ci)) {
  
  df <- merge(olink[,c("csid",overlap_cor_ci$olink_id[i])],somascan_normalised[,c("csid",overlap_cor_ci$somascan_id[i])],by="csid")
  overlap_cor_ci$rho_olink_soma_normal[i] <- cor.test(df[[2]], df[[3]], method = 'spearman')$estimate
  
  df <- merge(olink[,c("csid",overlap_cor_ci$olink_id[i])],somascan_non_normalised[,c("csid",overlap_cor_ci$somascan_id[i])],by="csid")
  overlap_cor_ci$rho_olink_soma_non_normal[i] <- cor.test(df[[2]], df[[3]], method = 'spearman')$estimate
  
}


## pearson

for (i in 1:nrow(overlap_cor_ci)) {
  
  df <- merge(olink[,c("csid",overlap_cor_ci$olink_id[i])],somascan_normalised[,c("csid",overlap_cor_ci$somascan_id[i])],by="csid")
  overlap_cor_ci$r_olink_soma_normal[i] <- cor.test(df[[2]], df[[3]], method = 'pearson')$estimate
  
  df <- merge(olink[,c("csid",overlap_cor_ci$olink_id[i])],somascan_normalised_log[,c("csid",overlap_cor_ci$somascan_id[i])],by="csid")
  overlap_cor_ci$r_olink_soma_normal_log[i] <- cor.test(df[[2]], df[[3]], method = 'pearson')$estimate
  
  df <- merge(olink[,c("csid",overlap_cor_ci$olink_id[i])],somascan_non_normalised[,c("csid",overlap_cor_ci$somascan_id[i])],by="csid")
  overlap_cor_ci$r_olink_soma_non_normal[i] <- cor.test(df[[2]], df[[3]], method = 'pearson')$estimate
  
  df <- merge(olink[,c("csid",overlap_cor_ci$olink_id[i])],somascan_non_normalised_log[,c("csid",overlap_cor_ci$somascan_id[i])],by="csid")
  overlap_cor_ci$r_olink_soma_non_normal_log[i] <- cor.test(df[[2]], df[[3]], method = 'pearson')$estimate
  
}

# overlap_cor <- read.csv("overlap_cor.csv")

# calculate median rho

median(overlap_cor$rho_olink_soma_normal)
median(overlap_cor$rho_olink_soma_non_normal)

# also try the weighted mean method

weighted <- aggregate(overlap_cor$rho_olink_soma_normal,by=list(overlap_cor$uniprot_id),data=overlap_cor,FUN=mean)
names(weighted) <- c("uniprot_id","rho_olink_soma_normal")
median(weighted$rho_olink_soma_normal)

# spearman olink vs soma

plot_rho_olink_soma_normal <- ggplot(overlap_cor, aes(x=rho_olink_soma_normal)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,450) +
  annotate("text", x = 0.5, y = 200, label = "Median rho = 0.23", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan-ANML", subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18), text = element_text(size = 14))  

plot_rho_olink_soma_normal

plot_rho_olink_soma_non_normal <- ggplot(overlap_cor, aes(x=rho_olink_soma_non_normal)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,450) +
  annotate("text", x = 0.55, y = 200, label = "Median rho = 0.29", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18),text = element_text(size = 14))  

plot_rho_olink_soma_non_normal

plot_rho_olink_soma <- ggarrange(plot_rho_olink_soma_normal,plot_rho_olink_soma_non_normal,ncol=2,nrow=1)

ggsave("",plot_rho_olink_soma,width=14,height=6)

# pearson olink vs soma

median(overlap_cor$r_olink_soma_normal_log)
median(overlap_cor$r_olink_soma_non_normal_log)

plot_r_olink_soma_normal <- ggplot(overlap_cor, aes(x=r_olink_soma_normal_log)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,500) +
  annotate("text", x = 0.4, y = 200, label = "Median r = 0.17", size = 12/.pt) +
  geom_vline(aes(xintercept=median(r_olink_soma_normal_log)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan-ANML", subtitle ="") + xlab("Pearson's r") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18), text = element_text(size = 14))  

plot_r_olink_soma_normal

plot_r_olink_soma_non_normal <- ggplot(overlap_cor, aes(x=r_olink_soma_non_normal_log)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.6, 1, 0.2), lim = c(-0.6, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,500) +
  annotate("text", x = 0.5, y = 200, label = "Median r = 0.23", size = 12/.pt) +
  geom_vline(aes(xintercept=median(r_olink_soma_non_normal_log)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") + xlab("Pearson's r") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18),text = element_text(size = 14))  

plot_r_olink_soma_non_normal

plot_r_olink_soma <- ggarrange(plot_r_olink_soma_normal,plot_r_olink_soma_non_normal,ncol=2,nrow=1)

ggsave("",plot_r_olink_soma,width=14,height=6)


# only keep 1_to_1

protein_dup <- unique(overlap_cor$uniprot_id[duplicated(overlap_cor$uniprot_id)]) # get proteins targeted by multiple aptamers

aptamer_dup <- unique(overlap_cor$somascan_id[duplicated(overlap_cor$somascan_id)]) # get aptamers targeted by multiple proteins

overlap_1_to_1_cor <- overlap_cor[-c(which(overlap_cor$uniprot_id %in% protein_dup),which(overlap_cor$somascan_id %in% aptamer_dup)), ]

write.csv(overlap_1_to_1_cor,"overlap_1_to_1_cor.csv", quote=F, row.names=F)

##################################### get median

names(overlap_1_to_1_cor)

coeff_sum <- data.frame(coeff = names(overlap_1_to_1_cor)[9:14],
                        median = numeric(6), max = numeric(6), min = numeric(6))

coeff_sum$median[1:6] <- lapply(overlap_1_to_1_cor[,c(9:14)],median)
coeff_sum$max[1:6] <- lapply(overlap_1_to_1_cor[,c(9:14)],max)
coeff_sum$min[1:6] <- lapply(overlap_1_to_1_cor[,c(9:14)],min)

print(coeff_sum, row.names=F)

##################################### plot

# spearman olink vs soma

plot_rho_olink_soma_normal <- ggplot(overlap_1_to_1_cor, aes(x=rho_olink_soma_normal)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  annotate("text", x = 0.4, y = 200, label = "Median rho = 0.20", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan-ANML", subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18), text = element_text(size = 14))  

plot_rho_olink_soma_normal

plot_rho_olink_soma_non_normal <- ggplot(overlap_1_to_1_cor, aes(x=rho_olink_soma_non_normal)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  annotate("text", x = 0.5, y = 200, label = "Median rho = 0.26", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18),text = element_text(size = 14))  

plot_rho_olink_soma_non_normal

plot_rho_olink_soma <- ggarrange(plot_rho_olink_soma_normal,plot_rho_olink_soma_non_normal,ncol=2,nrow=1)

ggsave("",plot_rho_olink_soma,width=14,height=6)

## pearson olink vs soma

plot_r_olink_soma_normal <- ggplot(overlap_1_to_1_cor, aes(x=r_olink_soma_normal_log)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.2, 1, 0.2), lim = c(-0.2, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,400) +
  annotate("text", x = 0.4, y = 200, label = "Median r = 0.15", size = 12/.pt) +
  geom_vline(aes(xintercept=median(r_olink_soma_normal_log)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan-ANML", subtitle ="") + xlab("Pearson's r") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18), text = element_text(size = 14))  

plot_r_olink_soma_normal

plot_r_olink_soma_non_normal <- ggplot(overlap_1_to_1_cor, aes(x=r_olink_soma_non_normal_log)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.2, 1, 0.2), lim = c(-0.2, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,400) +
  annotate("text", x = 0.4, y = 200, label = "Median r = 0.20", size = 12/.pt) +
  geom_vline(aes(xintercept=median(r_olink_soma_non_normal_log)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") + xlab("Pearson's r") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18),text = element_text(size = 14))  

plot_r_olink_soma_non_normal

plot_r_olink_soma <- ggarrange(plot_r_olink_soma_normal,plot_r_olink_soma_non_normal,ncol=2,nrow=1)

ggsave("",plot_r_olink_soma,width=14,height=6)


# overlay

names(overlap_1_to_1_cor)

# normal

r_olink_soma_normal <- pivot_longer(overlap_1_to_1_cor[,c(1,11,12)], cols=2:3)

plot_r_olink_soma_normal <- ggplot(r_olink_soma_normal, aes(x=value,fill=name)) +
  geom_histogram(binwidth = 0.05,alpha=0.5,boundary=0,position="identity") +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1)) +
  ylim(0,500) +
  ggtitle("Olink vs SomaScan-ANML", subtitle ="") + 
  xlab("Pearson's r") + 
  ylab("Frequency") +
  scale_fill_discrete(name = "Log-transformation", labels = c("Before", "After")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  theme(legend.position = c(.95, .95), legend.justification = c("right", "top"))

plot_r_olink_soma_normal

# non_normal

r_olink_soma_non_normal <- pivot_longer(overlap_1_to_1_cor[,c(1,13,14)], cols=2:3)

plot_r_olink_soma_non_normal <- ggplot(r_olink_soma_non_normal, aes(x=value,fill=name)) +
  geom_histogram(binwidth = 0.05,alpha=0.5,boundary=0,position="identity") +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1)) +
  ylim(0,500) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") + 
  xlab("Pearson's r") + 
  ylab("Frequency") +
  scale_fill_discrete(name = "Log-transformation", labels = c("Before", "After")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  theme(legend.position = c(.95, .95), legend.justification = c("right", "top"))

plot_r_olink_soma_non_normal

plot_r_olink_soma_overlay <- ggarrange(plot_r_olink_soma_normal,plot_r_olink_soma_non_normal,ncol=1,nrow=2)

ggsave("",plot_r_olink_soma_overlay,width=8,height=10)

## separate

olink_extra <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/proteomics/olink_extra_info.csv")
olink_extra <- olink_extra[!duplicated(olink_extra$UniProt),]

# dilution

olink_dilution <- olink_extra[,c("UniProt","Dilution.factor")]
names(olink_dilution) <- c("uniprot_id","dilution_olink")
olink_dilution$dilution_olink <- gsub(" ", "", olink_dilution$dilution_olink)
table(olink_dilution$dilution_olink)

# merge

overlap_1_to_1_cor <- merge(overlap_1_to_1_cor,olink_dilution,by="uniprot_id")

batch <- sort(unique(overlap_1_to_1_cor$batch))
dilution <- sort(unique(overlap_1_to_1_cor$dilution),decreasing=T)
dilution_label <- c("20%","0.5%","0.005%")
dilution_olink <- sort(unique(overlap_1_to_1_cor$dilution_olink))
dilution_olink_label <- c("1:1","1:10","1:100","1:1000","1:100000")

## table

table_1 <- overlap_1_to_1_cor %>% 
  select(rho_olink_soma_normal,rho_olink_soma_non_normal,batch) %>%
  tbl_summary(by = batch, 
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ 
                c("{median} ({p25}, {p75})","{min}, {max}")) 

table_1

table_2 <- overlap_1_to_1_cor %>% 
  select(rho_olink_soma_normal,rho_olink_soma_non_normal,dilution) %>%
  tbl_summary(by = dilution, 
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ 
                c("{median} ({p25}, {p75})","{min}, {max}")) 

table_2

table_3 <- overlap_1_to_1_cor %>% 
  select(rho_olink_soma_normal,rho_olink_soma_non_normal,dilution_olink) %>%
  tbl_summary(by = dilution_olink, 
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ 
                c("{median} ({p25}, {p75})","{min}, {max}")) 

table_3


# by batch

# normalised

olink_somascan_batch_normal_list <- list()

i <- 1

for (x in 1:2){
  olink_somascan_batch_normal_list[[i]] <- ggplot(overlap_1_to_1_cor[overlap_1_to_1_cor$batch==x,], aes(x=rho_olink_soma_normal)) + 
    geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
    scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
    ylim(0,300) +
    ggtitle(paste0("Olink (batch ",batch[x],") vs SomaScan-ANML"), subtitle ="")  + 
    xlab("Spearman's rho") + 
    ylab("Frequency") +  
    geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
    # annotate("text", x = 0.7, y = 300, label = paste("Median rho =", medians[x]), size = 12/.pt) +
    theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
    theme(text = element_text(size = 13))
  i <- i+1
}

olink_somascan_batch_normal <- ggarrange(olink_somascan_batch_normal_list[[1]],olink_somascan_batch_normal_list[[2]],ncol=2,nrow=1)

ggsave("",olink_somascan_batch_normal,width=10,height=4.5)


# non_normalised

# medians <- c(0.42, 0.13)

olink_somascan_batch_non_normal_list <- list()

i <- 1

for (x in 1:2){
  olink_somascan_batch_non_normal_list[[i]] <- ggplot(overlap_1_to_1_cor[overlap_1_to_1_cor$batch==x,], aes(x=rho_olink_soma_non_normal)) + 
    geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
    scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
    ylim(0,300) +
    ggtitle(paste0("Olink (batch ",batch[x],") vs SomaScan-non-ANML"), subtitle ="")  + 
    xlab("Spearman's rho") + 
    ylab("Frequency") +
    # annotate("text", x = 0.7, y = 300, label = paste("Median rho =", medians[x]), size = 12/.pt) +
    geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
    theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
    theme(text = element_text(size = 13))
  i <- i+1
}

olink_somascan_batch_non_normal <- ggarrange(olink_somascan_batch_non_normal_list[[1]],olink_somascan_batch_non_normal_list[[2]],ncol=2,nrow=1)

ggsave("",olink_somascan_batch_non_normal,width=10,height=4.5)



# by SomaScan dilution

# normalised

olink_somascan_dilution_normal_list <- list()

i <- 1

for (y in 1:3){
  olink_somascan_dilution_normal_list[[i]] <- ggplot(overlap_1_to_1_cor[overlap_1_to_1_cor$dilution==dilution[y],], aes(x=rho_olink_soma_normal)) + 
    geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
    scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
    ylim(0,250) +
    ggtitle(paste0("Olink vs SomaScan (",dilution_label[y], ", ANML)"), subtitle ="")  + 
    xlab("Spearman's rho") + 
    ylab("Frequency") +
    geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
    # annotate("text", x = 0.35, y = 300, label = paste("Median rho =", medians[y]), size = 12/.pt) +
    theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
    theme(text = element_text(size = 13))
  i <- i+1
}

olink_somascan_dilution_normal <- ggarrange(olink_somascan_dilution_normal_list[[1]],olink_somascan_dilution_normal_list[[2]],olink_somascan_dilution_normal_list[[3]],ncol=3,nrow=1)

ggsave("",olink_somascan_dilution_normal,width=15,height=4.5)


# non_normalised

olink_somascan_dilution_non_normal_list <- list()

i <- 1

for (y in 1:3){
  olink_somascan_dilution_non_normal_list[[i]] <- ggplot(overlap_1_to_1_cor[overlap_1_to_1_cor$dilution==dilution[y],], aes(x=rho_olink_soma_non_normal)) + 
    geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
    scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
    ylim(0,250) +
    ggtitle(paste0("Olink vs SomaScan (",dilution_label[y], ", non-ANML)"), subtitle ="")  + 
    xlab("Spearman's rho") + 
    ylab("Frequency") +
    geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
    # annotate("text", x = 0.42, y = 300, label = paste("Median rho =", medians[y]), size = 12/.pt) +
    theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
    theme(text = element_text(size = 13))
  i <- i+1
}

olink_somascan_dilution_non_normal <- ggarrange(olink_somascan_dilution_non_normal_list[[1]],olink_somascan_dilution_non_normal_list[[2]],olink_somascan_dilution_non_normal_list[[3]],ncol=3,nrow=1)

ggsave("",olink_somascan_dilution_non_normal,width=15,height=4.5)

olink_somascan_dilution <- ggarrange(olink_somascan_dilution_normal_list[[1]],olink_somascan_dilution_normal_list[[2]],olink_somascan_dilution_normal_list[[3]],
                                     olink_somascan_dilution_non_normal_list[[1]],olink_somascan_dilution_non_normal_list[[2]],olink_somascan_dilution_non_normal_list[[3]],
                                     ncol=3,nrow=2)

ggsave("",olink_somascan_dilution,width=15,height=9)



# by olink dilution

# normalised

olink_somascan_dilution_olink_normal_list <- list()

i <- 1

for (y in 1:5){
  olink_somascan_dilution_olink_normal_list[[i]] <- ggplot(overlap_1_to_1_cor[overlap_1_to_1_cor$dilution_olink==dilution_olink[y],], aes(x=rho_olink_soma_normal)) + 
    geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
    scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
    ylim(0,250) +
    ggtitle(paste0("Olink vs SomaScan (",dilution_olink_label[y], ", ANML)"), subtitle ="")  + 
    xlab("Spearman's rho") + 
    ylab("Frequency") +
    geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
    # annotate("text", x = 0.35, y = 300, label = paste("Median rho =", medians[y]), size = 12/.pt) +
    theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
    theme(text = element_text(size = 13))
  i <- i+1
}

olink_somascan_dilution_olink_normal <- ggarrange(olink_somascan_dilution_olink_normal_list[[1]],
                                                  olink_somascan_dilution_olink_normal_list[[2]],
                                                  olink_somascan_dilution_olink_normal_list[[3]],
                                                  olink_somascan_dilution_olink_normal_list[[4]],
                                                  olink_somascan_dilution_olink_normal_list[[5]],
                                                  ncol=5,nrow=1)

ggsave("",olink_somascan_dilution_olink_normal,width=25,height=4.5)


# non_normalised

# medians <- c("0.65", "0.68", "0.19")

olink_somascan_dilution_olink_non_normal_list <- list()

i <- 1

for (y in 1:5){
  olink_somascan_dilution_olink_non_normal_list[[i]] <- ggplot(overlap_1_to_1_cor[overlap_1_to_1_cor$dilution_olink==dilution_olink[y],], aes(x=rho_olink_soma_non_normal)) + 
    geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
    scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
    ylim(0,250) +
    ggtitle(paste0("Olink vs SomaScan (",dilution_olink_label[y], ", non-ANML)"), subtitle ="")  + 
    xlab("Spearman's rho") + 
    ylab("Frequency") +
    geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
    # annotate("text", x = 0.42, y = 300, label = paste("Median rho =", medians[y]), size = 12/.pt) +
    theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
    theme(text = element_text(size = 13))
  i <- i+1
}

olink_somascan_dilution_olink_non_normal <- ggarrange(olink_somascan_dilution_olink_non_normal_list[[1]],
                                                      olink_somascan_dilution_olink_non_normal_list[[2]],
                                                      olink_somascan_dilution_olink_non_normal_list[[3]],
                                                      olink_somascan_dilution_olink_non_normal_list[[4]],
                                                      olink_somascan_dilution_olink_non_normal_list[[5]],
                                                      ncol=5,nrow=1)

ggsave("",olink_somascan_dilution_olink_non_normal,width=25,height=4.5)

olink_somascan_dilution_olink <- ggarrange(olink_somascan_dilution_olink_normal_list[[1]],
                                           olink_somascan_dilution_olink_normal_list[[2]],
                                           olink_somascan_dilution_olink_normal_list[[3]],
                                           olink_somascan_dilution_olink_normal_list[[4]],
                                           olink_somascan_dilution_olink_normal_list[[5]],
                                           olink_somascan_dilution_olink_non_normal_list[[1]],
                                           olink_somascan_dilution_olink_non_normal_list[[2]],
                                           olink_somascan_dilution_olink_non_normal_list[[3]],
                                           olink_somascan_dilution_olink_non_normal_list[[4]],
                                           olink_somascan_dilution_olink_non_normal_list[[5]],
                                           ncol=5,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/olink_somascan_dilution_olink.png",olink_somascan_dilution_olink,width=25,height=9)




################### correlation between all possible pairs

olink_overlap <- olink[,names(olink) %in% overlap$olink_id]
olink_overlap <- olink_overlap[, unique(overlap$olink_id)]
olink_overlap <- cbind(olink$csid,olink_overlap)
names(olink_overlap)[1] <- "csid"

somascan_normalised_log_overlap <- somascan_normalised_log[,names(somascan_normalised_log) %in% overlap$somascan_id]
somascan_normalised_log_overlap <- somascan_normalised_log_overlap[, unique(overlap$somascan_id)]
somascan_normalised_log_overlap <- cbind(somascan_normalised_log$csid,somascan_normalised_log_overlap)
names(somascan_normalised_log_overlap)[1] <- "csid"

olink_overlap <- olink_overlap[olink_overlap$csid %in% somascan_normalised_log_overlap$csid,]

somascan_non_normalised_log_overlap <- somascan_non_normalised_log[,names(somascan_non_normalised_log) %in% overlap$somascan_id]
somascan_non_normalised_log_overlap <- somascan_non_normalised_log_overlap[, unique(overlap$somascan_id)]
somascan_non_normalised_log_overlap <- cbind(somascan_non_normalised_log$csid,somascan_non_normalised_log_overlap)
names(somascan_non_normalised_log_overlap)[1] <- "csid"

## normal

# correlation matrix
cor_mat_normal <- cor(olink_overlap[,-1], somascan_normalised_log_overlap[,-1], method="spearman")
# heatmap(cor_mat,Rowv=NA,Colv=NA,labRow=NA,labCol=NA, xlab="SomaScan", ylab="OLINK")

cordata_normal = melt(cor_mat_normal)

heatmap_normal <- ggplot(cordata_normal, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + xlab("Olink") + ylab("SomaScan-ANML") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(fill="Spearman's rho") +
  scale_fill_gradient2(low="blue", high="red" )

ggsave("",heatmap_normal,width=8,height=6)

# check the most correlated pair against the matching based on UniProt
cordata_normal_olink <- cordata_normal[order(cordata_normal$Var1,-cordata_normal$value),]
cordata_normal_olink <-cordata_normal_olink[!duplicated(cordata_normal_olink$Var1),]

cordata_normal_soma <- cordata_normal[order(cordata_normal$Var2,-cordata_normal$value),]
cordata_normal_soma <- cordata_normal_soma[!duplicated(cordata_normal_soma$Var2),]

cordata_normal_match <- merge(cordata_normal_olink,cordata_normal_soma,by=c("Var1","Var2","value"))
names(cordata_normal_match) <- c("olink_id","somascan_id","rho")
nrow(cordata_normal_match)
write.csv(cordata_normal_match,"mutual_best_hits_normal.csv", quote=F, row.names=F)

match_normal <- merge(cordata_normal_match, overlap, by=c("olink_id","somascan_id"))
nrow(match_normal)

## non-normal

# correlation matrix
cor_mat_non_normal <- cor(olink_overlap[,-1], somascan_non_normalised_log_overlap[,-1], method="spearman")
# heatmap(cor_mat,Rowv=NA,Colv=NA,labRow=NA,labCol=NA, xlab="SomaScan", ylab="OLINK")

cordata_non_normal = melt(cor_mat_non_normal)

heatmap_non_normal <- ggplot(cordata_non_normal, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + xlab("Olink") + ylab("SomaScan-non-ANML") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(fill="Spearman's rho") +
  scale_fill_gradient2(low="blue", high="red" )

ggsave("",heatmap_non_normal,width=8,height=6)

heatmap <- ggarrange(heatmap_normal,
                     heatmap_non_normal,
                     ncol=2,nrow=1, common.legend = TRUE, legend="right")

ggsave("",heatmap,width=14,height=6)


# check the most correlated pair against the matching based on UniProt
cordata_non_normal_olink <- cordata_non_normal[order(cordata_non_normal$Var1,-cordata_non_normal$value),]
cordata_non_normal_olink <-cordata_non_normal_olink[!duplicated(cordata_non_normal_olink$Var1),]

cordata_non_normal_soma <- cordata_non_normal[order(cordata_non_normal$Var2,-cordata_non_normal$value),]
cordata_non_normal_soma <- cordata_non_normal_soma[!duplicated(cordata_non_normal_soma$Var2),]

cordata_non_normal_match <- merge(cordata_non_normal_olink,cordata_non_normal_soma,by=c("Var1","Var2","value"))
names(cordata_non_normal_match) <- c("olink_id","somascan_id","rho")
nrow(cordata_non_normal_match)
write.csv(cordata_non_normal_match,"mutual_best_hits_non_normal.csv", quote=F, row.names=F)

match_non_normal <- merge(cordata_non_normal_match, overlap, by=c("olink_id","somascan_id"))
nrow(match_non_normal)



