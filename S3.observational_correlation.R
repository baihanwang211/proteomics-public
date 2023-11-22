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


# save

# overlap_cor <- read.csv("overlap_cor.csv")

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
  ggtitle("OLINK vs SomaScan-ANML", subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18), text = element_text(size = 14))  

plot_rho_olink_soma_normal

plot_rho_olink_soma_non_normal <- ggplot(overlap_1_to_1_cor, aes(x=rho_olink_soma_non_normal)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  annotate("text", x = 0.5, y = 200, label = "Median rho = 0.26", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("OLINK vs SomaScan-non-ANML", subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.title = element_text(size = 18),text = element_text(size = 14))  

plot_rho_olink_soma_non_normal

plot_rho_olink_soma <- ggarrange(plot_rho_olink_soma_normal,plot_rho_olink_soma_non_normal,ncol=2,nrow=1)

ggsave("plot_rho_olink_soma_1_to_1.png",plot_rho_olink_soma,width=14,height=6)

## pearson olink vs soma

# overlay

names(overlap_1_to_1_cor)

# normal

r_olink_soma_normal <- pivot_longer(overlap_1_to_1_cor[,c(1,11,12)], cols=2:3)

plot_r_olink_soma_normal <- ggplot(r_olink_soma_normal, aes(x=value,fill=name)) +
  geom_histogram(binwidth = 0.05,alpha=0.5,boundary=0,position="identity") +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1)) +
  ylim(0,500) +
  ggtitle("OLINK vs SomaScan-ANML", subtitle ="") + 
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
  ggtitle("OLINK vs SomaScan-non-ANML", subtitle ="") + 
  xlab("Pearson's r") + 
  ylab("Frequency") +
  scale_fill_discrete(name = "Log-transformation", labels = c("Before", "After")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  theme(legend.position = c(.95, .95), legend.justification = c("right", "top"))

plot_r_olink_soma_non_normal

plot_r_olink_soma_overlay <- ggarrange(plot_r_olink_soma_normal,plot_r_olink_soma_non_normal,ncol=1,nrow=2)

ggsave("plot_r_olink_soma_overlay_1_to_1.png",plot_r_olink_soma_overlay,width=8,height=10)



## separate

olink_extra <- read.csv("olink_extra_info.csv")
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


# batch

# normalised

olink_somascan_batch_normal_list <- list()

i <- 1

for (x in 1:2){
  olink_somascan_batch_normal_list[[i]] <- ggplot(overlap_1_to_1_cor[overlap_1_to_1_cor$batch==x,], aes(x=rho_olink_soma_normal)) + 
    geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
    scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
    ylim(0,300) +
    ggtitle(paste0("OLINK (batch ",batch[x],") vs SomaScan-ANML"), subtitle ="")  + 
    xlab("Spearman's rho") + 
    ylab("Frequency") +  
    geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
    theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
    theme(text = element_text(size = 13))
  i <- i+1
}

olink_somascan_batch_normal <- ggarrange(olink_somascan_batch_normal_list[[1]],olink_somascan_batch_normal_list[[2]],ncol=2,nrow=1)

ggsave("olink_somascan_batch_normal.png",olink_somascan_batch_normal,width=10,height=4.5)


# non_normalised

# medians <- c(0.42, 0.13)

olink_somascan_batch_non_normal_list <- list()

i <- 1

for (x in 1:2){
  olink_somascan_batch_non_normal_list[[i]] <- ggplot(overlap_1_to_1_cor[overlap_1_to_1_cor$batch==x,], aes(x=rho_olink_soma_non_normal)) + 
    geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
    scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
    ylim(0,300) +
    ggtitle(paste0("OLINK (batch ",batch[x],") vs SomaScan-non-ANML"), subtitle ="")  + 
    xlab("Spearman's rho") + 
    ylab("Frequency") +
    geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
    theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
    theme(text = element_text(size = 13))
  i <- i+1
}

olink_somascan_batch_non_normal <- ggarrange(olink_somascan_batch_non_normal_list[[1]],olink_somascan_batch_non_normal_list[[2]],ncol=2,nrow=1)

ggsave("olink_somascan_batch_non_normal.png",olink_somascan_batch_non_normal,width=10,height=4.5)



# dilution

# normalised

olink_somascan_dilution_normal_list <- list()

i <- 1

for (y in 1:3){
  olink_somascan_dilution_normal_list[[i]] <- ggplot(overlap_1_to_1_cor[overlap_1_to_1_cor$dilution==dilution[y],], aes(x=rho_olink_soma_normal)) + 
    geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
    scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
    ylim(0,250) +
    ggtitle(paste0("OLINK vs SomaScan (",dilution_label[y], ", ANML)"), subtitle ="")  + 
    xlab("Spearman's rho") + 
    ylab("Frequency") +
    geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
    theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
    theme(text = element_text(size = 13))
  i <- i+1
}

olink_somascan_dilution_normal <- ggarrange(olink_somascan_dilution_normal_list[[1]],olink_somascan_dilution_normal_list[[2]],olink_somascan_dilution_normal_list[[3]],ncol=3,nrow=1)

ggsave("olink_somascan_dilution_normal.png",olink_somascan_dilution_normal,width=15,height=4.5)


# non_normalised

olink_somascan_dilution_non_normal_list <- list()

i <- 1

for (y in 1:3){
  olink_somascan_dilution_non_normal_list[[i]] <- ggplot(overlap_1_to_1_cor[overlap_1_to_1_cor$dilution==dilution[y],], aes(x=rho_olink_soma_non_normal)) + 
    geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
    scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
    ylim(0,250) +
    ggtitle(paste0("OLINK vs SomaScan (",dilution_label[y], ", non-ANML)"), subtitle ="")  + 
    xlab("Spearman's rho") + 
    ylab("Frequency") +
    geom_vline(aes(xintercept=median(rho_olink_soma_non_normal)),color="black", linetype="dashed", linewidth=1) +
    theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
    theme(text = element_text(size = 13))
  i <- i+1
}

olink_somascan_dilution_non_normal <- ggarrange(olink_somascan_dilution_non_normal_list[[1]],olink_somascan_dilution_non_normal_list[[2]],olink_somascan_dilution_non_normal_list[[3]],ncol=3,nrow=1)

ggsave("olink_somascan_dilution_non_normal.png",olink_somascan_dilution_non_normal,width=15,height=4.5)

olink_somascan_dilution <- ggarrange(olink_somascan_dilution_normal_list[[1]],olink_somascan_dilution_normal_list[[2]],olink_somascan_dilution_normal_list[[3]],
                                     olink_somascan_dilution_non_normal_list[[1]],olink_somascan_dilution_non_normal_list[[2]],olink_somascan_dilution_non_normal_list[[3]],
                                     ncol=3,nrow=2)

ggsave("olink_somascan_dilution.png",olink_somascan_dilution,width=15,height=9)




# olink dilution

# normalised

olink_somascan_dilution_olink_normal_list <- list()

i <- 1

for (y in 1:5){
  olink_somascan_dilution_olink_normal_list[[i]] <- ggplot(overlap_1_to_1_cor[overlap_1_to_1_cor$dilution_olink==dilution_olink[y],], aes(x=rho_olink_soma_normal)) + 
    geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
    scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
    ylim(0,250) +
    ggtitle(paste0("OLINK vs SomaScan (",dilution_olink_label[y], ", ANML)"), subtitle ="")  + 
    xlab("Spearman's rho") + 
    ylab("Frequency") +
    geom_vline(aes(xintercept=median(rho_olink_soma_normal)),color="black", linetype="dashed", linewidth=1) +
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

ggsave("olink_somascan_dilution_olink_normal.png",olink_somascan_dilution_olink_normal,width=25,height=4.5)


# non_normalised

# medians <- c("0.65", "0.68", "0.19")

olink_somascan_dilution_olink_non_normal_list <- list()

i <- 1

for (y in 1:5){
  olink_somascan_dilution_olink_non_normal_list[[i]] <- ggplot(overlap_1_to_1_cor[overlap_1_to_1_cor$dilution_olink==dilution_olink[y],], aes(x=rho_olink_soma_non_normal)) + 
    geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
    scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
    ylim(0,250) +
    ggtitle(paste0("OLINK vs SomaScan (",dilution_olink_label[y], ", non-ANML)"), subtitle ="")  + 
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

ggsave("olink_somascan_dilution_olink_non_normal.png",olink_somascan_dilution_olink_non_normal,width=25,height=4.5)

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

ggsave("olink_somascan_dilution_olink.png",olink_somascan_dilution_olink,width=25,height=9)








################ extra: correlation between soma

for (i in 1:nrow(overlap_1_to_1_cor)) {
  df <- merge(somascan_normalised_log[,c("csid",overlap_1_to_1_cor$somascan_id[i])],somascan_non_normalised_log[,c("csid",overlap_1_to_1_cor$somascan_id[i])],by="csid")
  overlap_1_to_1_cor$rho_soma[i] <- cor.test(df[[2]], df[[3]], method = 'spearman')$estimate
}

min(overlap_1_to_1_cor$rho_soma)
max(overlap_1_to_1_cor$rho_soma)
median(overlap_1_to_1_cor$rho_soma)

plot_rho_soma <- ggplot(overlap_1_to_1_cor, aes(x=rho_soma)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black",fill="white") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), lim = c(0, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,350) +
  annotate("text", x = 0.7, y = 300, label = "Median rho = 0.84", size = 12/.pt) +
  geom_vline(aes(xintercept=median(rho_soma)),color="black", linetype="dashed", linewidth=1) +
  ggtitle("SomaScan ANML vs non-ANML", subtitle ="") + xlab("Spearman's rho") + ylab("Frequency") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 14))  

plot_rho_soma

ggsave("plot_rho_soma_1_to_1.png",plot_rho_soma,width=7,height=6)

write.csv(overlap_1_to_1_cor,"overlap_1_to_1_cor_soma.csv", quote=F, row.names=F)
