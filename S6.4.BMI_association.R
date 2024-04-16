### this script test associations between protein levels and bmi

rm(list = ls())

setwd("")

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(egg)
library(forcats)
library(ggpattern)
library(ckbplotr)
library(scales)
library(ggrepel)
library(VennDiagram)
library(RColorBrewer)
library(ckbplotr)

### run regression again to get SD

# load data

olink_pheno <- read.csv("olink_pheno.csv")
somascan_normalised_log_pheno <- read.csv("somascan_normalised_log_pheno.csv")
somascan_non_normalised_log_pheno <- read.csv("somascan_non_normalised_log_pheno.csv")

# change some variables to factor
olink_pheno$is_female <- factor(olink_pheno$is_female)
somascan_normalised_log_pheno$is_female <- factor(somascan_normalised_log_pheno$is_female)
somascan_non_normalised_log_pheno$is_female <- factor(somascan_non_normalised_log_pheno$is_female)

olink_pheno$region_code <- factor(olink_pheno$region_code)
somascan_normalised_log_pheno$region_code <- factor(somascan_normalised_log_pheno$region_code)
somascan_non_normalised_log_pheno$region_code <- factor(somascan_non_normalised_log_pheno$region_code)

olink_pheno$ascertainment <- factor(olink_pheno$ascertainment)
somascan_normalised_log_pheno$ascertainment <- factor(somascan_normalised_log_pheno$ascertainment)
somascan_non_normalised_log_pheno$ascertainment <- factor(somascan_non_normalised_log_pheno$ascertainment)

# keep same participants
olink_pheno <- olink_pheno[olink_pheno$csid %in% somascan_normalised_log_pheno$csid, ]

# run regression

bmi <- read.csv("overlap_1_to_1_cor.csv")

for (i in 1:nrow(bmi)){

  olink_id <- bmi$olink_id[i]
  somascan_id <- bmi$somascan_id[i]

  model <- lm(paste(olink_id, "~ bmi_calc + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                       olink_plt_id + ascertainment + age + age2 + is_female + region_code"), olink_pheno)
  bmi$olink_es_bmi_calc[i] <- summary(model)$coefficients[2,"Estimate"]
  bmi$olink_se_bmi_calc[i] <- summary(model)$coefficients[2,"Std. Error"]
  bmi$olink_lci_bmi_calc[i] <- confint(model,"bmi_calc")[1]
  bmi$olink_hci_bmi_calc[i] <- confint(model,"bmi_calc")[2]
  bmi$olink_p_bmi_calc[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]

  model <- lm(paste(somascan_id, "~ bmi_calc + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                       somascan_plt_id + ascertainment + age + age2 + is_female + region_code"), somascan_normalised_log_pheno)
  bmi$soma_normal_es_bmi_calc[i] <- summary(model)$coefficients[2,"Estimate"]
  bmi$soma_normal_se_bmi_calc[i] <- summary(model)$coefficients[2,"Std. Error"]
  bmi$soma_normal_lci_bmi_calc[i] <- confint(model,"bmi_calc")[1]
  bmi$soma_normal_hci_bmi_calc[i] <- confint(model,"bmi_calc")[2]
  bmi$soma_normal_p_bmi_calc[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]

  model <- lm(paste(somascan_id, "~ bmi_calc + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                       somascan_plt_id + ascertainment + age + age2 + is_female + region_code"), somascan_non_normalised_log_pheno)
  bmi$soma_non_normal_es_bmi_calc[i] <- summary(model)$coefficients[2,"Estimate"]
  bmi$soma_non_normal_se_bmi_calc[i] <- summary(model)$coefficients[2,"Std. Error"]
  bmi$soma_non_normal_lci_bmi_calc[i] <- confint(model,"bmi_calc")[1]
  bmi$soma_non_normal_hci_bmi_calc[i] <- confint(model,"bmi_calc")[2]
  bmi$soma_non_normal_p_bmi_calc[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]

  print(i)
}

# FDR

bmi$olink_p_fdr_bmi_calc <- p.adjust(bmi$olink_p_bmi_calc,method = "fdr")
bmi$soma_normal_p_fdr_bmi_calc <- p.adjust(bmi$soma_normal_p_bmi_calc,method = "fdr")
bmi$soma_non_normal_p_fdr_bmi_calc <- p.adjust(bmi$soma_non_normal_p_bmi_calc,method = "fdr")

bmi$olink_p_bonferroni_bmi_calc <- p.adjust(bmi$olink_p_bmi_calc,method = "bonferroni")
bmi$soma_normal_p_bonferroni_bmi_calc <- p.adjust(bmi$soma_normal_p_bmi_calc,method = "bonferroni")
bmi$soma_non_normal_p_bonferroni_bmi_calc <- p.adjust(bmi$soma_non_normal_p_bmi_calc,method = "bonferroni")

write.csv(bmi,"bmi.csv", quote=F, row.names=F)



############################### fdr

##### plot volcano

# load

bmi <- read.csv("bmi.csv")

hex <- brewer.pal(3,"Set1")[c(1,2)]

set.seed(4747)

# find top five in each

olink_5_pos <- bmi$olink_id[bmi$olink_es_bmi_calc>0][order(bmi$olink_p_bmi_calc[bmi$olink_es_bmi_calc>0])[1:5]]
soma_normal_5_pos <- bmi$olink_id[bmi$soma_normal_es_bmi_calc>0][order(bmi$soma_normal_p_bmi_calc[bmi$soma_normal_es_bmi_calc>0])[1:5]]
soma_non_normal_5_pos <- bmi$olink_id[bmi$soma_non_normal_es_bmi_calc>0][order(bmi$soma_non_normal_p_bmi_calc[bmi$soma_non_normal_es_bmi_calc>0])[1:5]]
top_pos <- unique(c(olink_5_pos,soma_normal_5_pos,soma_non_normal_5_pos))

olink_3_neg <- bmi$olink_id[bmi$olink_es_bmi_calc<0][order(bmi$olink_p_bmi_calc[bmi$olink_es_bmi_calc<0])[1:3]]
soma_normal_3_neg <- bmi$olink_id[bmi$soma_normal_es_bmi_calc<0][order(bmi$soma_normal_p_bmi_calc[bmi$soma_normal_es_bmi_calc<0])[1:3]]
soma_non_normal_3_neg <- bmi$olink_id[bmi$soma_non_normal_es_bmi_calc<0][order(bmi$soma_non_normal_p_bmi_calc[bmi$soma_non_normal_es_bmi_calc<0])[1:3]]
top_neg <- unique(c(olink_3_neg,soma_normal_3_neg,soma_non_normal_3_neg))

# combine

top <- unique(c(top_pos,top_neg))

# olink

bmi$olink_sig_fdr <- F
bmi$olink_sig_fdr[bmi$olink_p_fdr_bmi_calc<0.05 & sign(bmi$olink_es_bmi_calc)>0] <- "Positive"
bmi$olink_sig_fdr[bmi$olink_p_fdr_bmi_calc<0.05 & sign(bmi$olink_es_bmi_calc)<0] <- "Inverse"
bmi$olink_sig_fdr[bmi$olink_p_fdr_bmi_calc>=0.05] <- "n.s."
bmi$olink_sig_fdr <- factor(bmi$olink_sig_fdr,levels=c("Positive","Inverse","n.s."))
table(bmi$olink_sig_fdr)
n_sig <- sum(bmi$olink_sig_fdr!="n.s.")

# label top ten

olink_top <- bmi[bmi$olink_id %in% top,]
olink_top$olink_p_bmi_calc[-log(olink_top$olink_p_bmi_calc,10)>150] <- 10^(-150)

# data out of limit to present as arrows

olink_x_loc <- bmi$olink_es_bmi_calc[-log(bmi$olink_p_bmi_calc,10)>150]*bmi_sd

olink_threshold <- max(bmi$olink_p_bmi_calc[bmi$olink_p_fdr_bmi_calc<0.05])

volcano_olink <- ggplot(bmi, aes(x=olink_es_bmi_calc*bmi_sd, y=-log(olink_p_bmi_calc,10))) +
  geom_point(aes(color = olink_sig_fdr), size = 0.5) +
  ggtitle("Olink", subtitle = paste("\nn = ",n_sig,"/1694",sep="")) +
  xlab("Effect size per SD increase in BMI") + 
  ylab(expression("-log"[10]*"(p)")) +
  ylim(0,157) +
  xlim(-0.6,0.6) +
  labs(color='Association') +
  scale_color_manual(values=c(hex,"gray")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  geom_segment(aes(x=olink_x_loc, y=150, xend=olink_x_loc, yend=157), 
               arrow = arrow(length=unit(2, 'pt')), color=hex[1],lwd=0.5) +
  geom_text_repel(data = olink_top,
                  mapping = aes(x=olink_es_bmi_calc*bmi_sd, y=-log(olink_p_bmi_calc,10), label=olink_id),
                  size = 2.5) +
  theme(legend.position = "none",text = element_text(size = 14),plot.subtitle=element_text(hjust=0.5)) +
  geom_hline(aes(yintercept=-log(olink_threshold,10)),linetype="dashed")

volcano_olink

# ANML

bmi$soma_normal_sig_fdr <- F
bmi$soma_normal_sig_fdr[bmi$soma_normal_p_fdr_bmi_calc<0.05 & sign(bmi$soma_normal_es_bmi_calc)>0] <- "Positive"
bmi$soma_normal_sig_fdr[bmi$soma_normal_p_fdr_bmi_calc<0.05 & sign(bmi$soma_normal_es_bmi_calc)<0] <- "Inverse"
bmi$soma_normal_sig_fdr[bmi$soma_normal_p_fdr_bmi_calc>=0.05] <- "n.s."
bmi$soma_normal_sig_fdr <- factor(bmi$soma_normal_sig_fdr,levels=c("Positive","Inverse","n.s."))
table(bmi$soma_normal_sig_fdr)
n_sig <- sum(bmi$soma_normal_sig_fdr!="n.s.")

soma_normal_top <- bmi[bmi$olink_id %in% top,]
soma_normal_top$soma_normal_p_bmi_calc[-log(soma_normal_top$soma_normal_p_bmi_calc,10)>150] <- 10^(-150)

soma_normal_x_loc <- bmi$soma_normal_es_bmi_calc[-log(bmi$soma_normal_p_bmi_calc,10)>150]*bmi_sd

soma_normal_threshold <- max(bmi$soma_normal_p_bmi_calc[bmi$soma_normal_p_fdr_bmi_calc<0.05])

volcano_soma_norm <- ggplot(bmi, aes(x=soma_normal_es_bmi_calc*bmi_sd, y=-log(soma_normal_p_bmi_calc,10))) +
  geom_point(aes(color = soma_normal_sig_fdr), size = 0.5) +
  ggtitle("SomaScan-ANML", subtitle = paste("\nn = ",n_sig,"/1694",sep="")) +
  xlab("Effect size per SD increase in BMI") + 
  ylab(expression("-log"[10]*"(p)")) +
  ylim(0,157) +
  xlim(-0.6,0.6) +
  labs(color='Association') +
  scale_color_manual(values=c(hex,"gray")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  geom_segment(aes(x=soma_normal_x_loc, y=150, xend=soma_normal_x_loc, yend=157), 
               arrow = arrow(length=unit(2, 'pt')), color=hex[1],lwd=0.5) +
  geom_text_repel(data = soma_normal_top,
                  mapping = aes(x=soma_normal_es_bmi_calc*bmi_sd, y=-log(soma_normal_p_bmi_calc,10), label=olink_id),
                  size = 2.5) +
  theme(legend.position = "none",text = element_text(size = 14),plot.subtitle=element_text(hjust=0.5)) +
  geom_hline(aes(yintercept=-log(soma_normal_threshold,10)),linetype="dashed")

volcano_soma_norm

# non-ANML

bmi$soma_non_normal_sig_fdr <- F
bmi$soma_non_normal_sig_fdr[bmi$soma_non_normal_p_fdr_bmi_calc<0.05 & sign(bmi$soma_non_normal_es_bmi_calc)>0] <- "Positive"
bmi$soma_non_normal_sig_fdr[bmi$soma_non_normal_p_fdr_bmi_calc<0.05 & sign(bmi$soma_non_normal_es_bmi_calc)<0] <- "Inverse"
bmi$soma_non_normal_sig_fdr[bmi$soma_non_normal_p_fdr_bmi_calc>=0.05] <- "n.s."
bmi$soma_non_normal_sig_fdr <- factor(bmi$soma_non_normal_sig_fdr,levels=c("Positive","Inverse","n.s."))
table(bmi$soma_non_normal_sig_fdr)
n_sig <- sum(bmi$soma_non_normal_sig_fdr!="n.s.")
soma_non_normal_top <- bmi[bmi$olink_id %in% top,]
soma_non_normal_top$soma_non_normal_p_bmi_calc[-log(soma_non_normal_top$soma_non_normal_p_bmi_calc,10)>150] <- 10^(-150)

soma_non_normal_x_loc <- bmi$soma_non_normal_es_bmi_calc[-log(bmi$soma_non_normal_p_bmi_calc,10)>150]*bmi_sd

# remove very small p

bmi$soma_non_normal_p_bmi_calc[bmi$soma_non_normal_p_bmi_calc==0] <- NA

soma_non_normal_threshold <- max(bmi$soma_non_normal_p_bmi_calc[bmi$soma_non_normal_p_fdr_bmi_calc<0.05],na.rm = T)

volcano_soma_non_normal <- ggplot(bmi, aes(x=soma_non_normal_es_bmi_calc*bmi_sd, y=-log(soma_non_normal_p_bmi_calc,10))) +
  geom_point(aes(color = soma_non_normal_sig_fdr), size = 0.5) +
  ggtitle("SomaScan-non-ANML", subtitle = paste("\nn = ",n_sig,"/1694",sep="")) +
  xlab("Effect size per SD increase in BMI") + 
  ylab(expression("-log"[10]*"(p)")) +
  ylim(0,157) +
  xlim(-0.6,0.6) +
  labs(color='Association') +
  scale_color_manual(values=c(hex,"gray")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  geom_segment(aes(x=soma_non_normal_x_loc[1], y=150, xend=soma_non_normal_x_loc[1], yend=157), 
               arrow = arrow(length=unit(2, 'pt')), color=hex[1],lwd=0.5) +
  geom_segment(aes(x=soma_non_normal_x_loc[2], y=150, xend=soma_non_normal_x_loc[2], yend=157), 
               arrow = arrow(length=unit(2, 'pt')), color=hex[1],lwd=0.5) +
  geom_text_repel(data = soma_non_normal_top,
                  mapping = aes(x=soma_non_normal_es_bmi_calc*bmi_sd, y=-log(soma_non_normal_p_bmi_calc,10), label=olink_id),
                  size = 2.5) +
  theme(legend.position = "none",text = element_text(size = 14),plot.subtitle=element_text(hjust=0.5)) +
  geom_hline(aes(yintercept=-log(soma_non_normal_threshold,10)),linetype="dashed")

volcano_soma_non_normal

volcano_combined <- ggarrange(volcano_olink,volcano_soma_norm,volcano_soma_non_normal,ncol=3,nrow=1)

ggsave("",volcano_combined,width=14,height=6)




##### plot scatter plot

bmi <- read.csv("bmi.csv")
overlap_coloc <- read.csv("overlap_coloc_1_to_1.csv")
overlap_coloc$coloc_normal_cis_tf <- ifelse(overlap_coloc$coloc_normal_cis>0,T,F)
overlap_coloc$coloc_non_normal_cis_tf <- ifelse(overlap_coloc$coloc_non_normal_cis>0,T,F)
overlap_coloc <- overlap_coloc[,c("uniprot_id","coloc_normal_cis_tf","coloc_non_normal_cis_tf")]
n_coloc_normal <- sum(overlap_coloc$coloc_normal_cis_tf==T)
n_coloc_non_normal <- sum(overlap_coloc$coloc_non_normal_cis_tf==T)
bmi <- merge(bmi,overlap_coloc,by="uniprot_id")

## olink vs normal

## venn first

cols <- hue_pal()(3)

c1 <- as.numeric(sum(bmi$olink_p_fdr_bmi_calc<0.05))
c2 <- as.numeric(sum(bmi$soma_normal_p_fdr_bmi_calc<0.05))
c3 <- as.numeric(sum(bmi$olink_p_fdr_bmi_calc<0.05 & bmi$soma_normal_p_fdr_bmi_calc<0.05))
c4 <- as.numeric(sum(bmi$olink_p_fdr_bmi_calc<0.05 & bmi$soma_normal_p_fdr_bmi_calc<0.05 & sign(bmi$olink_es_bmi_calc)==sign(bmi$soma_normal_es_bmi_calc)))

venn_diagram_bmi_olink_soma_normal <- venn.diagram(
  list(olink = 1:c1,soma = (c1-c3+1):(c1-c3+c2)),
  filename = NULL,
  fontfamily = "sans",	cat.fontfamily = "sans",
  category.names = c("Olink" , "SomaScan-ANML"),
  height = 1600, 
  width = 1600, 
  cat.pos = c(-30, 10),
  cat.dist = c(0.04, 0.03),
  fill=cols[c(1,2)],
  margin=0.01,
  ext.text=T,
  disable.logging=T)

grid.newpage()
grid.draw(venn_diagram_bmi_olink_soma_normal)

venn_diagram_bmi_olink_soma_normal[[7]]$label <- paste(venn_diagram_bmi_olink_soma_normal[[7]]$label, "\n", c4 ," (", format(round(c4/c3*100,digits=1),nsmall=1), "%) shared", sep="")  

grid.newpage()
grid.draw(venn_diagram_bmi_olink_soma_normal)

png("",width = 1600, height = 1600, pointsize=70)
grid.draw(venn_diagram_bmi_olink_soma_normal)
dev.off()


bmi$olink_soma_normal_sig_fdr <- NA
bmi$olink_soma_normal_sig_fdr[bmi$olink_p_fdr_bmi_calc<0.05 & bmi$soma_normal_p_fdr_bmi_calc<0.05 &
                                       sign(bmi$olink_es_bmi_calc)==sign(bmi$soma_normal_es_bmi_calc)] <- "Shared"
bmi$olink_soma_normal_sig_fdr[bmi$olink_p_fdr_bmi_calc<0.05 & bmi$soma_normal_p_fdr_bmi_calc<0.05 &
                                       sign(bmi$olink_es_bmi_calc)!=sign(bmi$soma_normal_es_bmi_calc)] <- "Opposite direction"
bmi$olink_soma_normal_sig_fdr[bmi$olink_p_fdr_bmi_calc<0.05 & bmi$soma_normal_p_fdr_bmi_calc>=0.05] <- "Olink-specific"
bmi$olink_soma_normal_sig_fdr[bmi$olink_p_fdr_bmi_calc>=0.05 & bmi$soma_normal_p_fdr_bmi_calc<0.05] <- "SomaScan-specific"
bmi$olink_soma_normal_sig_fdr[bmi$olink_p_fdr_bmi_calc>=0.05 & bmi$soma_normal_p_fdr_bmi_calc>=0.05] <- "Neither"
bmi$olink_soma_normal_sig_fdr <- factor(bmi$olink_soma_normal_sig_fdr,levels = c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
table(bmi$olink_soma_normal_sig_fdr)

r_all_normal <- format(round(cor.test(bmi$olink_es_bmi_calc, bmi$soma_normal_es_bmi_calc)$estimate,2),nsmall=2)
r_shared_normal <- format(round(cor.test(bmi$olink_es_bmi_calc[bmi$olink_soma_normal_sig_fdr=="Shared"], bmi$soma_normal_es_bmi_calc[bmi$olink_soma_normal_sig_fdr=="Shared"])$estimate,2),nsmall=2)
r_coloc_normal <- format(round(cor.test(bmi$olink_es_bmi_calc[bmi$coloc_normal_cis_tf==T], bmi$soma_normal_es_bmi_calc[bmi$coloc_normal_cis_tf==T])$estimate,2),nsmall=2)

sum(bmi$coloc_normal_cis_tf==T)
sum(bmi$olink_soma_normal_sig_fdr=="Shared" & bmi$coloc_normal_cis_tf==T)
sum(bmi$olink_soma_normal_sig_fdr=="Shared" & bmi$coloc_normal_cis_tf==T)/sum(bmi$coloc_normal_cis_tf==T)

scatter_normal <- ggplot(bmi, aes(x=olink_es_bmi_calc*bmi_sd, y=soma_normal_es_bmi_calc*bmi_sd, color=olink_soma_normal_sig_fdr, stroke=NA)) +
  geom_errorbar(aes(xmin = olink_lci_bmi_calc*bmi_sd, xmax = olink_hci_bmi_calc*bmi_sd),color="gray80") +
  geom_errorbar(aes(ymin = soma_normal_lci_bmi_calc*bmi_sd, ymax = soma_normal_hci_bmi_calc*bmi_sd),color="gray80") +
  geom_point(size=1.5,alpha=0.5) +
  xlim(-0.6,0.6) +
  ylim(-0.6,0.6) +
  scale_color_manual(values=c("black","gray80","gray80","gray80","gray80")) +
  geom_abline(linetype = "dashed",color="gray60") +
  geom_vline(xintercept=0,linetype = "dashed",color="gray60") +
  geom_hline(yintercept=0,linetype = "dashed",color="gray60") +
  ggtitle("Olink vs SomaScan-ANML", subtitle ="") +
  xlab("Olink") +
  ylab("SomaScan-ANML") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  theme(legend.position = "none") +
  geom_text_repel(data = olink_top,
                  mapping = aes(x=olink_es_bmi_calc*bmi_sd, y=soma_normal_es_bmi_calc*bmi_sd, label=olink_id),
                  size = 2, color="black") +
  annotate("text", x = -0.6, y = 0.6, label = paste("r =", r_all_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.6, y = 0.55, label = paste0("r = ", r_shared_normal, " (shared, n = ", c4,")"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.6, y = 0.5, label = paste0("r = ", r_coloc_normal, " (coloc cis, n = ", n_coloc_normal, ")"), size = 12/.pt, hjust = 0)

scatter_normal

ggsave("",scatter_normal,width=6,height=6,bg = "white")

# olink vs non-normal

# venn

c1 <- as.numeric(sum(bmi$olink_p_fdr_bmi_calc<0.05))
c2 <- as.numeric(sum(bmi$soma_non_normal_p_fdr_bmi_calc<0.05))
c3 <- as.numeric(sum(bmi$olink_p_fdr_bmi_calc<0.05 & bmi$soma_non_normal_p_fdr_bmi_calc<0.05))
c4 <- as.numeric(sum(bmi$olink_p_fdr_bmi_calc<0.05 & bmi$soma_non_normal_p_fdr_bmi_calc<0.05 & sign(bmi$olink_es_bmi_calc)==sign(bmi$soma_non_normal_es_bmi_calc)))

venn_diagram_bmi_olink_soma_non_normal <- venn.diagram(
  list(olink = 1:c1,soma = (c1-c3+1):(c1-c3+c2)),
  filename = NULL,
  fontfamily = "sans",	cat.fontfamily = "sans",
  category.names = c("Olink" , "SomaScan-non-ANML"),
  height = 1600, 
  width = 1600, 
  cat.pos = c(10, -40),
  cat.dist = c(0.04, 0.03),
  fill=cols[c(1,3)],
  inverted=T,
  margin=0.01,
  ext.text=T,
  disable.logging=T)

grid.newpage()
grid.draw(venn_diagram_bmi_olink_soma_non_normal)

venn_diagram_bmi_olink_soma_non_normal[[7]]$label <- paste(venn_diagram_bmi_olink_soma_non_normal[[7]]$label, "\n", c4 ," (", format(round(c4/c3*100,digits=1),nsmall=1), "%) shared", sep="")  

grid.newpage()
grid.draw(venn_diagram_bmi_olink_soma_non_normal)

png("",width = 1600, height = 1600, pointsize=70)
grid.draw(venn_diagram_bmi_olink_soma_non_normal)
dev.off()


bmi$olink_soma_non_normal_sig_fdr <- NA
bmi$olink_soma_non_normal_sig_fdr[bmi$olink_p_fdr_bmi_calc<0.05 & bmi$soma_non_normal_p_fdr_bmi_calc<0.05 &
                                           sign(bmi$olink_es_bmi_calc)==sign(bmi$soma_non_normal_es_bmi_calc)] <- "Shared"
bmi$olink_soma_non_normal_sig_fdr[bmi$olink_p_fdr_bmi_calc<0.05 & bmi$soma_non_normal_p_fdr_bmi_calc<0.05 &
                                           sign(bmi$olink_es_bmi_calc)!=sign(bmi$soma_non_normal_es_bmi_calc)] <- "Opposite direction"
bmi$olink_soma_non_normal_sig_fdr[bmi$olink_p_fdr_bmi_calc<0.05 & bmi$soma_non_normal_p_fdr_bmi_calc>=0.05] <- "Olink-specific"
bmi$olink_soma_non_normal_sig_fdr[bmi$olink_p_fdr_bmi_calc>=0.05 & bmi$soma_non_normal_p_fdr_bmi_calc<0.05] <- "SomaScan-specific"
bmi$olink_soma_non_normal_sig_fdr[bmi$olink_p_fdr_bmi_calc>=0.05 & bmi$soma_non_normal_p_fdr_bmi_calc>=0.05] <- "Neither"
bmi$olink_soma_non_normal_sig_fdr <- factor(bmi$olink_soma_non_normal_sig_fdr,levels = c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
table(bmi$olink_soma_non_normal_sig_fdr)

r_all_non_normal <- format(round(cor.test(bmi$olink_es_bmi_calc, bmi$soma_non_normal_es_bmi_calc)$estimate,2),nsmall=2)
r_shared_non_normal <- format(round(cor.test(bmi$olink_es_bmi_calc[bmi$olink_soma_non_normal_sig_fdr=="Shared"], bmi$soma_non_normal_es_bmi_calc[bmi$olink_soma_non_normal_sig_fdr=="Shared"])$estimate,2),nsmall=2)
r_coloc_non_normal <- format(round(cor.test(bmi$olink_es_bmi_calc[bmi$coloc_non_normal_cis_tf==T], bmi$soma_non_normal_es_bmi_calc[bmi$coloc_non_normal_cis_tf==T])$estimate,2),nsmall=2)

sum(bmi$coloc_non_normal_cis_tf==T)
sum(bmi$olink_soma_non_normal_sig_fdr=="Shared" & bmi$coloc_non_normal_cis_tf==T)
sum(bmi$olink_soma_non_normal_sig_fdr=="Shared" & bmi$coloc_non_normal_cis_tf==T)/sum(bmi$coloc_non_normal_cis_tf==T)

scatter_non_normal <- ggplot(bmi, aes(x=olink_es_bmi_calc*bmi_sd, y=soma_non_normal_es_bmi_calc*bmi_sd, color=olink_soma_non_normal_sig_fdr, stroke=NA)) +
  geom_errorbar(aes(xmin = olink_lci_bmi_calc*bmi_sd, xmax = olink_hci_bmi_calc*bmi_sd),color="gray80") +
  geom_errorbar(aes(ymin = soma_non_normal_lci_bmi_calc*bmi_sd, ymax = soma_non_normal_hci_bmi_calc*bmi_sd),color="gray80") +
  geom_point(size=1.5,alpha=0.5) +
  xlim(-0.6,0.6) +
  ylim(-0.6,0.6) +
  scale_color_manual(values=c("black","gray80","gray80","gray80","gray80")) +
  geom_abline(linetype = "dashed",color="gray60") +
  geom_vline(xintercept=0,linetype = "dashed",color="gray60") +
  geom_hline(yintercept=0,linetype = "dashed",color="gray60") +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") +
  xlab("Olink") +
  ylab("SomaScan-non-ANML") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  theme(legend.position = "none") +
  geom_text_repel(data = olink_top,
                  mapping = aes(x=olink_es_bmi_calc*bmi_sd, y=soma_non_normal_es_bmi_calc*bmi_sd, label=olink_id),
                  size = 2, color="black")  +
  annotate("text", x = -0.6, y = 0.6, label = paste("r =", r_all_non_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.6, y = 0.55, label = paste0("r = ", r_shared_non_normal, " (shared, n = ", c4,")"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.6, y = 0.5, label = paste0("r = ", r_coloc_non_normal, " (coloc cis, n = ", n_coloc_non_normal, ")"), size = 12/.pt, hjust = 0)


scatter_non_normal

ggsave("",scatter_non_normal,width=6,height=6,bg = "white")

scatter_plot <- ggarrange(scatter_normal,scatter_non_normal,ncol=2,nrow=1)


## plot against rho

# all categories

hex <- hue_pal()(10)

show_col(hex)

bmi$olink_soma_normal_sig_fdr <- factor(bmi$olink_soma_normal_sig_fdr,levels = c("Olink-specific","Shared","SomaScan-specific","Opposite direction","Neither"))
bmi$olink_soma_non_normal_sig_fdr <- factor(bmi$olink_soma_non_normal_sig_fdr,levels = c("Olink-specific","Shared","SomaScan-specific","Opposite direction","Neither"))

median_rho_all_normal <- format(round(median(bmi$rho_olink_soma_normal),2),nsmall=2)
median_rho_shared_normal <- format(round(median(bmi$rho_olink_soma_normal[bmi$olink_soma_normal_sig_fdr=="Shared"]),2),nsmall=2)

hist_bmi_normal <- ggplot(bmi, aes(x=rho_olink_soma_normal, fill=olink_soma_normal_sig_fdr)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black") +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  ggtitle("Olink vs SomaScan-ANML", subtitle ="") +
  xlab("Spearman's rho") + 
  ylab("Frequency") +
  labs(fill="Concordance") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_fill_manual(values=c(hex[c(3,5,7,10)],"gray80")) +
  theme(legend.position = "none")  +
  annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_shared_normal, " (shared, n = ", sum(bmi$olink_soma_normal_sig_fdr=="Shared"), ")"), size = 12/.pt, hjust = 0)

hist_bmi_normal

median_rho_all_non_normal <- format(round(median(bmi$rho_olink_soma_non_normal),2),nsmall=2)
median_rho_shared_non_normal <- format(round(median(bmi$rho_olink_soma_non_normal[bmi$olink_soma_non_normal_sig_fdr=="Shared"]),2),nsmall=2)

hist_bmi_non_normal <- ggplot(bmi, aes(x=rho_olink_soma_non_normal, fill=olink_soma_non_normal_sig_fdr)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black") +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") +
  xlab("Spearman's rho") + 
  ylab("Frequency") +
  labs(fill="Concordance") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_fill_manual(values=c(hex[c(3,5,7,10)],"gray80"))  +
  annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_non_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_shared_non_normal, " (shared, n = ", sum(bmi$olink_soma_non_normal_sig_fdr=="Shared"), ")"), size = 12/.pt, hjust = 0)


hist_bmi_non_normal

hist_bmi <- ggarrange(hist_bmi_normal,hist_bmi_non_normal,ncol=2,nrow=1)

ggsave("",hist_bmi,width=14,height=6)



## plot shared only using shade

bmi$olink_soma_normal_shared <- "No"
bmi$olink_soma_normal_shared[bmi$olink_soma_normal_sig_fdr=="Shared"] <- "Yes"
table(bmi$olink_soma_normal_shared)

hist_bmi_shared_normal <- ggplot(bmi, aes(x=rho_olink_soma_normal,fill=olink_soma_normal_shared,pattern=olink_soma_normal_shared)) + 
  geom_histogram_pattern(binwidth = 0.05,
                         boundary=0,
                         pattern_fill = "black", 
                         fill    = 'white',
                         colour  = 'black',
                         pattern_spacing = 0.015,
                         pattern_frequency = 5, 
                         pattern_angle = 45,
                         pattern_density=0.01) +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  ggtitle("Olink vs SomaScan-ANML", subtitle ="") +
  xlab("Spearman's rho") + 
  ylab("Frequency") +
  labs(fill="Shared") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_pattern_manual(values=c("none","stripe")) +
  theme(legend.position = "none") +
  annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_shared_normal, " (shared, n = ", sum(bmi$olink_soma_normal_sig_fdr=="Shared"), ")"), size = 12/.pt, hjust = 0)

hist_bmi_shared_normal

bmi$olink_soma_non_normal_shared <- "No"
bmi$olink_soma_non_normal_shared[bmi$olink_soma_non_normal_sig_fdr=="Shared"] <- "Yes"

hist_bmi_shared_non_normal <- ggplot(bmi, aes(x=rho_olink_soma_non_normal,fill=olink_soma_non_normal_shared,pattern=olink_soma_non_normal_shared)) + 
  geom_histogram_pattern(binwidth = 0.05,
                         boundary=0,
                         pattern_fill = "black", 
                         fill    = 'white',
                         colour  = 'black',
                         pattern_spacing = 0.015,
                         pattern_frequency = 5, 
                         pattern_angle = 45,
                         pattern_density=0.01) +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") +
  xlab("Spearman's rho") + 
  ylab("Frequency") +
  labs(fill="Shared") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_pattern_manual(values=c("none","stripe")) +
  theme(legend.position = "none") +
  annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_non_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_shared_non_normal, " (shared, n = ", sum(bmi$olink_soma_non_normal_sig_fdr=="Shared"), ")"), size = 12/.pt, hjust = 0)

hist_bmi_shared_non_normal

hist_bmi_shared <- ggarrange(hist_bmi_shared_normal,hist_bmi_shared_non_normal,ncol=2,nrow=1)

ggsave("",hist_bmi_shared,width=14,height=6)
