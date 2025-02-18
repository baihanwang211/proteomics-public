##### this script loads proteomic associations with incident IHD and plots results

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
library(ggrepel)
library(VennDiagram)
library(RColorBrewer)
library(readxl)

##### load results

ihd <- read_excel("")

overlap_1_to_1_cor <- read.csv("overlap_1_to_1_cor.csv")

# keep only one-to-one matched
ihd <- ihd[ihd$uniprot_id %in% overlap_1_to_1_cor$uniprot_id,]



##### plot volcano

cols <- hue_pal()(3)

hex <- brewer.pal(3,"Set1")[c(1,2)]

## volcano

set.seed(4747)

olink_4_pos <- ihd$olink_id[ihd$olink_coef>0][order(ihd$olink_p[ihd$olink_coef>0])[1:4]]
soma_non_ANML_4_pos <- ihd$olink_id[ihd$soma_non_ANML_coef>0][order(ihd$soma_non_ANML_p[ihd$soma_non_ANML_coef>0])[1:4]]
top_pos <- unique(c(olink_4_pos,soma_non_ANML_4_pos))


olink_4_neg <- ihd$olink_id[ihd$olink_coef<0][order(ihd$olink_p[ihd$olink_coef<0])[1:4]]
soma_non_ANML_4_neg <- ihd$olink_id[ihd$soma_non_ANML_coef<0][order(ihd$soma_non_ANML_p[ihd$soma_non_ANML_coef<0])[1:4]]
top_neg <- unique(c(olink_4_neg,soma_non_ANML_4_neg))

top <- ihd[ihd$olink_id %in% unique(c(top_pos,top_neg)),]

# olink

ihd$olink_sig_fdr <- F
ihd$olink_sig_fdr[ihd$olink_p_fdr<0.05 & sign(ihd$olink_coef)>0] <- "Positive"
ihd$olink_sig_fdr[ihd$olink_p_fdr<0.05 & sign(ihd$olink_coef)<0] <- "Inverse"
ihd$olink_sig_fdr[ihd$olink_p_fdr>=0.05] <- "n.s."
ihd$olink_sig_fdr <- factor(ihd$olink_sig_fdr,levels=c("Positive","Inverse","n.s."))
table(ihd$olink_sig_fdr)
n_sig <- sum(ihd$olink_sig_fdr!="n.s.")

olink_threshold <- max(ihd$olink_p[ihd$olink_p_fdr<0.05])

volcano_olink <- ggplot(ihd, aes(x=olink_coef, y=-log(olink_p,10))) +
  geom_point(aes(color = olink_sig_fdr), size = 0.5) +
  ggtitle("Olink", subtitle = paste("\nn = ",n_sig,"/1694",sep="")) +
  xlab("Hazard ratio") + 
  ylab(expression("-log"[10]*"(p)")) +
  ylim(0,15) +
  scale_x_continuous(labels=c("0.75","1.00","1.25","1.50"), breaks=c(log(0.75),log(1),log(1.25),log(1.5)), limits=c(log(0.7),log(1.8))) +
  labs(color='Association') +
  scale_color_manual(values=c(hex,"gray")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  # geom_segment(aes(x=olink_x_loc, y=150, xend=olink_x_loc, yend=155), 
  #              arrow = arrow(length=unit(2, 'pt')), color=hex[1],lwd=1) +
  geom_text_repel(data = top,
                  mapping = aes(x=olink_coef, y=-log(olink_p,10), label=olink_id),
                  size = 3.5) +
  theme(legend.position = "none",text = element_text(size = 14),plot.subtitle=element_text(hjust=0.5)) +
  geom_hline(aes(yintercept=-log(olink_threshold,10)),linetype="dashed")

volcano_olink

# soma non normal

ihd$soma_non_ANML_sig_fdr <- F
ihd$soma_non_ANML_sig_fdr[ihd$soma_non_ANML_p_fdr<0.05 & sign(ihd$soma_non_ANML_coef)>0] <- "Positive"
ihd$soma_non_ANML_sig_fdr[ihd$soma_non_ANML_p_fdr<0.05 & sign(ihd$soma_non_ANML_coef)<0] <- "Inverse"
ihd$soma_non_ANML_sig_fdr[ihd$soma_non_ANML_p_fdr>=0.05] <- "n.s."
ihd$soma_non_ANML_sig_fdr <- factor(ihd$soma_non_ANML_sig_fdr,levels=c("Positive","Inverse","n.s."))
table(ihd$soma_non_ANML_sig_fdr)
n_sig <- sum(ihd$soma_non_ANML_sig_fdr!="n.s.")

soma_non_ANML_threshold <- max(ihd$soma_non_ANML_p[ihd$soma_non_ANML_p_fdr<0.05])

volcano_soma_non_ANML <- ggplot(ihd, aes(x=soma_non_ANML_coef, y=-log(soma_non_ANML_p,10))) +
  geom_point(aes(color = soma_non_ANML_sig_fdr), size = 0.5) +
  ggtitle("SomaScan", subtitle = paste("\nn = ",n_sig,"/1694",sep="")) +
  xlab("Hazard ratio") + 
  ylab(expression("-log"[10]*"(p)")) +
  ylim(0,15) +
  scale_x_continuous(labels=c("0.75","1.00","1.25","1.50"), breaks=c(log(0.75),log(1),log(1.25),log(1.5)), limits=c(log(0.7),log(1.8))) +
  labs(color='Association') +
  scale_color_manual(values=c(hex,"gray")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  # geom_segment(aes(x=soma_non_ANML_x_loc, y=150, xend=soma_non_ANML_x_loc, yend=155), 
  #              arrow = arrow(length=unit(2, 'pt')), color=hex[1],lwd=1) +
  geom_text_repel(data = top,
                  mapping = aes(x=soma_non_ANML_coef, y=-log(soma_non_ANML_p,10), label=olink_id),
                  size = 3.5) +
  theme(legend.position = "none",text = element_text(size = 14),plot.subtitle=element_text(hjust=0.5)) +
  geom_hline(aes(yintercept=-log(soma_non_ANML_threshold,10)),linetype="dashed")

volcano_soma_non_ANML

volcano_combined_figure <- ggarrange(volcano_olink,volcano_soma_non_ANML,ncol=2,nrow=1)



###### venn

overlap_pqtl <- read.csv("")
overlap_pqtl$coloc_non_ANML_cis_tf <- ifelse(overlap_pqtl$coloc_non_ANML_cis>0,T,F)
overlap_pqtl <- overlap_pqtl[,c("uniprot_id","coloc_non_ANML_cis_tf")]
n_coloc_non_ANML <- sum(overlap_pqtl$coloc_non_ANML_cis_tf==T)
ihd <- merge(ihd,overlap_pqtl,by="uniprot_id")

c1 <- as.numeric(sum(ihd$olink_p_fdr<0.05))
c2 <- as.numeric(sum(ihd$soma_non_ANML_p_fdr<0.05))
c3 <- as.numeric(sum(ihd$olink_p_fdr<0.05 & ihd$soma_non_ANML_p_fdr<0.05))
c4 <- as.numeric(sum(ihd$olink_p_fdr<0.05 & ihd$soma_non_ANML_p_fdr<0.05 & sign(ihd$olink_coef)==sign(ihd$soma_non_ANML_coef)))

venn_ihd_non_ANML <- venn.diagram(
  list(olink_sig=1:c1, soma_non_ANML_sig=(c1-c3+1):(c1-c3+c2)),
  filename = NULL,
  disable.logging =T,
  fontfamily = "sans",	cat.fontfamily = "sans",
  category.names = c("Olink", "SomaScan"),
  height = 1600, 
  width = 1600, 
  cat.pos = c(-25, 19),
  cat.dist = c(0.04, 0.04),
  fill=cols[c(1,3)],
  cat.cex = 0.8,
  cex = 0.8)

grid.newpage()
grid.draw(venn_ihd_non_ANML)

venn_ihd_non_ANML[[7]]$label <- paste(venn_ihd_non_ANML[[7]]$label, "\n(", c4, " shared)", sep="")  

grid.newpage()
grid.draw(venn_ihd_non_ANML)

##### scatter

ihd$olink_soma_non_ANML_sig_fdr <- NA
ihd$olink_soma_non_ANML_sig_fdr[ihd$olink_p_fdr<0.05 & ihd$soma_non_ANML_p_fdr<0.05 &
                                                          sign(ihd$olink_coef)==sign(ihd$soma_non_ANML_coef)] <- "Shared"
ihd$olink_soma_non_ANML_sig_fdr[ihd$olink_p_fdr<0.05 & ihd$soma_non_ANML_p_fdr<0.05 &
                                                          sign(ihd$olink_coef)!=sign(ihd$soma_non_ANML_coef)] <- "Opposite direction"
ihd$olink_soma_non_ANML_sig_fdr[ihd$olink_p_fdr<0.05 & ihd$soma_non_ANML_p_fdr>=0.05] <- "Olink-specific"
ihd$olink_soma_non_ANML_sig_fdr[ihd$olink_p_fdr>=0.05 & ihd$soma_non_ANML_p_fdr<0.05] <- "SomaScan-specific"
ihd$olink_soma_non_ANML_sig_fdr[ihd$olink_p_fdr>=0.05 & ihd$soma_non_ANML_p_fdr>=0.05] <- "Neither"
ihd$olink_soma_non_ANML_sig_fdr <- factor(ihd$olink_soma_non_ANML_sig_fdr,levels = c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
table(ihd$olink_soma_non_ANML_sig_fdr)

r_all_non_ANML <- format(round(cor.test(ihd$olink_coef, ihd$soma_non_ANML_coef)$estimate,2),nsmall=2)
r_shared_non_ANML <- format(round(cor.test(ihd$olink_coef[ihd$olink_soma_non_ANML_sig_fdr=="Shared"], 
                                             ihd$soma_non_ANML_coef[ihd$olink_soma_non_ANML_sig_fdr=="Shared"])$estimate,2),nsmall=2)
r_coloc_non_ANML <- format(round(cor.test(ihd$olink_coef[ihd$coloc_non_ANML_cis_tf==T], 
                                            ihd$soma_non_ANML_coef[ihd$coloc_non_ANML_cis_tf==T])$estimate,2),nsmall=2)

sum(ihd$coloc_non_ANML_cis_tf==T)
sum(ihd$olink_soma_non_ANML_sig_fdr=="Shared" & ihd$coloc_non_ANML_cis_tf==T)
sum(ihd$olink_soma_non_ANML_sig_fdr=="Shared" & ihd$coloc_non_ANML_cis_tf==T)/sum(ihd$coloc_non_ANML_cis_tf==T)

scatter_non_ANML <- ggplot(ihd, aes(x=olink_coef, y=soma_non_ANML_coef, color=olink_soma_non_ANML_sig_fdr, stroke=NA)) +
  geom_errorbar(aes(xmin = olink_lci, xmax = olink_hci),color="gray80") +
  geom_errorbar(aes(ymin = soma_non_ANML_lci, ymax = soma_non_ANML_hci),color="gray80") +
  geom_point(size=1.5,alpha=0.5) +
  xlim(-0.6,0.6) +
  ylim(-0.6,0.6) +
  scale_color_manual(values=c("black","gray80","gray80","gray80","gray80")) +
  geom_abline(linetype = "dashed",color="gray60") +
  geom_vline(xintercept=0,linetype = "dashed",color="gray60") +
  geom_hline(yintercept=0,linetype = "dashed",color="gray60") +
  # ggtitle("Olink vs SomaScan", subtitle ="") +
  xlab("Olink") +
  ylab("SomaScan") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  theme(legend.position = "none") +
  geom_text_repel(data = top,
                  mapping = aes(x=olink_coef, y=soma_non_ANML_coef, label=olink_id),
                  size = 3.5, color="black") +
  annotate("text", x = -0.6, y = 0.6, label = paste("r =", r_all_non_ANML, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.6, y = 0.55, label = paste0("r = ", r_shared_non_ANML, " (shared, n = ", c4,")"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.6, y = 0.5, label = paste0("r = ", r_coloc_non_ANML, " (coloc cis, n = ", n_coloc_non_ANML, ")"), size = 12/.pt, hjust = 0)

scatter_non_ANML





##### plot against rho

ihd$olink_soma_non_ANML_shared <- "No"
ihd$olink_soma_non_ANML_shared[ihd$olink_soma_non_ANML_sig_fdr=="Shared"] <- "Yes"

hist_ihd_shared_non_ANML <- ggplot(ihd, aes(x=rho_olink_soma_non_ANML,fill=olink_soma_non_ANML_shared,pattern=olink_soma_non_ANML_shared)) + 
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
  # ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") +
  xlab("Spearman's rho") + 
  ylab("Frequency") +
  labs(fill="Shared") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_pattern_manual(values=c("none","stripe")) +
  theme(legend.position = "none") +
  annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_non_ANML, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_shared_non_ANML, " (shared, n = ", sum(ihd$olink_soma_non_ANML_sig_fdr=="Shared"), ")"), size = 12/.pt, hjust = 0)

hist_ihd_shared_non_ANML
