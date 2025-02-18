##### this script loads proteomic associations with BMI and plots results

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
library(ckbplotr)

##### load results

assoc_pheno <- read_excel("")

overlap_1_to_1_cor <- read.csv("overlap_1_to_1_cor.csv")

bmi <- assoc_pheno[,c("uniprot_id", "olink_id", "somascan_id", "olink_es_bmi", "olink_se_bmi", "olink_p_bmi", "soma_non_ANML_es_bmi", "soma_non_ANML_se_bmi", "soma_non_ANML_p_bmi")]

# keep only one-to-one matched
bmi <- bmi[bmi$uniprot_id %in% overlap_1_to_1_cor$uniprot_id,]

##### plot volcano

# SD of BMI in CKB
bmi_sd <- 3.641146

hex <- brewer.pal(3,"Set1")[c(1,2)]

set.seed(4747)

# find top five in each

olink_5_pos <- bmi$olink_id[bmi$olink_es_bmi>0][order(bmi$olink_p_bmi[bmi$olink_es_bmi>0])[1:5]]
soma_non_ANML_5_pos <- bmi$olink_id[bmi$soma_non_ANML_es_bmi>0][order(bmi$soma_non_ANML_p_bmi[bmi$soma_non_ANML_es_bmi>0])[1:5]]
top_pos <- unique(c(olink_5_pos,soma_non_ANML_5_pos))

olink_3_neg <- bmi$olink_id[bmi$olink_es_bmi<0][order(bmi$olink_p_bmi[bmi$olink_es_bmi<0])[1:3]]
soma_non_ANML_3_neg <- bmi$olink_id[bmi$soma_non_ANML_es_bmi<0][order(bmi$soma_non_ANML_p_bmi[bmi$soma_non_ANML_es_bmi<0])[1:3]]
top_neg <- unique(c(olink_3_neg,soma_non_ANML_3_neg))

# combine

top <- unique(c(top_pos,top_neg))

# olink

bmi$olink_sig_fdr <- F
bmi$olink_sig_fdr[bmi$olink_p_fdr_bmi<0.05 & sign(bmi$olink_es_bmi)>0] <- "Positive"
bmi$olink_sig_fdr[bmi$olink_p_fdr_bmi<0.05 & sign(bmi$olink_es_bmi)<0] <- "Inverse"
bmi$olink_sig_fdr[bmi$olink_p_fdr_bmi>=0.05] <- "n.s."
bmi$olink_sig_fdr <- factor(bmi$olink_sig_fdr,levels=c("Positive","Inverse","n.s."))
table(bmi$olink_sig_fdr)
n_sig <- sum(bmi$olink_sig_fdr!="n.s.")

# label top ten

olink_top <- bmi[bmi$olink_id %in% top,]
olink_top$olink_p_bmi[-log(olink_top$olink_p_bmi,10)>150] <- 10^(-150)

# data out of limit to present as arrows

olink_x_loc <- bmi$olink_es_bmi[-log(bmi$olink_p_bmi,10)>150]*bmi_sd

olink_threshold <- max(bmi$olink_p_bmi[bmi$olink_p_fdr_bmi<0.05])

volcano_olink <- ggplot(bmi, aes(x=olink_es_bmi*bmi_sd, y=-log(olink_p_bmi,10))) +
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
                  mapping = aes(x=olink_es_bmi*bmi_sd, y=-log(olink_p_bmi,10), label=olink_id),
                  size = 3.5) +
  theme(legend.position = "none",text = element_text(size = 12),plot.subtitle=element_text(hjust=0.5)) +
  geom_hline(aes(yintercept=-log(olink_threshold,10)),linetype="dashed")

volcano_olink

# non-ANML

bmi$soma_non_ANML_sig_fdr <- F
bmi$soma_non_ANML_sig_fdr[bmi$soma_non_ANML_p_fdr_bmi<0.05 & sign(bmi$soma_non_ANML_es_bmi)>0] <- "Positive"
bmi$soma_non_ANML_sig_fdr[bmi$soma_non_ANML_p_fdr_bmi<0.05 & sign(bmi$soma_non_ANML_es_bmi)<0] <- "Inverse"
bmi$soma_non_ANML_sig_fdr[bmi$soma_non_ANML_p_fdr_bmi>=0.05] <- "n.s."
bmi$soma_non_ANML_sig_fdr <- factor(bmi$soma_non_ANML_sig_fdr,levels=c("Positive","Inverse","n.s."))
table(bmi$soma_non_ANML_sig_fdr)
n_sig <- sum(bmi$soma_non_ANML_sig_fdr!="n.s.")

soma_non_ANML_top <- bmi[bmi$olink_id %in% top,]
soma_non_ANML_top$soma_non_ANML_p_bmi[-log(soma_non_ANML_top$soma_non_ANML_p_bmi,10)>150] <- 10^(-150)

soma_non_ANML_x_loc <- bmi$soma_non_ANML_es_bmi[-log(bmi$soma_non_ANML_p_bmi,10)>150]*bmi_sd

# remove very small p

bmi$soma_non_ANML_p_bmi[bmi$soma_non_ANML_p_bmi==0] <- NA

soma_non_ANML_threshold <- max(bmi$soma_non_ANML_p_bmi[bmi$soma_non_ANML_p_fdr_bmi<0.05],na.rm = T)

volcano_soma_non_ANML_figure <- ggplot(bmi, aes(x=soma_non_ANML_es_bmi*bmi_sd, y=-log(soma_non_ANML_p_bmi,10))) +
  geom_point(aes(color = soma_non_ANML_sig_fdr), size = 0.5) +
  ggtitle("SomaScan", subtitle = paste("\nn = ",n_sig,"/1694",sep="")) +
  xlab("Effect size per SD increase in BMI") + 
  ylab(expression("-log"[10]*"(p)")) +
  ylim(0,157) +
  xlim(-0.6,0.6) +
  labs(color='Association') +
  scale_color_manual(values=c(hex,"gray")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  geom_segment(aes(x=soma_non_ANML_x_loc[1], y=150, xend=soma_non_ANML_x_loc[1], yend=157), 
               arrow = arrow(length=unit(2, 'pt')), color=hex[1],lwd=0.5) +
  geom_segment(aes(x=soma_non_ANML_x_loc[2], y=150, xend=soma_non_ANML_x_loc[2], yend=157), 
               arrow = arrow(length=unit(2, 'pt')), color=hex[1],lwd=0.5) +
  geom_text_repel(data = soma_non_ANML_top,
                  mapping = aes(x=soma_non_ANML_es_bmi*bmi_sd, y=-log(soma_non_ANML_p_bmi,10), label=olink_id),
                  size = 3.5) +
  theme(legend.position = "none",text = element_text(size = 12),plot.subtitle=element_text(hjust=0.5)) +
  geom_hline(aes(yintercept=-log(soma_non_ANML_threshold,10)),linetype="dashed")

volcano_soma_non_ANML_figure

volcano_combined_figure <- ggarrange(volcano_olink,volcano_soma_non_ANML_figure,ncol=2,nrow=1)




##### venn

overlap_pqtl <- read.csv("")
overlap_pqtl$coloc_non_ANML_cis_tf <- ifelse(overlap_pqtl$coloc_non_ANML_cis>0,T,F)
overlap_pqtl <- overlap_pqtl[,c("uniprot_id","coloc_non_ANML_cis_tf")]
n_coloc_non_ANML <- sum(overlap_pqtl$coloc_non_ANML_cis_tf==T)
bmi <- merge(bmi,overlap_pqtl,by="uniprot_id")

c1 <- as.numeric(sum(bmi$olink_p_fdr_bmi<0.05))
c2 <- as.numeric(sum(bmi$soma_non_ANML_p_fdr_bmi<0.05))
c3 <- as.numeric(sum(bmi$olink_p_fdr_bmi<0.05 & bmi$soma_non_ANML_p_fdr_bmi<0.05))
c4 <- as.numeric(sum(bmi$olink_p_fdr_bmi<0.05 & bmi$soma_non_ANML_p_fdr_bmi<0.05 & sign(bmi$olink_es_bmi)==sign(bmi$soma_non_ANML_es_bmi)))

venn_diagram_bmi_olink_soma_non_ANML <- venn.diagram(
  list(olink = 1:c1,soma = (c1-c3+1):(c1-c3+c2)),
  filename = NULL,
  fontfamily = "sans",	cat.fontfamily = "sans",
  category.names = c("Olink" , "SomaScan"),
  height = 1600, 
  width = 1600, 
  cat.pos = c(10, -40),
  cat.dist = c(0.04, 0.03),
  fill=cols[c(1,3)],
  inverted=T,
  margin=0.01,
  ext.text=T,
  disable.logging=T,
  cat.cex = 0.8,
  cex = 0.8)

grid.newpage()
grid.draw(venn_diagram_bmi_olink_soma_non_ANML)

venn_diagram_bmi_olink_soma_non_ANML[[7]]$label <- paste(venn_diagram_bmi_olink_soma_non_ANML[[7]]$label, "\n", c4 ," (", format(round(c4/c3*100,digits=1),nsmall=1), "%) shared", sep="")  

grid.newpage()
grid.draw(venn_diagram_bmi_olink_soma_non_ANML)

##### scatter

bmi$olink_soma_non_ANML_sig_fdr <- NA
bmi$olink_soma_non_ANML_sig_fdr[bmi$olink_p_fdr_bmi<0.05 & bmi$soma_non_ANML_p_fdr_bmi<0.05 &
                                           sign(bmi$olink_es_bmi)==sign(bmi$soma_non_ANML_es_bmi)] <- "Shared"
bmi$olink_soma_non_ANML_sig_fdr[bmi$olink_p_fdr_bmi<0.05 & bmi$soma_non_ANML_p_fdr_bmi<0.05 &
                                           sign(bmi$olink_es_bmi)!=sign(bmi$soma_non_ANML_es_bmi)] <- "Opposite direction"
bmi$olink_soma_non_ANML_sig_fdr[bmi$olink_p_fdr_bmi<0.05 & bmi$soma_non_ANML_p_fdr_bmi>=0.05] <- "Olink-specific"
bmi$olink_soma_non_ANML_sig_fdr[bmi$olink_p_fdr_bmi>=0.05 & bmi$soma_non_ANML_p_fdr_bmi<0.05] <- "SomaScan-specific"
bmi$olink_soma_non_ANML_sig_fdr[bmi$olink_p_fdr_bmi>=0.05 & bmi$soma_non_ANML_p_fdr_bmi>=0.05] <- "Neither"
bmi$olink_soma_non_ANML_sig_fdr <- factor(bmi$olink_soma_non_ANML_sig_fdr,levels = c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
table(bmi$olink_soma_non_ANML_sig_fdr)

r_all_non_ANML <- format(round(cor.test(bmi$olink_es_bmi, bmi$soma_non_ANML_es_bmi)$estimate,2),nsmall=2)
r_shared_non_ANML <- format(round(cor.test(bmi$olink_es_bmi[bmi$olink_soma_non_ANML_sig_fdr=="Shared"], bmi$soma_non_ANML_es_bmi[bmi$olink_soma_non_ANML_sig_fdr=="Shared"])$estimate,2),nsmall=2)
r_coloc_non_ANML <- format(round(cor.test(bmi$olink_es_bmi[bmi$coloc_non_ANML_cis_tf==T], bmi$soma_non_ANML_es_bmi[bmi$coloc_non_ANML_cis_tf==T])$estimate,2),nsmall=2)

sum(bmi$coloc_non_ANML_cis_tf==T)
sum(bmi$olink_soma_non_ANML_sig_fdr=="Shared" & bmi$coloc_non_ANML_cis_tf==T)
sum(bmi$olink_soma_non_ANML_sig_fdr=="Shared" & bmi$coloc_non_ANML_cis_tf==T)/sum(bmi$coloc_non_ANML_cis_tf==T)

scatter_non_ANML <- ggplot(bmi, aes(x=olink_es_bmi*bmi_sd, y=soma_non_ANML_es_bmi*bmi_sd, color=olink_soma_non_ANML_sig_fdr, stroke=NA)) +
  geom_errorbar(aes(xmin = olink_lci_bmi*bmi_sd, xmax = olink_hci_bmi*bmi_sd),color="gray80") +
  geom_errorbar(aes(ymin = soma_non_ANML_lci_bmi*bmi_sd, ymax = soma_non_ANML_hci_bmi*bmi_sd),color="gray80") +
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
  geom_text_repel(data = olink_top,
                  mapping = aes(x=olink_es_bmi*bmi_sd, y=soma_non_ANML_es_bmi*bmi_sd, label=olink_id),
                  size = 3.5, color="black")  +
  annotate("text", x = -0.6, y = 0.6, label = paste("r =", r_all_non_ANML, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.6, y = 0.55, label = paste0("r = ", r_shared_non_ANML, " (shared, n = ", c4,")"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.6, y = 0.5, label = paste0("r = ", r_coloc_non_ANML, " (coloc cis, n = ", n_coloc_non_ANML, ")"), size = 12/.pt, hjust = 0)

scatter_non_ANML




##### plot against rho

bmi$olink_soma_non_ANML_shared <- "No"
bmi$olink_soma_non_ANML_shared[bmi$olink_soma_non_ANML_sig_fdr=="Shared"] <- "Yes"

hist_bmi_shared_non_ANML <- ggplot(bmi, aes(x=rho_olink_soma_non_ANML,fill=olink_soma_non_ANML_shared,pattern=olink_soma_non_ANML_shared)) + 
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
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_shared_non_ANML, " (shared, n = ", sum(bmi$olink_soma_non_ANML_sig_fdr=="Shared"), ")"), size = 12/.pt, hjust = 0)

hist_bmi_shared_non_ANML
