rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(egg)
library(VennDiagram)
library(RColorBrewer)

overlap_assoc <- read.csv("overlap_assoc.csv")

### ihd

overlap_assoc$olink_id <- gsub("\\.","_",overlap_assoc$olink_id)

olink_result <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/Overlap, MM/Overlap/IHD/data used/OLINK_BothB_FinalModel.csv")
names(olink_result)
olink_result <- olink_result[,c(3,4,7)]
names(olink_result) <- c("olink_id","olink_es_ihd","olink_p_ihd")
olink_result <- olink_result[!duplicated(olink_result$olink_id),]
olink_result$olink_id <- gsub("-","_",olink_result$olink_id)
# test <- merge(overlap_assoc, olink_result, by="olink_id")
# setdiff(overlap_assoc$uniprot_id,test$uniprot_id)
# setdiff(test$uniprot_id,overlap_assoc$uniprot_id)
overlap_assoc <- merge(overlap_assoc, olink_result, by="olink_id")

soma_normal_result <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/Overlap, MM/Overlap/IHD/data used/Total_Normalized.csv")
names(soma_normal_result)
soma_normal_result <- soma_normal_result[,c(1,2,7)]
names(soma_normal_result) <- c("somascan_id","soma_normal_es_ihd","soma_normal_p_ihd")
overlap_assoc <- merge(overlap_assoc, soma_normal_result, by="somascan_id")

soma_non_normal_result <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/Overlap, MM/Overlap/IHD/data used/Total_UNnormalized.csv")
soma_non_normal_result <- soma_non_normal_result[,c(1,2,7)]
names(soma_non_normal_result) <- c("somascan_id","soma_non_normal_es_ihd","soma_non_normal_p_ihd")
overlap_assoc <- merge(overlap_assoc, soma_non_normal_result, by="somascan_id")

overlap_assoc$olink_sig_ihd <- F
overlap_assoc$olink_sig_ihd[overlap_assoc$olink_p_ihd<0.05] <- T

overlap_assoc$soma_normal_sig_ihd <- F
overlap_assoc$soma_normal_sig_ihd[overlap_assoc$soma_normal_p_ihd<0.05] <- T

overlap_assoc$soma_non_normal_sig_ihd <- F
overlap_assoc$soma_non_normal_sig_ihd[overlap_assoc$soma_non_normal_p_ihd<0.05] <- T

## restrict to 1-to-1

overlap_1_to_1_assoc <- overlap_assoc[!(duplicated(overlap_assoc$uniprot_id) | duplicated(overlap_assoc$uniprot_id, fromLast = TRUE) | !is.na(overlap_assoc$somascan_uniprot_2)), ]

names(overlap_1_to_1_assoc)

colSums(overlap_1_to_1_assoc[,111:113])


## get number of hits

overlap_1_to_1_assoc$sig_ihd_normal <- "n.s."
overlap_1_to_1_assoc$sig_ihd_normal[overlap_1_to_1_assoc$olink_sig_ihd==T] <- "sig. Olink"
overlap_1_to_1_assoc$sig_ihd_normal[overlap_1_to_1_assoc$soma_normal_sig_ihd==T] <- "sig. SomaScan"
overlap_1_to_1_assoc$sig_ihd_normal[overlap_1_to_1_assoc$olink_sig_ihd==T & overlap_1_to_1_assoc$soma_normal_sig_ihd==T] <- "sig. both"
overlap_1_to_1_assoc$sig_ihd_normal <- factor(overlap_1_to_1_assoc$sig_ihd_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
table(overlap_1_to_1_assoc$sig_ihd_normal)

cols <- brewer.pal(n = 3, name = "Set1")[c(1,2)]

venn.diagram(
  list(Olink = 1:519, SomaScan = 359:774),
  filename = 'K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/venn_diagram_IHD_olink_soma_normal.png',
  disable.logging = T,
  category.names = c("Olink" , "SomaScan (ANML)"),
  height = 1280, 
  width = 1280, 
  resolution = 300,
  cat.pos = c(-27, 20),
  cat.dist = c(0.03, 0.03),
  fontfamily = "sans",	cat.fontfamily = "sans",
  fill=cols)

overlap_1_to_1_assoc$sig_ihd_non_normal <- "n.s."
overlap_1_to_1_assoc$sig_ihd_non_normal[overlap_1_to_1_assoc$olink_sig_ihd==T] <- "sig. Olink"
overlap_1_to_1_assoc$sig_ihd_non_normal[overlap_1_to_1_assoc$soma_non_normal_sig_ihd==T] <- "sig. SomaScan"
overlap_1_to_1_assoc$sig_ihd_non_normal[overlap_1_to_1_assoc$olink_sig_ihd==T & overlap_1_to_1_assoc$soma_non_normal_sig_ihd==T] <- "sig. both"
overlap_1_to_1_assoc$sig_ihd_non_normal <- factor(overlap_1_to_1_assoc$sig_ihd_non_normal, levels=c("n.s.","sig. Olink","sig. SomaScan","sig. both"))
table(overlap_1_to_1_assoc$sig_ihd_non_normal)

venn.diagram(
  list(Olink = 1:519, SomaScan = 330:805),
  filename = 'K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/venn_diagram_IHD_olink_soma_non_normal.png',
  disable.logging = T,
  category.names = c("Olink" , "SomaScan (non-ANML)"),
  height = 1280, 
  width = 1280, 
  resolution = 300,
  cat.pos = c(-27, 19),
  cat.dist = c(0.04, 0.04),
  fontfamily = "sans",	cat.fontfamily = "sans",
  fill=cols)

overlap_1_to_1_assoc$sig_ihd_within_soma <- "n.s."
overlap_1_to_1_assoc$sig_ihd_within_soma[overlap_1_to_1_assoc$soma_normal_sig_ihd==T] <- "sig. SomaScan (ANML)"
overlap_1_to_1_assoc$sig_ihd_within_soma[overlap_1_to_1_assoc$soma_non_normal_sig_ihd==T] <- "sig. SomaScan (non-ANML)"
overlap_1_to_1_assoc$sig_ihd_within_soma[overlap_1_to_1_assoc$soma_normal_sig_ihd==T & overlap_1_to_1_assoc$soma_non_normal_sig_ihd==T] <- "sig. both"
overlap_1_to_1_assoc$sig_ihd_within_soma <- factor(overlap_1_to_1_assoc$sig_ihd_within_soma, levels=c("n.s.","sig. SomaScan (ANML)","sig. SomaScan (non-ANML)","sig. both"))
table(overlap_1_to_1_assoc$sig_ihd_within_soma)

venn.diagram(
  list(Soma1 = 1:416, Soma2 = 31:506),
  filename = 'K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/venn_diagram_IHD_within_soma.png',
  disable.logging = T,
  category.names = c("SomaScan (ANML)" , "SomaScan (non-ANML)"),
  inverted = T,
  height = 1280, 
  width = 1280, 
  resolution = 300,
  cat.pos = c(27, -27),
  cat.dist = c(0.04, 0.04),
  fontfamily = "sans",	cat.fontfamily = "sans",
  fill=cols)

## scatter

# beta

scatter_ihd_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_ihd, y=soma_normal_es_ihd,stroke=NA)) + 
  geom_point(size=3,alpha = 0.3) +
  xlim(-1.3,1.3) +
  ylim(-1.3,1.3) +
  geom_abline(linetype = "dashed") +
  ggtitle("ANML") +
  xlab("Olink beta") +
  ylab("SomaScan beta") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20))

scatter_ihd_normal

cor.test(overlap_1_to_1_assoc$olink_es_ihd,overlap_1_to_1_assoc$soma_normal_es_ihd)

# ihd non-normal

scatter_ihd_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=olink_es_ihd, y=soma_non_normal_es_ihd,stroke=NA)) + 
  geom_point(size=3,alpha = 0.3) +
  xlim(-1.3,1.3) +
  ylim(-1.3,1.3) +
  geom_abline(linetype = "dashed") +
  ggtitle("non-ANML") +
  xlab("Olink beta") +
  ylab("SomaScan beta") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20))

scatter_ihd_non_normal

cor.test(overlap_1_to_1_assoc$olink_es_ihd,overlap_1_to_1_assoc$soma_non_normal_es_ihd)

scatter_ihd <- ggarrange(scatter_ihd_normal,scatter_ihd_non_normal,ncol=2,nrow=1)

scatter_ihd <- annotate_figure(scatter_ihd, top = text_grob("Ischemic heart disease", face = "bold", size = 28))

scatter_ihd

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/scatter_ihd.png",scatter_ihd,width=16,height=8, bg = "white")

# histogram ihd

hist_ihd_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_normal,fill=sig_ihd_normal)) + 
  geom_histogram(color="black") +
  xlim(-0.4,1) +
  ylim(0,400) +
  ggtitle("ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  scale_fill_discrete(name = "Significance") +
  theme_few()

hist_ihd_normal

hist_ihd_non_normal <- ggplot(overlap_1_to_1_assoc, aes(x=rho_olink_soma_non_normal,fill=sig_ihd_non_normal)) + 
  geom_histogram(color="black") +
  xlim(-0.4,1) +
  ylim(0,400) +
  ggtitle("Non-ANML") + xlab("Spearman correlation coefficient") + ylab("Frequency") +
  scale_fill_discrete(name = "Significance") +
  theme_few()

hist_ihd_non_normal

hist_ihd <- ggarrange(hist_ihd_normal,hist_ihd_non_normal,ncol=2,nrow=1)

hist_ihd <- annotate_figure(hist_ihd, top = text_grob("Ischemic heart disease", face = "bold", size = 28))

hist_ihd

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_ihd.png",hist_ihd,width=16,height=6,bg="white")


