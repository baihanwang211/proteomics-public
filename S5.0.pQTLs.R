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
library(tidyr)

#### correlation of p values

p <- read.delim("logP.txt")
overlap <- read.csv("overlap_1_to_1_cor.csv")
names(p)[c(2,3)] <- c("olink_id","somascan_id")
setdiff(p$somascan_id,overlap$somascan_id)

# restrict to one-to-one

overlap_p <- merge(overlap,p,by="somascan_id")
names(overlap_p)
overlap_logp <- overlap_p
overlap_logp$olinklogp <- -log10(overlap_logp$OlinkP)
overlap_logp$anmllogp <- -log10(overlap_logp$anmlP)
overlap_logp$nanmllogp <- -log10(overlap_logp$nanmlP)

overlap_logp <- overlap_logp[is.finite(overlap_logp$olinklogp)&is.finite(overlap_logp$anmllogp)&is.finite(overlap_logp$nanmllogp),]

hex <- brewer.pal(3,"Set1")[c(1,2)]

scatter_p_anml <- ggplot(overlap_logp, aes(x=olinklogp, y=anmllogp, color=cis_trans, stroke=NA)) +
  geom_point(size=1.5,alpha=0.5) +
  xlim(0,330) +
  ylim(0,330) +
  scale_color_manual(values=c(hex)) +
  geom_abline(linetype = "dashed",color="gray60") +
  ggtitle("OLINK vs SomaScan-ANML", subtitle ="") +
  xlab(expression("-log"[10]*"(p) OLINK")) +
  ylab(expression("-log"[10]*"(p) SomaScan-ANML")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  theme(legend.position = "none")

scatter_p_non_anml <- ggplot(overlap_logp, aes(x=olinklogp, y=anmllogp, color=cis_trans, stroke=NA)) +
  geom_point(size=1.5,alpha=0.5) +
  xlim(0,330) +
  ylim(0,330) +
  scale_color_manual(values=c(hex)) +
  geom_abline(linetype = "dashed",color="gray60") +
  ggtitle("OLINK vs SomaScan-non-ANML", subtitle ="") +
  xlab(expression("-log"[10]*"(p) OLINK")) +
  ylab(expression("-log"[10]*"(p) SomaScan-non-ANML")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  theme(legend.position = "none")

scatter_p_soma <- ggplot(overlap_logp, aes(x=anmllogp, y=nanmllogp, color=cis_trans, stroke=NA)) +
  geom_point(size=1.5,alpha=0.5) +
  xlim(0,330) +
  ylim(0,330) +
  scale_color_manual(values=c(hex)) +
  geom_abline(linetype = "dashed",color="gray60") +
  ggtitle("SomaScan ANML vs non-ANML", subtitle ="") +
  xlab(expression("-log"[10]*"(p) SomaScan-ANML")) +
  ylab(expression("-log"[10]*"(p) SomaScan-non-ANML")) +
  labs(color='') +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12))

scatter_p <- ggarrange(scatter_p_anml,scatter_p_non_anml,scatter_p_soma,ncol=3,nrow=1)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/pqtl_scatter_p.png",scatter_p,width=18,height=6,bg = "white")


#### heritability

# olink_h2 <- read.delim("olink_h2.txt")
# names(olink_h2)[-1] <- paste('olink', names(olink_h2)[-1], sep = '_')
# names(olink_h2) <- tolower(names(olink_h2))
# names(olink_h2)[1] <- "olink_id"
# overlap_h2 <- merge(overlap,olink_h2,by="olink_id",all.x=T)
# 
# soma_normal_h2 <- read.delim("soma_anml_h2.txt")
# soma_normal_h2 <- soma_normal_h2[!is.na(soma_normal_h2$cis_trans),]
# dup <- soma_normal_h2[duplicated(soma_normal_h2[,1]),1]
# dup_eg <- soma_normal_h2[soma_normal_h2$Aptamer_ID %in% dup,]
# 
# soma_normal_h2 <- pivot_wider(soma_normal_h2,id_cols="Aptamer_ID",names_from="cis_trans",values_from = "h2")
# soma_normal_h2[soma_normal_h2=="NULL"] <- 0
# 
# soma_normal_h2 <- replace(soma_normal_h2,.=="NULL",0)
# 
# names(soma_normal_h2)[-1] <- paste('soma_normal', names(soma_normal_h2)[-1], sep = '_')
# names(soma_normal_h2) <- tolower(names(soma_normal_h2))
# names(soma_normal_h2)[1] <- "somascan_id"
# soma_normal_h2$soma_normal_total <- soma_normal_h2$soma_normal_trans + soma_normal_h2$soma_normal_cis
# 
# soma_non_normal_h2 <- read.delim("soma_nanml_h2.txt")
# soma_non_normal_h2 <- soma_non_normal_h2[!is.na(soma_non_normal_h2$cis_trans),]
# soma_non_normal_h2 <- pivot_wider(soma_non_normal_h2,id_cols="Aptamer_ID",names_from="cis_trans",values_from = "h2")
# 
# 
# 
# 
# 
# names(soma_normal_h2)[-1] <- paste('soma_normal', names(soma_normal_h2)[-1], sep = '_')
# names(soma_normal_h2) <- tolower(names(soma_normal_h2))
# names(soma_normal_h2)[1] <- "somascan_id"
# overlap_h2 <- merge(overlap_h2,soma_normal_h2,by="somascan_id",all.x=T)
# 
# 
# soma_non_normal_h2 <- read.delim("soma_nanml_h2.txt")




#### number of pQTLs

pqtl_normal <- read.delim("hit_count_anml (1).txt")
names(pqtl_normal)
pqtl_normal <- pqtl_normal[,c(2,6:9)]
names(pqtl_normal) <- c("somascan_id","olink_cis","olink_trans","soma_normal_cis","soma_normal_trans")

pqtl_non_normal <- read.delim("hit_count_non_anml (1).txt")
names(pqtl_non_normal)
pqtl_non_normal <- pqtl_non_normal[,c(2,6:9)]
names(pqtl_non_normal) <- c("somascan_id","olink_cis","olink_trans","soma_non_normal_cis","soma_non_normal_trans")

pqtl <- merge(pqtl_normal,pqtl_non_normal,by=c("somascan_id","olink_cis","olink_trans"))

# pqtl <- read.delim("overlap_pqtls.txt")
# names(pqtl)[1] <- "somascan_id"
# pqtl <- pqtl[,-c(2,3)]
# names(pqtl) <- c("somascan_id","olink_cis","olink_trans","soma_normal_cis","soma_normal_trans","soma_non_normal_cis","soma_non_normal_trans")

overlap_pqtl <- merge(overlap,pqtl,by="somascan_id")

# pqtl_normal <- read.delim("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/hit_count_anml.txt")
# names(pqtl_normal)
# pqtl_normal <- pqtl_normal[,c("soma","olink_cis","olink_trans","soma_cis","soma_trans")]
# names(pqtl_normal) <- c("somascan_id","olink_cis","olink_trans","soma_normal_cis","soma_normal_trans")
# 
# pqtl_non_normal <- read.delim("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/hit_count_non_anml.txt")
# names(pqtl_non_normal)
# pqtl_non_normal <- pqtl_non_normal[,c("soma","soma_cis","soma_trans")]
# names(pqtl_non_normal) <- c("somascan_id","soma_non_normal_cis","soma_non_normal_trans")
# 
# pqtl <- merge(pqtl_normal,pqtl_non_normal,by="somascan_id")
# names(pqtl)

## count total number of pQTLs

names(overlap_pqtl)
pqtl_count <- data.frame(colSums(overlap_pqtl[,c((ncol(overlap_pqtl)-5):ncol(overlap_pqtl))]))
names(pqtl_count) <- "Frequency"
pqtl_count$Platform <- c("OLINK","OLINK","SomaScan-ANML","SomaScan-ANML","SomaScan-non-ANML","SomaScan-non-ANML")
pqtl_count$pQTL <- c("Cis-pQTL","Trans-pQTL","Cis-pQTL","Trans-pQTL","Cis-pQTL","Trans-pQTL")

overlap_pqtl_tf <- overlap_pqtl
overlap_pqtl_tf[,c((ncol(overlap_pqtl_tf)-5):ncol(overlap_pqtl_tf))] <- T

for (i in (ncol(overlap_pqtl_tf)-5):ncol(overlap_pqtl_tf)) {
  overlap_pqtl_tf[overlap_pqtl[,i]==0,i] <- F
}

pqtl_tf_count <- data.frame(colSums(overlap_pqtl_tf[,c((ncol(overlap_pqtl_tf)-5):ncol(overlap_pqtl_tf))]))
names(pqtl_tf_count) <- "Frequency"
pqtl_tf_count$Platform <- c("OLINK","OLINK","SomaScan-ANML","SomaScan-ANML","SomaScan-non-ANML","SomaScan-non-ANML")
pqtl_tf_count$pQTL <- c("Cis-pQTL","Trans-pQTL","Cis-pQTL","Trans-pQTL","Cis-pQTL","Trans-pQTL")
pqtl_tf_count_add <- pqtl_tf_count
pqtl_tf_count_add$Frequency <- pqtl_count$Frequency - pqtl_tf_count$Frequency

pqtl_tf_count$cat<- "original"
pqtl_tf_count_add$cat<- "additional"

pqtl_count_all <- rbind(pqtl_tf_count,pqtl_tf_count_add)

pqtl_count_all$Frequency_all <- NA
pqtl_count_all$Frequency_all[pqtl_count_all$cat=="additional"] <- pqtl_count$Frequency

pqtl_count_all$Frequency_shade <- NA
pqtl_count_all$Frequency_shade[pqtl_count_all$cat=="original"] <- pqtl_count_all$Frequency[pqtl_count_all$cat=="original"] 

# plot stacked bar together

hex <- hue_pal()(3)

# plot_pqtl_hit_all <- ggplot(pqtl_count_all,aes(y=Frequency, x=Platform, fill=Platform, pattern=cat)) +
#   geom_bar_pattern(stat = "identity",
#                    position = "stack",
#                    colour = 'black',
#                    pattern_fill = "black",
#                    pattern_spacing = 0.015,
#                    pattern_frequency = 5, 
#                    pattern_angle = 45,
#                    pattern_density=0.01) +
#   facet_wrap( ~ pQTL) +
#   scale_fill_manual(values = hex) +
#   scale_pattern_manual(values=c("none","stripe")) +
#   xlab("") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
#   theme(text = element_text(size = 10)) +
#   theme(legend.position = "none") +
#   theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
# 
# plot_pqtl_hit_all

plot_pqtl_hit_all <- ggplot(pqtl_count_all,aes(y=Frequency, x=fct_rev(Platform), fill=fct_rev(Platform), pattern=cat)) +
  geom_bar_pattern(stat = "identity",
                   position = "stack",
                   colour = 'black',
                   pattern_fill = "black",
                   pattern_spacing = 0.03,
                   pattern_frequency = 5,
                   pattern_angle = 45,
                   pattern_density=0.01) +
  ylim(0,1600) +
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

plot_pqtl_hit_all

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_pqtl_hit_all.png",plot_pqtl_hit_all,width=9,height=6)



# ## count pQTLs per protein
# 
# cols <- brewer.pal(3,"Set1")
# 
# # olink
# 
# overlap_pqtl$olink_all <- overlap_pqtl$olink_cis + overlap_pqtl$olink_trans
# overlap_pqtl$olink_cat <- NA
# overlap_pqtl$olink_cat[overlap_pqtl$olink_cis!=0 & overlap_pqtl$olink_trans==0] <- "cis"
# overlap_pqtl$olink_cat[overlap_pqtl$olink_cis==0 & overlap_pqtl$olink_trans!=0] <- "trans"
# overlap_pqtl$olink_cat[overlap_pqtl$olink_cis!=0 & overlap_pqtl$olink_trans!=0] <- "both"
# overlap_pqtl$olink_cat[overlap_pqtl$olink_cis==0 & overlap_pqtl$olink_trans==0] <- "no pQTLs"
# overlap_pqtl$olink_cat <- factor(overlap_pqtl$olink_cat,levels=c("cis","both","trans","no pQTLs"))
# table(overlap_pqtl$olink_cat)
# 
# overlap_pqtl$soma_normal_all <- overlap_pqtl$soma_normal_cis + overlap_pqtl$soma_normal_trans
# overlap_pqtl$soma_normal_cat <- NA
# overlap_pqtl$soma_normal_cat[overlap_pqtl$soma_normal_cis!=0 & overlap_pqtl$soma_normal_trans==0] <- "cis"
# overlap_pqtl$soma_normal_cat[overlap_pqtl$soma_normal_cis==0 & overlap_pqtl$soma_normal_trans!=0] <- "trans"
# overlap_pqtl$soma_normal_cat[overlap_pqtl$soma_normal_cis!=0 & overlap_pqtl$soma_normal_trans!=0] <- "both"
# overlap_pqtl$soma_normal_cat[overlap_pqtl$soma_normal_cis==0 & overlap_pqtl$soma_normal_trans==0] <- "no pQTLs"
# overlap_pqtl$soma_normal_cat <- factor(overlap_pqtl$soma_normal_cat,levels=c("cis","both","trans","no pQTLs"))
# table(overlap_pqtl$soma_normal_cat)
# 
# overlap_pqtl$soma_non_normal_all <- overlap_pqtl$soma_non_normal_cis + overlap_pqtl$soma_non_normal_trans
# overlap_pqtl$soma_non_normal_cat <- NA
# overlap_pqtl$soma_non_normal_cat[overlap_pqtl$soma_non_normal_cis!=0 & overlap_pqtl$soma_non_normal_trans==0] <- "cis"
# overlap_pqtl$soma_non_normal_cat[overlap_pqtl$soma_non_normal_cis==0 & overlap_pqtl$soma_non_normal_trans!=0] <- "trans"
# overlap_pqtl$soma_non_normal_cat[overlap_pqtl$soma_non_normal_cis!=0 & overlap_pqtl$soma_non_normal_trans!=0] <- "both"
# overlap_pqtl$soma_non_normal_cat[overlap_pqtl$soma_non_normal_cis==0 & overlap_pqtl$soma_non_normal_trans==0] <- "no pQTLs"
# overlap_pqtl$soma_non_normal_cat <- factor(overlap_pqtl$soma_non_normal_cat,levels=c("cis","both","trans","no pQTLs"))
# table(overlap_pqtl$soma_non_normal_cat)
# 
# # olink
# 
# plot_pqtl_per_protein_olink <- ggplot(overlap_pqtl,aes(x=olink_all, fill=fct_rev(olink_cat))) +
#   geom_bar(stat = "count", position = "stack") +
#   ylim(0,800) +
#   xlab("Total number of pQTLs per protein") +
#   ylab("Number of proteins") + 
#   ggtitle("OLINK", subtitle ="") +
#   scale_fill_manual(values = c("gray",cols[c(2,3,1)])) +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
#   theme(text = element_text(size = 12)) +
#   theme(legend.position = "none")
# 
# # somascan anml
# 
# plot_pqtl_per_protein_soma_normal <- ggplot(overlap_pqtl,aes(x=soma_normal_all, fill=fct_rev(soma_normal_cat))) +
#   geom_bar(stat = "count", position = "stack") +
#   ylim(0,800) +
#   xlab("Total number of pQTLs per protein") +
#   ylab("Number of proteins") + 
#   ggtitle("SomaScan-ANML", subtitle ="") +
#   scale_fill_manual(values = c("gray",cols[c(2,3,1)])) +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
#   theme(text = element_text(size = 12)) +
#   theme(legend.position = "none")
# 
# # somascan non-anml
# 
# plot_pqtl_per_protein_soma_non_normal <- ggplot(overlap_pqtl,aes(x=soma_non_normal_all, fill=fct_rev(soma_non_normal_cat))) +
#   geom_bar(stat = "count", position = "stack") +
#   ylim(0,800) +
#   xlab("Total number of pQTLs per protein") +
#   ylab("Number of proteins") + 
#   ggtitle("SomaScan-non-ANML", subtitle ="") +
#   scale_fill_manual(name="",
#                     values = c("gray",cols[c(2,3,1)]),
#                     labels= c("no pQTLs","trans only", "both cis and trans", "cis only")) +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
#   theme(text = element_text(size = 12)) 
# 
# plot_pqtl_per_protein <- ggarrange(plot_pqtl_per_protein_olink,plot_pqtl_per_protein_soma_normal,plot_pqtl_per_protein_soma_non_normal,ncol=3,nrow=1)
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_pqtl_per_protein.png",plot_pqtl_per_protein,width=18,height=6,bg = "white")


# hex <- hue_pal()(3)
# 
# plot_pqtl_hit <- ggplot(pqtl_count, aes(y=Frequency, x=pQTL, fill=Platform)) + 
#   geom_bar(stat="identity",position="dodge",colour="black") +
#   xlab("") +
#   ylab("Number of pQTLs") +
#   labs(fill = "Platform") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
#   theme(text = element_text(size = 10)) +
#   theme(strip.placement = "outside") +
#   scale_fill_manual(values=hex)
# 
# plot_pqtl_hit
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_pqtl_hit.png",plot_pqtl_hit,width=7,height=6)

# # Whether or not having a pQTL
# 
# hex <- hue_pal()(3)
# 
# plot_pqtl_tf_hit <- ggplot(pqtl_tf_count, aes(y=Frequency, x=pQTL, fill=Platform)) + 
#   geom_bar(stat="identity",position="dodge",colour="black") +
#   xlab("") +
#   ylab("Number of protein targets with pQTL(s)") +
#   labs(fill = "Platform") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
#   theme(text = element_text(size = 10)) +
#   theme(strip.placement = "outside") +
#   scale_fill_manual(values=hex)
# 
# plot_pqtl_tf_hit
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_pqtl_tf_hit.png",plot_pqtl_tf_hit,width=7,height=6)




# load coloc results

coloc_normal <- read.delim("coloc_anml (1).txt")
names(coloc_normal)
coloc_normal <- coloc_normal[,c("soma","shared_cis","shared_trans")]
names(coloc_normal) <- c("somascan_id","coloc_normal_cis","coloc_normal_trans")
overlap_coloc <- merge(overlap_pqtl,coloc_normal,by="somascan_id",all.x=T)

coloc_non_normal <- read.delim("coloc_nanml (1).txt")
names(coloc_non_normal)
coloc_non_normal <- coloc_non_normal[,c("soma","shared_cis","shared_trans")]
names(coloc_non_normal) <- c("somascan_id","coloc_non_normal_cis","coloc_non_normal_trans")
overlap_coloc <- merge(overlap_coloc,coloc_non_normal,by="somascan_id",all.x=T)

coloc_soma <- read.delim("coloc_anml_naml (1).txt")
names(coloc_soma)
coloc_soma <- coloc_soma[,c("soma","shared_cis","shared_trans")]
names(coloc_soma) <- c("somascan_id","coloc_soma_cis","coloc_soma_trans")
overlap_coloc <- merge(overlap_coloc,coloc_soma,by="somascan_id",all.x=T)

# venn diagrams

# cols <- brewer.pal(n = 3, name = "Set1")[c(1,2)]

# cis normal

overlap_coloc$coloc_normal_cis[is.na(overlap_coloc$coloc_normal_cis)] <- 0
pqtl_tf_count
c1 <- pqtl_tf_count["olink_cis",]$Frequency
c2 <- pqtl_tf_count["soma_normal_cis",]$Frequency
c3 <- as.numeric(sum(overlap_coloc$olink_cis>0 & overlap_coloc$soma_normal_cis>0))
c4 <- as.numeric(sum(overlap_coloc$olink_cis>0 & overlap_coloc$soma_normal_cis>0 & overlap_coloc$coloc_normal_cis>0))
c5 <- as.numeric(sum(overlap_coloc$olink_cis>0 & overlap_coloc$soma_normal_cis==0 & overlap_coloc$coloc_normal_cis>0))
c6 <- as.numeric(sum(overlap_coloc$olink_cis==0 & overlap_coloc$soma_normal_cis>0 & overlap_coloc$coloc_normal_cis>0))
c7 <- as.numeric(sum(overlap_coloc$olink_cis==0 & overlap_coloc$soma_normal_cis==0 & overlap_coloc$coloc_normal_cis>0))
overlap_coloc$olink_id[overlap_coloc$olink_cis==0 & overlap_coloc$soma_normal_cis==0 & overlap_coloc$coloc_normal_cis>0]
# overlap_coloc$somascan_id[overlap_coloc$olink_cis==0 & overlap_coloc$soma_normal_cis==0 & overlap_coloc$coloc_normal_cis>0]

overlap_coloc$coloc_normal_cis_tf <- F
overlap_coloc$coloc_normal_cis_tf[overlap_coloc$coloc_normal_cis!=0] <- T
# overlap_coloc$coloc_normal_cis_tf[overlap_coloc$olink_cis==0 & overlap_coloc$soma_normal_cis==0 & overlap_coloc$coloc_normal_cis>0] <- F
table(overlap_coloc$coloc_normal_cis_tf)

venn_diagram_pqtl_olink_soma_normal <- venn.diagram(
      list(olink = 1:c1,soma = (c1-c3+1):(c1-c3+c2)),
      filename = NULL,
      fontfamily = "sans",	cat.fontfamily = "sans",
      category.names = c("OLINK" , "SomaScan-ANML"),
      height = 1600, 
      width = 1600, 
      cat.pos = c(-40, 28),
      cat.dist = c(0.03, 0.04),
      # margin=0.01,
      fill=hex[c(1,2)],
      ext.text=T,
      disable.logging=T)

grid.newpage()
grid.draw(venn_diagram_pqtl_olink_soma_normal)

venn_diagram_pqtl_olink_soma_normal[[5]]$label <- paste(venn_diagram_pqtl_olink_soma_normal[[5]]$label, 
                                                        "\n", c5, " (", format(round((c5/(c1-c3))*100,digits=1),nsmall=1), "%)\ncolocalise", sep="")  
venn_diagram_pqtl_olink_soma_normal[[6]]$label <- paste(venn_diagram_pqtl_olink_soma_normal[[6]]$label,
                                                        "\n", c6, " (", format(round((c6/(c2-c3))*100,digits=1),nsmall=1), "%)\ncolocalise", sep="")  
venn_diagram_pqtl_olink_soma_normal[[7]]$label <- paste(venn_diagram_pqtl_olink_soma_normal[[7]]$label,
                                                        "\n", c4, " (", format(round(c4/c3*100,digits=1),nsmall=1), "%) colocalise", sep="")  

grid.newpage()
grid.draw(venn_diagram_pqtl_olink_soma_normal)

png("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/venn_diagram_pqtl_olink_soma_normal.png",width = 1600, height = 1600, pointsize=45)
grid.draw(venn_diagram_pqtl_olink_soma_normal)
dev.off()


# cis non-normal

overlap_coloc$coloc_non_normal_cis[is.na(overlap_coloc$coloc_non_normal_cis)] <- 0
pqtl_tf_count
c1 <- pqtl_tf_count["olink_cis",]$Frequency
c2 <- pqtl_tf_count["soma_non_normal_cis",]$Frequency
c3 <- as.numeric(sum(overlap_coloc$olink_cis>0 & overlap_coloc$soma_non_normal_cis>0))
c4 <- as.numeric(sum(overlap_coloc$olink_cis>0 & overlap_coloc$soma_non_normal_cis>0 & overlap_coloc$coloc_non_normal_cis>0))
c5 <- as.numeric(sum(overlap_coloc$olink_cis>0 & overlap_coloc$soma_non_normal_cis==0 & overlap_coloc$coloc_non_normal_cis>0))
c6 <- as.numeric(sum(overlap_coloc$olink_cis==0 & overlap_coloc$soma_non_normal_cis>0 & overlap_coloc$coloc_non_normal_cis>0))
c7 <- as.numeric(sum(overlap_coloc$olink_cis==0 & overlap_coloc$soma_non_normal_cis==0 & overlap_coloc$coloc_non_normal_cis>0))
overlap_coloc$olink_id[overlap_coloc$olink_cis==0 & overlap_coloc$soma_non_normal_cis==0 & overlap_coloc$coloc_non_normal_cis>0]

overlap_coloc$coloc_non_normal_cis_tf <- F
overlap_coloc$coloc_non_normal_cis_tf[overlap_coloc$coloc_non_normal_cis!=0] <- T
# overlap_coloc$coloc_non_normal_cis_tf[overlap_coloc$olink_cis==0 & overlap_coloc$soma_non_normal_cis==0 & overlap_coloc$coloc_non_normal_cis>0] <- F
table(overlap_coloc$coloc_non_normal_cis_tf)

venn_diagram_pqtl_olink_soma_non_normal <- venn.diagram(
  list(olink = 1:c1,soma = (c1-c3+1):(c1-c3+c2)),
  filename = NULL,
  fontfamily = "sans",	cat.fontfamily = "sans",
  category.names = c("OLINK" , "SomaScan-non-ANML"),
  height = 1600, 
  width = 1600, 
  cat.pos = c(-40, 28),
  cat.dist = c(0.03, 0.04),
  # margin=0.01,
  fill=hex[c(1,3)],
  ext.text=T,
  disable.logging=T)

grid.newpage()
grid.draw(venn_diagram_pqtl_olink_soma_non_normal)

venn_diagram_pqtl_olink_soma_non_normal[[5]]$label <- paste(venn_diagram_pqtl_olink_soma_non_normal[[5]]$label, 
                                                        "\n", c5, " (", format(round((c5/(c1-c3))*100,digits=1),nsmall=1), "%)\ncolocalise", sep="")  
venn_diagram_pqtl_olink_soma_non_normal[[6]]$label <- paste(venn_diagram_pqtl_olink_soma_non_normal[[6]]$label,
                                                        "\n", c6, " (", format(round((c6/(c2-c3))*100,digits=1),nsmall=1), "%)\ncolocalise", sep="")  
venn_diagram_pqtl_olink_soma_non_normal[[7]]$label <- paste(venn_diagram_pqtl_olink_soma_non_normal[[7]]$label,
                                                        "\n", c4, " (", format(round(c4/c3*100,digits=1),nmsall=1), "%) colocalise", sep="")  

grid.newpage()
grid.draw(venn_diagram_pqtl_olink_soma_non_normal)

png("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/venn_diagram_pqtl_olink_soma_non_normal.png",width = 1600, height = 1600, pointsize=45)
grid.draw(venn_diagram_pqtl_olink_soma_non_normal)
dev.off()


# soma normal vs non-normal

overlap_coloc$coloc_soma_cis[is.na(overlap_coloc$coloc_soma_cis)] <- 0
pqtl_tf_count
c1 <- pqtl_tf_count["soma_normal_cis",]$Frequency
c2 <- pqtl_tf_count["soma_non_normal_cis",]$Frequency
c3 <- as.numeric(sum(overlap_coloc$soma_normal_cis>0 & overlap_coloc$soma_non_normal_cis>0))
c4 <- as.numeric(sum(overlap_coloc$soma_normal_cis>0 & overlap_coloc$soma_non_normal_cis>0 & overlap_coloc$coloc_soma_cis>0))
c5 <- as.numeric(sum(overlap_coloc$soma_normal_cis>0 & overlap_coloc$soma_non_normal_cis==0 & overlap_coloc$coloc_soma_cis>0))
c6 <- as.numeric(sum(overlap_coloc$soma_normal_cis==0 & overlap_coloc$soma_non_normal_cis>0 & overlap_coloc$coloc_soma_cis>0))
c7 <- as.numeric(sum(overlap_coloc$soma_normal_cis==0 & overlap_coloc$soma_non_normal_cis==0 & overlap_coloc$coloc_soma_cis>0))
overlap_coloc$somascan_id[overlap_coloc$soma_normal_cis==0 & overlap_coloc$soma_non_normal_cis==0 & overlap_coloc$coloc_soma_cis>0]

overlap_coloc$coloc_soma_cis_tf <- F
overlap_coloc$coloc_soma_cis_tf[overlap_coloc$coloc_soma_cis!=0] <- T
# overlap_coloc$coloc_soma_cis_tf[overlap_coloc$soma_normal_cis==0 & overlap_coloc$soma_non_normal_cis==0 & overlap_coloc$coloc_soma_cis>0] <- F
table(overlap_coloc$coloc_soma_cis_tf)

venn_diagram_pqtl_soma <- venn.diagram(
  list(normal = 1:c1, non_normal = (c1-c3+1):(c1-c3+c2)),
  filename = NULL,
  fontfamily = "sans",	cat.fontfamily = "sans",
  category.names = c("SomaScan\n(ANML)" , "SomaScan\n(non-ANML)"),
  height = 1600, 
  width = 1600, 
  cat.pos = c(-120, 120),
  cat.dist = c(0.07, 0.07),
  margin=0.05,
  fill=hex[c(2,3)],
  ext.text=T,
  ext.percent=c(0.5,0.5,0.5),
  ext.length=0.8,
  disable.logging=T)

grid.newpage()
grid.draw(venn_diagram_pqtl_soma)

venn_diagram_pqtl_soma[[5]]$label <- paste(venn_diagram_pqtl_soma[[5]]$label, 
                                                            "\n", c5, " (", format(round((c5/(c1-c3))*100,digits=1),nsmall=1), "%) colocalise", sep="")  
venn_diagram_pqtl_soma[[7]]$label <- paste(venn_diagram_pqtl_soma[[7]]$label,
                                                            "\n", c6, " (", format(round((c6/(c2-c3))*100,digits=1),nsmall=1), "%) colocalise", sep="")  
venn_diagram_pqtl_soma[[9]]$label <- paste(venn_diagram_pqtl_soma[[9]]$label,
                                                            "\n", c4, " (", format(round(c4/c3*100,digits=1),nsmall=1), "%) colocalise", sep="")  

grid.newpage()
grid.draw(venn_diagram_pqtl_soma)

png("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/venn_diagram_pqtl_soma.png",width = 1600, height = 1600, pointsize=45)
grid.draw(venn_diagram_pqtl_soma)
dev.off()

write.csv(overlap_coloc,"overlap_coloc.csv", quote=F, row.names=F)


# trans and save

# overlap_coloc$coloc_normal_trans[is.na(overlap_coloc$coloc_normal_trans)] <- 0
# overlap_coloc$coloc_normal_trans_tf <- F
# overlap_coloc$coloc_normal_trans_tf[overlap_coloc$coloc_normal_trans!=0] <- T
# table(overlap_coloc$coloc_normal_trans_tf)
# 
# overlap_coloc$coloc_non_normal_trans[is.na(overlap_coloc$coloc_non_normal_trans)] <- 0
# overlap_coloc$coloc_non_normal_trans_tf <- F
# overlap_coloc$coloc_non_normal_trans_tf[overlap_coloc$coloc_non_normal_trans!=0] <- T
# table(overlap_coloc$coloc_non_normal_trans_tf)

# write.csv(overlap_coloc,"overlap_coloc.csv", quote=F, row.names=F)
# overlap_coloc <- read.csv("overlap_coloc.csv")


# histogram

## plot shared only using shade

median_rho_all_normal <- format(round(median(overlap_coloc$rho_olink_soma_normal),2),nsmall=2)
median_rho_coloc_normal <- format(round(median(overlap_coloc$rho_olink_soma_normal[overlap_coloc$coloc_normal_cis_tf==T]),2),nsmall=2)

n_coloc_normal <- sum(overlap_coloc$coloc_normal_cis_tf==T)
n_coloc_non_normal <- sum(overlap_coloc$coloc_non_normal_cis_tf==T)

hist_pqtl_normal <- ggplot(overlap_coloc, aes(x=rho_olink_soma_normal,pattern=coloc_normal_cis_tf)) + 
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
  ggtitle("OLINK vs SomaScan-ANML", subtitle ="") +
  xlab("Spearman's rho") + 
  ylab("Frequency") +
  labs(fill="Shared") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_pattern_manual(values=c("none","stripe")) +
  theme(legend.position = "none") +
  annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_coloc_normal, " (coloc cis, n = ", n_coloc_normal, ")"), size = 12/.pt, hjust = 0)

hist_pqtl_normal

median_rho_all_non_normal <- format(round(median(overlap_coloc$rho_olink_soma_non_normal),2),nsmall=2)
median_rho_coloc_non_normal <- format(round(median(overlap_coloc$rho_olink_soma_non_normal[overlap_coloc$coloc_non_normal_cis_tf==T]),2),nsmall=2)

hist_pqtl_non_normal <- ggplot(overlap_coloc, aes(x=rho_olink_soma_non_normal,pattern=coloc_non_normal_cis_tf)) + 
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
  ggtitle("OLINK vs SomaScan-non-ANML", subtitle ="") +
  xlab("Spearman's rho") + 
  ylab("Frequency") +
  labs(fill="Shared") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_pattern_manual(values=c("none","stripe")) +
  theme(legend.position = "none") +
  annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_non_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_coloc_non_normal, " (coloc cis, n = ", n_coloc_non_normal, ")"), size = 12/.pt, hjust = 0)

hist_pqtl_non_normal

hist_pqtl <- ggarrange(hist_pqtl_normal,hist_pqtl_non_normal,ncol=2,nrow=1)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_pqtl.png",hist_pqtl,width=14,height=6)

write.csv(overlap_coloc,"overlap_coloc.csv", quote=F, row.names=F)
