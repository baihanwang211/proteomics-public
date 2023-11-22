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

ggsave("pqtl_scatter_p.png",scatter_p,width=18,height=6,bg = "white")



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

overlap_pqtl <- merge(overlap,pqtl,by="somascan_id")

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

ggsave("plot_pqtl_hit_all.png",plot_pqtl_hit_all,width=9,height=6)


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

overlap_coloc$coloc_normal_cis_tf <- F
overlap_coloc$coloc_normal_cis_tf[overlap_coloc$coloc_normal_cis!=0] <- T
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

png("venn_diagram_pqtl_olink_soma_normal.png",width = 1600, height = 1600, pointsize=45)
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

png("venn_diagram_pqtl_olink_soma_non_normal.png",width = 1600, height = 1600, pointsize=45)
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

png("venn_diagram_pqtl_soma.png",width = 1600, height = 1600, pointsize=45)
grid.draw(venn_diagram_pqtl_soma)
dev.off()

write.csv(overlap_coloc,"overlap_coloc.csv", quote=F, row.names=F)



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

ggsave("hist_pqtl.png",hist_pqtl,width=14,height=6)

write.csv(overlap_coloc,"overlap_coloc.csv", quote=F, row.names=F)
