### this script loads pQTL and colocalisation results from GWAS and plots results

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
library(stringr)

overlap_cor <- read.csv("overlap_cor.csv")
overlap_1_to_1_cor <- read.csv("overlap_1_to_1_cor.csv")

pqtl_normal <- read.delim("anml_overlap_out.txt")
names(pqtl_normal)[1:3] <- c("uniprot_id","olink_id","somascan_id")
pqtl_normal <- pqtl_normal[,-2]
overlap_pqtl_normal <- merge(overlap_cor,pqtl_normal,by=c("uniprot_id","somascan_id"),all.x=T)
overlap_pqtl_normal[is.na(overlap_pqtl_normal)] <- 0
overlap_pqtl_normal_1_to_1 <- overlap_pqtl_normal[overlap_pqtl_normal$uniprot_id %in% overlap_1_to_1_cor$uniprot_id,]

pqtl_non_normal <- read.delim("nanml_overlap_out.txt")
names(pqtl_non_normal)[1:3] <- c("uniprot_id","olink_id","somascan_id")
pqtl_non_normal <- pqtl_non_normal[,-2]
overlap_pqtl_non_normal <- merge(overlap_cor,pqtl_non_normal,by=c("uniprot_id","somascan_id"),all.x=T)
overlap_pqtl_non_normal[is.na(overlap_pqtl_non_normal)] <- 0
overlap_pqtl_non_normal_1_to_1 <- overlap_pqtl_non_normal[overlap_pqtl_non_normal$uniprot_id %in% overlap_1_to_1_cor$uniprot_id,]

pqtl_count_normal <- data.frame(colSums(overlap_pqtl_normal_1_to_1[,c("olink_cis","olink_trans","soma_cis","soma_trans")]))
names(pqtl_count_normal) <- "Frequency"
pqtl_count_non_normal <- data.frame(colSums(overlap_pqtl_non_normal_1_to_1[,c("soma_cis","soma_trans")]))
names(pqtl_count_non_normal) <- "Frequency"

pqtl_count <- rbind(pqtl_count_normal,pqtl_count_non_normal)

pqtl_count$Platform <- c("OLINK","OLINK","SomaScan-ANML","SomaScan-ANML","SomaScan-non-ANML","SomaScan-non-ANML")
pqtl_count$pQTL <- c("Cis-pQTL","Trans-pQTL","Cis-pQTL","Trans-pQTL","Cis-pQTL","Trans-pQTL")
rownames(pqtl_count) <- c("olink_cis","olink_trans","soma_normal_cis","soma_normal_trans","soma_non_normal_cis","soma_non_normal_trans")

# proteins with pqtls

overlap_pqtl_normal_1_to_1_tf <- overlap_pqtl_normal_1_to_1
overlap_pqtl_normal_1_to_1_tf[,c("olink_cis","olink_trans","soma_cis","soma_trans")] <- T

column_names <- c("olink_cis","olink_trans","soma_cis","soma_trans")

for (i in 1:length(column_names)) {
  overlap_pqtl_normal_1_to_1_tf[overlap_pqtl_normal_1_to_1[,column_names[i]]==0,column_names[i]] <- F
}

overlap_pqtl_non_normal_1_to_1_tf <- overlap_pqtl_non_normal_1_to_1
overlap_pqtl_non_normal_1_to_1_tf[,c("olink_cis","olink_trans","soma_cis","soma_trans")] <- T

for (i in 1:length(column_names)) {
  overlap_pqtl_non_normal_1_to_1_tf[overlap_pqtl_non_normal_1_to_1[,column_names[i]]==0,column_names[i]] <- F
}



pqtl_count_tf_normal <- data.frame(colSums(overlap_pqtl_normal_1_to_1_tf[,c("olink_cis","olink_trans","soma_cis","soma_trans")]))
names(pqtl_count_tf_normal) <- "Frequency"
pqtl_count_tf_non_normal <- data.frame(colSums(overlap_pqtl_non_normal_1_to_1_tf[,c("soma_cis","soma_trans")]))
names(pqtl_count_tf_non_normal) <- "Frequency"

pqtl_count_tf <- rbind(pqtl_count_tf_normal,pqtl_count_tf_non_normal)

pqtl_count_tf$Platform <- c("Olink","Olink","SomaScan-ANML","SomaScan-ANML","SomaScan-non-ANML","SomaScan-non-ANML")
pqtl_count_tf$pQTL <- c("Cis-pQTL","Trans-pQTL","Cis-pQTL","Trans-pQTL","Cis-pQTL","Trans-pQTL")
rownames(pqtl_count_tf) <- c("olink_cis","olink_trans","soma_normal_cis","soma_normal_trans","soma_non_normal_cis","soma_non_normal_trans")

pqtl_count_tf_add <- pqtl_count_tf
pqtl_count_tf_add$Frequency <- pqtl_count$Frequency - pqtl_count_tf$Frequency

pqtl_count_tf$cat<- "original"
pqtl_count_tf_add$cat<- "additional"

pqtl_count_all <- rbind(pqtl_count_tf,pqtl_count_tf_add)

pqtl_count_all$Frequency_all <- NA
pqtl_count_all$Frequency_all[pqtl_count_all$cat=="additional"] <- pqtl_count$Frequency

pqtl_count_all$Frequency_shade <- NA
pqtl_count_all$Frequency_shade[pqtl_count_all$cat=="original"] <- pqtl_count_all$Frequency[pqtl_count_all$cat=="original"] 

# plot stacked bar together

pqtl_count_all <- rbind(pqtl_count_tf,pqtl_count_tf_add)

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

ggsave("",plot_pqtl_hit_all,width=9,height=6)

## figure for main text

pqtl_count_all_figure <- pqtl_count_all[pqtl_count_all$Platform!="SomaScan-ANML",]

pqtl_count_all_figure$Platform <- str_remove(pqtl_count_all_figure$Platform, "-non-ANML")

plot_pqtl_hit_all_figure <- ggplot(pqtl_count_all_figure,aes(y=Frequency, x=fct_rev(Platform), fill=fct_rev(Platform), pattern=cat)) +
  geom_bar_pattern(stat = "identity",
                   position = "stack",
                   colour = 'black',
                   pattern_fill = "black",
                   pattern_spacing = 0.03,
                   pattern_frequency = 5,
                   pattern_angle = 45,
                   pattern_density=0.01) +
  ylim(0,1400) +
  geom_text(aes(label=Frequency_all),position="stack",hjust=-0.1) +
  geom_label(aes(label=Frequency_shade),position=position_stack(vjust = 0.5),label.size = NA) +
  scale_fill_manual(values = rev(hex[c(1,3)])) +
  scale_pattern_manual(values=c("none","stripe")) +
  xlab("") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(legend.position = "none") +
  coord_flip() +
  facet_grid(pQTL~.)

plot_pqtl_hit_all_figure

ggsave("",plot_pqtl_hit_all_figure,width=8,height=5)


library(tidyr)


# load coloc results

overlap_pqtl_normal <- overlap_pqtl_normal[,c("uniprot_id","olink_id","somascan_id","olink_cis","soma_cis","shared_cis","olink_trans","soma_trans","shared_trans")]
names(overlap_pqtl_normal) <- c("uniprot_id","olink_id","somascan_id","olink_cis","soma_normal_cis","coloc_normal_cis","olink_trans","soma_normal_trans","coloc_normal_trans")

overlap_pqtl_non_normal <- overlap_pqtl_non_normal[,c("uniprot_id","olink_id","somascan_id","olink_cis","soma_cis","shared_cis","olink_trans","soma_trans","shared_trans")]
names(overlap_pqtl_non_normal) <- c("uniprot_id","olink_id","somascan_id","olink_cis","soma_non_normal_cis","coloc_non_normal_cis","olink_trans","soma_non_normal_trans","coloc_non_normal_trans")

overlap_coloc <- merge(overlap_pqtl_normal,overlap_pqtl_non_normal,by=c("uniprot_id","olink_id","somascan_id","olink_cis","olink_trans"))

write.csv(overlap_coloc,"overlap_coloc.csv", quote=F, row.names=F)

overlap_coloc_1_to_1 <- overlap_coloc[overlap_coloc$uniprot_id %in% overlap_1_to_1_cor$uniprot_id,]

write.csv(overlap_coloc_1_to_1,"overlap_coloc_1_to_1.csv", quote=F, row.names=F)


# venn diagrams

# cis normal

overlap_coloc_1_to_1$coloc_normal_cis[is.na(overlap_coloc_1_to_1$coloc_normal_cis)] <- 0
pqtl_count_tf
c1 <- pqtl_count_tf["olink_cis",]$Frequency
c2 <- pqtl_count_tf["soma_normal_cis",]$Frequency
c3 <- as.numeric(sum(overlap_coloc_1_to_1$olink_cis>0 & overlap_coloc_1_to_1$soma_normal_cis>0))
c4 <- as.numeric(sum(overlap_coloc_1_to_1$olink_cis>0 & overlap_coloc_1_to_1$soma_normal_cis>0 & overlap_coloc_1_to_1$coloc_normal_cis>0))
c5 <- as.numeric(sum(overlap_coloc_1_to_1$olink_cis>0 & overlap_coloc_1_to_1$soma_normal_cis==0 & overlap_coloc_1_to_1$coloc_normal_cis>0))
c6 <- as.numeric(sum(overlap_coloc_1_to_1$olink_cis==0 & overlap_coloc_1_to_1$soma_normal_cis>0 & overlap_coloc_1_to_1$coloc_normal_cis>0))
c7 <- as.numeric(sum(overlap_coloc_1_to_1$olink_cis==0 & overlap_coloc_1_to_1$soma_normal_cis==0 & overlap_coloc_1_to_1$coloc_normal_cis>0))
overlap_coloc_1_to_1$olink_id[overlap_coloc_1_to_1$olink_cis==0 & overlap_coloc_1_to_1$soma_normal_cis==0 & overlap_coloc_1_to_1$coloc_normal_cis>0]

overlap_coloc_1_to_1$coloc_normal_cis_tf <- F
overlap_coloc_1_to_1$coloc_normal_cis_tf[overlap_coloc_1_to_1$coloc_normal_cis!=0] <- T
table(overlap_coloc_1_to_1$coloc_normal_cis_tf)

venn_diagram_pqtl_olink_soma_normal <- venn.diagram(
      list(olink = 1:c1,soma = (c1-c3+1):(c1-c3+c2)),
      filename = NULL,
      fontfamily = "sans",	cat.fontfamily = "sans",
      category.names = c("Olink" , "SomaScan-ANML"),
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

png("",width = 1600, height = 1600, pointsize=45)
grid.draw(venn_diagram_pqtl_olink_soma_normal)
dev.off()

# cis non-normal

overlap_coloc_1_to_1$coloc_non_normal_cis[is.na(overlap_coloc_1_to_1$coloc_non_normal_cis)] <- 0
pqtl_count_tf
c1 <- pqtl_count_tf["olink_cis",]$Frequency
c2 <- pqtl_count_tf["soma_non_normal_cis",]$Frequency
c3 <- as.numeric(sum(overlap_coloc_1_to_1$olink_cis>0 & overlap_coloc_1_to_1$soma_non_normal_cis>0))
c4 <- as.numeric(sum(overlap_coloc_1_to_1$olink_cis>0 & overlap_coloc_1_to_1$soma_non_normal_cis>0 & overlap_coloc_1_to_1$coloc_non_normal_cis>0))
c5 <- as.numeric(sum(overlap_coloc_1_to_1$olink_cis>0 & overlap_coloc_1_to_1$soma_non_normal_cis==0 & overlap_coloc_1_to_1$coloc_non_normal_cis>0))
c6 <- as.numeric(sum(overlap_coloc_1_to_1$olink_cis==0 & overlap_coloc_1_to_1$soma_non_normal_cis>0 & overlap_coloc_1_to_1$coloc_non_normal_cis>0))
c7 <- as.numeric(sum(overlap_coloc_1_to_1$olink_cis==0 & overlap_coloc_1_to_1$soma_non_normal_cis==0 & overlap_coloc_1_to_1$coloc_non_normal_cis>0))
overlap_coloc_1_to_1$olink_id[overlap_coloc_1_to_1$olink_cis==0 & overlap_coloc_1_to_1$soma_non_normal_cis==0 & overlap_coloc_1_to_1$coloc_non_normal_cis>0]

overlap_coloc_1_to_1$coloc_non_normal_cis_tf <- F
overlap_coloc_1_to_1$coloc_non_normal_cis_tf[overlap_coloc_1_to_1$coloc_non_normal_cis!=0] <- T
table(overlap_coloc_1_to_1$coloc_non_normal_cis_tf)

venn_diagram_pqtl_olink_soma_non_normal <- venn.diagram(
  list(olink = 1:c1,soma = (c1-c3+1):(c1-c3+c2)),
  filename = NULL,
  fontfamily = "sans",	cat.fontfamily = "sans",
  category.names = c("Olink" , "SomaScan-non-ANML"),
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

png("",width = 1600, height = 1600, pointsize=45)
grid.draw(venn_diagram_pqtl_olink_soma_non_normal)
dev.off()

## figure for main text

venn_diagram_pqtl_olink_soma_non_normal <- venn.diagram(
  list(olink = 1:c1,soma = (c1-c3+1):(c1-c3+c2)),
  filename = NULL,
  fontfamily = "sans",	cat.fontfamily = "sans",
  category.names = c("Olink" , "SomaScan"),
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

png("",width = 1600, height = 1600, pointsize=45)
grid.draw(venn_diagram_pqtl_olink_soma_non_normal)
dev.off()

write.csv(overlap_coloc_1_to_1,"overlap_coloc_1_to_1.csv", quote=F, row.names=F)


# histogram

## plot shared only using shade

overlap_coloc_1_to_1 <- merge(overlap_coloc_1_to_1, overlap_1_to_1_cor[,c("uniprot_id","rho_olink_soma_normal","rho_olink_soma_non_normal")],by="uniprot_id")

median_rho_all_normal <- format(round(median(overlap_coloc_1_to_1$rho_olink_soma_normal),2),nsmall=2)
median_rho_coloc_normal <- format(round(median(overlap_coloc_1_to_1$rho_olink_soma_normal[overlap_coloc_1_to_1$coloc_normal_cis_tf==T]),2),nsmall=2)

n_coloc_normal <- sum(overlap_coloc_1_to_1$coloc_normal_cis_tf==T)
n_coloc_non_normal <- sum(overlap_coloc_1_to_1$coloc_non_normal_cis_tf==T)

hist_pqtl_normal <- ggplot(overlap_coloc_1_to_1, aes(x=rho_olink_soma_normal,pattern=coloc_normal_cis_tf)) + 
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
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_coloc_normal, " (coloc cis, n = ", n_coloc_normal, ")"), size = 12/.pt, hjust = 0)

hist_pqtl_normal

median_rho_all_non_normal <- format(round(median(overlap_coloc_1_to_1$rho_olink_soma_non_normal),2),nsmall=2)
median_rho_coloc_non_normal <- format(round(median(overlap_coloc_1_to_1$rho_olink_soma_non_normal[overlap_coloc_1_to_1$coloc_non_normal_cis_tf==T]),2),nsmall=2)

hist_pqtl_non_normal <- ggplot(overlap_coloc_1_to_1, aes(x=rho_olink_soma_non_normal,pattern=coloc_non_normal_cis_tf)) + 
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
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_coloc_non_normal, " (coloc cis, n = ", n_coloc_non_normal, ")"), size = 12/.pt, hjust = 0)

hist_pqtl_non_normal

hist_pqtl <- ggarrange(hist_pqtl_normal,hist_pqtl_non_normal,ncol=2,nrow=1)

ggsave("",hist_pqtl,width=14,height=6)




### trans pqtls

overlap_coloc_1_to_1 <- read.csv("overlap_coloc_1_to_1.csv")

overlap_1_to_1_cor <- read.csv("overlap_1_to_1_cor.csv")

overlap_coloc_1_to_1 <- merge(overlap_coloc_1_to_1, overlap_1_to_1_cor[,c("uniprot_id","rho_olink_soma_normal","rho_olink_soma_non_normal")],by="uniprot_id")

sum(overlap_coloc_1_to_1$olink_trans>0 | overlap_coloc_1_to_1$soma_non_normal_trans>0)

sum(overlap_coloc_1_to_1$coloc_non_normal_trans>0)

median(overlap_coloc_1_to_1$rho_olink_soma_non_normal[overlap_coloc_1_to_1$coloc_non_normal_trans>0])





### read LD results

r2_anml <- read.delim("")
length(unique(r2_anml$u))
r2_anml <- r2_anml[!is.na(r2_anml$cis_trans),]
overlap_coloc_1_to_1 <- read.csv("overlap_coloc_1_to_1.csv")

cis_both <- overlap_coloc_1_to_1$uniprot_id[overlap_coloc_1_to_1$olink_cis>0 & overlap_coloc_1_to_1$soma_normal_cis>0]

r2_anml_cis <- r2_anml[r2_anml$cis_trans=="cis",]
r2_anml_cis_1_to_1 <- r2_anml_cis[r2_anml_cis$u %in% cis_both, ]
length(unique(r2_anml_cis_1_to_1$u))
mean(r2_anml_cis_1_to_1$r2)

overlap_coloc_1_to_1$r2_cis_anml_0.6=F

for (i in 1:nrow(overlap_coloc_1_to_1)) {
  df <- r2_anml_cis_1_to_1[r2_anml_cis_1_to_1$u==overlap_coloc_1_to_1$uniprot_id[i],]
  if (max(df$r2)>0.6) {
    overlap_coloc_1_to_1$r2_cis_anml_0.6[i] <- TRUE
  }
}
sum(overlap_coloc_1_to_1$r2_cis_anml_0.6)

r2_nanml <- read.delim("")
length(unique(r2_nanml$u))
r2_nanml <- r2_nanml[!is.na(r2_nanml$cis_trans),]

cis_both <- overlap_coloc_1_to_1$uniprot_id[overlap_coloc_1_to_1$olink_cis>0 & overlap_coloc_1_to_1$soma_non_normal_cis>0]

r2_nanml_cis <- r2_nanml[r2_nanml$cis_trans=="cis",]
r2_nanml_cis_1_to_1 <- r2_nanml_cis[r2_nanml_cis$u %in% cis_both, ]
length(unique(r2_nanml_cis_1_to_1$u))
mean(r2_nanml_cis_1_to_1$r2)

overlap_coloc_1_to_1$r2_cis_nanml_0.6=F

for (i in 1:nrow(overlap_coloc_1_to_1)) {
  df <- r2_nanml_cis_1_to_1[r2_nanml_cis_1_to_1$u==overlap_coloc_1_to_1$uniprot_id[i],]
  if (max(df$r2)>0.6) {
    overlap_coloc_1_to_1$r2_cis_nanml_0.6[i] <- TRUE
  }
}
sum(overlap_coloc_1_to_1$r2_cis_nanml_0.6)
