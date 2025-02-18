##### this script loads pQTL and colocalisation results from GWAS and plots results

rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data_test")

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
library(tidyr)
library(readxl)

##### load files

# load list of overlapping proteins and their correlation coefficients
overlap_cor <- read.csv("overlap_cor.csv")
overlap_1_to_1_cor <- read.csv("overlap_1_to_1_cor.csv")

# load pQTL counts for overlapping proteins
overlap_pqtl <- read_excel("")

overlap_pqtl_1_to_1 <- overlap_pqtl[overlap_pqtl$uniprot_id %in% overlap_1_to_1_cor$uniprot_id,]

pqtl_count <- data.frame(colSums(overlap_pqtl_1_to_1[,c("olink_cis","olink_trans","soma_non_ANML_cis","soma_non_ANML_trans")]))
names(pqtl_count) <- "Frequency"

pqtl_count$Platform <- c("OLINK","OLINK","SomaScan-non-ANML","SomaScan-non-ANML")
pqtl_count$pQTL <- c("Cis-pQTL","Trans-pQTL","Cis-pQTL","Trans-pQTL")


# proteins with pqtls
overlap_pqtl_1_to_1_tf <- overlap_pqtl_1_to_1
overlap_pqtl_1_to_1_tf[,c("olink_cis","olink_trans","soma_non_ANML_cis","soma_non_ANML_trans")] <- T

column_names <- c("olink_cis","olink_trans","soma_non_ANML_cis","soma_non_ANML_trans")

for (i in 1:length(column_names)) {
  overlap_pqtl_1_to_1_tf[overlap_pqtl_1_to_1[,column_names[i]]==0,column_names[i]] <- F
}

for (i in 1:length(column_names)) {
  overlap_pqtl_1_to_1_tf[overlap_pqtl_1_to_1[,column_names[i]]==0,column_names[i]] <- F
}

pqtl_count_tf <- data.frame(colSums(overlap_pqtl_1_to_1_tf[,c("olink_cis","olink_trans","soma_non_ANML_cis","soma_non_ANML_trans")]))
names(pqtl_count_tf) <- "Frequency"

pqtl_count_tf$Platform <- c("Olink","Olink","SomaScan-non-ANML","SomaScan-non-ANML")
pqtl_count_tf$pQTL <- c("Cis-pQTL","Trans-pQTL","Cis-pQTL","Trans-pQTL")
rownames(pqtl_count_tf) <- c("olink_cis","olink_trans","soma_non_ANML_cis","soma_non_ANML_trans")

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

plot_pqtl_hit <- ggplot(pqtl_count_all,aes(y=Frequency, x=fct_rev(Platform), fill=fct_rev(Platform), pattern=cat)) +
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
  theme(text = element_text(size = 12)) +
  theme(legend.position = "none") +
  coord_flip() +
  facet_grid(pQTL~.)

plot_pqtl_hit


# load coloc results

overlap_pqtl_1_to_1$coloc_non_ANML_cis[is.na(overlap_pqtl_1_to_1$coloc_non_ANML_cis)] <- 0

pqtl_count_tf
c1 <- pqtl_count_tf["olink_cis",]$Frequency
c2 <- pqtl_count_tf["soma_non_ANML_cis",]$Frequency
c3 <- as.numeric(sum(overlap_pqtl_1_to_1$olink_cis>0 & overlap_pqtl_1_to_1$soma_non_ANML_cis>0))
c4 <- as.numeric(sum(overlap_pqtl_1_to_1$olink_cis>0 & overlap_pqtl_1_to_1$soma_non_ANML_cis>0 & overlap_pqtl_1_to_1$coloc_non_ANML_cis>0))
c5 <- as.numeric(sum(overlap_pqtl_1_to_1$olink_cis>0 & overlap_pqtl_1_to_1$soma_non_ANML_cis==0 & overlap_pqtl_1_to_1$coloc_non_ANML_cis>0))
c6 <- as.numeric(sum(overlap_pqtl_1_to_1$olink_cis==0 & overlap_pqtl_1_to_1$soma_non_ANML_cis>0 & overlap_pqtl_1_to_1$coloc_non_ANML_cis>0))
c7 <- as.numeric(sum(overlap_pqtl_1_to_1$olink_cis==0 & overlap_pqtl_1_to_1$soma_non_ANML_cis==0 & overlap_pqtl_1_to_1$coloc_non_ANML_cis>0))
overlap_pqtl_1_to_1$olink_id[overlap_pqtl_1_to_1$olink_cis==0 & overlap_pqtl_1_to_1$soma_non_ANML_cis==0 & overlap_pqtl_1_to_1$coloc_non_ANML_cis>0]

overlap_pqtl_1_to_1$coloc_non_ANML_cis_tf <- F
overlap_pqtl_1_to_1$coloc_non_ANML_cis_tf[overlap_pqtl_1_to_1$coloc_non_ANML_cis!=0] <- T
table(overlap_pqtl_1_to_1$coloc_non_ANML_cis_tf)

# venn diagram

venn_diagram_pqtl_olink_soma <- venn.diagram(
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
grid.draw(venn_diagram_pqtl_olink_soma)

venn_diagram_pqtl_olink_soma[[5]]$label <- paste(venn_diagram_pqtl_olink_soma[[5]]$label, 
                                                            "\n", c5, " (", format(round((c5/(c1-c3))*100,digits=1),nsmall=1), "%)\ncolocalise", sep="")  
venn_diagram_pqtl_olink_soma[[6]]$label <- paste(venn_diagram_pqtl_olink_soma[[6]]$label,
                                                            "\n", c6, " (", format(round((c6/(c2-c3))*100,digits=1),nsmall=1), "%)\ncolocalise", sep="")  
venn_diagram_pqtl_olink_soma[[7]]$label <- paste(venn_diagram_pqtl_olink_soma[[7]]$label,
                                                            "\n", c4, " (", format(round(c4/c3*100,digits=1),nmsall=1), "%) colocalise", sep="")  

grid.newpage()
grid.draw(venn_diagram_pqtl_olink_soma)



# histogram 

## plot shared only using shade

overlap_pqtl_1_to_1 <- merge(overlap_pqtl_1_to_1, overlap_1_to_1_cor[,c("uniprot_id","rho_olink_soma_non_ANML")],by="uniprot_id")

median_rho_all <- format(round(median(overlap_pqtl_1_to_1$rho_olink_soma_non_ANML),2),nsmall=2)
median_rho_coloc <- format(round(median(overlap_pqtl_1_to_1$rho_olink_soma_non_ANML[overlap_pqtl_1_to_1$coloc_non_normal_cis_tf==T]),2),nsmall=2)

hist_pqtl <- ggplot(overlap_pqtl_1_to_1, aes(x=rho_olink_soma,pattern=coloc_non_ANML_cis_tf)) + 
  geom_histogram_pattern(binwidth = 0.05,
                         boundary=0,
                         pattern_fill = "black", 
                         fill    = 'white',
                         colour  = 'black',
                         pattern_spacing = 0.015,
                         pattern_frequency = 5, 
                         pattern_angle = 45,
                         pattern_density=0.01) +
  scale_x_continuous(breaks = seq(-0.3, 1, 0.1), lim = c(-0.3, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  # ggtitle("OLINK vs SomaScan", subtitle ="") +
  xlab("Spearman's rho") + 
  ylab("Frequency") +
  labs(fill="Shared") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_pattern_manual(values=c("none","stripe")) +
  theme(legend.position = "none") +
  annotate("text", x = -0.3, y = 300, label = paste("Median rho =", median_rho_all, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.3, y = 285, label = paste0("Median rho = ", median_rho_coloc, " (coloc cis, n = ", n_coloc, ")"), size = 12/.pt, hjust = 0)

hist_pqtl
