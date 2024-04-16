### this plot categorises proteins based on their consistencies across Olink and SomaScan

rm(list = ls())

library(VennDiagram)
library(ggvenn)
library(RColorBrewer)

setwd("")

overlap <- read.csv("overlap_1_to_1_cor.csv")

overlap_coloc <- read.csv("overlap_coloc_1_to_1.csv")

coloc_normal <- overlap_coloc[,c("somascan_id","coloc_normal_cis","coloc_normal_trans")]
overlap_coloc <- merge(overlap,coloc_normal,by="somascan_id",all.x=T)

overlap_coloc$coloc <- F
overlap_coloc$coloc[overlap_coloc$coloc_normal_cis!=0] <- T
table(overlap_coloc$coloc)

overlap_coloc$cor_0.4 <- F
overlap_coloc$cor_0.4[overlap_coloc$rho_olink_soma_normal>0.4] <- T
table(overlap_coloc$cor_0.4)

bmi <- read.csv("bmi.csv")
bmi$bmi_shared <- F
bmi$bmi_shared[bmi$olink_p_fdr_bmi_calc<0.05 & bmi$soma_normal_p_fdr_bmi_calc<0.05 &
                                sign(bmi$olink_es_bmi_calc)==sign(bmi$soma_normal_es_bmi_calc)] <- T
overlap_coloc_bmi <- merge(overlap_coloc,bmi,by=c("uniprot_id","olink_id","somascan_id"))
table(overlap_coloc_bmi$bmi_shared)

overlap_1_to_1_ihd <- read.csv("")
overlap_1_to_1_ihd$ihd_shared <- F
overlap_1_to_1_ihd$ihd_shared[overlap_1_to_1_ihd$olink_p_fdr<0.05 & overlap_1_to_1_ihd$soma_normal_p_fdr<0.05 &
                                               sign(overlap_1_to_1_ihd$olink_coef)==sign(overlap_1_to_1_ihd$soma_normal_coef)] <- T
overlap_coloc_bmi_ihd <- merge(overlap_coloc_bmi,overlap_1_to_1_ihd,by=c("uniprot_id"))
table(overlap_coloc_bmi_ihd$ihd_shared)

names(overlap_coloc_bmi_ihd)[1:3] <- c("uniprot_id","olink_id","somascan_id")

overlap_coloc_bmi_ihd <- overlap_coloc_bmi_ihd[,c("uniprot_id","olink_id","somascan_id","cor_0.4","coloc","bmi_shared","ihd_shared")]

# display.brewer.all()

venn <- ggplot(overlap_coloc_bmi_ihd, aes(A=cor_0.4, B=coloc, C=bmi_shared, D=ihd_shared)) +
  geom_venn(set_names = c("Rho > 0.4 \n(n = 648)", "Colocalisation \n(n = 404)", "Shared associations \nfor BMI (n = 510)", "Shared associations \nfor IHD (n = 78)"),
            fill_color = brewer.pal(4,"Set1")) + 
  xlim(-2.5,3) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("",venn,width=8,height=6)

write.csv(overlap_coloc_bmi_ihd, "summary_normal.csv", row.names = F, quote = F)


# non-normal

coloc_non_normal <- overlap_coloc[,c("somascan_id","coloc_non_normal_cis","coloc_non_normal_trans")]
overlap_coloc <- merge(overlap,coloc_non_normal,by="somascan_id",all.x=T)

overlap_coloc$coloc <- F
overlap_coloc$coloc[overlap_coloc$coloc_non_normal_cis!=0] <- T
table(overlap_coloc$coloc)

overlap_coloc$cor_0.4 <- F
overlap_coloc$cor_0.4[overlap_coloc$rho_olink_soma_non_normal>0.4] <- T
table(overlap_coloc$cor_0.4)

bmi <- read.csv("bmi.csv")
bmi$bmi_shared <- F
bmi$bmi_shared[bmi$olink_p_fdr_bmi_calc<0.05 & bmi$soma_non_normal_p_fdr_bmi_calc<0.05 &
                 sign(bmi$olink_es_bmi_calc)==sign(bmi$soma_non_normal_es_bmi_calc)] <- T
overlap_coloc_bmi <- merge(overlap_coloc,bmi,by=c("uniprot_id","olink_id","somascan_id"))
table(overlap_coloc_bmi$bmi_shared)

overlap_1_to_1_ihd <- read.csv("")
overlap_1_to_1_ihd$ihd_shared <- F
overlap_1_to_1_ihd$ihd_shared[overlap_1_to_1_ihd$olink_p_fdr<0.05 & overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05 &
                                sign(overlap_1_to_1_ihd$olink_coef)==sign(overlap_1_to_1_ihd$soma_non_normal_coef)] <- T
overlap_coloc_bmi_ihd <- merge(overlap_coloc_bmi,overlap_1_to_1_ihd,by=c("uniprot_id"))
table(overlap_coloc_bmi_ihd$ihd_shared)

names(overlap_coloc_bmi_ihd)[1:3] <- c("uniprot_id","olink_id","somascan_id")

overlap_coloc_bmi_ihd <- overlap_coloc_bmi_ihd[,c("uniprot_id","olink_id","somascan_id","cor_0.4","coloc","bmi_shared","ihd_shared")]

# display.brewer.all()

venn <- ggplot(overlap_coloc_bmi_ihd, aes(A=cor_0.4, B=coloc, C=bmi_shared, D=ihd_shared)) +
  geom_venn(set_names = c("Rho > 0.4 \n(n = 694)", "Colocalisation \n(n = 399)", "Shared associations \nfor BMI (n = 859)", "Shared associations \nfor IHD (n = 78)"),
            fill_color = brewer.pal(4,"Set1")) + 
  xlim(-2.5,3) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("",venn,width=8,height=6)

write.csv(overlap_coloc_bmi_ihd, "summary_non_normal.csv", row.names = F, quote = F)
