rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(UniProt.ws)
library(ckbplotr)
library(hausekeep)
library(stringr)
library(Boruta)
library(ggplot2)
library(reshape)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)
library(dataPreparation)
library(stringr)

olink_soma_normal <- read.csv("olink_soma_normal.csv")

olink_soma_non_normal <- read.csv("olink_soma_non_normal.csv")

overlap <- read.csv("overlap.csv")

# calculate percentage of outliers

# for (i in 1:nrow(overlap)) {
#   overlap$outlier_olink[i] <- sum(is.na(outliers_mad(olink_soma_normal[,overlap$olink_id[i]], threshold = 5)))/nrow(olink_soma_normal)
#   overlap$outlier_soma_normal[i] <- sum(is.na(outliers_mad(olink_soma_normal[,overlap$somascan_id[i]], threshold = 5)))/nrow(olink_soma_normal)
#   overlap$outlier_soma_non_normal[i] <- sum(is.na(outliers_mad(olink_soma_non_normal[,overlap$somascan_id[i]], threshold = 5)))/nrow(olink_soma_non_normal)
# }

for (i in 1:nrow(overlap)) {
  column <- pull(olink_soma_normal,overlap$olink_id[i])
  overlap$outlier_olink[i] <- sum((column < mean(column)-4*sd(column) | column > mean(column)+4*sd(column)))/nrow(olink_soma_normal)
  column <- pull(olink_soma_normal,overlap$somascan_id[i])
  overlap$outlier_soma_normal[i] <- sum((column < mean(column)-4*sd(column) | column > mean(column)+4*sd(column)))/nrow(olink_soma_normal)
  column <- pull(olink_soma_non_normal,overlap$somascan_id[i])
  overlap$outlier_soma_non_normal[i] <- sum((column < mean(column)-4*sd(column) | column > mean(column)+4*sd(column)))/nrow(olink_soma_normal)
}

# get uniprot annotations
# first define fields we need

field <- c("mass","length","cc_alternative_products", # Sequences
           "ft_init_met","ft_signal","ft_transit","ft_propep","ft_peptide", # Molecule processing
           "ft_topo_dom","ft_transmem","ft_intramem","ft_repeat","ft_zn_fing","ft_dna_bind","ft_coiled","ft_motif","ft_compbias", # Regions
           "ft_act_site", "ft_binding", # Sites
           "ft_non_std", "ft_mod_res", "ft_lipid", "ft_carbohyd", "ft_disulfid", "ft_crosslnk", # Amino acid modifications
           "ft_helix", "ft_turn", "ft_strand", #Secondary structure
           "go_p","go_c","go_f" # GO
           )

#### testing here ####

test <- queryUniProt(query = unique(overlap$uniprot_id)[1],fields=field)

label <- names(test)

#### testing here ####

annot <- as.data.frame(matrix(ncol = length(field)+1, nrow = length(unique(overlap$uniprot_id))))

names(annot) <- c("uniprot_id",label)

annot$uniprot_id <- unique(overlap$uniprot_id)

for (i in 1:nrow(annot)) {
  annot[i,c(2:ncol(annot))] <- queryUniProt(query = overlap$uniprot_id[i],fields=field)
}

saveRDS(annot, file = "annot.RDS") 

# annot <- readRDS("annot.RDS")

names(annot)

colSums(!is.na(annot))

# remove annotations with two few proteins

annot <- annot[, -which(names(annot) %in% c("Intramembrane","Non.standard.residue"))]

overlap_annot <- merge(overlap, annot, by="uniprot_id")

# define levels of batch and dilution

names(overlap_annot)

overlap_annot$batch <- factor(overlap_annot$batch, levels = c("1","2"))

overlap_annot$dilution <- factor(overlap_annot$dilution, levels = c("5e-05","0.005","0.2"), labels = c("0.005%","0.5%","20%"))

# get top ten annotations from GO

# GO biological process

annot_go_p <- separate_rows(overlap_annot[,c("uniprot_id","Gene.Ontology..biological.process.")],uniprot_id,Gene.Ontology..biological.process.,sep="; ")

annot_go_p_tb <- data.frame(count=sort(table(annot_go_p$Gene.Ontology..biological.process.), decreasing=TRUE))

annot_go_p_wide <- pivot_wider(annot_go_p,
                               id_cols = "uniprot_id",
                               names_from = "Gene.Ontology..biological.process.", 
                               values_from = 'Gene.Ontology..biological.process.', 
                               values_fill = F,
                               values_fn = function(x) T)

annot_go_p_10 <- annot_go_p_wide[, c("uniprot_id",as.vector(annot_go_p_tb$count.Var1[1:10]))]

# GO cellular component

annot_go_c <- separate_rows(overlap_annot[,c("uniprot_id","Gene.Ontology..cellular.component.")],uniprot_id,Gene.Ontology..cellular.component.,sep="; ")

annot_go_c_tb <- data.frame(count=sort(table(annot_go_c$Gene.Ontology..cellular.component.), decreasing=TRUE))

annot_go_c_wide <- pivot_wider(annot_go_c,
                               id_cols = "uniprot_id",
                               names_from = "Gene.Ontology..cellular.component.", 
                               values_from = 'Gene.Ontology..cellular.component.', 
                               values_fill = F,
                               values_fn = function(x) T)

annot_go_c_10 <- annot_go_c_wide[, c("uniprot_id",as.vector(annot_go_c_tb$count.Var1[1:10]))]

# GO molecular function

annot_go_f <- separate_rows(overlap_annot[,c("uniprot_id","Gene.Ontology..molecular.function.")],uniprot_id,Gene.Ontology..molecular.function.,sep="; ")

annot_go_f_tb <- data.frame(count=sort(table(annot_go_f$Gene.Ontology..molecular.function.), decreasing=TRUE))

annot_go_f_wide <- pivot_wider(annot_go_f,
                               id_cols = "uniprot_id",
                               names_from = "Gene.Ontology..molecular.function.", 
                               values_from = 'Gene.Ontology..molecular.function.', 
                               values_fill = F,
                               values_fn = function(x) T)

annot_go_f_10 <- annot_go_f_wide[, c("uniprot_id",as.vector(annot_go_f_tb$count.Var1[1:10]))]

# remove original columns and merge with overlap

overlap_annot <- overlap_annot[ , -which(names(overlap_annot) %in% c("Gene.Ontology..biological.process.","Gene.Ontology..cellular.component.","Gene.Ontology..molecular.function."))]

df_list <- list(overlap_annot, annot_go_p_10, annot_go_c_10, annot_go_f_10)

overlap_annot <- Reduce(function(x, y) merge(x, y, by = "uniprot_id"), df_list)

# make some factors binary

names(overlap_annot)

overlap_annot_tf <- overlap_annot[,20:42]

overlap_annot_tf[!is.na(overlap_annot_tf)] <- T

overlap_annot_tf[is.na(overlap_annot_tf)] <- F

overlap_annot[,20:42] <- overlap_annot_tf

# get number of isoforms

overlap_annot$n_isoform <- as.numeric(str_extract(overlap_annot$Alternative.products..isoforms., "(?<=Named isoforms=)\\d+"))

overlap_annot$n_isoform[is.na(overlap_annot$n_isoform)] <- 1

overlap_annot <- overlap_annot[ , -which(names(overlap_annot) %in% "Alternative.products..isoforms.")]

names(overlap_annot)

saveRDS(overlap_annot,"overlap_annot.RDS")

####### run Random-forest

# overlap_annot <- readRDS("overlap_annot.RDS")

overlap_annot_normal <- overlap_annot[,c(9,3,4,8,14,15,17:length(overlap_annot))]

overlap_annot_non_normal <- overlap_annot[,c(10,3,4,8,14,16,17:length(overlap_annot))]

## run machine learning

# normalised

set.seed(47)

boruta_normal <- Boruta(rho_olink_somascan_normal~., data = overlap_annot_normal, doTrace = 2, maxRuns = 1000)

saveRDS(boruta_normal,"boruta_normal.RDS")

# boruta_normal <- readRDS("boruta_normal.RDS")

print(boruta_normal)

boruta_normal_df <- attStats(boruta_normal)

print(boruta_normal_df)

# reshape data

boruta_normal_melt <- melt(as.data.frame(boruta_normal$ImpHistory))

boruta_normal_melt <- boruta_normal_melt[is.finite(boruta_normal_melt$value),]

# remove shadow variables

boruta_normal_melt <- boruta_normal_melt[-grep("shadow",boruta_normal_melt$variable),]

names(boruta_normal_melt) <- c("variable","importance")

# color by decision

boruta_normal_decision <- data.frame(variable=row.names(boruta_normal_df),decision=as.vector(boruta_normal_df$decision))

boruta_normal_melt <- merge(boruta_normal_melt,boruta_normal_decision,by="variable")

boruta_normal_melt$decision <- factor(boruta_normal_melt$decision,levels=c("Confirmed","Tentative","Rejected"))

# change display of names

levels(boruta_normal_melt$variable) <- str_to_sentence(levels(boruta_normal_melt$variable))

levels(boruta_normal_melt$variable)[1:11] <- c("Olink: panel","Olink: batch", "SomaScan: dilution", "Olink: mean NPX",
                                         "Olink: skewness", "SomaScan: skewness", "Olink: kurtosis", "SomaScan: kurtosis",
                                         "Olink: % below LOD", "Olink: % outliers", "SomaScan: % outliers")


levels(boruta_normal_melt$variable) <- gsub("\\.", " ", levels(boruta_normal_melt$variable))

levels(boruta_normal_melt$variable) <- gsub("Dna", "DNA", levels(boruta_normal_melt$variable))
levels(boruta_normal_melt$variable) <- gsub("dna", "DNA", levels(boruta_normal_melt$variable))

levels(boruta_normal_melt$variable) <- gsub("Rna", "RNA", levels(boruta_normal_melt$variable))
levels(boruta_normal_melt$variable) <- gsub("rna", "RNA", levels(boruta_normal_melt$variable))

levels(boruta_normal_melt$variable) <- gsub("Atp", "ATP", levels(boruta_normal_melt$variable))
levels(boruta_normal_melt$variable) <- gsub("atp", "ATP", levels(boruta_normal_melt$variable))

levels(boruta_normal_melt$variable) <- gsub("go:", "GO:", levels(boruta_normal_melt$variable))

levels(boruta_normal_melt$variable) <- gsub("ii", "II", levels(boruta_normal_melt$variable))

# plot

boruta_normal_plot <- ggplot(boruta_normal_melt, aes(x=reorder(variable, importance, FUN = median), y=importance)) + 
  geom_boxplot(aes(fill=decision)) + 
  coord_flip() +
  xlab("Variable") +
  ylab("Importance") +
  labs(fill="Decision") +
  theme_bw() +
  theme(text = element_text(size = 17))

boruta_normal_plot

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/boruta_normal_plot.png",boruta_normal_plot,width=19.2,height=10.8)

# non-normalised

set.seed(47)

boruta_non_normal <- Boruta(rho_olink_somascan_non_normal~., data = overlap_annot_non_normal, doTrace = 2, maxRuns = 1000)

saveRDS(boruta_non_normal,"boruta_non_normal.RDS")

# boruta_non_normal <- readRDS("boruta_non_normal.RDS")

print(boruta_non_normal)

boruta_non_normal_df <- attStats(boruta_non_normal)

print(boruta_non_normal_df)

# reshape data

boruta_non_normal_melt <- melt(as.data.frame(boruta_non_normal$ImpHistory))

boruta_non_normal_melt <- boruta_non_normal_melt[is.finite(boruta_non_normal_melt$value),]

# remove shadow variables

boruta_non_normal_melt <- boruta_non_normal_melt[-grep("shadow",boruta_non_normal_melt$variable),]

names(boruta_non_normal_melt) <- c("variable","importance")

# color by decision

boruta_non_normal_decision <- data.frame(variable=row.names(boruta_non_normal_df),decision=as.vector(boruta_non_normal_df$decision))

boruta_non_normal_melt <- merge(boruta_non_normal_melt,boruta_non_normal_decision,by="variable")

boruta_non_normal_melt$decision <- factor(boruta_non_normal_melt$decision,levels=c("Confirmed","Tentative","Rejected"))

# change display of names

levels(boruta_normal_melt$variable) <- str_to_sentence(levels(boruta_non_normal_melt$variable))

levels(boruta_non_normal_melt$variable)[1:11] <- c("Olink: panel","Olink: batch", "SomaScan: dilution", "Olink: mean NPX",
                                               "Olink: skewness", "SomaScan: skewness", "Olink: kurtosis", "SomaScan: kurtosis",
                                               "Olink: % below LOD", "Olink: % outliers", "SomaScan: % outliers")


levels(boruta_non_normal_melt$variable) <- gsub("\\.", " ", levels(boruta_non_normal_melt$variable))

levels(boruta_non_normal_melt$variable) <- gsub("Dna", "DNA", levels(boruta_non_normal_melt$variable))
levels(boruta_non_normal_melt$variable) <- gsub("dna", "DNA", levels(boruta_non_normal_melt$variable))

levels(boruta_non_normal_melt$variable) <- gsub("Rna", "RNA", levels(boruta_non_normal_melt$variable))
levels(boruta_non_normal_melt$variable) <- gsub("rna", "RNA", levels(boruta_non_normal_melt$variable))

levels(boruta_non_normal_melt$variable) <- gsub("Atp", "ATP", levels(boruta_non_normal_melt$variable))
levels(boruta_non_normal_melt$variable) <- gsub("atp", "ATP", levels(boruta_non_normal_melt$variable))

levels(boruta_non_normal_melt$variable) <- gsub("go:", "GO:", levels(boruta_non_normal_melt$variable))

levels(boruta_non_normal_melt$variable) <- gsub("ii", "II", levels(boruta_non_normal_melt$variable))

# plot

plot(boruta_non_normal)

boruta_non_normal_plot <- ggplot(boruta_non_normal_melt, aes(x=reorder(variable, importance, FUN = median), y=importance)) + 
  geom_boxplot(aes(fill=decision)) + 
  coord_flip() +
  xlab("Variable") +
  ylab("Importance") +
  labs(fill="Decision") +
  theme_bw() +
  theme(text = element_text(size = 17))

boruta_non_normal_plot

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/boruta_non_normal_plot.png",boruta_non_normal_plot,width=19.2,height=10.8)

# boruta_plot <- ggarrange(boruta_normal_plot,boruta_non_normal_plot,labels=c("A","B"),ncol=2,nrow=1)

# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/boruta_plot.png",boruta_plot,width=19.2*1.5,height=10.8*1.5)

# # run regression models to confirm the direction of associations

# normalised

names(overlap_annot_normal) <- gsub(" |:", "_", names(overlap_annot_normal))

names(overlap_annot_normal) <- gsub("\\[|\\]", "", names(overlap_annot_normal))

rownames(boruta_normal_df) <- gsub(" |:", "_", rownames(boruta_normal_df))

rownames(boruta_normal_df) <- gsub("\\[|\\]", "", rownames(boruta_normal_df))

confirmed <- rownames(boruta_normal_df)[boruta_normal_df$decision=="Confirmed"]

boruta_normal_regress <- data.frame(variable=confirmed)

boruta_normal_regress$coeff <- NA

for (i in 1:(length(confirmed)-2)) {
  lm_results <- lm(paste("rho_olink_somascan_normal ~", confirmed[i+2]), data = overlap_annot_normal)
  boruta_normal_regress$coeff[i+2] <- lm_results$coefficients[1]
}

# non_normal

names(overlap_annot_non_normal) <- gsub(" |:", "_", names(overlap_annot_non_normal))

names(overlap_annot_non_normal) <- gsub("\\[|\\]", "", names(overlap_annot_non_normal))

rownames(boruta_non_normal_df) <- gsub(" |:", "_", rownames(boruta_non_normal_df))

rownames(boruta_non_normal_df) <- gsub("\\[|\\]", "", rownames(boruta_non_normal_df))

confirmed <- rownames(boruta_non_normal_df)[boruta_non_normal_df$decision=="Confirmed"]

boruta_non_normal_regress <- data.frame(variable=confirmed)

boruta_non_normal_regress$coeff <- NA

for (i in 1:(length(confirmed)-2)) {
  lm_results <- lm(paste("rho_olink_somascan_non_normal ~", confirmed[i+2]), data = overlap_annot_non_normal)
  boruta_non_normal_regress$coeff[i+2] <- lm_results$coefficients[1]
}
