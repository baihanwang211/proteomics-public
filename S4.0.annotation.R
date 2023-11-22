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
library(moments)
library(datawizard)
library(gtsummary)

## load proteomics data

olink <- read.csv("olink.csv")
somascan_normalised <- read.csv("somascan_normalised.csv")
somascan_non_normalised <- read.csv("somascan_non_normalised.csv")
somascan_normalised_log <- read.csv("somascan_normalised_log.csv")
somascan_non_normalised_log <- read.csv("somascan_non_normalised_log.csv")

## load overlapping proteins

overlap_cor <- read.csv("overlap_cor.csv")

## create dummy variables

table(overlap_cor$panel)

overlap_cor$panel_Cardiometabolic <- F
overlap_cor$panel_Cardiometabolic[overlap_cor$panel=="Cardiometabolic"] <- T

overlap_cor$panel_Cardiometabolic_II <- F
overlap_cor$panel_Cardiometabolic_II[overlap_cor$panel=="Cardiometabolic_II"] <- T

overlap_cor$panel_Inflammation <- F
overlap_cor$panel_Inflammation[overlap_cor$panel=="Inflammation"] <- T

overlap_cor$panel_Inflammation_II <- F
overlap_cor$panel_Inflammation_II[overlap_cor$panel=="Inflammation_II"] <- T

overlap_cor$panel_Neurology <- F
overlap_cor$panel_Neurology[overlap_cor$panel=="Neurology"] <- T

overlap_cor$panel_Neurology_II <- F
overlap_cor$panel_Neurology_II[overlap_cor$panel=="Neurology_II"] <- T

overlap_cor$panel_Oncology <- F
overlap_cor$panel_Oncology[overlap_cor$panel=="Oncology"] <- T

overlap_cor$panel_Oncology_II <- F
overlap_cor$panel_Oncology_II[overlap_cor$panel=="Oncology_II"] <- T

table(overlap_cor$dilution)

overlap_cor$dilution_5e05 <- F
overlap_cor$dilution_5e05[overlap_cor$dilution==5e-05] <- T

overlap_cor$dilution_0.005 <- F
overlap_cor$dilution_0.005[overlap_cor$dilution==0.005] <- T

overlap_cor$dilution_0.2 <- F
overlap_cor$dilution_0.2[overlap_cor$dilution==0.2] <- T

## calculate mean npx olink

mean_npx_olink <- data.frame(colMeans(olink[, 2:ncol(olink)]))
names(mean_npx_olink) <- "mean_npx_olink"
mean_npx_olink$olink_id <- row.names(mean_npx_olink)
overlap_cor <- merge(overlap_cor,mean_npx_olink,by="olink_id")

## calculate mean log soma

mean_log_soma_normal <- data.frame(colMeans(somascan_normalised_log[, 2:ncol(somascan_normalised_log)]))
names(mean_log_soma_normal) <- "mean_log_soma_normal"
mean_log_soma_normal$somascan_id <- row.names(mean_log_soma_normal)
overlap_cor <- merge(overlap_cor,mean_log_soma_normal,by="somascan_id")
mean_log_soma_non_normal <- data.frame(colMeans(somascan_non_normalised_log[, 2:ncol(somascan_non_normalised_log)]))
names(mean_log_soma_non_normal) <- "mean_log_soma_non_normal"
mean_log_soma_non_normal$somascan_id <- row.names(mean_log_soma_non_normal)
overlap_cor <- merge(overlap_cor,mean_log_soma_non_normal,by="somascan_id")

## calculate skewness 

skewness_olink <- data.frame(skewness(olink[,2:ncol(olink)]))
names(skewness_olink) <- c("skewness_olink")
skewness_olink$olink_id <- row.names(skewness_olink)
overlap_cor <- merge(overlap_cor,skewness_olink,by="olink_id")

skewness_soma_normal <- data.frame(skewness(somascan_normalised_log[,2:ncol(somascan_normalised_log)]))
names(skewness_soma_normal) <- c("skewness_soma_normal")
skewness_soma_normal$somascan_id <- row.names(skewness_soma_normal)
overlap_cor <- merge(overlap_cor,skewness_soma_normal,by="somascan_id")

skewness_soma_non_normal <- data.frame(skewness(somascan_non_normalised_log[,2:ncol(somascan_non_normalised_log)]))
names(skewness_soma_non_normal) <- c("skewness_soma_non_normal")
skewness_soma_non_normal$somascan_id <- row.names(skewness_soma_non_normal)
overlap_cor <- merge(overlap_cor,skewness_soma_non_normal,by="somascan_id")

## calculate kurtosis 

kurtosis_olink <- data.frame(kurtosis(olink[,2:ncol(olink)]))
names(kurtosis_olink) <- c("kurtosis_olink")
kurtosis_olink$olink_id <- row.names(kurtosis_olink)
overlap_cor <- merge(overlap_cor,kurtosis_olink,by="olink_id")

kurtosis_soma_normal <- data.frame(kurtosis(somascan_normalised_log[,2:ncol(somascan_normalised_log)]))
names(kurtosis_soma_normal) <- c("kurtosis_soma_normal")
kurtosis_soma_normal$somascan_id <- row.names(kurtosis_soma_normal)
overlap_cor <- merge(overlap_cor,kurtosis_soma_normal,by="somascan_id")

kurtosis_soma_non_normal <- data.frame(kurtosis(somascan_non_normalised_log[,2:ncol(somascan_non_normalised_log)]))
names(kurtosis_soma_non_normal) <- c("kurtosis_soma_non_normal")
kurtosis_soma_non_normal$somascan_id <- row.names(kurtosis_soma_non_normal)
overlap_cor <- merge(overlap_cor,kurtosis_soma_non_normal,by="somascan_id")

## calculate percentage below LOD and QC/assay warnings for olink

olink_all <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/data_baseline_olink_explore.csv")
olink_unique <- olink_all[!duplicated(olink_all[c("csid","assay")]),]
olink_unique$below_lod <- F
olink_unique$below_lod[olink_unique$npx < olink_unique$lod] <- T
olink_unique <- olink_unique[,c("csid","assay","below_lod", "qc_warning", "assay_warning")]

# lod

olink_lod <- pivot_wider(olink_unique, id_cols = csid, names_from = assay, values_from = below_lod)
olink_lod_count <- data.frame(colSums(olink_lod[2:ncol(olink_lod)])/nrow(olink_lod))
names(olink_lod_count) <- "below_lod_olink"
olink_lod_count$olink_id <- row.names(olink_lod_count)
olink_lod_count$olink_id <- gsub("-",".",olink_lod_count$olink_id)
overlap_cor <- merge(overlap_cor, olink_lod_count, by = "olink_id")
hist(overlap_cor$below_lod_olink)

# qc warning

olink_qc <- pivot_wider(olink_unique, id_cols = csid, names_from = assay, values_from = qc_warning)
olink_qc_count <- data.frame(colSums(olink_qc[2:ncol(olink_qc)])/nrow(olink_qc))
names(olink_qc_count) <- "qc_warning_olink"
olink_qc_count$olink_id <- row.names(olink_qc_count)
olink_qc_count$olink_id <- gsub("-",".",olink_qc_count$olink_id)
overlap_cor <- merge(overlap_cor, olink_qc_count, by = "olink_id")
hist(overlap_cor$qc_warning_olink)

# assay warning

olink_assay <- pivot_wider(olink_unique, id_cols = csid, names_from = assay, values_from = assay_warning)
olink_assay_count <- data.frame(colSums(olink_assay[2:ncol(olink_assay)])/nrow(olink_assay))
names(olink_assay_count) <- "assay_warning_olink"
olink_assay_count$olink_id <- row.names(olink_assay_count)
olink_assay_count$olink_id <- gsub("-",".",olink_assay_count$olink_id)
overlap_cor <- merge(overlap_cor, olink_assay_count, by = "olink_id")
hist(overlap_cor$assay_warning_olink)

# QC check for soma

somascan_qc <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00244-V1/somalogic_meta.csv")[,c("aptname","col_check")]
names(somascan_qc) <- c("somascan_id","qc_check_soma")
overlap_cor <- merge(overlap_cor, somascan_qc, by = "somascan_id")

# calculate percentage outlier

for (i in 1:nrow(overlap_cor)) {
  column <- pull(olink,overlap_cor$olink_id[i])
  overlap_cor$outlier_olink[i] <- sum((column < mean(column)-4*sd(column) | column > mean(column)+4*sd(column)))/nrow(olink)
  column <- pull(somascan_normalised_log,overlap_cor$somascan_id[i])
  overlap_cor$outlier_soma_normal[i] <- sum((column < mean(column)-4*sd(column) | column > mean(column)+4*sd(column)))/nrow(somascan_normalised_log)
  column <- pull(somascan_non_normalised_log,overlap_cor$somascan_id[i])
  overlap_cor$outlier_soma_non_normal[i] <- sum((column < mean(column)-4*sd(column) | column > mean(column)+4*sd(column)))/nrow(somascan_non_normalised_log)
}

##### other variables for soma

# lod

soma_extra <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/proteomics/somascan_extra_info.csv")
soma_extra$SeqId <- gsub("-",".",soma_extra$SeqId)
soma_extra$SeqId <- paste0("seq.",soma_extra$SeqId)
soma_lod <-soma_extra[,c("SeqId","LoDB..RFU.")]
names(soma_lod) <- c("somascan_id","lod")

overlap_cor$below_lod_soma_normal <- NA
overlap_cor$below_lod_soma_non_normal <- NA

for (i in 1:nrow(overlap_cor)) {
  overlap_cor$below_lod_soma_normal[i] <- sum(somascan_normalised[[overlap_cor$somascan_id[i]]]<soma_lod$lod[soma_lod$somascan_id==overlap_cor$somascan_id[i]])/nrow(somascan_normalised)
  overlap_cor$below_lod_soma_non_normal[i] <- sum(somascan_non_normalised[[overlap_cor$somascan_id[i]]]<soma_lod$lod[soma_lod$somascan_id==overlap_cor$somascan_id[i]])/nrow(somascan_non_normalised)
}

# saveRDS(overlap_cor,"overlap_cor_annot_mid.RDS")
# overlap_cor <- readRDS("overlap_cor_annot_mid.RDS")

# kdm

soma_kd <-soma_extra[,c("SeqId","Apparent.Kd..M.")]
names(soma_lod) <- c("somascan_id","kd_soma")
overlap_cor <- merge(overlap_cor, soma_lod, by = "somascan_id")

##### other variables for olink

olink_extra <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/proteomics/olink_extra_info.csv")
olink_extra <- olink_extra[!duplicated(olink_extra$UniProt),]

# dilution

olink_dilution <- olink_extra[,c("UniProt","Dilution.factor")]
names(olink_dilution) <- c("uniprot_id","dilution_olink")
olink_dilution$dilution_olink <- gsub(" ", "", olink_dilution$dilution_olink)
table(olink_dilution$dilution_olink)

olink_dilution$dilution_1_1 <- F
olink_dilution$dilution_1_1[olink_dilution$dilution_olink=="1:1"] <- T
table(olink_dilution$dilution_1_1)

olink_dilution$dilution_1_10 <- F
olink_dilution$dilution_1_10[olink_dilution$dilution_olink=="1:10"] <- T
table(olink_dilution$dilution_1_10)

olink_dilution$dilution_1_100 <- F
olink_dilution$dilution_1_100[olink_dilution$dilution_olink=="1:100"] <- T
table(olink_dilution$dilution_1_100)

olink_dilution$dilution_1_1000 <- F
olink_dilution$dilution_1_1000[olink_dilution$dilution_olink=="1:1000"] <- T
table(olink_dilution$dilution_1_1000)

olink_dilution$dilution_1_100000 <- F
olink_dilution$dilution_1_100000[olink_dilution$dilution_olink=="1:100000"] <- T
table(olink_dilution$dilution_1_100000)

olink_dilution <- olink_dilution[,-which(names(olink_dilution) %in% "dilution_olink")]

overlap_cor <- merge(overlap_cor,olink_dilution,by="uniprot_id")

saveRDS(overlap_cor,"overlap_cor_annot_mid.RDS")
# overlap_cor <- readRDS("overlap_cor_annot_mid.RDS")

## calculate CV

# cv_olink <- data.frame(colMeans(olink[, 2:ncol(olink)]))
# names(cv_olink) <- "mean"
# cv_olink$olink_id <- row.names(cv_olink)
# cv_olink$sd <- sapply(olink[, 2:ncol(olink)], sd)
# cv_olink$cv_olink <- cv_olink$sd/cv_olink$mean
# overlap_cor <- merge(overlap_cor, cv_olink[,c("olink_id","cv_olink")], by="olink_id")
# 
# cv_soma_normal <- data.frame(colMeans(somascan_normalised_log[, 2:ncol(somascan_normalised_log)]))
# names(cv_soma_normal) <- "mean"
# cv_soma_normal$somascan_id <- row.names(cv_soma_normal)
# cv_soma_normal$sd <- sapply(somascan_normalised_log[, 2:ncol(somascan_normalised_log)], sd)
# cv_soma_normal$cv_soma_normal <- cv_soma_normal$sd/cv_soma_normal$mean
# overlap_cor <- merge(overlap_cor, cv_soma_normal[,c("somascan_id","cv_soma_normal")], by="somascan_id")
# 
# cv_soma_non_normal <- data.frame(colMeans(somascan_non_normalised_log[, 2:ncol(somascan_non_normalised_log)]))
# names(cv_soma_non_normal) <- "mean"
# cv_soma_non_normal$somascan_id <- row.names(cv_soma_non_normal)
# cv_soma_non_normal$sd <- sapply(somascan_non_normalised_log[, 2:ncol(somascan_non_normalised_log)], sd)
# cv_soma_non_normal$cv_soma_non_normal <- cv_soma_non_normal$sd/cv_soma_non_normal$mean
# overlap_cor <- merge(overlap_cor, cv_soma_non_normal[,c("somascan_id","cv_soma_non_normal")], by="somascan_id")

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

test <- queryUniProt(query = unique(overlap_cor$uniprot_id)[1],fields=field)

label <- names(test)

#### testing here ####

annot <- as.data.frame(matrix(ncol = length(field)+1, nrow = length(unique(overlap_cor$uniprot_id))))

names(annot) <- c("uniprot_id",label)

annot$uniprot_id <- unique(overlap_cor$uniprot_id)

for (i in 1:nrow(annot)) {
  annot[i,c(2:ncol(annot))] <- queryUniProt(query = overlap_cor$uniprot_id[i],fields=field)
}

# saveRDS(annot, file = "annot.RDS") 

# annot <- readRDS("annot.RDS")

overlap_annot <- merge(overlap_cor, annot, by="uniprot_id")

saveRDS(overlap_annot,"overlap_annot.RDS")

# overlap_annot <- readRDS("overlap_annot.RDS")

# only keep 1_to_1

protein_dup <- unique(overlap_annot$uniprot_id[duplicated(overlap_annot$uniprot_id)]) # get proteins targeted by multiple aptamers

aptamer_dup <- unique(overlap_annot$somascan_id[duplicated(overlap_annot$somascan_id)]) # get aptamers targeted by multiple proteins

overlap_1_to_1_annot <- overlap_annot[-c(which(overlap_annot$uniprot_id %in% protein_dup),which(overlap_annot$somascan_id %in% aptamer_dup)), ]

# get top ten annotations from GO

# GO biological process

annot_go_p <- separate_rows(overlap_1_to_1_annot[,c("uniprot_id","Gene.Ontology..biological.process.")],uniprot_id,Gene.Ontology..biological.process.,sep="; ")

annot_go_p_tb <- data.frame(count=sort(table(annot_go_p$Gene.Ontology..biological.process.), decreasing=TRUE))

annot_go_p_wide <- pivot_wider(annot_go_p,
                               id_cols = "uniprot_id",
                               names_from = "Gene.Ontology..biological.process.", 
                               values_from = 'Gene.Ontology..biological.process.', 
                               values_fill = F,
                               values_fn = function(x) T)

annot_go_p_10 <- annot_go_p_wide[, c("uniprot_id",as.vector(annot_go_p_tb$count.Var1[1:10]))]

# GO cellular component

annot_go_c <- separate_rows(overlap_1_to_1_annot[,c("uniprot_id","Gene.Ontology..cellular.component.")],uniprot_id,Gene.Ontology..cellular.component.,sep="; ")

annot_go_c_tb <- data.frame(count=sort(table(annot_go_c$Gene.Ontology..cellular.component.), decreasing=TRUE))

annot_go_c_wide <- pivot_wider(annot_go_c,
                               id_cols = "uniprot_id",
                               names_from = "Gene.Ontology..cellular.component.", 
                               values_from = 'Gene.Ontology..cellular.component.', 
                               values_fill = F,
                               values_fn = function(x) T)

annot_go_c_10 <- annot_go_c_wide[, c("uniprot_id",as.vector(annot_go_c_tb$count.Var1[1:10]))]

# GO molecular function

annot_go_f <- separate_rows(overlap_1_to_1_annot[,c("uniprot_id","Gene.Ontology..molecular.function.")],uniprot_id,Gene.Ontology..molecular.function.,sep="; ")

annot_go_f_tb <- data.frame(count=sort(table(annot_go_f$Gene.Ontology..molecular.function.), decreasing=TRUE))

annot_go_f_wide <- pivot_wider(annot_go_f,
                               id_cols = "uniprot_id",
                               names_from = "Gene.Ontology..molecular.function.", 
                               values_from = 'Gene.Ontology..molecular.function.', 
                               values_fill = F,
                               values_fn = function(x) T)

annot_go_f_10 <- annot_go_f_wide[, c("uniprot_id",as.vector(annot_go_f_tb$count.Var1[1:10]))]

# remove original columns and merge with overlap

overlap_1_to_1_annot <- overlap_1_to_1_annot[ , -which(names(overlap_1_to_1_annot) %in% c("Gene.Ontology..biological.process.","Gene.Ontology..cellular.component.","Gene.Ontology..molecular.function."))]

df_list <- list(overlap_1_to_1_annot, annot_go_p_10, annot_go_c_10, annot_go_f_10)

overlap_1_to_1_annot <- Reduce(function(x, y) merge(x, y, by = "uniprot_id"), df_list)

# make some factors binary

names(overlap_1_to_1_annot)

overlap_1_to_1_annot_tf <- overlap_1_to_1_annot[,53:77]

overlap_1_to_1_annot_tf[!is.na(overlap_1_to_1_annot_tf)] <- T

overlap_1_to_1_annot_tf[is.na(overlap_1_to_1_annot_tf)] <- F

overlap_1_to_1_annot[,53:77] <- overlap_1_to_1_annot_tf

# get number of isoforms

overlap_1_to_1_annot$n_isoform <- as.numeric(str_extract(overlap_1_to_1_annot$Alternative.products..isoforms., "(?<=Named isoforms=)\\d+"))

overlap_1_to_1_annot$n_isoform[is.na(overlap_1_to_1_annot$n_isoform)] <- 1

overlap_1_to_1_annot <- overlap_1_to_1_annot[ , -which(names(overlap_1_to_1_annot) %in% "Alternative.products..isoforms.")]

# merge with shared cis-pQTL
overlap_coloc <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/overlap_coloc.csv")
overlap_coloc <- overlap_coloc[,c("uniprot_id","coloc_normal_cis_tf","coloc_non_normal_cis_tf","coloc_normal_trans_tf","coloc_non_normal_trans_tf")]
overlap_1_to_1_annot <- merge(overlap_1_to_1_annot,overlap_coloc,by="uniprot_id")

# remove columns with too few observations

names(overlap_1_to_1_annot)
sapply(overlap_1_to_1_annot[52:106],table)

# remove annotations with two few proteins (n<15)

overlap_1_to_1_annot <- overlap_1_to_1_annot[, -which(names(overlap_1_to_1_annot) %in% c("Intramembrane","DNA.binding","non_standard_residue"))]

names(overlap_1_to_1_annot)

names(overlap_1_to_1_annot) <- tolower(names(overlap_1_to_1_annot))

names(overlap_1_to_1_annot) <- gsub(" |:", "_", names(overlap_1_to_1_annot))

names(overlap_1_to_1_annot) <- gsub("\\[|\\]", "", names(overlap_1_to_1_annot))

names(overlap_1_to_1_annot) <- gsub("\\.", "_", names(overlap_1_to_1_annot))

names(overlap_1_to_1_annot)

# define levels of batch and dilution

overlap_1_to_1_annot$batch <- factor(overlap_1_to_1_annot$batch, levels = c("1","2"))

saveRDS(overlap_1_to_1_annot,"overlap_1_to_1_annot.RDS")

# overlap_1_to_1_annot <- readRDS("overlap_1_to_1_annot.RDS")

############# prepare file for running on cluster

####### select only relevant variables

names(overlap_1_to_1_annot)

overlap_1_to_1_annot_normal <- overlap_1_to_1_annot[,c(9,5,15:25,26,27,29,30,32,33,35:38,39,40,42,44:49,50:105)]
names(overlap_1_to_1_annot_normal)
saveRDS(overlap_1_to_1_annot_normal,"overlap_1_to_1_annot_normal.RDS")

overlap_1_to_1_annot_non_normal <- overlap_1_to_1_annot[,c(10,5,15:25,26,28,29,31,32,34,35:38,39,41,43,44:49,50:105)]
names(overlap_1_to_1_annot_non_normal)
saveRDS(overlap_1_to_1_annot_non_normal,"overlap_1_to_1_annot_non_normal.RDS")






## check
overlap_1_to_1_annot_normal <- readRDS("overlap_1_to_1_annot_normal.RDS")
plot(overlap_1_to_1_annot_normal$kurtosis_olink,overlap_1_to_1_annot_normal$rho_olink_soma_normal)
plot(overlap_1_to_1_annot_normal$skewness_olink,overlap_1_to_1_annot_normal$rho_olink_soma_normal)
plot(overlap_1_to_1_annot_normal$kurtosis_soma_normal,overlap_1_to_1_annot_normal$rho_olink_soma_normal)
plot(overlap_1_to_1_annot_normal$skewness_soma_normal,overlap_1_to_1_annot_normal$rho_olink_soma_normal)





## count number of proteins per go category

names(overlap_1_to_1_annot)
count_go <- data.frame(sapply(overlap_1_to_1_annot[75:104],sum))
names(count_go)[1] <- "freq"
# count_go$term <- row.names(count_go)
count_go$percent <- count_go$freq/1694
count_go



overlap_1_to_1_annot[,50:105] %>%
  tbl_summary(statistic = list(
                all_continuous() ~ "{mean} {sd}",
                all_categorical() ~ "{n} {p}"),
              digits = everything() ~ 1)