rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(tidyverse)
library(datawizard)
library(RNOmni)

## load olink data

# firstly load linkage file for uniprot id

olink_protein <- read_excel("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00211-V1/olink_uniprot.xlsx")[c(1,3,4)]

names(olink_protein) <- c("uniprot_id","olink_id","panel")

olink_protein <- olink_protein[!duplicated(olink_protein$olink_id),]

olink_protein$olink_id <- tolower(olink_protein$olink_id)

olink_protein$olink_id <- paste0("ol_", olink_protein$olink_id)

olink_protein$olink_id <- gsub("-","_",olink_protein$olink_id)

olink_protein$batch <- 1

olink_protein[grep("II",olink_protein$panel),]$batch <-2

olink_protein$panel <- gsub("_II","",olink_protein$panel)

olink_protein_batch1 <- olink_protein$olink_id[olink_protein$batch==1]

olink_protein_batch2 <- olink_protein$olink_id[olink_protein$batch==2]

## load olink ckb data

olink_cardiometabolic <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00211-V1/data_baseline_olink_cardiometabolic.csv")

olink_inflammation <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00211-V1/data_baseline_olink_inflammation.csv")

olink_neurology <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00211-V1/data_baseline_olink_neurology.csv")

olink_oncology <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00211-V1/data_baseline_olink_oncology.csv")

olink_all <- olink_cardiometabolic %>% inner_join(olink_inflammation, by='csid', suffix=c("","_dup")) %>% inner_join(olink_neurology, by='csid', suffix=c("","_dup")) %>% inner_join(olink_oncology, by='csid', suffix=c("","_dup"))

# find difference between proteins in the assay and ckb data

setdiff(olink_protein$olink_id,names(olink_all)[-1])

setdiff(names(olink_all)[-1],olink_protein$olink_id)

# only keep those proteins that are also included in assay list

olink_unique <- olink_all[,-grep("dup",names(olink_all))]

olink <- olink_all[, c(1, which(names(olink_all) %in% olink_protein$olink_id))]

# load removed individuals

olink_batch1_rm <- c("G1585511", "G1007262", "G1690469")

olink_batch2_rm <- c("G1097667", "G2076351", "G1464504", "G1677005", "G1585511", "G1824411", "G1277351", "G1007262", "G2037800")

# load gwas linkage

gwas <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00211-V1/data_gwas_genetics.csv")

olink_batch1_rm <- gwas$csid[which(gwas$ccvid %in% olink_batch1_rm)]

olink_batch2_rm <- gwas$csid[which(gwas$ccvid %in% olink_batch2_rm)]

## load somascan data

# load linkage file first

somascan_protein <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00211-V1/somalogic_meta.csv")[,c("aptname","uniprot","dilution2")]

names(somascan_protein) <- c("somascan_id","uniprot_id","dilution")

# split uniprot column into 2 because some  aptamers have more than 1 uniprot id

somascan_protein <- separate(data = somascan_protein, col = uniprot_id, into = c("uniprot_id_1", "uniprot_id_2","uniprot_id_3"), sep = "\\|")

# load normalised somascan data and only keep overlapping proteins

somascan_normalised <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00211-V1/data_baseline_somalogic.csv")

# somascan_normalised <- somascan_normalised[which(somascan_normalised$variable %in% somascan_protein$somascan_id),]

# change to wide format

somascan_normalised <- pivot_wider(somascan_normalised, names_from = variable, values_from = value)

somascan_normalised <- somascan_normalised[somascan_normalised$qc==0,-c(2,3)]

## load non-normalised somascan data and only keep overlapping proteins

somascan_non_normalised <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00211-V1/data_baseline_somalogic_qc.csv")

# somascan_non_normalised <- somascan_non_normalised[which(somascan_non_normalised$variable %in% somascan_protein$somascan_id),]

# change to wide format

somascan_non_normalised <- pivot_wider(somascan_non_normalised, names_from = variable, values_from = value)

somascan_non_normalised <- somascan_non_normalised[somascan_non_normalised$qc==0,-c(2,3)]

## log

somascan_normalised_log <- somascan_normalised

somascan_normalised_log[2:length(somascan_normalised_log)] <- lapply(somascan_normalised_log[2:length(somascan_normalised_log)],log)


somascan_non_normalised_log <- somascan_non_normalised

somascan_non_normalised_log[2:length(somascan_non_normalised_log)] <- lapply(somascan_non_normalised_log[2:length(somascan_non_normalised_log)],log)


### pca for log

# normalised

comp_normal_log = svd(somascan_normalised_log[,2:ncol(somascan_normalised_log)])

means=c(1:10)

for (i in c(1:10)) {
  means[i]=mean(comp_normal_log$u[,i])
}

stdevs=c(1:10)

for (i in c(1:10)) {
  stdevs[i]=sd(comp_normal_log$u[,i])
}

Zs=comp_normal_log$u[,c(1:10)]

for (i in c(1:length(Zs[,1]))) {
  Zs[i,]=(Zs[i,]-means)/stdevs
}

exclude=c(1:nrow(somascan_normalised_log))

for (i in c(1:nrow(somascan_normalised_log))) {
  exclude[i]=max(abs(Zs[i,c(1:10)]))>5 # first 10 PCs based on scree plot
}

list_ids_rm_normal_log=somascan_normalised_log$csid[exclude==T]
  
write.table(list_ids_rm_normal_log,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/list_ids_rm_normal_log.txt", quote=F, row.names=F, col.names=F)

# non-normalised

comp_non_normal_log = svd(somascan_non_normalised_log[,2:ncol(somascan_non_normalised_log)])

means=c(1:10)

for (i in c(1:10)) {
  means[i]=mean(comp_non_normal_log$u[,i])
}

stdevs=c(1:10)

for (i in c(1:10)) {
  stdevs[i]=sd(comp_non_normal_log$u[,i])
}

Zs=comp_non_normal_log$u[,c(1:10)]

for (i in c(1:length(Zs[,1]))) {
  Zs[i,]=(Zs[i,]-means)/stdevs
}

exclude=c(1:nrow(somascan_non_normalised_log))

for (i in c(1:nrow(somascan_non_normalised_log))) {
  exclude[i]=max(abs(Zs[i,c(1:10)]))>5 # first 10 PCs based on scree plot
}

list_ids_rm_non_normal_log=somascan_non_normalised_log$csid[exclude==T]

write.table(list_ids_rm_non_normal_log,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/list_ids_rm_non_normal_log.txt", quote=F, row.names=F, col.names=F)



## 1/x

somascan_normalised_1x <- somascan_normalised

somascan_normalised_1x[2:length(somascan_normalised_1x)] <- 1/somascan_normalised_1x[2:length(somascan_normalised_1x)]


somascan_non_normalised_1x <- somascan_non_normalised

somascan_non_normalised_1x[2:length(somascan_non_normalised_1x)] <- 1/somascan_non_normalised_1x[2:length(somascan_non_normalised_1x)]

### pca for 1/x

# normalised

comp_normal_1x = svd(somascan_normalised_1x[,2:ncol(somascan_normalised_1x)])

means=c(1:10)

for (i in c(1:10)) {
  means[i]=mean(comp_normal_1x$u[,i])
}

stdevs=c(1:10)

for (i in c(1:10)) {
  stdevs[i]=sd(comp_normal_1x$u[,i])
}

Zs=comp_normal_1x$u[,c(1:10)]

for (i in c(1:length(Zs[,1]))) {
  Zs[i,]=(Zs[i,]-means)/stdevs
}

exclude=c(1:nrow(somascan_normalised_1x))

for (i in c(1:nrow(somascan_normalised_1x))) {
  exclude[i]=max(abs(Zs[i,c(1:10)]))>5 # first 10 PCs based on scree plot
}

list_ids_rm_normal_1x=somascan_normalised_1x$csid[exclude==T]

write.table(list_ids_rm_normal_1x,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/list_ids_rm_normal_1x.txt", quote=F, row.names=F, col.names=F)

# non-normalised

comp_non_normal_1x = svd(somascan_non_normalised_1x[,2:ncol(somascan_non_normalised_1x)])

means=c(1:10)

for (i in c(1:10)) {
  means[i]=mean(comp_non_normal_1x$u[,i])
}

stdevs=c(1:10)

for (i in c(1:10)) {
  stdevs[i]=sd(comp_non_normal_1x$u[,i])
}

Zs=comp_non_normal_1x$u[,c(1:10)]

for (i in c(1:length(Zs[,1]))) {
  Zs[i,]=(Zs[i,]-means)/stdevs
}

exclude=c(1:nrow(somascan_non_normalised_1x))

for (i in c(1:nrow(somascan_non_normalised_1x))) {
  exclude[i]=max(abs(Zs[i,c(1:10)]))>5 # first 10 PCs based on scree plot
}

list_ids_rm_non_normal_1x=somascan_non_normalised_1x$csid[exclude==T]

write.table(list_ids_rm_non_normal_1x,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/list_ids_rm_non_normal_1x.txt", quote=F, row.names=F, col.names=F)


# remove individuals and save dataset

olink_clean <- olink
olink_clean[which(olink_clean$csid %in% olink_batch1_rm),which(names(olink_clean) %in% olink_protein_batch1)] <- NA
olink_clean[which(olink_clean$csid %in% olink_batch2_rm),which(names(olink_clean) %in% olink_protein_batch2)] <- NA
colSums(is.na(olink_clean))
write.csv(olink_clean,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/olink_clean.csv", quote=F, row.names=F)


somascan_normalised_log_clean <- somascan_normalised_log
somascan_normalised_log_clean[somascan_normalised_log_clean$csid %in% list_ids_rm_normal_log, 2:ncol(somascan_normalised_log_clean)] <- NA
colSums(is.na(somascan_normalised_log_clean))
write.csv(somascan_normalised_log_clean,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/somascan_normalised_log_clean.csv", quote=F, row.names=F)


somascan_non_normalised_log_clean <- somascan_non_normalised_log
somascan_non_normalised_log_clean[somascan_non_normalised_log_clean$csid %in% list_ids_rm_non_normal_log, 2:ncol(somascan_non_normalised_log_clean)] <- NA
colSums(is.na(somascan_non_normalised_log_clean))
write.csv(somascan_non_normalised_log_clean,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/somascan_non_normalised_log_clean.csv", quote=F, row.names=F)


somascan_normalised_1x_clean <- somascan_normalised_1x
somascan_normalised_1x_clean[somascan_normalised_1x_clean$csid %in% list_ids_rm_normal_1x, 2:ncol(somascan_normalised_1x_clean)] <- NA
colSums(is.na(somascan_normalised_1x_clean))
write.csv(somascan_normalised_1x_clean,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/somascan_normalised_1x_clean.csv", quote=F, row.names=F)


somascan_non_normalised_1x_clean <- somascan_non_normalised_1x
somascan_non_normalised_1x_clean[somascan_non_normalised_1x_clean$csid %in% list_ids_rm_non_normal_1x, 2:ncol(somascan_non_normalised_1x_clean)] <- NA
colSums(is.na(somascan_non_normalised_1x_clean))
write.csv(somascan_non_normalised_1x_clean,"K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/somascan_non_normalised_1x_clean.csv", quote=F, row.names=F)


##### reload

olink_clean <- read.csv("olink_clean")
somascan_normalised_log_clean <- read.csv("somascan_normalised_log_clean.csv")
somascan_non_normalised_log_clean <- read.csv("somascan_non_normalised_log_clean.csv")
somascan_normalised_1x_clean <- read.csv("somascan_normalised_1x_clean.csv")
somascan_non_normalised_1x_clean <- read.csv("somascan_non_normalised_1x_clean.csv")

##### adjustment and RINT

# load demographic variables

baseline <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00211-V1/data_baseline_questionnaires.csv")

baseline_var <- baseline[,c("csid","age_at_study_date_x100","is_female", "region_code")]

baseline_var$age_at_study_date_x100 <- baseline_var$age_at_study_date_x100/100

names(baseline_var)[2] <- "age"

baseline_var$age2 <- baseline_var$age * baseline_var$age

baseline_var$age <- as.numeric(baseline_var$age)
baseline_var$age2 <- as.numeric(baseline_var$age2)
baseline_var$is_female <- as.factor(baseline_var$is_female)
baseline_var$region_code <- as.factor(baseline_var$region_code)

# load pcs

gwas_var <- gwas[c(1,grep("national",names(gwas)))]

covar <- merge(baseline_var,gwas_var,by="csid")

# merge 

olink_clean_adj <- merge(covar,olink_clean,by="csid")
somascan_normalised_log_clean_adj <- merge(covar,somascan_normalised_log_clean,by="csid")
somascan_non_normalised_log_clean_adj <- merge(covar,somascan_non_normalised_log_clean,by="csid")
somascan_normalised_1x_clean_adj <- merge(covar,somascan_normalised_1x_clean,by="csid")
somascan_non_normalised_1x_clean_adj <- merge(covar,somascan_non_normalised_1x_clean,by="csid")

#adjust

olink_clean_adj <- data_adjust(olink_clean_adj,effect=names(covar)[-1],select=17:ncol(olink_clean_adj))

somascan_normalised_log_clean_adj <- data_adjust(somascan_normalised_log_clean_adj,effect=names(covar)[-1],select=17:ncol(somascan_normalised_log_clean_adj))
somascan_non_normalised_log_clean_adj <- data_adjust(somascan_non_normalised_log_clean_adj,effect=names(covar)[-1],select=17:ncol(somascan_non_normalised_log_clean_adj))
somascan_normalised_1x_clean_adj <- data_adjust(somascan_normalised_1x_clean_adj,effect=names(covar)[-1],select=17:ncol(somascan_normalised_1x_clean_adj))
somascan_non_normalised_1x_clean_adj <- data_adjust(somascan_non_normalised_1x_clean_adj,effect=names(covar)[-1],select=17:ncol(somascan_non_normalised_1x_clean_adj))

###### rint

# olink

olink_clean_adj_rint_batch1 <- olink_clean_adj[,c(1,which(names(olink_clean_adj) %in% olink_protein_batch1))]
olink_clean_adj_rint_batch1 <- olink_clean_adj_rint_batch1[complete.cases(olink_clean_adj_rint_batch1),]

olink_clean_adj_rint_batch2 <- olink_clean_adj[,c(1,which(names(olink_clean_adj) %in% olink_protein_batch2))]
olink_clean_adj_rint_batch2 <- olink_clean_adj_rint_batch2[complete.cases(olink_clean_adj_rint_batch2),]

for (i in c(2:ncol(olink_clean_adj_rint_batch1))) {
  olink_clean_adj_rint_batch1[,i] <- RankNorm(olink_clean_adj_rint_batch1[,i])
}

for (i in c(2:ncol(olink_clean_adj_rint_batch2))) {
  olink_clean_adj_rint_batch2[,i] <- RankNorm(olink_clean_adj_rint_batch2[,i])
}

olink_clean_adj_rint <- as.data.frame(olink_clean_adj$csid)
names(olink_clean_adj_rint) <- "csid"

olink_clean_adj_rint <- merge(olink_clean_adj_rint,olink_clean_adj_rint_batch1,by="csid",all.x = T)
olink_clean_adj_rint <- merge(olink_clean_adj_rint,olink_clean_adj_rint_batch2,by="csid",all.x = T)

##### somascan

somascan_normalised_log_clean_adj_rint <- somascan_normalised_log_clean_adj[complete.cases(somascan_normalised_log_clean_adj),]
somascan_normalised_log_clean_adj_rint <- somascan_normalised_log_clean_adj_rint[,-c(2:16)]

somascan_non_normalised_log_clean_adj_rint <- somascan_non_normalised_log_clean_adj[complete.cases(somascan_non_normalised_log_clean_adj),]
somascan_non_normalised_log_clean_adj_rint <- somascan_non_normalised_log_clean_adj_rint[,-c(2:16)]

somascan_normalised_1x_clean_adj_rint <- somascan_normalised_1x_clean_adj[complete.cases(somascan_normalised_1x_clean_adj),]
somascan_normalised_1x_clean_adj_rint <- somascan_normalised_1x_clean_adj_rint[,-c(2:16)]

somascan_non_normalised_1x_clean_adj_rint <- somascan_non_normalised_1x_clean_adj[complete.cases(somascan_non_normalised_1x_clean_adj),]
somascan_non_normalised_1x_clean_adj_rint <- somascan_non_normalised_1x_clean_adj_rint[,-c(2:16)]

names(somascan_normalised_log_clean_adj_rint)

for (i in c(2:ncol(somascan_normalised_log_clean_adj_rint))) {
  somascan_normalised_log_clean_adj_rint[,i] <- RankNorm(somascan_normalised_log_clean_adj_rint[,i])
  somascan_non_normalised_log_clean_adj_rint[,i] <- RankNorm(somascan_non_normalised_log_clean_adj_rint[,i])
  somascan_normalised_1x_clean_adj_rint[,i] <- RankNorm(somascan_normalised_1x_clean_adj_rint[,i])
  somascan_non_normalised_1x_clean_adj_rint[,i] <- RankNorm(somascan_non_normalised_1x_clean_adj_rint[,i])
}

somascan_csid <- as.data.frame(somascan_normalised_log_clean_adj$csid)
names(somascan_csid) <- "csid"

somascan_normalised_log_clean_adj_rint <- merge(somascan_csid,somascan_normalised_log_clean_adj_rint,by="csid")
somascan_non_normalised_log_clean_adj_rint <- merge(somascan_csid,somascan_non_normalised_log_clean_adj_rint,by="csid")
somascan_normalised_1x_clean_adj_rint <- merge(somascan_csid,somascan_normalised_1x_clean_adj_rint,by="csid")
somascan_non_normalised_1x_clean_adj_rint <- merge(somascan_csid,somascan_non_normalised_1x_clean_adj_rint,by="csid")

# save

write.csv(olink_clean_adj_rint,"olink_clean_adj_rint.csv", quote=F, row.names=F)

write.csv(somascan_normalised,"somascan_normalised.csv", quote=F, row.names=F)
write.csv(somascan_normalised_log_clean_adj_rint,"somascan_normalised_log_clean_adj_rint.csv", quote=F, row.names=F)
write.csv(somascan_non_normalised_log_clean_adj_rint,"somascan_non_normalised_log_clean_adj_rint.csv", quote=F, row.names=F)

write.csv(somascan_non_normalised,"somascan_non_normalised.csv", quote=F, row.names=F)
write.csv(somascan_normalised_1x_clean_adj_rint,"somascan_normalised_1x_clean_adj_rint.csv", quote=F, row.names=F)
write.csv(somascan_non_normalised_1x_clean_adj_rint,"somascan_non_normalised_1x_clean_adj_rint.csv", quote=F, row.names=F)

