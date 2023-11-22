rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(tidyverse)
library(datawizard)
library(RNOmni)

# load somascan data

somascan_normalised_log <- read.csv("somascan_normalised_log.csv")

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

gwas <- read.csv("K:/kadoorie/Staff_Folders/BaihanW/data/DAR-2023-00211-V1/data_gwas_genetics.csv")

gwas_var <- gwas[c(1,grep("national",names(gwas)))]

covar <- merge(baseline_var,gwas_var,by="csid")

# merge 

somascan_normalised_log_adj <- merge(covar,somascan_normalised_log,by="csid")

# adjust

somascan_normalised_log_adj <- data_adjust(somascan_normalised_log_adj,effect=names(covar)[-1],select=17:ncol(somascan_normalised_log_adj))

# adjsut method 2

somascan_normalised_log_adj2 <- merge(covar,somascan_normalised_log,by="csid")

for(i in 17:ncol(somascan_normalised_log_adj2)){
  somascan_normalised_log_adj2[,i] = lm(paste(names(somascan_normalised_log_adj2)[i], "~", paste(names(covar)[-1],collapse=" + ")), somascan_normalised_log_adj2)$resid
  print(i)
  }

##################################### 13/06/23

## RINT

somascan_normalised_log_adj_rint <- somascan_normalised_log_adj[,-c(2:16)]

names(somascan_normalised_log_adj_rint)

for (i in c(2:ncol(somascan_normalised_log_adj_rint))) {
  somascan_normalised_log_adj_rint[,i] <- RankNorm(somascan_normalised_log_adj_rint[,i],)
}


somascan_normalised_log_adj2_rint <- somascan_normalised_log_adj2[,-c(2:16)]

rntransform = function(data){
  out = rank(data) - 0.5
  out[is.na(data)] = NA
  out = out/(max(out, na.rm = T) + 0.5)
  out = qnorm(out)
  out
}

for(i in 2:ncol(somascan_normalised_log_adj2_rint)){
  somascan_normalised_log_adj2_rint[,i] = rntransform(somascan_normalised_log_adj2_rint[,i])
  print(i)
}


somascan_normalised_log_adj_rint2 <- somascan_normalised_log_adj2[,-c(2:16)]

names(somascan_normalised_log_adj_rint2)

for (i in c(2:ncol(somascan_normalised_log_adj_rint2))) {
  somascan_normalised_log_adj_rint2[,i] <- RankNorm(somascan_normalised_log_adj_rint2[,i],)
}


### pca for log

# normalised

comp_normal_log_adj_rint=svd(somascan_normalised_log_adj_rint[,2:ncol(somascan_normalised_log_adj_rint)])

means=c(1:10)

for (i in c(1:10)) {
  means[i]=mean(comp_normal_log_adj_rint$u[,i])
}

stdevs=c(1:10)

for (i in c(1:10)) {
  stdevs[i]=sd(comp_normal_log_adj_rint$u[,i])
}

Zs=comp_normal_log_adj_rint$u[,c(1:10)]

for (i in c(1:length(Zs[,1]))) {
  Zs[i,]=(Zs[i,]-means)/stdevs
}

exclude=c(1:nrow(somascan_normalised_log_adj_rint))

for (i in c(1:nrow(somascan_normalised_log_adj_rint))) {
  exclude[i]=max(abs(Zs[i,c(1:10)]))>5 # first 10 PCs based on scree plot
}

list_ids_rm_normal_log_adj_rint=somascan_normalised_log_adj_rint$csid[exclude==T]


# non-normalised

comp_non_normal_log_adj_rint=svd(somascan_non_normalised_log_adj_rint[,2:ncol(somascan_non_normalised_log_adj_rint)])

means=c(1:10)

for (i in c(1:10)) {
  means[i]=mean(comp_non_normal_log_adj_rint$u[,i])
}

stdevs=c(1:10)

for (i in c(1:10)) {
  stdevs[i]=sd(comp_non_normal_log_adj_rint$u[,i])
}

Zs=comp_non_normal_log_adj_rint$u[,c(1:10)]

for (i in c(1:length(Zs[,1]))) {
  Zs[i,]=(Zs[i,]-means)/stdevs
}

exclude=c(1:nrow(somascan_non_normalised_log_adj_rint))

for (i in c(1:nrow(somascan_non_normalised_log_adj_rint))) {
  exclude[i]=max(abs(Zs[i,c(1:10)]))>5 # first 10 PCs based on scree plot
}

list_ids_rm_non_normal_log_adj_rint=somascan_non_normalised_log_adj_rint$csid[exclude==T]


## exclude individuals

somascan_normalised_log_clean <- somascan_normalised_log
somascan_normalised_log_clean[somascan_normalised_log_clean$csid %in% list_ids_rm_normal_log_adj_rint, 2:ncol(somascan_normalised_log_clean)] <- NA
colSums(is.na(somascan_normalised_log_clean))

somascan_non_normalised_log_clean <- somascan_normalised_log
somascan_non_normalised_log_clean[somascan_non_normalised_log_clean$csid %in% list_ids_rm_non_normal_log_adj_rint, 2:ncol(somascan_non_normalised_log_clean)] <- NA
colSums(is.na(somascan_non_normalised_log_clean))



# merge 

somascan_normalised_log_clean_adj <- merge(covar,somascan_normalised_log_clean,by="csid")
somascan_non_normalised_log_clean_adj <- merge(covar,somascan_non_normalised_log_clean,by="csid")

# adjust

somascan_normalised_log_clean_adj <- data_adjust(somascan_normalised_log_clean_adj,effect=names(covar)[-1],select=17:ncol(somascan_normalised_log_clean_adj))
somascan_non_normalised_log_clean_adj <- data_adjust(somascan_non_normalised_log_clean_adj,effect=names(covar)[-1],select=17:ncol(somascan_non_normalised_log_clean_adj))



## RINT

somascan_normalised_log_clean_adj_rint <- somascan_normalised_log_clean_adj[complete.cases(somascan_normalised_log_clean_adj),]
somascan_normalised_log_clean_adj_rint <- somascan_normalised_log_clean_adj_rint[,-c(2:16)]

somascan_non_normalised_log_clean_adj_rint <- somascan_non_normalised_log_clean_adj[complete.cases(somascan_non_normalised_log_clean_adj),]
somascan_non_normalised_log_clean_adj_rint <- somascan_non_normalised_log_clean_adj_rint[,-c(2:16)]

names(somascan_normalised_log_clean_adj_rint)

for (i in c(2:ncol(somascan_normalised_log_clean_adj_rint))) {
  somascan_normalised_log_clean_adj_rint[,i] <- RankNorm(somascan_normalised_log_clean_adj_rint[,i],)
  somascan_non_normalised_log_clean_adj_rint[,i] <- RankNorm(somascan_non_normalised_log_clean_adj_rint[,i])
}


# link to ccvid

linkage <- gwas[,c("csid","ccvid")]

somascan_normalised_log_clean_adj_rint_link <- merge(linkage,somascan_normalised_log_clean_adj_rint)

somascan_non_normalised_log_clean_adj_rint_link <- merge(linkage,somascan_non_normalised_log_clean_adj_rint)

write.csv(somascan_normalised_log_clean_adj_rint_link[,1:17],"somascan_normalised_log_clean_adj_rint_link.csv", quote=F, row.names=F)

write.csv(somascan_non_normalised_log_clean_adj_rint_link[,1:17],"somascan_non_normalised_log_clean_adj_rint_link.csv", quote=F, row.names=F)
