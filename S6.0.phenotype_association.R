### this script tests associations between protein levels and baseline characteristics in ckb

rm(list = ls())

setwd("")

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(egg)

## load proteomics data

olink <- read.csv("olink.csv")
somascan_normalised_log <- read.csv("somascan_normalised_log.csv")
somascan_non_normalised_log <- read.csv("somascan_non_normalised_log.csv")

## load overlapping proteins

overlap_cor <- read.csv("overlap_cor.csv")

## only keep overlapping proteins

olink <- olink[,c(1, which(names(olink) %in% overlap_cor$olink_id))]
somascan_normalised_log <- somascan_normalised_log[,c(1, which(names(somascan_normalised_log) %in% overlap_cor$somascan_id))]
somascan_non_normalised_log <- somascan_non_normalised_log[,c(1, which(names(somascan_non_normalised_log) %in% overlap_cor$somascan_id))]

## standardise

olink[,c(2:ncol(olink))] <- scale(olink[,c(2:ncol(olink))])
somascan_normalised_log[,c(2:ncol(somascan_normalised_log))] <- scale(somascan_normalised_log[,c(2:ncol(somascan_normalised_log))])
somascan_non_normalised_log[,c(2:ncol(somascan_non_normalised_log))] <- scale(somascan_non_normalised_log[,c(2:ncol(somascan_non_normalised_log))])

## load plate id

olink_plate <- read.csv("")
olink_plate <- olink_plate[olink_plate$panel_full=="Cardiometabolic",c(1,5)]
names(olink_plate)[2] <- "olink_plt_id"

somascan_plate <- read.csv("")[c("csid","qc","plateid")]
somascan_plate <- somascan_plate[somascan_plate$qc==0,]
somascan_plate <- somascan_plate[c(1,3)]
names(somascan_plate)[2] <- "somascan_plt_id"

## load other variables

baseline <- read.csv("")

names(baseline)

baseline_var <- baseline[,c("csid",
                            "region_mean_temp", "hours_since_last_ate_x10", 
                            "age_at_study_date_x100","is_female", "region_code", "region_is_urban", 
                            "marital_status", "highest_education", "household_income",
                            "bmi_calc", "dbp_mean", "sbp_mean", "heart_rate_mean", "random_glucose_x10",
                            "met", "self_rated_health", "diabetes_diag", "kidney_dis_diag", "cancer_diag",
                            "smoking_category", "alcohol_category")]

names(baseline_var)

# mean temp2
baseline_var$region_mean_temp2 <- baseline_var$region_mean_temp * baseline_var$region_mean_temp
summary(baseline_var$region_mean_temp2)

# hours since last ate
baseline_var$hours_since_last_ate <- baseline_var$hours_since_last_ate_x10/10
summary(baseline_var$hours_since_last_ate)

# hours since last ate2
baseline_var$hours_since_last_ate2 <- baseline_var$hours_since_last_ate * baseline_var$hours_since_last_ate
summary(baseline_var$hours_since_last_ate2)

# age
baseline_var$age <- baseline_var$age_at_study_date_x100/100
summary(baseline_var$age)

# age2
baseline_var$age2 <- baseline_var$age * baseline_var$age
summary(baseline_var$age2)

# sex
baseline_var$is_female <- factor(baseline_var$is_female)
summary(baseline_var$is_female)

# region
baseline_var$region_code <- factor(baseline_var$region_code)
summary(baseline_var$region_code)

# urban
baseline_var$region_is_urban <- factor(baseline_var$region_is_urban)
summary(baseline_var$region_is_urban)

# married
table(baseline_var$marital_status)
sum(is.na(baseline_var$marital_status))
baseline_var$married <- 0
baseline_var$married[baseline_var$marital_status==0] <- 1
baseline_var$married <- factor(baseline_var$married)
summary(baseline_var$married)

# went to school >= 6 years
table(baseline_var$highest_education)
sum(is.na(baseline_var$highest_education))
baseline_var$school <- 1
baseline_var$school[baseline_var$highest_education==0 | baseline_var$highest_education==1] <- 0
baseline_var$school <- factor(baseline_var$school)
summary(baseline_var$school)

############## household income?

# random glucose
baseline_var$random_glucose <- baseline_var$random_glucose_x10/10
summary(baseline_var$random_glucose)

# poor self-rated health
table(baseline_var$self_rated_health)
sum(is.na(baseline_var$self_rated_health))
baseline_var$poor_health <- 0
baseline_var$poor_health[baseline_var$self_rated_health==3] <- 1
baseline_var$poor_health <- factor(baseline_var$poor_health)
summary(baseline_var$poor_health)

# diabetes
baseline_var$diabetes_diag <- factor(baseline_var$diabetes_diag)
summary(baseline_var$diabetes_diag)

# kidney disease
baseline_var$kidney_dis_diag <- factor(baseline_var$kidney_dis_diag)
summary(baseline_var$kidney_dis_diag)

# cancer
baseline_var$cancer_diag <- factor(baseline_var$cancer_diag)
summary(baseline_var$cancer_diag)

# smoking ever regular
table(baseline_var$smoking_category)
sum(is.na(baseline_var$smoking_category))
baseline_var$smoking_ever_regular <- 0
baseline_var$smoking_ever_regular[baseline_var$smoking_category==3 | baseline_var$smoking_category==4] <- 1
baseline_var$smoking_ever_regular[baseline_var$is_female==1] <- NA
baseline_var$smoking_ever_regular <- factor(baseline_var$smoking_ever_regular)
summary(baseline_var$smoking_ever_regular)

# alcohol current regular vs occasion
table(baseline_var$alcohol_category)
sum(is.na(baseline_var$alcohol_category))
baseline_var$alcohol_regular_vs_occasion <- NA
baseline_var$alcohol_regular_vs_occasion[baseline_var$alcohol_category==3 | baseline_var$alcohol_category==4] <- 0
baseline_var$alcohol_regular_vs_occasion[baseline_var$alcohol_category==6] <- 1
baseline_var$alcohol_regular_vs_occasion[baseline_var$is_female==1] <- NA
baseline_var$alcohol_regular_vs_occasion <- factor(baseline_var$alcohol_regular_vs_occasion)
summary(baseline_var$alcohol_regular_vs_occasion)

## load ascertainment

ascertainment <- read.csv("")
table(ascertainment$olinkexp1536_chd_b1_subcohort)
ascertainment <- ascertainment[,c("csid","olinkexp1536_chd_b1_subcohort")]
names(ascertainment) <- c("csid","ascertainment")
ascertainment <- ascertainment[!is.na(ascertainment$ascertainment),]
ascertainment$ascertainment <- factor(ascertainment$ascertainment,levels=c("1","0"))
table(ascertainment$ascertainment)
# subcohort_id <- ascertainment$csid[ascertainment$olinkexp1536_chd_b1_subcohort==1]

## load blood sample processing time

## merge

olink_pheno <- merge(olink, baseline_var, by="csid", all.x=T)
olink_pheno <- merge(olink_pheno, olink_plate, by="csid", all.x=T)
olink_pheno <- merge(olink_pheno, ascertainment, by="csid", all.x=T)

somascan_normalised_log_pheno <- merge(somascan_normalised_log, baseline_var, by="csid", all.x=T)
somascan_normalised_log_pheno <- merge(somascan_normalised_log_pheno, somascan_plate, by="csid", all.x=T)
somascan_normalised_log_pheno <- merge(somascan_normalised_log_pheno, ascertainment, by="csid", all.x=T)

somascan_non_normalised_log_pheno <- merge(somascan_non_normalised_log, baseline_var, by="csid", all.x=T)
somascan_non_normalised_log_pheno <- merge(somascan_non_normalised_log_pheno, somascan_plate, by="csid", all.x=T)
somascan_non_normalised_log_pheno <- merge(somascan_non_normalised_log_pheno, ascertainment, by="csid", all.x=T)

# keep same participants
olink_pheno <- olink_pheno[olink_pheno$csid %in% somascan_normalised_log_pheno$csid, ]

# save
write.csv(olink_pheno, "olink_pheno.csv", row.names = F, quote = F)
write.csv(somascan_normalised_log_pheno, "somascan_normalised_log_pheno.csv", row.names = F, quote = F)
write.csv(somascan_non_normalised_log_pheno, "somascan_non_normalised_log_pheno.csv", row.names = F, quote = F)

############################################################## run regression

overlap_assoc <- overlap_cor

## firstly, run variables that need to remove certain covariates

# mean temp

for (i in 1:nrow(overlap_assoc)){
  
  olink_id <- overlap_assoc$olink_id[i]
  somascan_id <- overlap_assoc$somascan_id[i]
  
  model <- lm(paste(olink_id, "~", "region_mean_temp + hours_since_last_ate + hours_since_last_ate2 + 
                     olink_plt_id + ascertainment + age + age2 + is_female + region_code"), olink_pheno)
  overlap_assoc$olink_es[i] <- model$coefficients[2]
  overlap_assoc$olink_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$olink_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]  
  
  model <- lm(paste(somascan_id, "~", "region_mean_temp + hours_since_last_ate + hours_since_last_ate2 + 
                     somascan_plt_id + ascertainment + age + age2 + is_female + region_code"), somascan_normalised_log_pheno)
  overlap_assoc$soma_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  model <- lm(paste(somascan_id, "~", "region_mean_temp + hours_since_last_ate + hours_since_last_ate2 + 
                     somascan_plt_id + ascertainment + age + age2 + is_female + region_code"), somascan_non_normalised_log_pheno)
  overlap_assoc$soma_non_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_non_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_non_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  print(i)
}

names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)] <- 
  paste(names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)],"region_mean_temp",sep="_")

# hours since last ate

for (i in 1:nrow(overlap_assoc)){
  
  olink_id <- overlap_assoc$olink_id[i]
  somascan_id <- overlap_assoc$somascan_id[i]
  
  model <- lm(paste(olink_id, "~", "hours_since_last_ate + region_mean_temp + region_mean_temp2 +
                     olink_plt_id + ascertainment + age + age2 + is_female + region_code"), olink_pheno)
  overlap_assoc$olink_es[i] <- model$coefficients[2]
  overlap_assoc$olink_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$olink_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]  
  
  model <- lm(paste(somascan_id, "~", "hours_since_last_ate + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + age + age2 + is_female + region_code"), somascan_normalised_log_pheno)
  overlap_assoc$soma_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  model <- lm(paste(somascan_id, "~", "hours_since_last_ate + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + age + age2 + is_female + region_code"), somascan_non_normalised_log_pheno)
  overlap_assoc$soma_non_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_non_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_non_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  print(i)
}

names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)] <- 
  paste(names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)],"hours_since_last_ate",sep="_")

# age

for (i in 1:nrow(overlap_assoc)){
  
  olink_id <- overlap_assoc$olink_id[i]
  somascan_id <- overlap_assoc$somascan_id[i]
  
  model <- lm(paste(olink_id, "~", "age + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     olink_plt_id + ascertainment + is_female + region_code"), olink_pheno)
  overlap_assoc$olink_es[i] <- model$coefficients[2]
  overlap_assoc$olink_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$olink_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]  
  
  model <- lm(paste(somascan_id, "~", "age + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + is_female + region_code"), somascan_normalised_log_pheno)
  overlap_assoc$soma_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  model <- lm(paste(somascan_id, "~", "age + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + is_female + region_code"), somascan_non_normalised_log_pheno)
  overlap_assoc$soma_non_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_non_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_non_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  print(i)
}

names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)] <- 
  paste(names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)],"age",sep="_")

# sex

for (i in 1:nrow(overlap_assoc)){
  
  olink_id <- overlap_assoc$olink_id[i]
  somascan_id <- overlap_assoc$somascan_id[i]
  
  model <- lm(paste(olink_id, "~", "is_female + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     olink_plt_id + ascertainment + age + age2 + region_code"), olink_pheno)
  overlap_assoc$olink_es[i] <- model$coefficients[2]
  overlap_assoc$olink_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$olink_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]  
  
  model <- lm(paste(somascan_id, "~", "is_female + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + age + age2 + region_code"), somascan_normalised_log_pheno)
  overlap_assoc$soma_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  model <- lm(paste(somascan_id, "~", "is_female + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + age + age2 + region_code"), somascan_non_normalised_log_pheno)
  overlap_assoc$soma_non_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_non_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_non_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  print(i)
}

names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)] <- 
  paste(names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)],"is_female",sep="_")

# urban

for (i in 1:nrow(overlap_assoc)){
  
  olink_id <- overlap_assoc$olink_id[i]
  somascan_id <- overlap_assoc$somascan_id[i]
  
  model <- lm(paste(olink_id, "~", "region_is_urban + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     olink_plt_id + ascertainment + age + age2 + is_female"), olink_pheno)
  overlap_assoc$olink_es[i] <- model$coefficients[2]
  overlap_assoc$olink_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$olink_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]  
  
  model <- lm(paste(somascan_id, "~", "region_is_urban + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + age + age2 + is_female"), somascan_normalised_log_pheno)
  overlap_assoc$soma_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  model <- lm(paste(somascan_id, "~", "region_is_urban + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + age + age2 + is_female"), somascan_non_normalised_log_pheno)
  overlap_assoc$soma_non_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_non_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_non_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  print(i)
}

names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)] <- 
  paste(names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)],"region_is_urban",sep="_")

# smoking

for (i in 1:nrow(overlap_assoc)){
  
  olink_id <- overlap_assoc$olink_id[i]
  somascan_id <- overlap_assoc$somascan_id[i]
  
  model <- lm(paste(olink_id, "~", "smoking_ever_regular + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     olink_plt_id + ascertainment + age + age2 + region_code"), olink_pheno)
  overlap_assoc$olink_es[i] <- model$coefficients[2]
  overlap_assoc$olink_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$olink_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]  
  
  model <- lm(paste(somascan_id, "~", "smoking_ever_regular + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + age + age2 + region_code"), somascan_normalised_log_pheno)
  overlap_assoc$soma_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  model <- lm(paste(somascan_id, "~", "smoking_ever_regular + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + age + age2 + region_code"), somascan_non_normalised_log_pheno)
  overlap_assoc$soma_non_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_non_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_non_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  print(i)
}

names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)] <- 
  paste(names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)],"smoking_ever_regular",sep="_")

# alcohol drinking

for (i in 1:nrow(overlap_assoc)){
  
  olink_id <- overlap_assoc$olink_id[i]
  somascan_id <- overlap_assoc$somascan_id[i]
  
  model <- lm(paste(olink_id, "~", "alcohol_regular_vs_occasion + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     olink_plt_id + ascertainment + age + age2 + region_code"), olink_pheno)
  overlap_assoc$olink_es[i] <- model$coefficients[2]
  overlap_assoc$olink_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$olink_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]  
  
  model <- lm(paste(somascan_id, "~", "alcohol_regular_vs_occasion + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + age + age2 + region_code"), somascan_normalised_log_pheno)
  overlap_assoc$soma_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  model <- lm(paste(somascan_id, "~", "alcohol_regular_vs_occasion + hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + age + age2 + region_code"), somascan_non_normalised_log_pheno)
  overlap_assoc$soma_non_normal_es[i] <- model$coefficients[2]
  overlap_assoc$soma_non_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
  overlap_assoc$soma_non_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
  
  print(i)
}

names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)] <- 
  paste(names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)],"alcohol_regular_vs_occasion",sep="_")

## finally, all other variables

variable <- c("married", "school", "bmi_calc", "dbp_mean", "sbp_mean", "heart_rate_mean",
              "random_glucose", "met", "poor_health", "diabetes_diag", "kidney_dis_diag", "cancer_diag")

for (x in 1:length(variable)) {
  
  print(variable[x])
  
  for (i in 1:nrow(overlap_assoc)){
    
    olink_id <- overlap_assoc$olink_id[i]
    somascan_id <- overlap_assoc$somascan_id[i]
    
    model <- lm(paste(olink_id, "~", variable[x], "+ hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     olink_plt_id + ascertainment + age + age2 + is_female + region_code"), olink_pheno)
    overlap_assoc$olink_es[i] <- model$coefficients[2]
    overlap_assoc$olink_se[i] <- summary(model)$coefficients[2,"Std. Error"]
    overlap_assoc$olink_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]  
    
    model <- lm(paste(somascan_id, "~", variable[x], "+ hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + age + age2 + is_female + region_code"), somascan_normalised_log_pheno)
    overlap_assoc$soma_normal_es[i] <- model$coefficients[2]
    overlap_assoc$soma_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
    overlap_assoc$soma_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
    
    model <- lm(paste(somascan_id, "~", variable[x], "+ hours_since_last_ate + hours_since_last_ate2 + region_mean_temp + region_mean_temp2 +
                     somascan_plt_id + ascertainment + age + age2 + is_female + region_code"), somascan_non_normalised_log_pheno)
    overlap_assoc$soma_non_normal_es[i] <- model$coefficients[2]
    overlap_assoc$soma_non_normal_se[i] <- summary(model)$coefficients[2,"Std. Error"]
    overlap_assoc$soma_non_normal_p[i] <- summary(model)$coefficients[2,"Pr(>|t|)"]   
    
    print(i)
  }
  
  names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)] <- 
    paste(names(overlap_assoc)[(ncol(overlap_assoc)-8):ncol(overlap_assoc)],variable[x],sep="_")
  
}

write.csv(overlap_assoc, "overlap_assoc.csv", row.names = F, quote = F)


# keep 1 to 1

protein_dup <- unique(overlap_assoc$uniprot_id[duplicated(overlap_assoc$uniprot_id)]) # get proteins targeted by multiple aptamers

aptamer_dup <- unique(overlap_assoc$somascan_id[duplicated(overlap_assoc$somascan_id)]) # get aptamers targeted by multiple proteins

overlap_1_to_1_assoc <- overlap_assoc[-c(which(overlap_assoc$uniprot_id %in% protein_dup),which(overlap_assoc$somascan_id %in% aptamer_dup)), ]

# save

write.csv(overlap_1_to_1_assoc,"overlap_1_to_1_assoc.csv", quote=F, row.names=F)
