rm(list = ls())

library(datawizard)
library(gtsummary)

setwd("")

olink_pheno <- read.csv("olink_pheno.csv")

pheno <- olink_pheno[,c(1,match("region_mean_temp",names(olink_pheno)):ncol(olink_pheno))]

names(pheno)

pheno <- pheno[,c("csid","region_mean_temp","hours_since_last_ate",
                  "age","is_female","region_is_urban","region_code","married", "school", 
                  "bmi_calc", "dbp_mean", "sbp_mean", "heart_rate_mean","random_glucose", 
                  "poor_health", "diabetes_diag", "kidney_dis_diag", "cancer_diag",
                  "met","smoking_category","alcohol_category","smoking_ever_regular","alcohol_regular_vs_occasion",
                  "ascertainment")]

pheno$is_female <- factor(pheno$is_female)
pheno$region_is_urban <- factor(pheno$region_is_urban)
pheno$region_code <- factor(pheno$region_code)
pheno$married <- factor(pheno$married)
pheno$school <- factor(pheno$school)
pheno$poor_health <- factor(pheno$poor_health)
pheno$diabetes_diag <- factor(pheno$diabetes_diag)
pheno$kidney_dis_diag <- factor(pheno$kidney_dis_diag)
pheno$cancer_diag <- factor(pheno$cancer_diag)
pheno$ascertainment <- factor(pheno$ascertainment,levels=c(0,1))

names(pheno)

pheno_adj <- data_adjust(pheno, effect=c("age","is_female","region_code"),
                         select=c("region_mean_temp","hours_since_last_ate","bmi_calc", "dbp_mean", "sbp_mean", "heart_rate_mean","random_glucose","met"),
                         keep_intercept=T)

pheno_adj <- data_adjust(pheno,effect=c("is_female","region_code"),select="age",keep_intercept=T)

pheno_adj_1 <- pheno_adj[,-which(names(pheno_adj) %in% c("csid","region_code","smoking_category","alcohol_category","smoking_ever_regular","alcohol_regular_vs_occasion"))]

pheno_adj_1 %>%
  tbl_summary(by = ascertainment,
              statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    all_categorical() ~ "{p}"),
  digits = everything() ~ 1) %>%
  add_overall(last=T)

pheno_adj_2 <- pheno_adj[,c("is_female","smoking_category","alcohol_category","smoking_ever_regular","alcohol_regular_vs_occasion","ascertainment")]

pheno_adj_2$smoking_ever_regular <- 0
pheno_adj_2$smoking_ever_regular[pheno_adj_2$smoking_category==3 | pheno_adj_2$smoking_category==4] <- 1

pheno_adj_2$alcohol_current_regular <- 0
pheno_adj_2$alcohol_current_regular[pheno_adj_2$alcohol_category==6] <- 1

pheno_adj_2[pheno_adj_2$is_female==0,] %>%
  tbl_summary(by = ascertainment,
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{p}"
              ),
              digits = everything() ~ 1) %>%
  add_overall(last=T)

pheno_adj_2[pheno_adj_2$is_female==1,] %>%
  tbl_summary(by = ascertainment,
              statistic = list(
                all_continuous() ~ "{mean} ({sd})",
                all_categorical() ~ "{p}"
              ),
              digits = everything() ~ 1) %>%
  add_overall(last=T)
