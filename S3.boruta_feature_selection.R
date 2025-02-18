### this script runs Boruta feature selection on HPC cluster

rm(list = ls())

setwd("")

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
library(ggthemes)
library(egg)
library(scales)
library(grid)

overlap_1_to_1_annot_non_ANML <- readRDS("overlap_1_to_1_annot_non_ANML.RDS")

##### define levels of batch

overlap_1_to_1_annot_non_ANML$batch <- factor(overlap_1_to_1_annot_non_ANML$batch, levels = c("1","2"))

##### run boruta

set.seed(47)

boruta_non_ANML <- Boruta(rho_olink_soma_non_ANML~., data = overlap_1_to_1_annot_non_ANML, doTrace = 2, maxRuns = 50000)

##### visualise results

# reshape data

boruta_non_ANML_melt <- melt(as.data.frame(boruta_non_ANML$ImpHistory))

boruta_non_ANML_melt <- boruta_non_ANML_melt[is.finite(boruta_non_ANML_melt$value),]

# remove shadow variables

boruta_non_ANML_melt <- boruta_non_ANML_melt[-grep("shadow",boruta_non_ANML_melt$variable),]

names(boruta_non_ANML_melt) <- c("variable","importance")

# Classify by decision

boruta_non_ANML_decision <- data.frame(variable=row.names(boruta_non_ANML_df),decision=as.vector(boruta_non_ANML_df$decision))

boruta_non_ANML_melt <- merge(boruta_non_ANML_melt,boruta_non_ANML_decision,by="variable")

boruta_non_ANML_melt$decision <- factor(boruta_non_ANML_melt$decision,levels=c("Confirmed","Tentative","Rejected"))

# plot only significant associations and classify by direction

boruta_non_ANML_melt_confirmed <- boruta_non_ANML_melt[boruta_non_ANML_melt$decision=="Confirmed",]

# ANML

variable <- unique(boruta_ANML_melt_confirmed$variable)

# plot only significant associations and classify by direction

boruta_non_ANML_melt_confirmed <- boruta_non_ANML_melt[boruta_non_ANML_melt$decision=="Confirmed",]

variable <- unique(boruta_non_ANML_melt_confirmed$variable)

boruta_non_ANML_regress <- data.frame(variable=variable)

boruta_non_ANML_regress$coeff <- NA

boruta_non_ANML_regress$p <- NA

overlap_1_to_1_annot_non_ANML$qc_check_soma <- factor(overlap_1_to_1_annot_non_ANML$qc_check_soma, levels=c("PASS","FLAG"))

for (i in 1:(length(variable))) {
  lm_results <- lm(paste("rho_olink_soma_non_ANML ~", variable[i]), data = overlap_1_to_1_annot_non_ANML)
  boruta_non_ANML_regress$coeff[i] <- lm_results$coefficients[2]
  boruta_non_ANML_regress$p[i] <- summary(lm_results)$coefficients[2,"Pr(>|t|)"]
}

boruta_non_ANML_regress$direction <- "Inverse"

boruta_non_ANML_regress$direction[boruta_non_ANML_regress$coeff>0] <- "Positive"

boruta_non_ANML_regress

boruta_non_ANML_melt_confirmed <- merge(boruta_non_ANML_melt,boruta_non_ANML_regress,by="variable",all.y=T)

xlevel_non_ANML <- str_replace(xlevel, "ANML", "non_ANML")

boruta_non_ANML_plot_confirmed <- ggplot(boruta_non_ANML_melt_confirmed, aes(x=factor(variable, levels=xlevel_non_ANML), y=importance)) + 
  geom_boxplot(aes(fill=factor(direction,levels=c("Positive","Inverse","Rejected"))),lwd=0.2,outlier.size=0.2) + 
  coord_flip() +
  ylim(-5,45) +
  ggtitle("Olink vs SomaScan-non-ANML") +
  xlab("") +
  ylab("Importance") +
  labs(fill="Direction of association") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12),axis.text.y = element_blank()) +
  scale_fill_manual(values=c(hex,"#808080"),
                    labels=c("Positive","Inverse","Unconfirmed")) 

boruta_non_ANML_plot_confirmed
