### this script plots the results of Boruta

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

# read data and define levels of batch and dilution

overlap_1_to_1_annot_normal <- readRDS("overlap_1_to_1_annot_normal.RDS")
overlap_1_to_1_annot_normal$batch <- factor(overlap_1_to_1_annot_normal$batch, levels = c("1","2"))

overlap_1_to_1_annot_non_normal <- readRDS("overlap_1_to_1_annot_non_normal.RDS")
overlap_1_to_1_annot_non_normal$batch <- factor(overlap_1_to_1_annot_non_normal$batch, levels = c("1","2"))


## load data

boruta_normal <- readRDS("boruta_normal.RDS")

boruta_non_normal <- readRDS("boruta_non_normal.RDS")

boruta_normal_df <- attStats(boruta_normal)

boruta_non_normal_df <- attStats(boruta_non_normal)

# reshape data

boruta_normal_melt <- melt(as.data.frame(boruta_normal$ImpHistory))

boruta_normal_melt <- boruta_normal_melt[is.finite(boruta_normal_melt$value),]

boruta_non_normal_melt <- melt(as.data.frame(boruta_non_normal$ImpHistory))

boruta_non_normal_melt <- boruta_non_normal_melt[is.finite(boruta_non_normal_melt$value),]

# remove shadow variables

boruta_normal_melt <- boruta_normal_melt[-grep("shadow",boruta_normal_melt$variable),]

names(boruta_normal_melt) <- c("variable","importance")

boruta_non_normal_melt <- boruta_non_normal_melt[-grep("shadow",boruta_non_normal_melt$variable),]

names(boruta_non_normal_melt) <- c("variable","importance")

# Classify by decision

boruta_normal_decision <- data.frame(variable=row.names(boruta_normal_df),decision=as.vector(boruta_normal_df$decision))

boruta_normal_melt <- merge(boruta_normal_melt,boruta_normal_decision,by="variable")

boruta_normal_melt$decision <- factor(boruta_normal_melt$decision,levels=c("Confirmed","Tentative","Rejected"))

boruta_non_normal_decision <- data.frame(variable=row.names(boruta_non_normal_df),decision=as.vector(boruta_non_normal_df$decision))

boruta_non_normal_melt <- merge(boruta_non_normal_melt,boruta_non_normal_decision,by="variable")

boruta_non_normal_melt$decision <- factor(boruta_non_normal_melt$decision,levels=c("Confirmed","Tentative","Rejected"))

# plot only significant associations and classify by direction

boruta_normal_melt_confirmed <- boruta_normal_melt[boruta_normal_melt$decision=="Confirmed",]

boruta_non_normal_melt_confirmed <- boruta_non_normal_melt[boruta_non_normal_melt$decision=="Confirmed",]

# normal

variable <- unique(boruta_normal_melt_confirmed$variable)

# run regression

boruta_normal_regress <- data.frame(variable=variable)

boruta_normal_regress$coeff <- NA

boruta_normal_regress$p <- NA

overlap_1_to_1_annot_normal$qc_check_soma <- factor(overlap_1_to_1_annot_normal$qc_check_soma, levels=c("PASS","FLAG"))

for (i in 1:(length(variable))) {
  lm_results <- lm(paste("rho_olink_soma_normal ~", variable[i]), data = overlap_1_to_1_annot_normal)
  boruta_normal_regress$coeff[i] <- lm_results$coefficients[2]
  boruta_normal_regress$p[i] <- summary(lm_results)$coefficients[2,"Pr(>|t|)"]
}

boruta_normal_regress$direction <- "Inverse"

boruta_normal_regress$direction[boruta_normal_regress$coeff>0] <- "Positive"

boruta_normal_regress

# add variables that are not in the ANML to the dataframe

setdiff(unique(boruta_non_normal_melt_confirmed$variable),variable)

variable_extra <- c("mass","panel_oncology")

boruta_normal_extra <- data.frame(variable=variable_extra,
                                  coeff=numeric(length(variable_extra)),
                                  p=numeric(length(variable_extra)),
                                  direction="Rejected")

boruta_normal_regress <- rbind(boruta_normal_regress,boruta_normal_extra)

boruta_normal_melt_confirmed <- merge(boruta_normal_melt,boruta_normal_regress,by="variable",all.y = T)


## non_normal

# plot only significant associations and classify by direction

boruta_non_normal_melt_confirmed <- boruta_non_normal_melt[boruta_non_normal_melt$decision=="Confirmed",]

variable <- unique(boruta_non_normal_melt_confirmed$variable)

boruta_non_normal_regress <- data.frame(variable=variable)

boruta_non_normal_regress$coeff <- NA

boruta_non_normal_regress$p <- NA

overlap_1_to_1_annot_non_normal$qc_check_soma <- factor(overlap_1_to_1_annot_non_normal$qc_check_soma, levels=c("PASS","FLAG"))

for (i in 1:(length(variable))) {
  lm_results <- lm(paste("rho_olink_soma_non_normal ~", variable[i]), data = overlap_1_to_1_annot_non_normal)
  boruta_non_normal_regress$coeff[i] <- lm_results$coefficients[2]
  boruta_non_normal_regress$p[i] <- summary(lm_results)$coefficients[2,"Pr(>|t|)"]
}

boruta_non_normal_regress$direction <- "Inverse"

boruta_non_normal_regress$direction[boruta_non_normal_regress$coeff>0] <- "Positive"

boruta_non_normal_regress

boruta_non_normal_melt_confirmed <- merge(boruta_non_normal_melt,boruta_non_normal_regress,by="variable",all.y=T)

xlevel_non_normal <- str_replace(xlevel, "normal", "non_normal")

boruta_non_normal_plot_confirmed <- ggplot(boruta_non_normal_melt_confirmed, aes(x=factor(variable, levels=xlevel_non_normal), y=importance)) + 
  geom_boxplot(aes(fill=factor(direction,levels=c("Positive","Inverse","Rejected"))),lwd=0.2,outlier.size=0.2) + 
  coord_flip() +
  ylim(-5,45) +
  ggtitle("Olink vs SomaScan-non-ANML") +
  xlab("") +
  ylab("Importance") +
  labs(fill="Direction of association") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12),axis.text.y = element_blank()) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c(hex,"#808080"),
                    labels=c("Positive","Inverse","Unconfirmed")) 

boruta_non_normal_plot_confirmed

legend <- as_ggplot(get_legend( ggplot(boruta_normal_melt_confirmed, aes(x=factor(variable,levels = xlevel), y=importance)) + 
                                  geom_boxplot(aes(fill=factor(direction,levels=c("Positive","Inverse","Rejected"))),lwd=0.2,outlier.size=0.2) + 
                                  coord_flip() +
                                  ylim(-5,45) +
                                  labs(fill="Direction of association") +
                                  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
                                  scale_fill_manual(values=c(hex,"#808080"),labels=c("Positive","Inverse","Rejected"))))

g1 <- ggplotGrob(boruta_normal_plot_confirmed)
g2 <- ggplotGrob(boruta_non_normal_plot_confirmed)
g3 <- ggplotGrob(legend)

fg12 <-
  gtable_frame(gtable_cbind(g1, g2),
               width = unit(6, "null"),
               height = unit(1, "null"))

fg3 <-
  gtable_frame(
    g3,
    width = unit(1, "null"),
    height = unit(1, "null"))

grid.newpage()
combined <- gtable_cbind(fg12, fg3)
grid.draw(combined)

ggsave("",combined,width=14,height=6)
