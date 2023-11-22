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

# plot

xlevel <- levels(reorder(droplevels(boruta_normal_melt_confirmed$variable), boruta_normal_melt_confirmed$importance, FUN = median))

hex <- brewer.pal(3,"Set1")[c(1,2)]

boruta_normal_plot_confirmed <- ggplot(boruta_normal_melt_confirmed, aes(x=factor(variable,levels = xlevel), y=importance)) + 
  geom_boxplot(aes(fill=factor(direction,levels=c("Positive","Inverse","Unconfirmed"))),lwd=0.2,outlier.size=0.2) + 
  coord_flip() +
  ylim(-5,45) +
  ggtitle("OLINK vs SomaScan-ANML") +
  xlab("Feature") +
  ylab("Importance") +
  labs(fill="Direction of association") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_x_discrete(labels=c("OLINK panel: Oncology I",
                            "Protein mass",
                            "OLINK panel: Inflammation II",
                            "OLINK panel: Cardiometabolic II",
                            "OLINK dilution: 1:100000",
                            "OLINK panel: Neurology II",
                            "OLINK dilution: 1:1000",
                            "OLINK panel: Oncology II",
                            "OLINK panel: Cardiometabolic I",
                            "OLINK dilution: 1:100",
                            "SomaScan dilution: 0.005%",
                            "SomaScan dilution: 0.5%",
                            "SomaScan Apparent Kd (M)",
                            "SomaScan QC flag",
                            "OLINK batch 2",
                            "SomaScan dilution 20%",
                            "OLINK dilution: 1:10",
                            "OLINK mean on NPX scale",
                            "OLINK % outliers",
                            "OLINK skewness",
                            "OLINK dilution: 1:1",
                            "SomaScan mean on log scale",
                            "OLINK QC warning",
                            "OLINK kurtosis",
                            "SomaScan % outliers",
                            "SomaScan skewness",
                            "SomaScan kurtosis",
                            "OLINK % below LOD")) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c(hex,"#808080"),
                    labels=c("Positive","Inverse","Unconfirmed")) 

boruta_normal_plot_confirmed


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

# add variables that are not in the ANML to the dataframe

setdiff(unique(boruta_normal_melt_confirmed$variable),variable)

# variable_extra <- c("lipidation")
# 
# boruta_non_normal_extra <- data.frame(variable=variable_extra,
#                                   coeff=numeric(1),
#                                   p=numeric(1),
#                                   direction="Rejected")

# boruta_non_normal_regress <- rbind(boruta_non_normal_regress,boruta_non_normal_extra)

boruta_non_normal_melt_confirmed <- merge(boruta_non_normal_melt,boruta_non_normal_regress,by="variable",all.y=T)

# plot

xlevel_non_normal <- str_replace(xlevel, "normal", "non_normal")

# boruta_non_normal_plot_confirmed <- ggplot(boruta_non_normal_melt_confirmed, aes(x=reorder(variable, importance, FUN = median), y=importance)) + 
#   geom_boxplot(aes(fill=factor(direction,levels=c("Positive","Inverse","Rejected"))),lwd=0.2,outlier.size=0.2) + 
#   coord_flip() +
#   ylim(-5,45) +
#   ggtitle("Boruta feature selection (non-ANML)") +
#   ylab("Importance") +
#   labs(fill="Direction of association") +
#   theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
#   theme(text = element_text(size = 12))
# 
# boruta_non_normal_plot_confirmed

boruta_non_normal_plot_confirmed <- ggplot(boruta_non_normal_melt_confirmed, aes(x=factor(variable, levels=xlevel_non_normal), y=importance)) + 
  geom_boxplot(aes(fill=factor(direction,levels=c("Positive","Inverse","Rejected"))),lwd=0.2,outlier.size=0.2) + 
  coord_flip() +
  ylim(-5,45) +
  ggtitle("OLINK vs SomaScan-non-ANML") +
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

# boruta_plot_confirmed <- ggarrange(boruta_normal_plot_confirmed,boruta_non_normal_plot_confirmed,legend,ncol=3,nrow=1)

# boruta_plot_confirmed <- grid.arrange(boruta_normal_plot_confirmed,boruta_non_normal_plot_confirmed,legend,widths=c(3,3,1))

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

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/boruta_plot_confirmed.png",combined,width=14,height=6)

# combine with histogram

# load("plot_rho_olink_soma.RData")

# figure2 <- ggarrange(plot_rho_olink_soma_normal,plot_rho_olink_soma_non_normal,boruta_normal_plot_confirmed,boruta_non_normal_plot_confirmed,
#                      ncol=2,nrow=2)
# 
# ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/figure2.png",figure2,width=16,height=9)
