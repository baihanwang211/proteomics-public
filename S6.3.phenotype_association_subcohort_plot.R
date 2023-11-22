rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(egg)
library(forcats)
library(ggpattern)
library(ckbplotr)
library(scales)
library(RColorBrewer)

overlap_1_to_1_assoc <- read.csv("overlap_assoc_subcohort.csv")

## number of hits

# define variables 

variable <- c("region_mean_temp","hours_since_last_ate","age","is_female","region_is_urban",
              "smoking_ever_regular","alcohol_regular_vs_occasion",
              "married", "school", "bmi_calc", "dbp_mean", "sbp_mean", "heart_rate_mean",
              "random_glucose", "met", "poor_health", "diabetes_diag", "kidney_dis_diag", "cancer_diag")

for (i in 1:length(variable)) {
  overlap_1_to_1_assoc$olink_sig_raw <-F
  overlap_1_to_1_assoc$olink_sig_raw[overlap_1_to_1_assoc[8+i*9]<0.05] <- T
  
  overlap_1_to_1_assoc$soma_normal_sig_raw <- F
  overlap_1_to_1_assoc$soma_normal_sig_raw[overlap_1_to_1_assoc[11+i*9]<0.05] <- T
  
  overlap_1_to_1_assoc$soma_non_normal_sig_raw <- F
  overlap_1_to_1_assoc$soma_non_normal_sig_raw[overlap_1_to_1_assoc[14+i*9]<0.05] <- T
  
  names(overlap_1_to_1_assoc)[(ncol(overlap_1_to_1_assoc)-2):ncol(overlap_1_to_1_assoc)] <-
    paste(names(overlap_1_to_1_assoc)[(ncol(overlap_1_to_1_assoc)-2):ncol(overlap_1_to_1_assoc)],variable[i],sep="_")
}







# bonferroni

names(overlap_1_to_1_assoc)

for (i in 1:length(variable)) {
  overlap_1_to_1_assoc$olink_p_bonferroni <- F
  overlap_1_to_1_assoc$olink_p_bonferroni <- p.adjust(overlap_1_to_1_assoc[[8+i*9]],method = "bonferroni")
  overlap_1_to_1_assoc$olink_sig_bonferroni <-F
  overlap_1_to_1_assoc$olink_sig_bonferroni[overlap_1_to_1_assoc$olink_p_bonferroni<0.05] <- T
  
  overlap_1_to_1_assoc$soma_normal_p_bonferroni <- F
  overlap_1_to_1_assoc$soma_normal_p_bonferroni <- p.adjust(overlap_1_to_1_assoc[[11+i*9]],method = "bonferroni")
  overlap_1_to_1_assoc$soma_normal_sig_bonferroni <-F
  overlap_1_to_1_assoc$soma_normal_sig_bonferroni[overlap_1_to_1_assoc$soma_normal_p_bonferroni<0.05] <- T
  
  overlap_1_to_1_assoc$soma_non_normal_p_bonferroni <- F
  overlap_1_to_1_assoc$soma_non_normal_p_bonferroni <- p.adjust(overlap_1_to_1_assoc[[14+i*9]],method = "bonferroni")
  overlap_1_to_1_assoc$soma_non_normal_sig_bonferroni <-F
  overlap_1_to_1_assoc$soma_non_normal_sig_bonferroni[overlap_1_to_1_assoc$soma_non_normal_p_bonferroni<0.05] <- T
  
  names(overlap_1_to_1_assoc)[(ncol(overlap_1_to_1_assoc)-5):ncol(overlap_1_to_1_assoc)] <-
    paste(names(overlap_1_to_1_assoc)[(ncol(overlap_1_to_1_assoc)-5):ncol(overlap_1_to_1_assoc)],variable[i],sep="_")
}

# ## count uncorrected
# 
# n_hit <- data.frame(colSums(overlap_1_to_1_assoc[grep("sig_raw",names(overlap_1_to_1_assoc))]))
# n_hit$variable <- row.names(n_hit)
# names(n_hit)[1] <- "n"
# 
# put them into a table

label <- c("Ambient temperature", "Hours since last meal", "Age", "Female", "Urban area",
           "Smoking", "Alcohol drinking", "Married", "≥ 6 years education", "BMI", "DBP", "SBP", "Heart rate",
           "Random glucose", "Physical activity", "Poor self-rated health", "Diabetes", "Kidney disease", "Cancer")

level <- c("Age", "Female", "Urban area", "Married", "≥ 6 years education",
           "Ambient temperature", "Hours since last meal",
           "BMI", "DBP", "SBP", "Heart rate", "Random glucose",
           "Poor self-rated health", "Diabetes", "Kidney disease", "Cancer",
           "Physical activity", "Smoking", "Alcohol drinking")

var_demo <- c("Age","Female","Urban area","Married", "≥ 6 years education")
var_phy <- c("Ambient temperature","Hours since last meal", "BMI", "DBP", "SBP", "Heart rate", "Random glucose")
var_med <- c("Poor self-rated health", "Diabetes", "Kidney disease", "Cancer")
var_life <- c("Physical activity","Smoking","Alcohol drinking")


## count bonferroni

names(overlap_1_to_1_assoc)
n_hit_bonferroni <- data.frame(colSums(overlap_1_to_1_assoc[grep("sig_bonferroni",names(overlap_1_to_1_assoc))]))
n_hit_bonferroni$variable <- row.names(n_hit_bonferroni)
names(n_hit_bonferroni)[1] <- "n"

for (i in 1:length(variable)) {
  n_hit_bonferroni$label[grep(variable[i],n_hit_bonferroni$variable)] <- label[i]
}

n_hit_bonferroni$label <- factor(n_hit_bonferroni$label, levels = level)

n_hit_bonferroni$platform <- NA
n_hit_bonferroni$platform[grep("olink",n_hit_bonferroni$variable)] <- "OLINK"
n_hit_bonferroni$platform[grep("soma_normal",n_hit_bonferroni$variable)] <- "SomaScan-ANML"
n_hit_bonferroni$platform[grep("soma_non_normal",n_hit_bonferroni$variable)] <- "SomaScan-non-ANML"

n_hit_bonferroni$platform <- factor(n_hit_bonferroni$platform, c("OLINK","SomaScan-ANML","SomaScan-non-ANML"))

n_hit_bonferroni$cat <- NA
n_hit_bonferroni$cat[n_hit_bonferroni$label %in% var_demo] <- "Socio-\ndemographics"
n_hit_bonferroni$cat[n_hit_bonferroni$label %in% var_phy] <- "Clinical\nmeasurements"
n_hit_bonferroni$cat[n_hit_bonferroni$label %in% var_med] <- "Medical\nhistory"
n_hit_bonferroni$cat[n_hit_bonferroni$label %in% var_life ] <- "Lifestyle\nfactors"

n_hit_bonferroni$cat <- factor(n_hit_bonferroni$cat, levels=c("Socio-\ndemographics","Clinical\nmeasurements","Medical\nhistory","Lifestyle\nfactors"))

hex <- rev(hue_pal()(3))

plot_n_hit_bonferroni <- ggplot(n_hit_bonferroni, aes(y=n, x=fct_rev(label), fill=fct_rev(platform))) + 
  geom_bar(stat="identity",position="dodge",colour="black") +
  ylim (0, 1694) +
  xlab("") +
  ylab("Number of significant associations") +
  labs(fill = "Platform") +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(strip.placement = "outside") +
  facet_grid(cat~.,scales = "free_y",space = "free",switch = "y") +
  scale_fill_manual(values=hex) +
  theme(legend.position = c(.95, .05), legend.justification = c("right", "bottom"))

plot_n_hit_bonferroni

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_n_hit_bonferroni_subcohort.png",plot_n_hit_bonferroni,width=9,height=9)









## whether hits are shared for bonferroni

trans_normal_bonferroni <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(variable)) {
  df <- overlap_1_to_1_assoc[grep(variable[i],names(overlap_1_to_1_assoc))]
  
  df$concord[df[13]<0.05 & df[15]<0.05 & sign(df[1])==sign(df[4])] <- "Shared"
  df$concord[df[13]<0.05 & df[15]>=0.05] <- "Olink-specific"
  df$concord[df[13]>=0.05 & df[15]<0.05] <- "SomaScan-specific"
  df$concord[df[13]<0.05 & df[15]<0.05 & sign(df[1])!=sign(df[4])] <- "Opposite direction"
  df$concord[df[13]>=0.05 & df[15]>=0.05] <- "Neither"
  
  df$concord <- factor(df$concord, levels=c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
  # table(df$concord)
  
  df_sum <- data.frame(table(df$concord))
  names(df_sum) <- c("concord","freq")
  df_sum$variable <- variable[i]
  df_sum$label <- label[i]
  df_sum$soma <- "OLINK vs SomaScan-ANML"
  
  trans_normal_bonferroni <- rbind(trans_normal_bonferroni,df_sum)
}

trans_non_normal_bonferroni <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(variable)) {
  df <- overlap_1_to_1_assoc[grep(variable[i],names(overlap_1_to_1_assoc))]
  
  df$concord[df[13]<0.05 & df[17]<0.05 & sign(df[1])==sign(df[7])] <- "Shared"
  df$concord[df[13]<0.05 & df[17]>=0.05] <- "Olink-specific"
  df$concord[df[13]>=0.05 & df[17]<0.05] <- "SomaScan-specific"
  df$concord[df[13]<0.05 & df[17]<0.05 & sign(df[1])!=sign(df[7])] <- "Opposite direction"
  df$concord[df[13]>=0.05 & df[17]>=0.05] <- "Neither"
  
  df$concord <- factor(df$concord, levels=c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
  # table(df$concord)
  
  df_sum <- data.frame(table(df$concord))
  names(df_sum) <- c("concord","freq")
  df_sum$variable <- variable[i]
  df_sum$label <- label[i]
  df_sum$soma <- "OLINK vs SomaScan-non-ANML"
  
  trans_non_normal_bonferroni <- rbind(trans_non_normal_bonferroni,df_sum)
}

trans_bonferroni <- rbind(trans_normal_bonferroni,trans_non_normal_bonferroni)

# plot

trans_bonferroni$label <- factor(trans_bonferroni$label, levels = level)

trans_bonferroni$concord <- factor(trans_bonferroni$concord, levels=c("Olink-specific","Shared","SomaScan-specific","Opposite direction","Neither"))

trans_bonferroni$percentage <- trans_bonferroni$freq/1694*100

trans_bonferroni$cat <- NA
trans_bonferroni$cat[trans_bonferroni$label %in% var_demo] <- "Socio-\ndemographics"
trans_bonferroni$cat[trans_bonferroni$label %in% var_phy] <- "Clinical\nmeasurements"
trans_bonferroni$cat[trans_bonferroni$label %in% var_med] <- "Medical\nhistory"
trans_bonferroni$cat[trans_bonferroni$label %in% var_life ] <- "Lifestyle\nfactors"

trans_bonferroni$cat <- factor(trans_bonferroni$cat, levels=c("Sample-\nrelated","Socio-\ndemographics","Clinical\nmeasurements","Medical\nhistory","Lifestyle\nfactors"))

hex <- rev(hue_pal()(9))

plot_hit_concordance_bonferroni <- ggplot(trans_bonferroni, aes(y=percentage, x=fct_rev(label), fill=fct_rev(concord))) + 
  geom_bar(stat="identity",position="stack",colour="black") +
  facet_grid(~ soma) +
  xlab("") +
  ylab("Percentage (%)") +
  labs(fill = "Concordance") +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(text = element_text(size = 10)) +
  theme(strip.placement = "outside") +
  facet_grid(cat~soma,scales = "free_y",space = "free",switch = "y") +
  scale_fill_manual(values=c("gray80",hex[c(1,5,7,9)]))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_hit_concordance_bonferroni_subcohort.png",plot_hit_concordance_bonferroni,width=9,height=6)


## whether hits are shared between ANML and non-ANML (bonferroni)

trans_soma_bonferroni <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(variable)) {
  df <- overlap_1_to_1_assoc[grep(variable[i],names(overlap_1_to_1_assoc))]
  
  df$concord[df[15]<0.05 & df[17]<0.05 & sign(df[4])==sign(df[7])] <- "Shared"
  df$concord[df[15]<0.05 & df[17]>=0.05] <- "ANML-specific"
  df$concord[df[15]>=0.05 & df[17]<0.05] <- "Non-ANML-specific"
  df$concord[df[15]<0.05 & df[17]<0.05 & sign(df[4])!=sign(df[7])] <- "Opposite direction"
  df$concord[df[15]>=0.05 & df[17]>=0.05] <- "Neither"
  
  df$concord <- factor(df$concord, levels=c("Shared","Opposite direction","ANML-specific","Non-ANML-specific","Neither"))
  # table(df$concord)
  
  df_sum <- data.frame(table(df$concord))
  names(df_sum) <- c("concord","freq")
  df_sum$variable <- variable[i]
  df_sum$label <- label[i]
  
  trans_soma_bonferroni <- rbind(trans_soma_bonferroni,df_sum)
}

# plot

trans_soma_bonferroni$label <- factor(trans_soma_bonferroni$label, levels = level)

trans_soma_bonferroni$concord <- factor(trans_soma_bonferroni$concord, levels=c("ANML-specific","Shared","Non-ANML-specific","Opposite direction","Neither"))

trans_soma_bonferroni$percentage <- trans_soma_bonferroni$freq/1694*100

trans_soma_bonferroni$cat <- NA
trans_soma_bonferroni$cat[trans_soma_bonferroni$label %in% var_demo] <- "Socio-\ndemographics"
trans_soma_bonferroni$cat[trans_soma_bonferroni$label %in% var_phy] <- "Clinical\nmeasurements"
trans_soma_bonferroni$cat[trans_soma_bonferroni$label %in% var_med] <- "Medical\nhistory"
trans_soma_bonferroni$cat[trans_soma_bonferroni$label %in% var_life ] <- "Lifestyle\nfactors"

trans_soma_bonferroni$cat <- factor(trans_soma_bonferroni$cat, levels=c("Sample-\nrelated","Socio-\ndemographics","Clinical\nmeasurements","Medical\nhistory","Lifestyle\nfactors"))

plot_hit_concordance_soma_bonferroni <- ggplot(trans_soma_bonferroni, aes(y=percentage, x=fct_rev(label), fill=fct_rev(concord))) + 
  geom_bar(stat="identity",position="stack",colour="black") +
  xlab("Variable") +
  ylab("Percentage (%)") +
  ggtitle("SomaScan ANML vs non-ANML") +
  labs(fill = "Concordance") +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(strip.placement = "outside") +
  facet_grid(cat~.,scales = "free_y",space = "free",switch = "y") +
  scale_fill_manual(values=c("gray80",hex[c(1,5,7,9)]))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_hit_concordance_soma_bonferroni_subcohort.png",plot_hit_concordance_soma_bonferroni,width=6,height=6)








## correlation coefficient between beta

cor_es_normal <- data.frame(matrix(ncol = 5, nrow = 0))
cor_es_non_normal <- data.frame(matrix(ncol = 5, nrow = 0))
cor_es_soma <- data.frame(matrix(ncol = 5, nrow = 0))

for (i in 1:length(variable)) {
  
  df <- overlap_1_to_1_assoc[grep(paste0("es_",variable[i]),names(overlap_1_to_1_assoc))]
  
  df_1 <- data.frame(variable=variable[i],label=label[i],
                     r=cor(df[[1]],df[[2]]),
                     lci=cor.test(df[[1]],df[[2]])$conf.int[1],
                     uci=cor.test(df[[1]],df[[2]])$conf.int[2])
  
  cor_es_normal <- rbind(cor_es_normal,df_1)
  
  df_2 <- data.frame(variable=variable[i],label=label[i],
                     r=cor(df[[1]],df[[3]]),
                     lci=cor.test(df[[1]],df[[3]])$conf.int[1],
                     uci=cor.test(df[[1]],df[[3]])$conf.int[2])
  
  cor_es_non_normal <- rbind(cor_es_non_normal,df_2)
  
  df_3 <- data.frame(variable=variable[i],label=label[i],
                     r=cor(df[[2]],df[[3]]),
                     lci=cor.test(df[[2]],df[[3]])$conf.int[1],
                     uci=cor.test(df[[2]],df[[3]])$conf.int[2])
  
  cor_es_soma <- rbind(cor_es_soma,df_3)
}

cor_es_normal$label <- factor(cor_es_normal$label, levels = level)
shared_normal <- trans_normal_bonferroni[trans_normal_bonferroni$concord=="Shared",c("label","concord","freq")]
cor_es_normal <- merge(cor_es_normal,shared_normal,by="label")

cor_es_non_normal$label <- factor(cor_es_non_normal$label, levels = level)
shared_non_normal <- trans_non_normal_bonferroni[trans_non_normal_bonferroni$concord=="Shared",c("label","concord","freq")]
cor_es_non_normal <- merge(cor_es_non_normal,shared_non_normal,by="label")

cor_es_soma$label <- factor(cor_es_soma$label, levels = level)
shared_soma <- trans_soma_bonferroni[trans_soma_bonferroni$concord=="Shared",c("label","concord","freq")]
cor_es_soma <- merge(cor_es_soma,shared_soma,by="label")

cor_es <- list(cor_es_normal,cor_es_non_normal)

row_labels <- data.frame(
  subgroup = level,
  group    = c("Socio-demographics","Socio-demographics","Socio-demographics","Socio-demographics","Socio-demographics",
               "Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements",
               "Medical history","Medical history","Medical history","Medical history",
               "Lifestyle factors","Lifestyle factors","Lifestyle factors"),
  label = level)

plot_cor_es <- forest_plot(panels = cor_es,
                           exponentiate = F,
                           panel.headings = c("OLINK vs SomaScan-ANML","OLINK vs SomaScan-non-ANML"),
                           row.labels = row_labels,
                           row.labels.levels = c("group","subgroup"),
                           rows = unique(row_labels$group),
                           col.key = "label",
                           col.estimate = "r",
                           col.lci = "lci",
                           col.uci = "uci",
                           col.left = "freq",
                           col.left.heading = "No. shared hits",
                           xlab = "Pearson's r (95% CI)",
                           xlim = c(0.2,0.7),
                           nullval = 0.5,
                           # estcolumn = F,
                           pointsize = 2,
                           base_size = 12,
                           col.right.heading = "Pearson's r (95% CI)"
)


# add dashed line

# Prepare data to be plotted using ckbplotr::forest_data()
datatoplot <- ckbplotr::forest_data(row.labels = row_labels,
                                    row.labels.levels = c("group", "subgroup"),
                                    rows = c("Socio-demographics", "Clinical measurements", "Medical history", "Lifestyle factors"),
                                    panels = cor_es,
                                    panel.names = c("1", "2"),
                                    col.key = "label",
                                    col.estimate = "r",
                                    col.lci = "lci",
                                    col.uci = "uci",
                                    col.left = "freq",
                                    exponentiate = FALSE)

# Get a character vector of the row labels, so these can be used in the plot
rowlabels <- datatoplot %>%
  dplyr::group_by(row) %>%
  dplyr::summarise(row.label = dplyr::first(row.label),
                   bold = all(is.na(estimate_transformed) | all(key %in% character(0))),
                   .groups = "drop") %>%
  dplyr::mutate(row.label = dplyr::if_else(bold & row.label != "",
                                           paste0("**", row.label, "**"),
                                           as.character(row.label))) %>% 
  dplyr::arrange(row) %>%
  dplyr::pull(row.label)

# Identify any CIs that extend outside axis limits
datatoplot <- datatoplot %>%
  dplyr::mutate(cioverright  = (uci_transformed > 0.7),
                uci_transformed = pmin(uci_transformed, 0.7),
                lci_transformed = pmin(lci_transformed, 0.7),
                cioverleft  = (lci_transformed < 0.2),
                lci_transformed = pmax(lci_transformed, 0.2),
                uci_transformed = pmax(uci_transformed, 0.2))

# Create the ggplot
plot_cor_es <- ggplot(datatoplot, aes(y=-row, x=estimate_transformed)) +
  
  # Put the different panels in side-by-side plots using facets
  facet_wrap(~panel, nrow = 1) +
  
  # Add a line at null effect
  annotate(geom = "segment",
           y = -0.7, yend = -Inf,
           x = 0.5, xend = 0.5,
           linewidth = 0.545454545454545,
           linetype = 3,
           colour = "black") +
  
  # Plot points at the transformed estimates
  ## Scale by inverse of the SE
  geom_point(aes(size = size),
             data = ~ dplyr::filter(.x, estimate_transformed > 0.2, estimate_transformed < 0.7),
             shape = 15,
             colour = "black",
             fill = "black",
             stroke = 0,
             na.rm = TRUE) +
  
  # Scale the size of points by their side length
  # and make the scale range from zero upwards
  scale_radius(limits = c(0, 1),
               range = c(0, 2)) +
  
  # Plot the CIs
  geom_errorbar(aes(xmin = lci_transformed,
                    xmax = uci_transformed),
                data = ~ dplyr::filter(.x, !is.na(estimate_transformed)),
                colour = "black",
                width = 0,
                linewidth = 0.545454545454545,
                na.rm = TRUE) +
  
  # Add tiny segments with arrows when the CIs go outside axis limits
  geom_segment(aes(y = -row,
                   yend = -row,
                   x = uci_transformed-0.000001,
                   xend = uci_transformed),
               data = ~ dplyr::filter(.x, cioverright == TRUE),
               colour = "black",
               linewidth = 0.545454545454545,
               arrow = arrow(type = "closed", length = unit(4.36363636363636, "pt")),
               na.rm = TRUE) +
  geom_segment(aes(y = -row,
                   yend = -row,
                   x = lci_transformed+0.000001,
                   xend = lci_transformed),
               data = ~ dplyr::filter(.x, cioverleft == TRUE),
               colour = "black",
               linewidth = 0.545454545454545,
               arrow = arrow(type = "closed", length = unit(4.36363636363636, "pt")),
               na.rm = TRUE) +
  
  # Use identity for aesthetic scales
  scale_shape_identity() +
  scale_fill_identity() +
  scale_colour_identity() +
  
  # Set coordinate system
  coord_cartesian(clip = "off",
                  xlim = c(0.2, 0.7)) +
  
  # Add columns to right side of plots
  ## column auto_estcolumn
  ckbplotr::geom_text_move(aes(y = -row,
                               x = 0.7,
                               label = `auto_estcolumn`),
                           move_x = unit(0.6, "mm"),
                           hjust = 0,
                           size  = 3.27272727272727,
                           colour  = "black",
                           na.rm = TRUE,
                           parse = FALSE) +
  ckbplotr::geom_text_move(aes(y = 0,
                               x = 0.7,
                               label = title),
                           move_x = unit(0.6, "mm"),
                           hjust    = 0,
                           size     = 3.27272727272727,
                           colour  = "black",
                           fontface = "bold",
                           data = dplyr::tibble(panel = factor(c("1", "2"),
                                                               levels = c("1", "2"),
                                                               ordered = TRUE),
                                                title = "Pearson\'s r (95% CI)")) +
  
  # Add columns to left side of plots
  ## column freq
  ckbplotr::geom_text_move(aes(y = -row,
                               x = 0.2,
                               label = freq,
                               fontface = "plain"),
                           move_x = unit(-0.6, "mm"),
                           hjust = 1,
                           size  = 3.27272727272727,
                           colour  = "black",
                           na.rm = TRUE) +
  ckbplotr::geom_text_move(aes(y = 0,
                               x = 0.2,
                               label = title),
                           move_x = unit(-0.6, "mm"),
                           hjust    = 1,
                           size     = 3.27272727272727,
                           colour  = "black",
                           fontface = "bold",
                           data = dplyr::tibble(panel = factor(c("1", "2"),
                                                               levels = c("1", "2"),
                                                               ordered = TRUE),
                                                title = "No. shared hits")) +
  
  # Add xlab below each axis
  geom_text(aes(y = -Inf, x = 0.45, label = xlab),
            hjust = 0.5,
            size  = 3.27272727272727,
            colour  = "black",
            vjust = 4.4,
            fontface = "bold",
            data = dplyr::tibble(panel = factor(c("1", "2"),
                                                levels = c("1", "2"),
                                                ordered = TRUE),
                                 xlab = "Pearson\'s r (95% CI)")) +
  
  # Add panel name above each panel
  geom_text(aes(y = 0, x = 0.45, label = title),
            hjust = 0.5,
            nudge_y = 2,
            size  = 3.27272727272727,
            colour  = "black",
            fontface = "bold",
            data = dplyr::tibble(panel = factor(c("1", "2"),
                                                levels = c("1", "2"),
                                                ordered = TRUE),
                                 title = c("OLINK vs SomaScan-ANML", "OLINK vs SomaScan-non-ANML"))) +
  
  # Set the scale for the x axis (the estimates and CIs)
  scale_x_continuous(trans  = "identity",
                     breaks = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
                     expand = c(0,0)) +
  
  # Set the scale for the y axis (the rows)
  scale_y_continuous(breaks = -1:-max(datatoplot$row),
                     labels = rowlabels,
                     limits = c(-max(datatoplot$row) - 0.7, NA),
                     expand = c(0,0)) +
  
  # Control the overall look of the plots
  theme(text             = element_text(size = 12, colour = "black"),
        line             = element_line(linewidth = 0.545454545454545),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title       = element_blank(),
        axis.line.x      = element_line(colour = "black", linewidth = 0.545454545454545, lineend = "round"),
        axis.title       = element_blank(),
        axis.ticks.x     = element_line(colour = "black"),
        axis.text.x      = element_text(colour = "black",,
                                        margin = margin(t = 4.8),
                                        vjust  = 1),
        axis.ticks.y     = element_blank(),
        axis.line.y      = element_blank(),
        axis.text.y      = ggtext::element_markdown(hjust  = 0,
                                                    colour = "black",
                                                    margin = margin(r = 25.8, unit = "mm")),
        panel.border     = element_blank(),
        panel.spacing    = unit(31.1, "mm") + unit(5, "mm") + unit(25.8, "mm"),
        strip.background = element_blank(),
        strip.placement  = "outside",
        strip.text       = element_blank(),
        legend.position  = "none",
        plot.background  = element_blank(),
        plot.margin      = margin(8, 8, 8, 8, "mm") + unit(c(0, 31.1, 0, 0), "mm"))

plot_cor_es_bg <- plot_cor_es + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_cor_es_bonferroni_subcohort.png",plot_cor_es_bg,width=14,height=6)

plot_cor_es_soma <- forest_plot(panels = list(cor_es_soma),
                                exponentiate = F,
                                panel.headings = c("SomaScan ANML vs non-ANML"),
                                row.labels = row_labels,
                                row.labels.levels = c("group","subgroup"),
                                rows = unique(row_labels$group),
                                col.key = "label",
                                col.estimate = "r",
                                col.lci = "lci",
                                col.uci = "uci",
                                col.left = "freq",
                                col.left.heading = "No. shared hits",
                                xlab = "Pearson's r (95% CI)",
                                # xlim = c(0.2,0.7),
                                # estcolumn = F,
                                pointsize = 2,
                                base_size = 12,
                                col.right.heading = "Pearson's r (95% CI)"
)

plot_cor_es_soma_bg <- plot_cor_es_soma$plot + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_cor_es_soma_bonferroni_subcohort.png",plot_cor_es_soma_bg,width=7,height=6)



## correlation coefficient between beta, only for shared associations

cor_es_normal_sig <- data.frame(matrix(ncol = 5, nrow = 0))
cor_es_non_normal_sig <- data.frame(matrix(ncol = 5, nrow = 0))
cor_es_soma_sig <- data.frame(matrix(ncol = 5, nrow = 0))

# remove variables with two few shared hits

variable_sig <- variable[! variable %in% c("married","school","poor_health","kidney_dis_diag","cancer_diag","met")]
label_sig <- label[! label %in% c("Married","≥ 6 years education","Poor self-rated health","Kidney disease","Cancer","Physical activity")]
level_sig <- level[! level %in% c("Married","≥ 6 years education","Poor self-rated health","Kidney disease","Cancer","Physical activity")]

for (i in 1:length(variable_sig)) {
  
  df <- overlap_1_to_1_assoc[c(grep(paste0("es_",variable_sig[i]),names(overlap_1_to_1_assoc)),
                               grep(paste0("p_bonferroni_",variable_sig[i]),names(overlap_1_to_1_assoc)))]
  
  df_1 <- df[df[4]<0.05 & df[5]<0.05 & sign(df[1])==sign(df[2]),]
  
  df_1_result <- data.frame(variable=variable_sig[i],label=label_sig[i],
                            r=cor(df_1[[1]],df_1[[2]]),
                            lci=cor.test(df_1[[1]],df_1[[2]])$conf.int[1],
                            uci=cor.test(df_1[[1]],df_1[[2]])$conf.int[2])
  
  cor_es_normal_sig <- rbind(cor_es_normal_sig,df_1_result)
  
  df_2 <- df[df[4]<0.05 & df[6]<0.05 & sign(df[1])==sign(df[3]),]
  
  df_2_result <- data.frame(variable=variable_sig[i],label=label_sig[i],
                            r=cor(df_2[[1]],df_2[[3]]),
                            lci=cor.test(df_2[[1]],df_2[[3]])$conf.int[1],
                            uci=cor.test(df_2[[1]],df_2[[3]])$conf.int[2])
  
  cor_es_non_normal_sig <- rbind(cor_es_non_normal_sig,df_2_result)
  
  df_3 <- df[df[5]<0.05 & df[6]<0.05 & sign(df[2])==sign(df[3]),]
  
  df_3_result <- data.frame(variable=variable_sig[i],label=label_sig[i],
                            r=cor(df_3[[2]],df_3[[3]]),
                            lci=cor.test(df_3[[2]],df_3[[3]])$conf.int[1],
                            uci=cor.test(df_3[[2]],df_3[[3]])$conf.int[2])
  
  cor_es_soma_sig <- rbind(cor_es_soma_sig,df_3_result)
}

cor_es_normal_sig$label <- factor(cor_es_normal_sig$label, levels = level_sig)
shared_normal_sig <- trans_normal_bonferroni[trans_normal_bonferroni$concord=="Shared",c("label","concord","freq")]
cor_es_normal_sig <- merge(cor_es_normal_sig,shared_normal,by="label")

cor_es_non_normal_sig$label <- factor(cor_es_non_normal_sig$label, levels = level_sig)
shared_non_normal_sig <- trans_non_normal_bonferroni[trans_non_normal_bonferroni$concord=="Shared",c("label","concord","freq")]
cor_es_non_normal_sig <- merge(cor_es_non_normal_sig,shared_non_normal,by="label")

cor_es_soma_sig$label <- factor(cor_es_soma_sig$label, levels = level_sig)
shared_soma_sig <- trans_soma_bonferroni[trans_soma_bonferroni$concord=="Shared",c("label","concord","freq")]
cor_es_soma_sig <- merge(cor_es_soma_sig,shared_soma,by="label")

cor_es_sig <- list(cor_es_normal_sig,cor_es_non_normal_sig)

row_labels_sig <- data.frame(
  subgroup = level_sig,
  group    = c("Socio-demographics","Socio-demographics","Socio-demographics",
               "Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements",
               "Medical history",
               "Lifestyle factors","Lifestyle factors"),
  label = level_sig)

plot_cor_es_sig <- forest_plot(panels = cor_es_sig,
                               exponentiate = F,
                               panel.headings = c("OLINK vs SomaScan-ANML","OLINK vs SomaScan-non-ANML"),
                               row.labels = row_labels_sig,
                               row.labels.levels = c("group","subgroup"),
                               rows = unique(row_labels$group),
                               col.key = "label",
                               col.estimate = "r",
                               col.lci = "lci",
                               col.uci = "uci",
                               col.left = "freq",
                               col.left.heading = "No. shared hits",
                               xlab = "Pearson's r (95% CI)",
                               # xlim = c(0.2,0.7),
                               # estcolumn = F,
                               pointsize = 2,
                               base_size = 12,
                               col.right.heading = "Pearson's r (95% CI)"
)

plot_cor_es_sig_bg <- plot_cor_es_sig$plot + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_cor_es_sig_bonferroni_subcohort.png",plot_cor_es_sig_bg,width=14,height=6)

plot_cor_es_soma_sig <- forest_plot(panels = list(cor_es_soma_sig),
                                    exponentiate = F,
                                    panel.headings = c("SomaScan ANML vs non-ANML"),
                                    row.labels = row_labels_sig,
                                    row.labels.levels = c("group","subgroup"),
                                    rows = unique(row_labels$group),
                                    col.key = "label",
                                    col.estimate = "r",
                                    col.lci = "lci",
                                    col.uci = "uci",
                                    col.left = "freq",
                                    col.left.heading = "No. shared hits",
                                    xlab = "Pearson's r (95% CI)",
                                    # xlim = c(0.2,0.7),
                                    # estcolumn = F,
                                    pointsize = 2,
                                    base_size = 12,
                                    col.right.heading = "Pearson's r (95% CI)"
)

plot_cor_es_soma_sig_bg <- plot_cor_es_soma_sig$plot + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_cor_es_soma_sig_bonferroni_subcohort.png",plot_cor_es_soma_sig_bg,width=7,height=6)


# write.csv(overlap_1_to_1_assoc,"overlap_1_to_1_assoc_sig.csv", quote=F, row.names=F)









######################## fdr


names(overlap_1_to_1_assoc)

for (i in 1:length(variable)) {
  overlap_1_to_1_assoc$olink_p_fdr <- F
  overlap_1_to_1_assoc$olink_p_fdr <- p.adjust(overlap_1_to_1_assoc[[8+i*9]],method = "fdr")
  overlap_1_to_1_assoc$olink_sig_fdr <-F
  overlap_1_to_1_assoc$olink_sig_fdr[overlap_1_to_1_assoc$olink_p_fdr<0.05] <- T
  
  overlap_1_to_1_assoc$soma_normal_p_fdr <- F
  overlap_1_to_1_assoc$soma_normal_p_fdr <- p.adjust(overlap_1_to_1_assoc[[11+i*9]],method = "fdr")
  overlap_1_to_1_assoc$soma_normal_sig_fdr <-F
  overlap_1_to_1_assoc$soma_normal_sig_fdr[overlap_1_to_1_assoc$soma_normal_p_fdr<0.05] <- T
  
  overlap_1_to_1_assoc$soma_non_normal_p_fdr <- F
  overlap_1_to_1_assoc$soma_non_normal_p_fdr <- p.adjust(overlap_1_to_1_assoc[[14+i*9]],method = "fdr")
  overlap_1_to_1_assoc$soma_non_normal_sig_fdr <-F
  overlap_1_to_1_assoc$soma_non_normal_sig_fdr[overlap_1_to_1_assoc$soma_non_normal_p_fdr<0.05] <- T
  
  names(overlap_1_to_1_assoc)[(ncol(overlap_1_to_1_assoc)-5):ncol(overlap_1_to_1_assoc)] <-
    paste(names(overlap_1_to_1_assoc)[(ncol(overlap_1_to_1_assoc)-5):ncol(overlap_1_to_1_assoc)],variable[i],sep="_")
}


## count fdr

names(overlap_1_to_1_assoc)
n_hit_fdr <- data.frame(colSums(overlap_1_to_1_assoc[grep("sig_fdr",names(overlap_1_to_1_assoc))]))
n_hit_fdr$variable <- row.names(n_hit_fdr)
names(n_hit_fdr)[1] <- "n"

for (i in 1:length(variable)) {
  n_hit_fdr$label[grep(variable[i],n_hit_fdr$variable)] <- label[i]
}

n_hit_fdr$label <- factor(n_hit_fdr$label, levels = level)

n_hit_fdr$platform <- NA
n_hit_fdr$platform[grep("olink",n_hit_fdr$variable)] <- "OLINK"
n_hit_fdr$platform[grep("soma_normal",n_hit_fdr$variable)] <- "SomaScan-ANML"
n_hit_fdr$platform[grep("soma_non_normal",n_hit_fdr$variable)] <- "SomaScan-non-ANML"

n_hit_fdr$platform <- factor(n_hit_fdr$platform, c("OLINK","SomaScan-ANML","SomaScan-non-ANML"))

n_hit_fdr$cat <- NA
n_hit_fdr$cat[n_hit_fdr$label %in% var_demo] <- "Socio-\ndemographics"
n_hit_fdr$cat[n_hit_fdr$label %in% var_phy] <- "Clinical\nmeasurements"
n_hit_fdr$cat[n_hit_fdr$label %in% var_med] <- "Medical\nhistory"
n_hit_fdr$cat[n_hit_fdr$label %in% var_life ] <- "Lifestyle\nfactors"

n_hit_fdr$cat <- factor(n_hit_fdr$cat, levels=c("Sample-\nrelated","Socio-\ndemographics","Clinical\nmeasurements","Medical\nhistory","Lifestyle\nfactors"))

hex <- rev(hue_pal()(3))

plot_n_hit_fdr <- ggplot(n_hit_fdr, aes(y=n, x=fct_rev(label), fill=fct_rev(platform))) + 
  geom_bar(stat="identity",position="dodge",colour="black") +
  ylim (0, 1694) +
  xlab("") +
  ylab("Number of significant associations") +
  labs(fill = "Platform") +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(strip.placement = "outside") +
  facet_grid(cat~.,scales = "free_y",space = "free",switch = "y") +
  scale_fill_manual(values=hex) +
  theme(legend.position = c(.95, .05), legend.justification = c("right", "bottom"))

plot_n_hit_fdr

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_n_hit_fdr_subcohort.png",plot_n_hit_fdr,width=9,height=9)









## whether hits are shared for fdr

trans_normal_fdr <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(variable)) {
  df <- overlap_1_to_1_assoc[grep(variable[i],names(overlap_1_to_1_assoc))]
  
  df$concord[df[19]<0.05 & df[21]<0.05 & sign(df[1])==sign(df[4])] <- "Shared"
  df$concord[df[19]<0.05 & df[21]>=0.05] <- "Olink-specific"
  df$concord[df[19]>=0.05 & df[21]<0.05] <- "SomaScan-specific"
  df$concord[df[19]<0.05 & df[21]<0.05 & sign(df[1])!=sign(df[4])] <- "Opposite direction"
  df$concord[df[19]>=0.05 & df[21]>=0.05] <- "Neither"
  
  df$concord <- factor(df$concord, levels=c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
  # table(df$concord)
  
  df_sum <- data.frame(table(df$concord))
  names(df_sum) <- c("concord","freq")
  df_sum$variable <- variable[i]
  df_sum$label <- label[i]
  df_sum$soma <- "OLINK vs SomaScan-ANML"
  
  trans_normal_fdr <- rbind(trans_normal_fdr,df_sum)
}

trans_non_normal_fdr <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(variable)) {
  df <- overlap_1_to_1_assoc[grep(variable[i],names(overlap_1_to_1_assoc))]
  
  df$concord[df[19]<0.05 & df[23]<0.05 & sign(df[1])==sign(df[7])] <- "Shared"
  df$concord[df[19]<0.05 & df[23]>=0.05] <- "Olink-specific"
  df$concord[df[19]>=0.05 & df[23]<0.05] <- "SomaScan-specific"
  df$concord[df[19]<0.05 & df[23]<0.05 & sign(df[1])!=sign(df[7])] <- "Opposite direction"
  df$concord[df[19]>=0.05 & df[23]>=0.05] <- "Neither"
  
  df$concord <- factor(df$concord, levels=c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
  # table(df$concord)
  
  df_sum <- data.frame(table(df$concord))
  names(df_sum) <- c("concord","freq")
  df_sum$variable <- variable[i]
  df_sum$label <- label[i]
  df_sum$soma <- "OLINK vs SomaScan-non-ANML"
  
  trans_non_normal_fdr <- rbind(trans_non_normal_fdr,df_sum)
}

trans_fdr <- rbind(trans_normal_fdr,trans_non_normal_fdr)

# plot

trans_fdr$label <- factor(trans_fdr$label, levels = level)

trans_fdr$concord <- factor(trans_fdr$concord, levels=c("Olink-specific","Shared","SomaScan-specific","Opposite direction","Neither"))

trans_fdr$percentage <- trans_fdr$freq/1694*100

trans_fdr$cat <- NA
trans_fdr$cat[trans_fdr$label %in% var_demo] <- "Socio-\ndemographics"
trans_fdr$cat[trans_fdr$label %in% var_phy] <- "Clinical\nmeasurements"
trans_fdr$cat[trans_fdr$label %in% var_med] <- "Medical\nhistory"
trans_fdr$cat[trans_fdr$label %in% var_life ] <- "Lifestyle\nfactors"

trans_fdr$cat <- factor(trans_fdr$cat, levels=c("Sample-\nrelated","Socio-\ndemographics","Clinical\nmeasurements","Medical\nhistory","Lifestyle\nfactors"))

hex <- rev(hue_pal()(9))

plot_hit_concordance_fdr <- ggplot(trans_fdr, aes(y=percentage, x=fct_rev(label), fill=fct_rev(concord))) + 
  geom_bar(stat="identity",position="stack",colour="black") +
  facet_grid(~ soma) +
  xlab("") +
  ylab("Percentage (%)") +
  labs(fill = "Concordance") +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(text = element_text(size = 10)) +
  theme(strip.placement = "outside") +
  facet_grid(cat~soma,scales = "free_y",space = "free",switch = "y") +
  scale_fill_manual(values=c("gray80",hex[c(1,5,7,9)]))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_hit_concordance_fdr_subcohort.png",plot_hit_concordance_fdr,width=9,height=6)


## whether hits are shared between ANML and non-ANML (fdr)

trans_soma_fdr <- data.frame(matrix(ncol = 4, nrow = 0))

for (i in 1:length(variable)) {
  df <- overlap_1_to_1_assoc[grep(variable[i],names(overlap_1_to_1_assoc))]
  
  df$concord[df[21]<0.05 & df[23]<0.05 & sign(df[4])==sign(df[7])] <- "Shared"
  df$concord[df[21]<0.05 & df[23]>=0.05] <- "ANML-specific"
  df$concord[df[21]>=0.05 & df[23]<0.05] <- "Non-ANML-specific"
  df$concord[df[21]<0.05 & df[23]<0.05 & sign(df[4])!=sign(df[7])] <- "Opposite direction"
  df$concord[df[21]>=0.05 & df[23]>=0.05] <- "Neither"
  
  df$concord <- factor(df$concord, levels=c("Shared","Opposite direction","ANML-specific","Non-ANML-specific","Neither"))
  # table(df$concord)
  
  df_sum <- data.frame(table(df$concord))
  names(df_sum) <- c("concord","freq")
  df_sum$variable <- variable[i]
  df_sum$label <- label[i]
  
  trans_soma_fdr <- rbind(trans_soma_fdr,df_sum)
}

# plot

trans_soma_fdr$label <- factor(trans_soma_fdr$label, levels = level)

trans_soma_fdr$concord <- factor(trans_soma_fdr$concord, levels=c("ANML-specific","Shared","Non-ANML-specific","Opposite direction","Neither"))

trans_soma_fdr$percentage <- trans_soma_fdr$freq/1694*100

trans_soma_fdr$cat <- NA
trans_soma_fdr$cat[trans_soma_fdr$label %in% var_demo] <- "Socio-\ndemographics"
trans_soma_fdr$cat[trans_soma_fdr$label %in% var_phy] <- "Clinical\nmeasurements"
trans_soma_fdr$cat[trans_soma_fdr$label %in% var_med] <- "Medical\nhistory"
trans_soma_fdr$cat[trans_soma_fdr$label %in% var_life ] <- "Lifestyle\nfactors"

trans_soma_fdr$cat <- factor(trans_soma_fdr$cat, levels=c("Sample-\nrelated","Socio-\ndemographics","Clinical\nmeasurements","Medical\nhistory","Lifestyle\nfactors"))

plot_hit_concordance_soma_fdr <- ggplot(trans_soma_fdr, aes(y=percentage, x=fct_rev(label), fill=fct_rev(concord))) + 
  geom_bar(stat="identity",position="stack",colour="black") +
  xlab("Variable") +
  ylab("Percentage (%)") +
  ggtitle("SomaScan ANML vs non-ANML") +
  labs(fill = "Concordance") +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip() +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 10)) +
  theme(strip.placement = "outside") +
  facet_grid(cat~.,scales = "free_y",space = "free",switch = "y") +
  scale_fill_manual(values=c("gray80",hex[c(1,5,7,9)]))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_hit_concordance_soma_fdr_subcohort.png",plot_hit_concordance_soma_fdr,width=6,height=6)








## correlation coefficient between beta

cor_es_normal <- data.frame(matrix(ncol = 5, nrow = 0))
cor_es_non_normal <- data.frame(matrix(ncol = 5, nrow = 0))
cor_es_soma <- data.frame(matrix(ncol = 5, nrow = 0))

for (i in 1:length(variable)) {
  
  df <- overlap_1_to_1_assoc[grep(paste0("es_",variable[i]),names(overlap_1_to_1_assoc))]
  
  df_1 <- data.frame(variable=variable[i],label=label[i],
                     r=cor(df[[1]],df[[2]]),
                     lci=cor.test(df[[1]],df[[2]])$conf.int[1],
                     uci=cor.test(df[[1]],df[[2]])$conf.int[2])
  
  cor_es_normal <- rbind(cor_es_normal,df_1)
  
  df_2 <- data.frame(variable=variable[i],label=label[i],
                     r=cor(df[[1]],df[[3]]),
                     lci=cor.test(df[[1]],df[[3]])$conf.int[1],
                     uci=cor.test(df[[1]],df[[3]])$conf.int[2])
  
  cor_es_non_normal <- rbind(cor_es_non_normal,df_2)
  
  df_3 <- data.frame(variable=variable[i],label=label[i],
                     r=cor(df[[2]],df[[3]]),
                     lci=cor.test(df[[2]],df[[3]])$conf.int[1],
                     uci=cor.test(df[[2]],df[[3]])$conf.int[2])
  
  cor_es_soma <- rbind(cor_es_soma,df_3)
}

cor_es_normal$label <- factor(cor_es_normal$label, levels = level)
shared_normal <- trans_normal_fdr[trans_normal_fdr$concord=="Shared",c("label","concord","freq")]
cor_es_normal <- merge(cor_es_normal,shared_normal,by="label")

cor_es_non_normal$label <- factor(cor_es_non_normal$label, levels = level)
shared_non_normal <- trans_non_normal_fdr[trans_non_normal_fdr$concord=="Shared",c("label","concord","freq")]
cor_es_non_normal <- merge(cor_es_non_normal,shared_non_normal,by="label")

cor_es_soma$label <- factor(cor_es_soma$label, levels = level)
shared_soma <- trans_soma_fdr[trans_soma_fdr$concord=="Shared",c("label","concord","freq")]
cor_es_soma <- merge(cor_es_soma,shared_soma,by="label")

cor_es <- list(cor_es_normal,cor_es_non_normal)

row_labels <- data.frame(
  subgroup = level,
  group    = c("Socio-demographics","Socio-demographics","Socio-demographics","Socio-demographics","Socio-demographics",
               "Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements",
               "Medical history","Medical history","Medical history","Medical history",
               "Lifestyle factors","Lifestyle factors","Lifestyle factors"),
  label = level)

plot_cor_es <- forest_plot(panels = cor_es,
                           exponentiate = F,
                           panel.headings = c("OLINK vs SomaScan-ANML","OLINK vs SomaScan-non-ANML"),
                           row.labels = row_labels,
                           row.labels.levels = c("group","subgroup"),
                           rows = unique(row_labels$group),
                           col.key = "label",
                           col.estimate = "r",
                           col.lci = "lci",
                           col.uci = "uci",
                           col.left = "freq",
                           col.left.heading = "No. shared hits",
                           xlab = "Pearson's r (95% CI)",
                           xlim = c(0.2,0.7),
                           nullval = 0.5,
                           # estcolumn = F,
                           pointsize = 2,
                           base_size = 12,
                           col.right.heading = "Pearson's r (95% CI)"
)

# add dashed line

# Prepare data to be plotted using ckbplotr::forest_data()
datatoplot <- ckbplotr::forest_data(row.labels = row_labels,
                                    row.labels.levels = c("group", "subgroup"),
                                    rows = c("Socio-demographics", "Clinical measurements", "Medical history", "Lifestyle factors"),
                                    panels = cor_es,
                                    panel.names = c("1", "2"),
                                    col.key = "label",
                                    col.estimate = "r",
                                    col.lci = "lci",
                                    col.uci = "uci",
                                    col.left = "freq",
                                    exponentiate = FALSE)

# Get a character vector of the row labels, so these can be used in the plot
rowlabels <- datatoplot %>%
  dplyr::group_by(row) %>%
  dplyr::summarise(row.label = dplyr::first(row.label),
                   bold = all(is.na(estimate_transformed) | all(key %in% character(0))),
                   .groups = "drop") %>%
  dplyr::mutate(row.label = dplyr::if_else(bold & row.label != "",
                                           paste0("**", row.label, "**"),
                                           as.character(row.label))) %>% 
  dplyr::arrange(row) %>%
  dplyr::pull(row.label)

# Identify any CIs that extend outside axis limits
datatoplot <- datatoplot %>%
  dplyr::mutate(cioverright  = (uci_transformed > 0.7),
                uci_transformed = pmin(uci_transformed, 0.7),
                lci_transformed = pmin(lci_transformed, 0.7),
                cioverleft  = (lci_transformed < 0.2),
                lci_transformed = pmax(lci_transformed, 0.2),
                uci_transformed = pmax(uci_transformed, 0.2))

# Create the ggplot
plot_cor_es <- ggplot(datatoplot, aes(y=-row, x=estimate_transformed)) +
  
  # Put the different panels in side-by-side plots using facets
  facet_wrap(~panel, nrow = 1) +
  
  # Add a line at null effect
  annotate(geom = "segment",
           y = -0.7, yend = -Inf,
           x = 0.5, xend = 0.5,
           linewidth = 0.545454545454545,
           linetype = 3,
           colour = "black") +
  
  # Plot points at the transformed estimates
  ## Scale by inverse of the SE
  geom_point(aes(size = size),
             data = ~ dplyr::filter(.x, estimate_transformed > 0.2, estimate_transformed < 0.7),
             shape = 15,
             colour = "black",
             fill = "black",
             stroke = 0,
             na.rm = TRUE) +
  
  # Scale the size of points by their side length
  # and make the scale range from zero upwards
  scale_radius(limits = c(0, 1),
               range = c(0, 2)) +
  
  # Plot the CIs
  geom_errorbar(aes(xmin = lci_transformed,
                    xmax = uci_transformed),
                data = ~ dplyr::filter(.x, !is.na(estimate_transformed)),
                colour = "black",
                width = 0,
                linewidth = 0.545454545454545,
                na.rm = TRUE) +
  
  # Add tiny segments with arrows when the CIs go outside axis limits
  geom_segment(aes(y = -row,
                   yend = -row,
                   x = uci_transformed-0.000001,
                   xend = uci_transformed),
               data = ~ dplyr::filter(.x, cioverright == TRUE),
               colour = "black",
               linewidth = 0.545454545454545,
               arrow = arrow(type = "closed", length = unit(4.36363636363636, "pt")),
               na.rm = TRUE) +
  geom_segment(aes(y = -row,
                   yend = -row,
                   x = lci_transformed+0.000001,
                   xend = lci_transformed),
               data = ~ dplyr::filter(.x, cioverleft == TRUE),
               colour = "black",
               linewidth = 0.545454545454545,
               arrow = arrow(type = "closed", length = unit(4.36363636363636, "pt")),
               na.rm = TRUE) +
  
  # Use identity for aesthetic scales
  scale_shape_identity() +
  scale_fill_identity() +
  scale_colour_identity() +
  
  # Set coordinate system
  coord_cartesian(clip = "off",
                  xlim = c(0.2, 0.7)) +
  
  # Add columns to right side of plots
  ## column auto_estcolumn
  ckbplotr::geom_text_move(aes(y = -row,
                               x = 0.7,
                               label = `auto_estcolumn`),
                           move_x = unit(0.6, "mm"),
                           hjust = 0,
                           size  = 3.27272727272727,
                           colour  = "black",
                           na.rm = TRUE,
                           parse = FALSE) +
  ckbplotr::geom_text_move(aes(y = 0,
                               x = 0.7,
                               label = title),
                           move_x = unit(0.6, "mm"),
                           hjust    = 0,
                           size     = 3.27272727272727,
                           colour  = "black",
                           fontface = "bold",
                           data = dplyr::tibble(panel = factor(c("1", "2"),
                                                               levels = c("1", "2"),
                                                               ordered = TRUE),
                                                title = "Pearson\'s r (95% CI)")) +
  
  # Add columns to left side of plots
  ## column freq
  ckbplotr::geom_text_move(aes(y = -row,
                               x = 0.2,
                               label = freq,
                               fontface = "plain"),
                           move_x = unit(-0.6, "mm"),
                           hjust = 1,
                           size  = 3.27272727272727,
                           colour  = "black",
                           na.rm = TRUE) +
  ckbplotr::geom_text_move(aes(y = 0,
                               x = 0.2,
                               label = title),
                           move_x = unit(-0.6, "mm"),
                           hjust    = 1,
                           size     = 3.27272727272727,
                           colour  = "black",
                           fontface = "bold",
                           data = dplyr::tibble(panel = factor(c("1", "2"),
                                                               levels = c("1", "2"),
                                                               ordered = TRUE),
                                                title = "No. shared hits")) +
  
  # Add xlab below each axis
  geom_text(aes(y = -Inf, x = 0.45, label = xlab),
            hjust = 0.5,
            size  = 3.27272727272727,
            colour  = "black",
            vjust = 4.4,
            fontface = "bold",
            data = dplyr::tibble(panel = factor(c("1", "2"),
                                                levels = c("1", "2"),
                                                ordered = TRUE),
                                 xlab = "Pearson\'s r (95% CI)")) +
  
  # Add panel name above each panel
  geom_text(aes(y = 0, x = 0.45, label = title),
            hjust = 0.5,
            nudge_y = 2,
            size  = 3.27272727272727,
            colour  = "black",
            fontface = "bold",
            data = dplyr::tibble(panel = factor(c("1", "2"),
                                                levels = c("1", "2"),
                                                ordered = TRUE),
                                 title = c("OLINK vs SomaScan-ANML", "OLINK vs SomaScan-non-ANML"))) +
  
  # Set the scale for the x axis (the estimates and CIs)
  scale_x_continuous(trans  = "identity",
                     breaks = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
                     expand = c(0,0)) +
  
  # Set the scale for the y axis (the rows)
  scale_y_continuous(breaks = -1:-max(datatoplot$row),
                     labels = rowlabels,
                     limits = c(-max(datatoplot$row) - 0.7, NA),
                     expand = c(0,0)) +
  
  # Control the overall look of the plots
  theme(text             = element_text(size = 12, colour = "black"),
        line             = element_line(linewidth = 0.545454545454545),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title       = element_blank(),
        axis.line.x      = element_line(colour = "black", linewidth = 0.545454545454545, lineend = "round"),
        axis.title       = element_blank(),
        axis.ticks.x     = element_line(colour = "black"),
        axis.text.x      = element_text(colour = "black",,
                                        margin = margin(t = 4.8),
                                        vjust  = 1),
        axis.ticks.y     = element_blank(),
        axis.line.y      = element_blank(),
        axis.text.y      = ggtext::element_markdown(hjust  = 0,
                                                    colour = "black",
                                                    margin = margin(r = 25.8, unit = "mm")),
        panel.border     = element_blank(),
        panel.spacing    = unit(31.1, "mm") + unit(5, "mm") + unit(25.8, "mm"),
        strip.background = element_blank(),
        strip.placement  = "outside",
        strip.text       = element_blank(),
        legend.position  = "none",
        plot.background  = element_blank(),
        plot.margin      = margin(8, 8, 8, 8, "mm") + unit(c(0, 31.1, 0, 0), "mm"))

plot_cor_es_bg <- plot_cor_es + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_cor_es_fdr_subcohort.png",plot_cor_es_bg,width=14,height=6)

plot_cor_es_soma <- forest_plot(panels = list(cor_es_soma),
                                exponentiate = F,
                                panel.headings = c("SomaScan ANML vs non-ANML"),
                                row.labels = row_labels,
                                row.labels.levels = c("group","subgroup"),
                                rows = unique(row_labels$group),
                                col.key = "label",
                                col.estimate = "r",
                                col.lci = "lci",
                                col.uci = "uci",
                                col.left = "freq",
                                col.left.heading = "No. shared hits",
                                xlab = "Pearson's r (95% CI)",
                                # xlim = c(0.2,0.7),
                                # estcolumn = F,
                                pointsize = 2,
                                base_size = 12,
                                col.right.heading = "Pearson's r (95% CI)"
)

plot_cor_es_soma_bg <- plot_cor_es_soma$plot + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_cor_es_soma_fdr_subcohort.png",plot_cor_es_soma_bg,width=7,height=6)



## correlation coefficient between beta, only for shared associations

cor_es_normal_sig <- data.frame(matrix(ncol = 5, nrow = 0))
cor_es_non_normal_sig <- data.frame(matrix(ncol = 5, nrow = 0))
cor_es_soma_sig <- data.frame(matrix(ncol = 5, nrow = 0))

# remove variables with two few shared hits

variable_sig <- variable[! variable %in% c("married","school","poor_health","kidney_dis_diag","cancer_diag","met")]
label_sig <- label[! label %in% c("Married","≥ 6 years education","Poor self-rated health","Kidney disease","Cancer","Physical activity")]
level_sig <- level[! level %in% c("Married","≥ 6 years education","Poor self-rated health","Kidney disease","Cancer","Physical activity")]

for (i in 1:length(variable_sig)) {
  
  df <- overlap_1_to_1_assoc[c(grep(paste0("es_",variable_sig[i]),names(overlap_1_to_1_assoc)),
                               grep(paste0("p_fdr_",variable_sig[i]),names(overlap_1_to_1_assoc)))]
  
  df_1 <- df[df[4]<0.05 & df[5]<0.05 & sign(df[1])==sign(df[2]),]
  
  df_1_result <- data.frame(variable=variable_sig[i],label=label_sig[i],
                            r=cor(df_1[[1]],df_1[[2]]),
                            lci=cor.test(df_1[[1]],df_1[[2]])$conf.int[1],
                            uci=cor.test(df_1[[1]],df_1[[2]])$conf.int[2])
  
  cor_es_normal_sig <- rbind(cor_es_normal_sig,df_1_result)
  
  df_2 <- df[df[4]<0.05 & df[6]<0.05 & sign(df[1])==sign(df[3]),]
  
  df_2_result <- data.frame(variable=variable_sig[i],label=label_sig[i],
                            r=cor(df_2[[1]],df_2[[3]]),
                            lci=cor.test(df_2[[1]],df_2[[3]])$conf.int[1],
                            uci=cor.test(df_2[[1]],df_2[[3]])$conf.int[2])
  
  cor_es_non_normal_sig <- rbind(cor_es_non_normal_sig,df_2_result)
  
  df_3 <- df[df[5]<0.05 & df[6]<0.05 & sign(df[2])==sign(df[3]),]
  
  df_3_result <- data.frame(variable=variable_sig[i],label=label_sig[i],
                            r=cor(df_3[[2]],df_3[[3]]),
                            lci=cor.test(df_3[[2]],df_3[[3]])$conf.int[1],
                            uci=cor.test(df_3[[2]],df_3[[3]])$conf.int[2])
  
  cor_es_soma_sig <- rbind(cor_es_soma_sig,df_3_result)
}

cor_es_normal_sig$label <- factor(cor_es_normal_sig$label, levels = level_sig)
shared_normal_sig <- trans_normal_fdr[trans_normal_fdr$concord=="Shared",c("label","concord","freq")]
cor_es_normal_sig <- merge(cor_es_normal_sig,shared_normal,by="label")

cor_es_non_normal_sig$label <- factor(cor_es_non_normal_sig$label, levels = level_sig)
shared_non_normal_sig <- trans_non_normal_fdr[trans_non_normal_fdr$concord=="Shared",c("label","concord","freq")]
cor_es_non_normal_sig <- merge(cor_es_non_normal_sig,shared_non_normal,by="label")

cor_es_soma_sig$label <- factor(cor_es_soma_sig$label, levels = level_sig)
shared_soma_sig <- trans_soma_fdr[trans_soma_fdr$concord=="Shared",c("label","concord","freq")]
cor_es_soma_sig <- merge(cor_es_soma_sig,shared_soma,by="label")

cor_es_sig <- list(cor_es_normal_sig,cor_es_non_normal_sig)

row_labels_sig <- data.frame(
  subgroup = level_sig,
  group    = c("Socio-demographics","Socio-demographics","Socio-demographics",
               "Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements","Clinical measurements",
               "Medical history",
               "Lifestyle factors","Lifestyle factors"),
  label = level_sig)

plot_cor_es_sig <- forest_plot(panels = cor_es_sig,
                               exponentiate = F,
                               panel.headings = c("OLINK vs SomaScan-ANML","OLINK vs SomaScan-non-ANML"),
                               row.labels = row_labels_sig,
                               row.labels.levels = c("group","subgroup"),
                               rows = unique(row_labels$group),
                               col.key = "label",
                               col.estimate = "r",
                               col.lci = "lci",
                               col.uci = "uci",
                               col.left = "freq",
                               col.left.heading = "No. shared hits",
                               xlab = "Pearson's r (95% CI)",
                               # xlim = c(0.2,0.7),
                               # estcolumn = F,
                               pointsize = 2,
                               base_size = 12,
                               col.right.heading = "Pearson's r (95% CI)"
)

plot_cor_es_sig_bg <- plot_cor_es_sig$plot + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_cor_es_sig_fdr_subcohort.png",plot_cor_es_sig_bg,width=14,height=6)

plot_cor_es_soma_sig <- forest_plot(panels = list(cor_es_soma_sig),
                                    exponentiate = F,
                                    panel.headings = c("SomaScan ANML vs non-ANML"),
                                    row.labels = row_labels_sig,
                                    row.labels.levels = c("group","subgroup"),
                                    rows = unique(row_labels$group),
                                    col.key = "label",
                                    col.estimate = "r",
                                    col.lci = "lci",
                                    col.uci = "uci",
                                    col.left = "freq",
                                    col.left.heading = "No. shared hits",
                                    xlab = "Pearson's r (95% CI)",
                                    # xlim = c(0.2,0.7),
                                    # estcolumn = F,
                                    pointsize = 2,
                                    base_size = 12,
                                    col.right.heading = "Pearson's r (95% CI)"
)

plot_cor_es_soma_sig_bg <- plot_cor_es_soma_sig$plot + theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/plot_cor_es_soma_sig_fdr_subcohort.png",plot_cor_es_soma_sig_bg,width=7,height=6)


# write.csv(overlap_1_to_1_assoc,"overlap_1_to_1_assoc_sig.csv", quote=F, row.names=F)