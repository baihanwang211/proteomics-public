### this script plots the results of IHD risk prediction

rm(list = ls())

setwd("")

library(ckbplotr)

c_stat <- read.csv("c_statistics.csv")

nri<- read.csv("nri.csv")

# get NRI

nri$uci <- formatC(nri$uci * 100,digits=1,format="f")
nri$lci <- formatC(nri$lci * 100,digits=1,format="f")
nri$value <- formatC(nri$value * 100,digits=1,format="f")

nri$nri <- paste(nri$value,"% (",nri$lci,"%, ",nri$uci,"%)",sep="")

nri <- nri[,c("model","type","nri")]

c_stat_nri <- merge(c_stat,nri,by="model")

c_stat_nri$nri[c_stat_nri$name == "All overlapping" | c_stat_nri$name == "Conventional" | c_stat_nri$name == "All overlapping (both)" |
               c_stat_nri$name == "FDR significant" |c_stat_nri$name == "Shared (ANML)" | c_stat_nri$name == "Shared (non-ANML)"] <- "-"

# order 

c_stat_nri$model <- factor(c_stat_nri$model,levels = c("conventional",
                                                       "all_overlap_both",
                                                       "all_overlap_both_noanml",
                                               "all_overlap_olink",
                                               "all_overlap_somascan",
                                               "all_overlap_somascan_noanml",
                                               "significant_olink",
                                               "significant_somascan",
                                               "significant_somascan_noanml",
                                               "shared_significant_olink",
                                               "shared_significant_somascan",
                                               "shared_significant_olink_noanml",
                                               "shared_significant_somascan_noanml",
                                               "conventional_plus_all_overlap_olink",
                                               "conventional_plus_all_overlap_somascan",
                                               "conventional_plus_all_overlap_somascan_noanml",
                                               "conventional_plus_all_overlap_both",
                                               "conventional_plus_all_overlap_both_noanml",
                                               "conventional_plus_significant_olink",
                                               "conventional_plus_significant_somascan",
                                               "conventional_plus_significant_somascan_noanml",
                                               "conventional_plus_shared_significant_olink",
                                               "conventional_plus_shared_significant_somascan",
                                               "conventional_plus_shared_significant_olink_noanml",
                                               "conventional_plus_shared_significant_somascan_noanml"))

c_stat_nri <- c_stat_nri[order(c_stat_nri$model),]

# add number of proteins

c_stat_nri$n <- "-"


c_stat_nri$n[c_stat_nri$name=="All overlapping (both)"] <- '1694 * 2'
c_stat_nri$n[c_stat_nri$name=="All overlapping (both)"] <- '1694 * 2'

c_stat_nri$n[c_stat_nri$name=="Conventional + all overlapping (both)"] <- '1694 * 2'
c_stat_nri$n[c_stat_nri$name=="Conventional + all overlapping (both)"] <- '1694 * 2'

c_stat_nri$n[c_stat_nri$name=="All overlapping"] <- 1694
c_stat_nri$n[c_stat_nri$name=="Conventional + all overlapping"] <- 1694

c_stat_nri$n[c_stat_nri$model=="significant_olink"] <- 279
c_stat_nri$n[c_stat_nri$model=="conventional_plus_significant_olink"] <- 279

c_stat_nri$n[c_stat_nri$model=="significant_somascan"] <- 165
c_stat_nri$n[c_stat_nri$model=="conventional_plus_significant_somascan"] <- 165

c_stat_nri$n[c_stat_nri$model=="significant_somascan_noanml"] <- 154
c_stat_nri$n[c_stat_nri$model=="conventional_plus_significant_somascan_noanml"] <- 154

c_stat_nri$n[c_stat_nri$name=="Shared (ANML)"] <- 78
c_stat_nri$n[c_stat_nri$name=="Conventional + shared (ANML)"] <- 78

c_stat_nri$n[c_stat_nri$name=="Shared (non-ANML)"] <- 78
c_stat_nri$n[c_stat_nri$name=="Conventional + shared (non-ANML)"] <- 78

c_stat_nri$n <- as.character(c_stat_nri$n)

# cat model

c_stat_nri_cat <- c_stat_nri[c_stat_nri$type=="percentile_nri",]

row_labels <- data.frame(
  model = as.character(c_stat_nri_cat$model),
  group    = c_stat_nri_cat$name,
  label    = c_stat_nri_cat$platform)

## plot using ckbplotr
plot_prediction <- forest_plot(panels = list(c_stat_nri_cat),
            exponentiate = F,
            col.key = "model",
            row.labels = row_labels,
            row.labels.levels = c("group","label"),
            rows = unique(row_labels$group),
            col.estimate = "value",
            col.stderr = "se",
            xlab = "C-statistic",
            xlim = c(0.75,0.9),
            digits=3,
            pointsize = 2,
            base_size = 12,
            nullval = c_stat_nri_cat$value[1],
            col.right.heading = c("C-statistic (95% CI)","NRI (95% CI)"),
            col.right="nri",
            col.left = "n",
            col.left.heading = "No. proteins")

# Prepare data to be plotted using ckbplotr::forest_data()
datatoplot <- ckbplotr::forest_data(row.labels = row_labels,
                                    row.labels.levels = c("group", "label"),
                                    rows = c("Conventional", "All overlapping (both)", "All overlapping", "FDR significant", "Shared (ANML)", "Shared (non-ANML)", "Conventional + all overlapping (both)", "Conventional + all overlapping", "Conventional + FDR significant", "Conventional + shared (ANML)", "Conventional + shared (non-ANML)"),
                                    panels = list(c_stat_nri_cat),
                                    panel.names = "1",
                                    col.key = "model",
                                    col.estimate = "value",
                                    col.stderr = "se",
                                    col.left = "n",
                                    col.right = "nri",
                                    digits = 3,
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
  dplyr::mutate(cioverright  = (uci_transformed > 0.9),
                uci_transformed = pmin(uci_transformed, 0.9),
                lci_transformed = pmin(lci_transformed, 0.9),
                cioverleft  = (lci_transformed < 0.75),
                lci_transformed = pmax(lci_transformed, 0.75),
                uci_transformed = pmax(uci_transformed, 0.75))

# Create the ggplot
plot_prediction <- ggplot(datatoplot, aes(y=-row, x=estimate_transformed)) +
  
  # Put the different panels in side-by-side plots using facets
  facet_wrap(~panel, nrow = 1) +
  
  # Add a line at null effect
  annotate(geom = "segment",
           y = -0.7, yend = -Inf,
           x = 0.844965482, xend = 0.844965482,
           linewidth = 0.545454545454545,
           linetype = 3,
           colour = "black") +
  
  # Plot points at the transformed estimates
  ## Scale by inverse of the SE
  geom_point(aes(size = size),
             data = ~ dplyr::filter(.x, estimate_transformed > 0.75, estimate_transformed < 0.9),
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
                  xlim = c(0.75, 0.9)) +
  
  # Add columns to right side of plots
  ## column auto_estcolumn
  ckbplotr::geom_text_move(aes(y = -row,
                               x = 0.9,
                               label = `auto_estcolumn`),
                           move_x = unit(0.8, "mm"),
                           hjust = 0,
                           size  = 3.27272727272727,
                           colour  = "black",
                           na.rm = TRUE,
                           parse = FALSE) +
  ckbplotr::geom_text_move(aes(y = 0,
                               x = 0.9,
                               label = title),
                           move_x = unit(0.8, "mm"),
                           hjust    = 0,
                           size     = 3.27272727272727,
                           colour  = "black",
                           fontface = "bold",
                           data = dplyr::tibble(panel = factor("1",
                                                               levels = "1",
                                                               ordered = TRUE),
                                                title = "C-statistic (95% CI)")) +
  
  ## column nri
  ckbplotr::geom_text_move(aes(y = -row,
                               x = 0.9,
                               label = nri),
                           move_x = unit(35.4, "mm"),
                           hjust = 0,
                           size  = 3.27272727272727,
                           colour  = "black",
                           na.rm = TRUE,
                           parse = FALSE) +
  ckbplotr::geom_text_move(aes(y = 0,
                               x = 0.9,
                               label = title),
                           move_x = unit(35.4, "mm"),
                           hjust    = 0,
                           size     = 3.27272727272727,
                           colour  = "black",
                           fontface = "bold",
                           data = dplyr::tibble(panel = factor("1",
                                                               levels = "1",
                                                               ordered = TRUE),
                                                title = "NRI (95% CI)")) +
  
  # Add columns to left side of plots
  ## column n
  ckbplotr::geom_text_move(aes(y = -row,
                               x = 0.75,
                               label = n,
                               fontface = "plain"),
                           move_x = unit(-0.8, "mm"),
                           hjust = 1,
                           size  = 3.27272727272727,
                           colour  = "black",
                           na.rm = TRUE) +
  ckbplotr::geom_text_move(aes(y = 0,
                               x = 0.75,
                               label = title),
                           move_x = unit(-0.8, "mm"),
                           hjust    = 1,
                           size     = 3.27272727272727,
                           colour  = "black",
                           fontface = "bold",
                           data = dplyr::tibble(panel = factor("1",
                                                               levels = "1",
                                                               ordered = TRUE),
                                                title = "No. proteins")) +
  
  # Add xlab below each axis
  geom_text(aes(y = -Inf, x = 0.825, label = xlab),
            hjust = 0.5,
            size  = 3.27272727272727,
            colour  = "black",
            vjust = 4.4,
            fontface = "bold",
            data = dplyr::tibble(panel = factor("1",
                                                levels = "1",
                                                ordered = TRUE),
                                 xlab = "C-statistic")) +
  
  # Set the scale for the x axis (the estimates and CIs)
  scale_x_continuous(trans  = "identity",
                     breaks = c(0.75, 0.8, 0.85, 0.9),
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
                                                    margin = margin(r = 22.2, unit = "mm")),
        panel.border     = element_blank(),
        panel.spacing    = unit(72.1, "mm") + unit(5, "mm") + unit(22.2, "mm"),
        strip.background = element_blank(),
        strip.placement  = "outside",
        strip.text       = element_blank(),
        legend.position  = "none",
        plot.background  = element_blank(),
        plot.margin      = margin(8, 8, 8, 8, "mm") + unit(c(0, 72.1, 0, 0), "mm"))

plot_prediction_bg <- plot_prediction + theme(plot.background = element_rect(fill = "white", colour = NA))

plot_prediction_bg

ggsave("",plot_prediction_bg,width=10,height=8)

######## categorical-free

c_stat_nri_cat_free <- c_stat_nri[c_stat_nri$type=="nri",]

row_labels <- data.frame(
  model = as.character(c_stat_nri_cat_free$model),
  group    = c_stat_nri_cat_free$name,
  label    = c_stat_nri_cat_free$platform)

## plot using ckbplotr
plot_prediction <- forest_plot(panels = list(c_stat_nri_cat_free),
                               exponentiate = F,
                               col.key = "model",
                               row.labels = row_labels,
                               row.labels.levels = c("group","label"),
                               rows = unique(row_labels$group),
                               col.estimate = "value",
                               col.stderr = "se",
                               xlab = "C-statistic",
                               xlim = c(0.75,0.9),
                               digits=3,
                               pointsize = 2,
                               base_size = 12,
                               nullval = c_stat_nri_cat_free$value[1],
                               col.right.heading = c("C-statistic (95% CI)","NRI (95% CI)"),
                               col.right="nri",
                               col.left = "n",
                               col.left.heading = "No. proteins")

# Prepare data to be plotted using ckbplotr::forest_data()
datatoplot <- ckbplotr::forest_data(row.labels = row_labels,
                                    row.labels.levels = c("group", "label"),
                                    rows = c("Conventional", "All overlapping (both)", "All overlapping", "FDR significant", "Shared (ANML)", "Shared (non-ANML)", "Conventional + all overlapping (both)", "Conventional + all overlapping", "Conventional + FDR significant", "Conventional + shared (ANML)", "Conventional + shared (non-ANML)"),
                                    panels = list(c_stat_nri_cat_free),
                                    panel.names = "1",
                                    col.key = "model",
                                    col.estimate = "value",
                                    col.stderr = "se",
                                    col.left = "n",
                                    col.right = "nri",
                                    digits = 3,
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
  dplyr::mutate(cioverright  = (uci_transformed > 0.9),
                uci_transformed = pmin(uci_transformed, 0.9),
                lci_transformed = pmin(lci_transformed, 0.9),
                cioverleft  = (lci_transformed < 0.75),
                lci_transformed = pmax(lci_transformed, 0.75),
                uci_transformed = pmax(uci_transformed, 0.75))

# Create the ggplot
plot_prediction <- ggplot(datatoplot, aes(y=-row, x=estimate_transformed)) +
  
  # Put the different panels in side-by-side plots using facets
  facet_wrap(~panel, nrow = 1) +
  
  # Add a line at null effect
  annotate(geom = "segment",
           y = -0.7, yend = -Inf,
           x = 0.844965482, xend = 0.844965482,
           linewidth = 0.545454545454545,
           linetype = 3,
           colour = "black") +
  
  # Plot points at the transformed estimates
  ## Scale by inverse of the SE
  geom_point(aes(size = size),
             data = ~ dplyr::filter(.x, estimate_transformed > 0.75, estimate_transformed < 0.9),
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
                  xlim = c(0.75, 0.9)) +
  
  # Add columns to right side of plots
  ## column auto_estcolumn
  ckbplotr::geom_text_move(aes(y = -row,
                               x = 0.9,
                               label = `auto_estcolumn`),
                           move_x = unit(0.8, "mm"),
                           hjust = 0,
                           size  = 3.27272727272727,
                           colour  = "black",
                           na.rm = TRUE,
                           parse = FALSE) +
  ckbplotr::geom_text_move(aes(y = 0,
                               x = 0.9,
                               label = title),
                           move_x = unit(0.8, "mm"),
                           hjust    = 0,
                           size     = 3.27272727272727,
                           colour  = "black",
                           fontface = "bold",
                           data = dplyr::tibble(panel = factor("1",
                                                               levels = "1",
                                                               ordered = TRUE),
                                                title = "C-statistic (95% CI)")) +
  
  ## column nri
  ckbplotr::geom_text_move(aes(y = -row,
                               x = 0.9,
                               label = nri),
                           move_x = unit(35.4, "mm"),
                           hjust = 0,
                           size  = 3.27272727272727,
                           colour  = "black",
                           na.rm = TRUE,
                           parse = FALSE) +
  ckbplotr::geom_text_move(aes(y = 0,
                               x = 0.9,
                               label = title),
                           move_x = unit(35.4, "mm"),
                           hjust    = 0,
                           size     = 3.27272727272727,
                           colour  = "black",
                           fontface = "bold",
                           data = dplyr::tibble(panel = factor("1",
                                                               levels = "1",
                                                               ordered = TRUE),
                                                title = "NRI (95% CI)")) +
  
  # Add columns to left side of plots
  ## column n
  ckbplotr::geom_text_move(aes(y = -row,
                               x = 0.75,
                               label = n,
                               fontface = "plain"),
                           move_x = unit(-0.8, "mm"),
                           hjust = 1,
                           size  = 3.27272727272727,
                           colour  = "black",
                           na.rm = TRUE) +
  ckbplotr::geom_text_move(aes(y = 0,
                               x = 0.75,
                               label = title),
                           move_x = unit(-0.8, "mm"),
                           hjust    = 1,
                           size     = 3.27272727272727,
                           colour  = "black",
                           fontface = "bold",
                           data = dplyr::tibble(panel = factor("1",
                                                               levels = "1",
                                                               ordered = TRUE),
                                                title = "No. proteins")) +
  
  # Add xlab below each axis
  geom_text(aes(y = -Inf, x = 0.825, label = xlab),
            hjust = 0.5,
            size  = 3.27272727272727,
            colour  = "black",
            vjust = 4.4,
            fontface = "bold",
            data = dplyr::tibble(panel = factor("1",
                                                levels = "1",
                                                ordered = TRUE),
                                 xlab = "C-statistic")) +
  
  # Set the scale for the x axis (the estimates and CIs)
  scale_x_continuous(trans  = "identity",
                     breaks = c(0.75, 0.8, 0.85, 0.9),
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
                                                    margin = margin(r = 22.2, unit = "mm")),
        panel.border     = element_blank(),
        panel.spacing    = unit(72.1, "mm") + unit(5, "mm") + unit(22.2, "mm"),
        strip.background = element_blank(),
        strip.placement  = "outside",
        strip.text       = element_blank(),
        legend.position  = "none",
        plot.background  = element_blank(),
        plot.margin      = margin(8, 8, 8, 8, "mm") + unit(c(0, 72.1, 0, 0), "mm"))

plot_prediction_bg <- plot_prediction + theme(plot.background = element_rect(fill = "white", colour = NA))

plot_prediction_bg

ggsave("",plot_prediction_bg,width=10,height=8)
