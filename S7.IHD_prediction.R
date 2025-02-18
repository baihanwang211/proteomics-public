##### this script plots the results of IHD risk prediction

rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data/ihd")

library(ckbplotr)

# load files with results of c statistics and nri

c_stat <- read.csv("c_statistics.csv")

nri<- read.csv("nri.csv")

# get NRI with CI

nri$uci <- formatC(nri$uci * 100,digits=1,format="f")
nri$lci <- formatC(nri$lci * 100,digits=1,format="f")
nri$value <- formatC(nri$value * 100,digits=1,format="f")

nri$nri <- paste(nri$value,"% (",nri$lci,"%, ",nri$uci,"%)",sep="")

nri <- nri[,c("model","type","nri")]

c_stat_nri <- merge(c_stat,nri,by="model")

row_labels <- data.frame(
  model = as.character(c_stat_nri$model),
  group    = c_stat_nri$name_figure,
  label    = c_stat_nri$platform)

## plot using ckbplotr
plot_prediction_figure <- forest_plot(panels = list(c_stat_nri),
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
                               nullval = c_stat_nri$value[1],
                               col.right.heading = c("C-statistic (95% CI)","NRI (95% CI)"),
                               col.right="nri",
                               col.left = "n",
                               col.left.heading = "No. proteins")


## add dotted line as null

# Prepare data to be plotted using ckbplotr::forest_data()
datatoplot <- ckbplotr::forest_data(panels = list(c_stat_nri),
                                    row.labels = row_labels,
                                    row.labels.levels = c("group", "label"),
                                    rows = c("Conventional", "Proteins only", "Conventional + proteins"),
                                    panel.names = "1",
                                    col.key = "model",
                                    col.estimate = "value",
                                    col.stderr = "se",
                                    col.left = "n",
                                    col.right = "nri",
                                    digits = 3,
                                    exponentiate = FALSE)

# Create the ggplot
plot_prediction_figure <- ggplot(datatoplot, aes(y = row, x = estimate_transformed)) +
  
  # Put the different panels in side-by-side plots using facets
  facet_wrap(vars(panel), nrow = 1) +
  
  # Add a line at null effect
  annotate(geom      = "segment",
           y         = 0.7,
           yend      = Inf,
           x         = 0.844965482,
           xend      = 0.844965482,
           linewidth = 0.545454545454545,
           colour    = "black",
           linetype = 3) +
  
  # Plot points at the transformed estimates
  geom_point(data   = ~ dplyr::filter(.x,
                                      estimate_transformed > 0.75,
                                      estimate_transformed < 0.9,
                                      !as_diamond),
             shape  = 15,
             size   = 2,
             colour = "black",
             fill   = "black",
             stroke = 0,
             na.rm  = TRUE) +
  
  # Plot the CIs
  geom_errorbar(aes(xmin = pmin(pmax(lci_transformed, 0.75), 0.9),
                    xmax = pmin(pmax(uci_transformed, 0.75), 0.9)),
                data = ~ dplyr::filter(.x, !is.na(estimate_transformed), !as_diamond),
                colour    = "black",
                width     = 0,
                linewidth = 0.545454545454545,
                na.rm     = TRUE) +
  
  # Add columns to right side of panels
  ## column auto_estcolumn
  ckbplotr::geom_text_move(aes(y = row,
                               x = 0.9,
                               label = `auto_estcolumn`),
                           move_x  = unit(0.8, "mm"),
                           hjust   = 0,
                           size    = 3.374014,
                           colour  = "black",
                           na.rm   = TRUE,
                           parse   = FALSE) +
  ckbplotr::geom_text_move(aes(y     = - 0,
                               x     = 0.9,
                               label = title),
                           move_x  = unit(0.8, "mm"),
                           hjust    = 0,
                           size     = 3.374014,
                           colour   = "black",
                           fontface = "bold",
                           lineheight = 1,
                           data = ~ dplyr::tibble(panel = sort(unique(.[["panel"]])),
                                                  title = "C-statistic (95% CI)")) +
  
  ## column nri
  ckbplotr::geom_text_move(aes(y = row,
                               x = 0.9,
                               label = nri),
                           move_x  = unit(36.1, "mm"),
                           hjust   = 0,
                           size    = 3.374014,
                           colour  = "black",
                           na.rm   = TRUE,
                           parse   = FALSE) +
  ckbplotr::geom_text_move(aes(y     = - 0,
                               x     = 0.9,
                               label = title),
                           move_x  = unit(36.1, "mm"),
                           hjust    = 0,
                           size     = 3.374014,
                           colour   = "black",
                           fontface = "bold",
                           lineheight = 1,
                           data = ~ dplyr::tibble(panel = sort(unique(.[["panel"]])),
                                                  title = "NRI (95% CI)")) +
  
  # Add columns to left side of panel
  ## column n
  ckbplotr::geom_text_move(aes(y = row,
                               x = 0.75,
                               label = n),
                           move_x  = unit(-0.8, "mm"),
                           hjust   = 1,
                           size    = 3.374014,
                           colour  = "black",
                           na.rm   = TRUE,
                           parse   = FALSE) +
  ckbplotr::geom_text_move(aes(y     = - 0,
                               x     = 0.75,
                               label = title),
                           move_x  = unit(-0.8, "mm"),
                           hjust    = 1,
                           size     = 3.374014,
                           colour   = "black",
                           fontface = "bold",
                           lineheight = 1,
                           data = ~ dplyr::tibble(panel = sort(unique(.[["panel"]])),
                                                  title = "No. proteins")) +
  
  # Add xlab below each axis
  geom_text(aes(y = Inf,
                x = 0.825,,
                label = xlab),
            hjust    = 0.5,
            size     = 3.374014,
            colour   = "black",
            vjust    = 4.4,
            fontface = "bold",
            data = ~ dplyr::tibble(panel = sort(unique(.[["panel"]])),
                                   xlab = "C-statistic")) +
  
  # Set coordinate system
  coord_cartesian(clip = "off",
                  xlim = c(0.75, 0.9)) +
  
  # Set the scale for the x axis (the estimates and CIs)
  scale_x_continuous(trans  = "identity",
                     limits = c(0.75, 0.9),
                     breaks = c(0.75, 0.8, 0.85, 0.9),
                     expand = c(0,0)) +
  
  # Set the scale for the y axis (the rows)
  scale_y_continuous(trans = "reverse",
                     breaks = attr(datatoplot, "rowlabels")$row,
                     labels = attr(datatoplot, "rowlabels")$row.label,
                     limits = c(max(attr(datatoplot, "rowlabels")$row) + 0.7, NA),
                     expand = c(0,0)) +
  
  # Control the overall look of the plot
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
        axis.ticks.length.y = unit(0, "pt"),
        axis.line.y      = element_blank(),
        axis.text.y      = ggtext::element_markdown(hjust  = 0,
                                                    colour = "black",
                                                    margin = margin(r = 23, unit = "mm")),
        panel.border     = element_blank(),
        panel.spacing    = unit(71.1, "mm") + unit(5, "mm") + unit(23, "mm"),
        strip.background = element_blank(),
        strip.placement  = "outside",
        strip.text       = element_blank(),
        legend.position  = "none",
        plot.background  = element_blank(),
        plot.margin      = margin(8, 8, 8, 8, "mm") + unit(c(0, 71.1, 0, 0), "mm"))

plot_prediction_bg_figure <- plot_prediction_figure + theme(plot.background = element_rect(fill = "white", colour = NA))

plot_prediction_bg_figure