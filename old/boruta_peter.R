library(Boruta)
library(ggplot2)

## run machine learning

# normalised

set.seed(47)

boruta_normal <- Boruta(rho_olink_soma_normal~., data = overlap_1_to_1_annot_normal, doTrace = 2, maxRuns = 5000)

saveRDS(boruta_normal,"boruta_normal.RDS")

# boruta_normal <- readRDS("boruta_normal.RDS")

print(boruta_normal)

boruta_normal_df <- attStats(boruta_normal)

print(boruta_normal_df)

# reshape data

boruta_normal_melt <- melt(as.data.frame(boruta_normal$ImpHistory))

boruta_normal_melt <- boruta_normal_melt[is.finite(boruta_normal_melt$value),]

# remove shadow variables

boruta_normal_melt <- boruta_normal_melt[-grep("shadow",boruta_normal_melt$variable),]

names(boruta_normal_melt) <- c("variable","importance")

# group by decision

boruta_normal_decision <- data.frame(variable=row.names(boruta_normal_df),decision=as.vector(boruta_normal_df$decision))

boruta_normal_melt <- merge(boruta_normal_melt,boruta_normal_decision,by="variable")

boruta_normal_melt$decision <- factor(boruta_normal_melt$decision,levels=c("Confirmed","Tentative","Rejected"))

# plot

boruta_normal_plot <- ggplot(boruta_normal_melt, aes(x=reorder(variable, importance, FUN = median), y=importance)) + 
  geom_boxplot(aes(fill=decision)) + 
  coord_flip() +
  ggtitle("Results of Boruta feature selection (ANML)") +
  xlab("Variable") +
  ylab("Importance") +
  labs(fill="Decision") +
  theme_few() +
  theme(text = element_text(size = 15))

boruta_normal_plot
