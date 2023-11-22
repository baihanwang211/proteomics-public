rm(list = ls())

library(ggplot2)
library(ggpubr)

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

olink_soma_normal <- read.csv("olink_soma_normal.csv")

olink_soma_non_normal <- read.csv("olink_soma_non_normal.csv")

overlap <- read.csv("overlap.csv")

# standardise protein levels

olink_soma_normal[c(15:ncol(olink_soma_normal))] <- scale(olink_soma_normal[c(15:ncol(olink_soma_normal))])

olink_soma_non_normal[c(15:ncol(olink_soma_non_normal))] <- scale(olink_soma_non_normal[c(15:ncol(olink_soma_non_normal))])

# get age2

olink_soma_normal$age2 <- olink_soma_normal$age * olink_soma_normal$age

olink_soma_non_normal$age2 <- olink_soma_non_normal$age * olink_soma_non_normal$age

# demographic results

results_demo_olink <- list()
results_demo_soma_normal <- list()
results_demo_soma_non_normal <- list()

# run regression

for (i in 1:nrow(overlap)){
  results_demo_olink[[i]] <- lm(paste(overlap$olink_id[i], "~", "age + is_female + age2 + region_code + somascan_plt_id + olink_plt_Neurology + olink_plt_Cardiometabolic + olink_plt_Inflammation + olink_plt_Oncology + olink_plt_Neurology_II + olink_plt_Cardiometabolic_II + olink_plt_Oncology_II + olink_plt_Inflammation_II"),olink_soma_normal)
  results_demo_soma_normal[[i]] <- lm(paste(overlap$somascan_id[i], "~", "age + is_female + age2 + region_code + somascan_plt_id + olink_plt_Neurology + olink_plt_Cardiometabolic + olink_plt_Inflammation + olink_plt_Oncology + olink_plt_Neurology_II + olink_plt_Cardiometabolic_II + olink_plt_Oncology_II + olink_plt_Inflammation_II"),olink_soma_normal)
  results_demo_soma_non_normal[[i]] <- lm(paste(overlap$somascan_id[i], "~", "age + is_female + age2 + region_code + somascan_plt_id + olink_plt_Neurology + olink_plt_Cardiometabolic + olink_plt_Inflammation + olink_plt_Oncology + olink_plt_Neurology_II + olink_plt_Cardiometabolic_II + olink_plt_Oncology_II + olink_plt_Inflammation_II"),olink_soma_non_normal)
  }

# age

model_age <- overlap

for (i in 1:nrow(overlap)) {
  model_age$coeff_olink[i] <- results_demo_olink[[i]]$coefficients["age"]
  model_age$coeff_somascan_normal[i] <- results_demo_soma_normal[[i]]$coefficients["age"]
  model_age$coeff_somascan_non_normal[i] <- results_demo_soma_non_normal[[i]]$coefficients["age"]
  model_age$p_olink[i] <- summary(results_demo_olink[[i]])$coefficients["age","Pr(>|t|)"]
  model_age$p_somascan_normal[i] <- summary(results_demo_soma_normal[[i]])$coefficients["age","Pr(>|t|)"]
  model_age$p_somascan_non_normal[i] <- summary(results_demo_soma_non_normal[[i]])$coefficients["age","Pr(>|t|)"]
}

model_age$significance_normal <- "None"
model_age$significance_normal[model_age$p_olink < 0.05 & model_age$p_somascan_normal < 0.05 ] <- "Both"
model_age$significance_normal[model_age$p_olink < 0.05 & model_age$p_somascan_normal >= 0.05 ] <- "Olink only"
model_age$significance_normal[model_age$p_olink >= 0.05 & model_age$p_somascan_normal < 0.05 ] <- "SomaScan only"
table(model_age$significance_normal)

model_age$significance_non_normal <- "None"
model_age$significance_non_normal[model_age$p_olink < 0.05 & model_age$p_somascan_non_normal < 0.05 ] <- "Both"
model_age$significance_non_normal[model_age$p_olink < 0.05 & model_age$p_somascan_non_normal >= 0.05 ] <- "Olink only"
model_age$significance_non_normal[model_age$p_olink >= 0.05 & model_age$p_somascan_non_normal < 0.05 ] <- "SomaScan only"
table(model_age$significance_non_normal)

# scatter plot

model_age_normal_plot <- ggplot(model_age, aes(x=coeff_olink, y=coeff_somascan_normal,col=rho_olink_somascan_normal)) + 
  geom_point(size=3) +
  xlim(-0.3,0.3) +
  ylim(-0.3,0.3) +
  scale_colour_viridis_c() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20))

cor.test(model_age$coeff_olink,model_age$coeff_somascan_normal, method="spearman")

model_age_non_normal_plot <- ggplot(model_age, aes(x=coeff_olink, y=coeff_somascan_non_normal,col=rho_olink_somascan_non_normal)) + 
  geom_point(size=3) +
  xlim(-0.3,0.3) +
  ylim(-0.3,0.3) +
  scale_colour_viridis_c() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20))

cor.test(model_age$coeff_olink,model_age$coeff_somascan_non_normal, method="spearman")

model_age_hist <- ggarrange(model_age_normal_plot,model_age_non_normal_plot,labels=c("A","B"),ncol=2,nrow=1)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/correlation_age.png",model_age_hist,width=20,height=6)

# histogram

# normal

model_age_normal_hist <- ggplot(model_age, aes(x=rho_olink_somascan_normal, fill=factor(significance_normal,levels = c("None","Olink only","SomaScan only","Both")))) + 
  geom_histogram(colour="black") +
  xlim(-0.5,1) +
  ylim(0,550) +
  ggtitle("Number of significant hits for age") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Significance") +
  theme_few()

# non_normal

model_age_non_normal_hist <- ggplot(model_age, aes(x=rho_olink_somascan_non_normal, fill=factor(significance_non_normal,levels = c("None","Olink only","SomaScan only","Both")))) + 
  geom_histogram(colour="black") +
  xlim(-0.5,1) +
  ylim(0,550) +
  ggtitle("Number of significant hits for age") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Significance") +
  theme_few()

model_age_hist <- ggarrange(model_age_normal_hist,model_age_non_normal_hist,labels=c("A","B"),ncol=1,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_age.png",model_age_hist,width=9,height=9)

# sex

model_sex <- overlap

for (i in 1:nrow(overlap)) {
  model_sex$coeff_olink[i] <- results_demo_olink[[i]]$coefficients["is_female"]
  model_sex$coeff_somascan_normal[i] <- results_demo_soma_normal[[i]]$coefficients["is_female"]
  model_sex$coeff_somascan_non_normal[i] <- results_demo_soma_non_normal[[i]]$coefficients["is_female"]
  model_sex$p_olink[i] <- summary(results_demo_olink[[i]])$coefficients["is_female","Pr(>|t|)"]
  model_sex$p_somascan_normal[i] <- summary(results_demo_soma_normal[[i]])$coefficients["is_female","Pr(>|t|)"]
  model_sex$p_somascan_non_normal[i] <- summary(results_demo_soma_non_normal[[i]])$coefficients["is_female","Pr(>|t|)"]
}

model_sex$significance_normal <- "None"
model_sex$significance_normal[model_sex$p_olink < 0.05 & model_sex$p_somascan_normal < 0.05 ] <- "Both"
model_sex$significance_normal[model_sex$p_olink < 0.05 & model_sex$p_somascan_normal >= 0.05 ] <- "Olink only"
model_sex$significance_normal[model_sex$p_olink >= 0.05 & model_sex$p_somascan_normal < 0.05 ] <- "SomaScan only"
table(model_sex$significance_normal)

model_sex$significance_non_normal <- "None"
model_sex$significance_non_normal[model_sex$p_olink < 0.05 & model_sex$p_somascan_non_normal < 0.05 ] <- "Both"
model_sex$significance_non_normal[model_sex$p_olink < 0.05 & model_sex$p_somascan_non_normal >= 0.05 ] <- "Olink only"
model_sex$significance_non_normal[model_sex$p_olink >= 0.05 & model_sex$p_somascan_non_normal < 0.05 ] <- "SomaScan only"
table(model_sex$significance_non_normal)

# scatter plot

model_sex_normal_plot <- ggplot(model_sex, aes(x=coeff_olink, y=coeff_somascan_normal,col=rho_olink_somascan_normal)) + 
  geom_point(size=3) +
  xlim(-1.5,1.5) +
  ylim(-1.5,1.5) +
  scale_colour_viridis_c() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/correlation_sex_normal.png",model_sex_normal_plot,width=15,height=10.8)

cor.test(model_sex$coeff_olink,model_sex$coeff_somascan_normal, method="spearman")

model_sex_non_normal_plot <- ggplot(model_sex, aes(x=coeff_olink, y=coeff_somascan_non_normal,col=rho_olink_somascan_non_normal)) + 
  geom_point(size=3) +
  xlim(-1.5,1.5) +
  ylim(-1.5,1.5) +
  scale_colour_viridis_c() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20))

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/correlation_sex_non_normal.png",model_sex_non_normal_plot,width=15,height=10.8)

cor.test(model_sex$coeff_olink,model_sex$coeff_somascan_non_normal, method="spearman")

# histogram

# normal

model_sex_normal_hist <- ggplot(model_sex, aes(x=rho_olink_somascan_normal, fill=factor(significance_normal,levels = c("None","Olink only","SomaScan only","Both")))) + 
  geom_histogram(colour="black") +
  xlim(-0.5,1) +
  ylim(0,550) +
  ggtitle("Number of significant hits for sex") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Significance") +
  theme_few()

# non_normal

model_sex_non_normal_hist <- ggplot(model_sex, aes(x=rho_olink_somascan_non_normal, fill=factor(significance_non_normal,levels = c("None","Olink only","SomaScan only","Both")))) + 
  geom_histogram(colour="black") +
  xlim(-0.5,1) +
  ylim(0,550) +
  ggtitle("Number of significant hits for sex") + 
  xlab("Spearman correlation coefficient") + 
  ylab("Frequency") +
  labs(fill="Significance") +
  theme_few()

model_sex_hist <- ggarrange(model_sex_normal_hist,model_sex_non_normal_hist,labels=c("A","B"),ncol=1,nrow=2)

ggsave("K:/kadoorie/Staff_Folders/BaihanW/proteomics/results/hist_sex.png",model_sex_hist,width=9,height=9)