rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(ggplot2)
### this script loads results from protein & IHD associations and plots the results

library(ggpubr)
library(ggthemes)
library(egg)
library(forcats)
library(ggpattern)
library(ckbplotr)
library(scales)
library(ggrepel)
library(VennDiagram)
library(RColorBrewer)
library(readxl)

# load results
olink_1 <- read.csv("")
olink_2 <- read.csv("")
olink <- rbind(olink_1,olink_2)

soma_normal <- read_excel("",sheet = "Log transformed")
soma_non_normal <- read_excel("",sheet = "Log transformed, No ANML")

# edit columns
names(olink)
olink <- olink[,-1]
names(olink) <- c("id","coef","hr","se","p")
olink$lci <- olink$coef-olink$se*qnorm(0.975) 
olink$hci <- olink$coef+olink$se*qnorm(0.975) 
names(olink) <- paste0("olink_",names(olink))
olink <- olink[!duplicated(olink$olink_id),]
olink$olink_id <- gsub("-",".",olink$olink_id)
olink$olink_id <- gsub("_",".",olink$olink_id)

names(soma_normal)
soma_normal <- soma_normal[,c("aptname","estimate","std.error","p.value.orig")]
names(soma_normal) <- c("id","coef","se","p")
soma_normal$lci <- soma_normal$coef-soma_normal$se*qnorm(0.975) 
soma_normal$hci <- soma_normal$coef+soma_normal$se*qnorm(0.975) 
soma_normal$hr <- exp(soma_normal$coef)
names(soma_normal) <- paste0("soma_normal_",names(soma_normal))
names(soma_normal)[1] <- "somascan_id"

names(soma_non_normal)
soma_non_normal <- soma_non_normal[,c("aptname","estimate","std.error","p.value.orig")]
names(soma_non_normal) <- c("id","coef","se","p")
soma_non_normal$lci <- soma_non_normal$coef-soma_non_normal$se*qnorm(0.975) 
soma_non_normal$hci <- soma_non_normal$coef+soma_non_normal$se*qnorm(0.975) 
soma_non_normal$hr <- exp(soma_non_normal$coef)
names(soma_non_normal) <- paste0("soma_non_normal_",names(soma_non_normal))
names(soma_non_normal)[1] <- "somascan_id"

# overlapping list all
overlap_cor <- read.csv("overlap_cor.csv")
overlap_cor$olink_id <- gsub("_",".",overlap_cor$olink_id)

# merge
overlap_ihd <- merge(overlap_cor,olink,by="olink_id")
overlap_ihd <- merge(overlap_ihd,soma_normal,by="somascan_id")
overlap_ihd <- merge(overlap_ihd,soma_non_normal,by="somascan_id")

write.csv(overlap_ihd,"", quote=F, row.names=F)

# load overlapping list
overlap_1_to_1_cor <- read.csv("")
overlap_1_to_1_cor$olink_id <- gsub("_",".",overlap_1_to_1_cor$olink_id)

# merge
overlap_1_to_1_ihd <- merge(overlap_1_to_1_cor,olink,by="olink_id")
overlap_1_to_1_ihd <- merge(overlap_1_to_1_ihd,soma_normal,by="somascan_id")
overlap_1_to_1_ihd <- merge(overlap_1_to_1_ihd,soma_non_normal,by="somascan_id")

# adjust p values
overlap_1_to_1_ihd$olink_p_bonferroni <- p.adjust(overlap_1_to_1_ihd$olink_p,method="bonferroni")
overlap_1_to_1_ihd$soma_normal_p_bonferroni <- p.adjust(overlap_1_to_1_ihd$soma_normal_p,method="bonferroni")
overlap_1_to_1_ihd$soma_non_normal_p_bonferroni <- p.adjust(overlap_1_to_1_ihd$soma_non_normal_p,method="bonferroni")

overlap_1_to_1_ihd$olink_p_fdr <- p.adjust(overlap_1_to_1_ihd$olink_p,method="fdr")
overlap_1_to_1_ihd$soma_normal_p_fdr <- p.adjust(overlap_1_to_1_ihd$soma_normal_p,method="fdr")
overlap_1_to_1_ihd$soma_non_normal_p_fdr <- p.adjust(overlap_1_to_1_ihd$soma_non_normal_p,method="fdr")

write.csv(overlap_1_to_1_ihd,"", quote=F, row.names=F)

ihd_all_overlap <- overlap_1_to_1_ihd[,c(1:3)]
write.csv(ihd_all_overlap,"", quote=F, row.names=F)

ihd_olink_sig <- overlap_1_to_1_ihd[overlap_1_to_1_ihd$olink_p_bonferroni<0.05,c(1:3)]
write.csv(ihd_olink_sig,"", quote=F, row.names=F)

ihd_soma_normal_sig <- overlap_1_to_1_ihd[overlap_1_to_1_ihd$soma_normal_p_bonferroni<0.05,c(1:3)]
write.csv(ihd_soma_normal_sig,"", quote=F, row.names=F)

ihd_soma_non_normal_sig <- overlap_1_to_1_ihd[overlap_1_to_1_ihd$soma_non_normal_p_bonferroni<0.05,c(1:3)]
write.csv(ihd_soma_non_normal_sig,"", quote=F, row.names=F)

ihd_olink_soma_normal_shared <- overlap_1_to_1_ihd[overlap_1_to_1_ihd$olink_p_bonferroni<0.05 & 
                                                 overlap_1_to_1_ihd$soma_normal_p_bonferroni<0.05 &
                                                 sign(overlap_1_to_1_ihd$olink_coef)==sign(overlap_1_to_1_ihd$soma_normal_coef),
                                                 c(1:3)]
write.csv(ihd_olink_soma_normal_shared,"", quote=F, row.names=F)

ihd_olink_soma_non_normal_shared <- overlap_1_to_1_ihd[overlap_1_to_1_ihd$olink_p_bonferroni<0.05 & 
                                                 overlap_1_to_1_ihd$soma_non_normal_p_bonferroni<0.05 &
                                                 sign(overlap_1_to_1_ihd$olink_coef)==sign(overlap_1_to_1_ihd$soma_non_normal_coef),
                                               c(1:3)]
write.csv(ihd_olink_soma_non_normal_shared,"", quote=F, row.names=F)



ihd_olink_sig <- overlap_1_to_1_ihd[overlap_1_to_1_ihd$olink_p_fdr<0.05,c(1:3)]
write.csv(ihd_olink_sig,"", quote=F, row.names=F)

ihd_soma_normal_sig <- overlap_1_to_1_ihd[overlap_1_to_1_ihd$soma_normal_p_fdr<0.05,c(1:3)]
write.csv(ihd_soma_normal_sig,"", quote=F, row.names=F)

ihd_soma_non_normal_sig <- overlap_1_to_1_ihd[overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05,c(1:3)]
write.csv(ihd_soma_non_normal_sig,"", quote=F, row.names=F)

ihd_olink_soma_normal_shared <- overlap_1_to_1_ihd[overlap_1_to_1_ihd$olink_p_fdr<0.05 & 
                                                     overlap_1_to_1_ihd$soma_normal_p_fdr<0.05 &
                                                     sign(overlap_1_to_1_ihd$olink_coef)==sign(overlap_1_to_1_ihd$soma_normal_coef),
                                                   c(1:3)]
write.csv(ihd_olink_soma_normal_shared,"", quote=F, row.names=F)

ihd_olink_soma_non_normal_shared <- overlap_1_to_1_ihd[overlap_1_to_1_ihd$olink_p_fdr<0.05 & 
                                                         overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05 &
                                                         sign(overlap_1_to_1_ihd$olink_coef)==sign(overlap_1_to_1_ihd$soma_non_normal_coef),
                                                       c(1:3)]
write.csv(ihd_olink_soma_non_normal_shared,"", quote=F, row.names=F)




# merge with coloc
overlap_coloc <- read.csv("overlap_coloc_1_to_1.csv")
overlap_coloc <- overlap_coloc[,c("uniprot_id","coloc_normal_cis_tf","coloc_non_normal_cis_tf")]
n_coloc_normal <- sum(overlap_coloc$coloc_normal_cis_tf==T)
n_coloc_non_normal <- sum(overlap_coloc$coloc_non_normal_cis_tf==T)
overlap_1_to_1_ihd <- merge(overlap_1_to_1_ihd,overlap_coloc,by="uniprot_id")





####### fdr

# count

sum(overlap_1_to_1_ihd$olink_p_fdr<0.05)
sum(overlap_1_to_1_ihd$soma_normal_p_fdr<0.05)
sum(overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05)

cols <- hue_pal()(3)

hex <- brewer.pal(3,"Set1")[c(1,2)]

## volcano

set.seed(4747)

olink_4_pos <- overlap_1_to_1_ihd$olink_id[overlap_1_to_1_ihd$olink_coef>0][order(overlap_1_to_1_ihd$olink_p[overlap_1_to_1_ihd$olink_coef>0])[1:4]]
soma_normal_4_pos <- overlap_1_to_1_ihd$olink_id[overlap_1_to_1_ihd$soma_normal_coef>0][order(overlap_1_to_1_ihd$soma_normal_p[overlap_1_to_1_ihd$soma_normal_coef>0])[1:4]]
soma_non_normal_4_pos <- overlap_1_to_1_ihd$olink_id[overlap_1_to_1_ihd$soma_non_normal_coef>0][order(overlap_1_to_1_ihd$soma_non_normal_p[overlap_1_to_1_ihd$soma_non_normal_coef>0])[1:4]]
top_pos <- unique(c(olink_4_pos,soma_normal_4_pos,soma_non_normal_4_pos))


olink_4_neg <- overlap_1_to_1_ihd$olink_id[overlap_1_to_1_ihd$olink_coef<0][order(overlap_1_to_1_ihd$olink_p[overlap_1_to_1_ihd$olink_coef<0])[1:4]]
soma_normal_4_neg <- overlap_1_to_1_ihd$olink_id[overlap_1_to_1_ihd$soma_normal_coef<0][order(overlap_1_to_1_ihd$soma_normal_p[overlap_1_to_1_ihd$soma_normal_coef<0])[1:4]]
soma_non_normal_4_neg <- overlap_1_to_1_ihd$olink_id[overlap_1_to_1_ihd$soma_non_normal_coef<0][order(overlap_1_to_1_ihd$soma_non_normal_p[overlap_1_to_1_ihd$soma_non_normal_coef<0])[1:4]]
top_neg <- unique(c(olink_4_neg,soma_normal_4_neg,soma_non_normal_4_neg))

top <- overlap_1_to_1_ihd[overlap_1_to_1_ihd$olink_id %in% unique(c(top_pos,top_neg)),]

# olink

overlap_1_to_1_ihd$olink_sig_fdr <- F
overlap_1_to_1_ihd$olink_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr<0.05 & sign(overlap_1_to_1_ihd$olink_coef)>0] <- "Positive"
overlap_1_to_1_ihd$olink_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr<0.05 & sign(overlap_1_to_1_ihd$olink_coef)<0] <- "Inverse"
overlap_1_to_1_ihd$olink_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr>=0.05] <- "n.s."
overlap_1_to_1_ihd$olink_sig_fdr <- factor(overlap_1_to_1_ihd$olink_sig_fdr,levels=c("Positive","Inverse","n.s."))
table(overlap_1_to_1_ihd$olink_sig_fdr)
n_sig <- sum(overlap_1_to_1_ihd$olink_sig_fdr!="n.s.")

olink_threshold <- max(overlap_1_to_1_ihd$olink_p[overlap_1_to_1_ihd$olink_p_fdr<0.05])

volcano_olink <- ggplot(overlap_1_to_1_ihd, aes(x=olink_coef, y=-log(olink_p,10))) +
  geom_point(aes(color = olink_sig_fdr), size = 0.5) +
  ggtitle("Olink", subtitle = paste("\nn = ",n_sig,"/1694",sep="")) +
  xlab("Hazard ratio") + 
  ylab(expression("-log"[10]*"(p)")) +
  ylim(0,15) +
  scale_x_continuous(labels=c("0.75","1.00","1.25","1.50"), breaks=c(log(0.75),log(1),log(1.25),log(1.5)), limits=c(log(0.7),log(1.8))) +
  labs(color='Association') +
  scale_color_manual(values=c(hex,"gray")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  # geom_segment(aes(x=olink_x_loc, y=150, xend=olink_x_loc, yend=155), 
  #              arrow = arrow(length=unit(2, 'pt')), color=hex[1],lwd=1) +
  geom_text_repel(data = top,
                  mapping = aes(x=olink_coef, y=-log(olink_p,10), label=olink_id),
                  size = 2.5) +
  theme(legend.position = "none",text = element_text(size = 14),plot.subtitle=element_text(hjust=0.5)) +
  geom_hline(aes(yintercept=-log(olink_threshold,10)),linetype="dashed")

volcano_olink

# soma normal

overlap_1_to_1_ihd$soma_normal_sig_fdr <- F
overlap_1_to_1_ihd$soma_normal_sig_fdr[overlap_1_to_1_ihd$soma_normal_p_fdr<0.05 & sign(overlap_1_to_1_ihd$soma_normal_coef)>0] <- "Positive"
overlap_1_to_1_ihd$soma_normal_sig_fdr[overlap_1_to_1_ihd$soma_normal_p_fdr<0.05 & sign(overlap_1_to_1_ihd$soma_normal_coef)<0] <- "Inverse"
overlap_1_to_1_ihd$soma_normal_sig_fdr[overlap_1_to_1_ihd$soma_normal_p_fdr>=0.05] <- "n.s."
overlap_1_to_1_ihd$soma_normal_sig_fdr <- factor(overlap_1_to_1_ihd$soma_normal_sig_fdr,levels=c("Positive","Inverse","n.s."))
table(overlap_1_to_1_ihd$soma_normal_sig_fdr)
n_sig <- sum(overlap_1_to_1_ihd$soma_normal_sig_fdr!="n.s.")

soma_normal_threshold <- max(overlap_1_to_1_ihd$soma_normal_p[overlap_1_to_1_ihd$soma_normal_p_fdr<0.05])

volcano_soma_normal <- ggplot(overlap_1_to_1_ihd, aes(x=soma_normal_coef, y=-log(soma_normal_p,10))) +
  geom_point(aes(color = soma_normal_sig_fdr), size = 0.5) +
  ggtitle("SomaScan-ANML", subtitle = paste("\nn = ",n_sig,"/1694",sep="")) +
  xlab("Hazard ratio") + 
  ylab(expression("-log"[10]*"(p)")) +
  ylim(0,15) +
  scale_x_continuous(labels=c("0.75","1.00","1.25","1.50"), breaks=c(log(0.75),log(1),log(1.25),log(1.5)), limits=c(log(0.7),log(1.8))) +
  labs(color='Association') +
  scale_color_manual(values=c(hex,"gray")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  # geom_segment(aes(x=soma_normal_x_loc, y=150, xend=soma_normal_x_loc, yend=155), 
  #              arrow = arrow(length=unit(2, 'pt')), color=hex[1],lwd=1) +
  geom_text_repel(data = top,
                  mapping = aes(x=soma_normal_coef, y=-log(soma_normal_p,10), label=olink_id),
                  size = 2.5) +
  theme(legend.position = "none",text = element_text(size = 14),plot.subtitle=element_text(hjust=0.5)) +
  geom_hline(aes(yintercept=-log(soma_normal_threshold,10)),linetype="dashed")

volcano_soma_normal

# soma non normal

overlap_1_to_1_ihd$soma_non_normal_sig_fdr <- F
overlap_1_to_1_ihd$soma_non_normal_sig_fdr[overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05 & sign(overlap_1_to_1_ihd$soma_non_normal_coef)>0] <- "Positive"
overlap_1_to_1_ihd$soma_non_normal_sig_fdr[overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05 & sign(overlap_1_to_1_ihd$soma_non_normal_coef)<0] <- "Inverse"
overlap_1_to_1_ihd$soma_non_normal_sig_fdr[overlap_1_to_1_ihd$soma_non_normal_p_fdr>=0.05] <- "n.s."
overlap_1_to_1_ihd$soma_non_normal_sig_fdr <- factor(overlap_1_to_1_ihd$soma_non_normal_sig_fdr,levels=c("Positive","Inverse","n.s."))
table(overlap_1_to_1_ihd$soma_non_normal_sig_fdr)
n_sig <- sum(overlap_1_to_1_ihd$soma_non_normal_sig_fdr!="n.s.")

soma_non_normal_threshold <- max(overlap_1_to_1_ihd$soma_non_normal_p[overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05])

volcano_soma_non_normal <- ggplot(overlap_1_to_1_ihd, aes(x=soma_non_normal_coef, y=-log(soma_non_normal_p,10))) +
  geom_point(aes(color = soma_non_normal_sig_fdr), size = 0.5) +
  ggtitle("SomaScan-non-ANML", subtitle = paste("\nn = ",n_sig,"/1694",sep="")) +
  xlab("Hazard ratio") + 
  ylab(expression("-log"[10]*"(p)")) +
  ylim(0,15) +
  scale_x_continuous(labels=c("0.75","1.00","1.25","1.50"), breaks=c(log(0.75),log(1),log(1.25),log(1.5)), limits=c(log(0.7),log(1.8))) +
  labs(color='Association') +
  scale_color_manual(values=c(hex,"gray")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  # geom_segment(aes(x=soma_non_normal_x_loc, y=150, xend=soma_non_normal_x_loc, yend=155), 
  #              arrow = arrow(length=unit(2, 'pt')), color=hex[1],lwd=1) +
  geom_text_repel(data = top,
                  mapping = aes(x=soma_non_normal_coef, y=-log(soma_non_normal_p,10), label=olink_id),
                  size = 2.5) +
  theme(legend.position = "none",text = element_text(size = 14),plot.subtitle=element_text(hjust=0.5)) +
  geom_hline(aes(yintercept=-log(soma_non_normal_threshold,10)),linetype="dashed")

volcano_soma_non_normal

volcano_combined <- ggarrange(volcano_olink,volcano_soma_normal,volcano_soma_non_normal,ncol=3,nrow=1)

ggsave("",volcano_combined,width=14,height=6)

# for text 

volcano_soma_non_normal <- ggplot(overlap_1_to_1_ihd, aes(x=soma_non_normal_coef, y=-log(soma_non_normal_p,10))) +
  geom_point(aes(color = soma_non_normal_sig_fdr), size = 0.5) +
  ggtitle("SomaScan", subtitle = paste("\nn = ",n_sig,"/1694",sep="")) +
  xlab("Hazard ratio") + 
  ylab(expression("-log"[10]*"(p)")) +
  ylim(0,15) +
  scale_x_continuous(labels=c("0.75","1.00","1.25","1.50"), breaks=c(log(0.75),log(1),log(1.25),log(1.5)), limits=c(log(0.7),log(1.8))) +
  labs(color='Association') +
  scale_color_manual(values=c(hex,"gray")) +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  # geom_segment(aes(x=soma_non_normal_x_loc, y=150, xend=soma_non_normal_x_loc, yend=155), 
  #              arrow = arrow(length=unit(2, 'pt')), color=hex[1],lwd=1) +
  geom_text_repel(data = top,
                  mapping = aes(x=soma_non_normal_coef, y=-log(soma_non_normal_p,10), label=olink_id),
                  size = 2.5) +
  theme(legend.position = "none",text = element_text(size = 14),plot.subtitle=element_text(hjust=0.5)) +
  geom_hline(aes(yintercept=-log(soma_non_normal_threshold,10)),linetype="dashed")

volcano_soma_non_normal

volcano_combined_figure <- ggarrange(volcano_olink,volcano_soma_non_normal,ncol=2,nrow=1)

ggsave("",volcano_combined_figure,width=9,height=6)

# olink vs soma normal

c1 <- as.numeric(sum(overlap_1_to_1_ihd$olink_p_fdr<0.05))
c2 <- as.numeric(sum(overlap_1_to_1_ihd$soma_normal_p_fdr<0.05))
c3 <- as.numeric(sum(overlap_1_to_1_ihd$olink_p_fdr<0.05 & overlap_1_to_1_ihd$soma_normal_p_fdr<0.05))
c4 <- as.numeric(sum(overlap_1_to_1_ihd$olink_p_fdr<0.05 & overlap_1_to_1_ihd$soma_normal_p_fdr<0.05 & sign(overlap_1_to_1_ihd$olink_coef)==sign(overlap_1_to_1_ihd$soma_normal_coef)))

venn_ihd_normal <- venn.diagram(
  list(olink_sig=1:c1, soma_normal_sig=(c1-c3+1):(c1-c3+c2)),
  filename = NULL,
  disable.logging =T,
  fontfamily = "sans",	cat.fontfamily = "sans",
  category.names = c("Olink", "SomaScan-\nANML"),
  height = 1600, 
  width = 1600, 
  cat.pos = c(-25, 25),
  cat.dist = c(0.04, 0.06),
  fill=cols[c(1,2)])

grid.newpage()
grid.draw(venn_ihd_normal)

venn_ihd_normal[[7]]$label <- paste(venn_ihd_normal[[7]]$label, "\n(", c4, " shared)", sep="")  

grid.newpage()
grid.draw(venn_ihd_normal)

png("",width = 1600, height = 1600, pointsize=70)
grid.draw(venn_ihd_normal)
dev.off()

overlap_1_to_1_ihd$olink_soma_normal_sig_fdr <- NA
overlap_1_to_1_ihd$olink_soma_normal_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr<0.05 & overlap_1_to_1_ihd$soma_normal_p_fdr<0.05 &
                                                      sign(overlap_1_to_1_ihd$olink_coef)==sign(overlap_1_to_1_ihd$soma_normal_coef)] <- "Shared"
overlap_1_to_1_ihd$olink_soma_normal_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr<0.05 & overlap_1_to_1_ihd$soma_normal_p_fdr<0.05 &
                                                      sign(overlap_1_to_1_ihd$olink_coef)!=sign(overlap_1_to_1_ihd$soma_normal_coef)] <- "Opposite direction"
overlap_1_to_1_ihd$olink_soma_normal_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr<0.05 & overlap_1_to_1_ihd$soma_normal_p_fdr>=0.05] <- "Olink-specific"
overlap_1_to_1_ihd$olink_soma_normal_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr>=0.05 & overlap_1_to_1_ihd$soma_normal_p_fdr<0.05] <- "SomaScan-specific"
overlap_1_to_1_ihd$olink_soma_normal_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr>=0.05 & overlap_1_to_1_ihd$soma_normal_p_fdr>=0.05] <- "Neither"
overlap_1_to_1_ihd$olink_soma_normal_sig_fdr <- factor(overlap_1_to_1_ihd$olink_soma_normal_sig_fdr,levels = c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
table(overlap_1_to_1_ihd$olink_soma_normal_sig_fdr)

cor.test(overlap_1_to_1_ihd$olink_coef, overlap_1_to_1_ihd$soma_normal_coef)

cor.test(overlap_1_to_1_ihd$olink_coef[overlap_1_to_1_ihd$olink_soma_normal_sig_fdr=="Shared"], 
         overlap_1_to_1_ihd$soma_normal_coef[overlap_1_to_1_ihd$olink_soma_normal_sig_fdr=="Shared"])

r_all_normal <- format(round(cor.test(overlap_1_to_1_ihd$olink_coef, overlap_1_to_1_ihd$soma_normal_coef)$estimate,2),nsmall=2)
r_shared_normal <- format(round(cor.test(overlap_1_to_1_ihd$olink_coef[overlap_1_to_1_ihd$olink_soma_normal_sig_fdr=="Shared"], 
                                         overlap_1_to_1_ihd$soma_normal_coef[overlap_1_to_1_ihd$olink_soma_normal_sig_fdr=="Shared"])$estimate,2),nsmall=2)
r_coloc_normal <- format(round(cor.test(overlap_1_to_1_ihd$olink_coef[overlap_1_to_1_ihd$coloc_normal_cis_tf==T], 
                                        overlap_1_to_1_ihd$soma_normal_coef[overlap_1_to_1_ihd$coloc_normal_cis_tf==T])$estimate,2),nsmall=2)

sum(overlap_1_to_1_ihd$coloc_normal_cis_tf==T)
sum(overlap_1_to_1_ihd$olink_soma_normal_sig_fdr=="Shared" & overlap_1_to_1_ihd$coloc_normal_cis_tf==T)
sum(overlap_1_to_1_ihd$olink_soma_normal_sig_fdr=="Shared" & overlap_1_to_1_ihd$coloc_normal_cis_tf==T)/sum(overlap_1_to_1_ihd$coloc_normal_cis_tf==T)

scatter_normal <- ggplot(overlap_1_to_1_ihd, aes(x=olink_coef, y=soma_normal_coef, color=olink_soma_normal_sig_fdr, stroke=NA)) +
  geom_errorbar(aes(xmin = olink_lci, xmax = olink_hci),color="gray80") +
  geom_errorbar(aes(ymin = soma_normal_lci, ymax = soma_normal_hci),color="gray80") +
  geom_point(size=1.5,alpha=0.5) +
  xlim(-0.6,0.6) +
  ylim(-0.6,0.6) +
  scale_color_manual(values=c("black","gray80","gray80","gray80","gray80")) +
  geom_abline(linetype = "dashed",color="gray60") +
  geom_vline(xintercept=0,linetype = "dashed",color="gray60") +
  geom_hline(yintercept=0,linetype = "dashed",color="gray60") +
  ggtitle("Olink vs SomaScan-ANML", subtitle ="") +
  xlab("Olink") +
  ylab("SomaScan-ANML") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  theme(legend.position = "none") +
  geom_text_repel(data = top,
                  mapping = aes(x=olink_coef, y=soma_normal_coef, label=olink_id),
                  size = 2, color="black") +
  annotate("text", x = -0.6, y = 0.6, label = paste("r =", r_all_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.6, y = 0.55, label = paste0("r = ", r_shared_normal, " (shared, n = ", c4,")"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.6, y = 0.5, label = paste0("r = ", r_coloc_normal, " (coloc cis, n = ", n_coloc_normal, ")"), size = 12/.pt, hjust = 0)

scatter_normal

ggsave("",scatter_normal,width=6,height=6,bg = "white")


# olink vs soma non-normal

c1 <- as.numeric(sum(overlap_1_to_1_ihd$olink_p_fdr<0.05))
c2 <- as.numeric(sum(overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05))
c3 <- as.numeric(sum(overlap_1_to_1_ihd$olink_p_fdr<0.05 & overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05))
c4 <- as.numeric(sum(overlap_1_to_1_ihd$olink_p_fdr<0.05 & overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05 & sign(overlap_1_to_1_ihd$olink_coef)==sign(overlap_1_to_1_ihd$soma_non_normal_coef)))

venn_ihd_non_normal <- venn.diagram(
  list(olink_sig=1:c1, soma_non_normal_sig=(c1-c3+1):(c1-c3+c2)),
  filename = NULL,
  disable.logging =T,
  fontfamily = "sans",	cat.fontfamily = "sans",
  category.names = c("Olink", "SomaScan-\nnon-ANML"),
  height = 1600, 
  width = 1600, 
  cat.pos = c(-25, 19),
  cat.dist = c(0.04, 0.06),
  fill=cols[c(1,3)])

grid.newpage()
grid.draw(venn_ihd_non_normal)

venn_ihd_non_normal[[7]]$label <- paste(venn_ihd_non_normal[[7]]$label, "\n(", c4, " shared)", sep="")  

grid.newpage()
grid.draw(venn_ihd_non_normal)

png("",width = 1600, height = 1600, pointsize=70)
grid.draw(venn_ihd_non_normal)
dev.off()


overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr <- NA
overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr<0.05 & overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05 &
                                                          sign(overlap_1_to_1_ihd$olink_coef)==sign(overlap_1_to_1_ihd$soma_non_normal_coef)] <- "Shared"
overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr<0.05 & overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05 &
                                                          sign(overlap_1_to_1_ihd$olink_coef)!=sign(overlap_1_to_1_ihd$soma_non_normal_coef)] <- "Opposite direction"
overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr<0.05 & overlap_1_to_1_ihd$soma_non_normal_p_fdr>=0.05] <- "Olink-specific"
overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr>=0.05 & overlap_1_to_1_ihd$soma_non_normal_p_fdr<0.05] <- "SomaScan-specific"
overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr[overlap_1_to_1_ihd$olink_p_fdr>=0.05 & overlap_1_to_1_ihd$soma_non_normal_p_fdr>=0.05] <- "Neither"
overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr <- factor(overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr,levels = c("Shared","Opposite direction","Olink-specific","SomaScan-specific","Neither"))
table(overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr)

r_all_non_normal <- format(round(cor.test(overlap_1_to_1_ihd$olink_coef, overlap_1_to_1_ihd$soma_non_normal_coef)$estimate,2),nsmall=2)
r_shared_non_normal <- format(round(cor.test(overlap_1_to_1_ihd$olink_coef[overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr=="Shared"], 
                                             overlap_1_to_1_ihd$soma_non_normal_coef[overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr=="Shared"])$estimate,2),nsmall=2)
r_coloc_non_normal <- format(round(cor.test(overlap_1_to_1_ihd$olink_coef[overlap_1_to_1_ihd$coloc_non_normal_cis_tf==T], 
                                            overlap_1_to_1_ihd$soma_non_normal_coef[overlap_1_to_1_ihd$coloc_non_normal_cis_tf==T])$estimate,2),nsmall=2)

sum(overlap_1_to_1_ihd$coloc_non_normal_cis_tf==T)
sum(overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr=="Shared" & overlap_1_to_1_ihd$coloc_non_normal_cis_tf==T)
sum(overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr=="Shared" & overlap_1_to_1_ihd$coloc_non_normal_cis_tf==T)/sum(overlap_1_to_1_ihd$coloc_non_normal_cis_tf==T)

scatter_non_normal <- ggplot(overlap_1_to_1_ihd, aes(x=olink_coef, y=soma_non_normal_coef, color=olink_soma_non_normal_sig_fdr, stroke=NA)) +
  geom_errorbar(aes(xmin = olink_lci, xmax = olink_hci),color="gray80") +
  geom_errorbar(aes(ymin = soma_non_normal_lci, ymax = soma_non_normal_hci),color="gray80") +
  geom_point(size=1.5,alpha=0.5) +
  xlim(-0.6,0.6) +
  ylim(-0.6,0.6) +
  scale_color_manual(values=c("black","gray80","gray80","gray80","gray80")) +
  geom_abline(linetype = "dashed",color="gray60") +
  geom_vline(xintercept=0,linetype = "dashed",color="gray60") +
  geom_hline(yintercept=0,linetype = "dashed",color="gray60") +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") +
  xlab("Olink") +
  ylab("SomaScan-non-ANML") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  theme(legend.position = "none") +
  geom_text_repel(data = top,
                  mapping = aes(x=olink_coef, y=soma_non_normal_coef, label=olink_id),
                  size = 2, color="black") +
  annotate("text", x = -0.6, y = 0.6, label = paste("r =", r_all_non_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.6, y = 0.55, label = paste0("r = ", r_shared_non_normal, " (shared, n = ", c4,")"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.6, y = 0.5, label = paste0("r = ", r_coloc_non_normal, " (coloc cis, n = ", n_coloc_non_normal, ")"), size = 12/.pt, hjust = 0)

scatter_non_normal

ggsave("",scatter_non_normal,width=6,height=6,bg = "white")



### histogram plot against rho

# all categories

hex <- hue_pal()(10)

show_col(hex)

overlap_1_to_1_ihd$olink_soma_normal_sig_fdr <- factor(overlap_1_to_1_ihd$olink_soma_normal_sig_fdr,levels = c("Olink-specific","SomaScan-specific","Shared","Opposite direction","Neither"))
overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr <- factor(overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr,levels = c("Olink-specific","Shared","SomaScan-specific","Opposite direction","Neither"))

median_rho_all_normal <- format(round(median(overlap_1_to_1_ihd$rho_olink_soma_normal),2),nsmall=2)
median_rho_shared_normal <- format(round(median(overlap_1_to_1_ihd$rho_olink_soma_normal[overlap_1_to_1_ihd$olink_soma_normal_sig_fdr=="Shared"]),2),nsmall=2)

hist_ihd_normal <- ggplot(overlap_1_to_1_ihd, aes(x=rho_olink_soma_normal, fill=olink_soma_normal_sig_fdr)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black") +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  ggtitle("Olink vs SomaScan-ANML", subtitle ="") +
  xlab("Spearman's rho") + 
  ylab("Frequency") +
  labs(fill="Concordance") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_fill_manual(values=c(hex[c(3,5,7,10)],"gray80")) +
  theme(legend.position = "none") +
  annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_shared_normal, " (shared, n = ", sum(overlap_1_to_1_ihd$olink_soma_normal_sig_fdr=="Shared"), ")"), size = 12/.pt, hjust = 0)

hist_ihd_normal

median_rho_all_non_normal <- format(round(median(overlap_1_to_1_ihd$rho_olink_soma_non_normal),2),nsmall=2)
median_rho_shared_non_normal <- format(round(median(overlap_1_to_1_ihd$rho_olink_soma_non_normal[overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr=="Shared"]),2),nsmall=2)

hist_ihd_non_normal <- ggplot(overlap_1_to_1_ihd, aes(x=rho_olink_soma_non_normal, fill=olink_soma_non_normal_sig_fdr)) + 
  geom_histogram(binwidth = 0.05,boundary=0,color="black") +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") +
  xlab("Spearman's rho") + 
  ylab("Frequency") +
  labs(fill="Concordance") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_fill_manual(values=c(hex[c(3,5,7)],"gray80")) +
  theme(legend.position = "none") +
  annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_non_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_shared_non_normal, " (shared, n = ", sum(overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr=="Shared"), ")"), size = 12/.pt, hjust = 0)

hist_ihd_non_normal

legend <- as_ggplot(get_legend(ggplot(overlap_1_to_1_ihd, aes(x=rho_olink_soma_normal, fill=olink_soma_normal_sig_fdr)) + 
                                 geom_histogram(boundary=0,color="black") +
                                 labs(fill="Concordance") +
                                 scale_fill_manual(values=c(hex[c(3,5,7,10)],"gray80"))))

hist_ihd <- grid.arrange(hist_ihd_normal,hist_ihd_non_normal,legend,widths=c(3,3,1))

ggsave("",hist_ihd,width=14,height=6,bg = "white")




## plot shared with shade

overlap_1_to_1_ihd$olink_soma_normal_shared <- "No"
overlap_1_to_1_ihd$olink_soma_normal_shared[overlap_1_to_1_ihd$olink_soma_normal_sig_fdr=="Shared"] <- "Yes"
table(overlap_1_to_1_ihd$olink_soma_normal_shared)

hist_ihd_shared_normal <- ggplot(overlap_1_to_1_ihd, aes(x=rho_olink_soma_normal,fill=olink_soma_normal_shared,pattern=olink_soma_normal_shared)) + 
  geom_histogram_pattern(binwidth = 0.05,
                         boundary=0,
                         pattern_fill = "black", 
                         fill    = 'white',
                         colour  = 'black',
                         pattern_spacing = 0.015,
                         pattern_frequency = 5, 
                         pattern_angle = 45,
                         pattern_density=0.01) +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  ggtitle("Olink vs SomaScan-ANML", subtitle ="") +
  xlab("Spearman's rho") + 
  ylab("Frequency") +
  labs(fill="Shared") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_pattern_manual(values=c("none","stripe")) +
  theme(legend.position = "none") +
  annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_shared_normal, " (shared, n = ", sum(overlap_1_to_1_ihd$olink_soma_normal_sig_fdr=="Shared"), ")"), size = 12/.pt, hjust = 0)

hist_ihd_shared_normal

overlap_1_to_1_ihd$olink_soma_non_normal_shared <- "No"
overlap_1_to_1_ihd$olink_soma_non_normal_shared[overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr=="Shared"] <- "Yes"

hist_ihd_shared_non_normal <- ggplot(overlap_1_to_1_ihd, aes(x=rho_olink_soma_non_normal,fill=olink_soma_non_normal_shared,pattern=olink_soma_non_normal_shared)) + 
  geom_histogram_pattern(binwidth = 0.05,
                         boundary=0,
                         pattern_fill = "black", 
                         fill    = 'white',
                         colour  = 'black',
                         pattern_spacing = 0.015,
                         pattern_frequency = 5, 
                         pattern_angle = 45,
                         pattern_density=0.01) +
  scale_x_continuous(breaks = seq(-0.4, 1, 0.2), lim = c(-0.4, 1), labels = label_number(accuracy = 0.1)) +
  ylim(0,300) +
  ggtitle("Olink vs SomaScan-non-ANML", subtitle ="") +
  xlab("Spearman's rho") + 
  ylab("Frequency") +
  labs(fill="Shared") +
  theme_ckb() + theme(axis.line = element_line(),plot.background = element_rect(fill = "white", colour = NA)) +
  theme(text = element_text(size = 12)) +
  scale_pattern_manual(values=c("none","stripe")) +
  theme(legend.position = "none") +
  annotate("text", x = -0.4, y = 300, label = paste("Median rho =", median_rho_all_non_normal, "(all, n = 1694)"), size = 12/.pt, hjust = 0) +
  annotate("text", x = -0.4, y = 285, label = paste0("Median rho = ", median_rho_shared_non_normal, " (shared, n = ", sum(overlap_1_to_1_ihd$olink_soma_non_normal_sig_fdr=="Shared"), ")"), size = 12/.pt, hjust = 0)

hist_ihd_shared_non_normal

hist_ihd_shared <- ggarrange(hist_ihd_shared_normal,hist_ihd_shared_non_normal,ncol=2,nrow=1)

ggsave("",hist_ihd_shared,width=14,height=6,bg = "white")

