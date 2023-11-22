rm(list = ls())

setwd("K:/kadoorie/Staff_Folders/BaihanW/proteomics/data")

library(UniProt.ws)

overlap <- read.csv("overlap.csv")

# define fields

field <- c("gene_names","gene_primary","gene_synonym","gene_oln","gene_orf")

#### testing here ####

test <- queryUniProt(query = unique(overlap$uniprot_id)[1],fields=field)

label <- names(test)

#### testing here ####

annot <- as.data.frame(matrix(ncol = length(field)+1, nrow = length(unique(overlap$uniprot_id))))

names(annot) <- c("uniprot_id",label)

annot$uniprot_id <- unique(overlap$uniprot_id)

for (i in 1:nrow(annot)) {
  annot[i,c(2:ncol(annot))] <- queryUniProt(query = overlap$uniprot_id[i],fields=field)
}

saveRDS(annot, file = "annot_gene_names.RDS") 

write.csv(annot,"annot_gene_names.csv", quote=F, row.names=F)
