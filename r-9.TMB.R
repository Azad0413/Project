rm(list = ls())
setwd("/data/nas1/luchunlin/project/HZ0301-3/")
if (! dir.exists("./09_TMB")){
  dir.create("./09_TMB")
}
setwd("./09_TMB")

## PART A TMB------
library(TCGAmutations)
maf <- TCGAmutations::tcga_load(study = "STAD")
library(maftools)
maf_sample <- data.frame(barcode = maf@clinical.data$Tumor_Sample_Barcode)
maf_sample$sample <- stringr::str_sub(maf_sample$barcode, 1, 16)
sample<-subset(maf_sample)$barcode
maf<-subsetMaf(maf,tsb = sample)
#maf<-subsetMaf(maf,tsb = sample,genes = model_gene$V1)

pdf(file = "01.oncoplot.pdf", height = 8, width = 10)
oncoplot(maf = maf, top = 20)
dev.off()

png(file = "01.oncoplot.png", family = "Times", height = 8, width = 10, units = "in", res = 600)
oncoplot(maf = maf, top = 20)
dev.off()

lasso.gene <- read.delim2('../06_Lasso/lasso_genes.csv',header = F)
##预后基因突变率--------
maf.mod<-subsetMaf(maf,tsb = sample,genes = lasso.gene$V1)

pdf(file = "02.mut_modgene.pdf", family = "Times", height = 5, width = 6)
plotmafSummary(maf = maf.mod, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

png(file = "02.mut_modgene.png", family = "Times", height = 5, width = 6, units = "in", res = 600)
plotmafSummary(maf = maf.mod, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()
