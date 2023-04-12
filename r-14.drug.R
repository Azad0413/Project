rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./13_drug")){
  dir.create("./13_drug")
}
setwd("./13_drug")

lassogene <- read.delim2('../08_Lasso/lasso_genes.csv',header = F)
