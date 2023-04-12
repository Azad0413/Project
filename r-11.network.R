rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/BJTC-356/")
if (! dir.exists("./10_network")){
  dir.create("./10_network")
}
setwd("./10_network")

hubgene <- read.delim2('../05_quadrant/final.gene.xls')
