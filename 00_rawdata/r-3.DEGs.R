rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/SJZZK-428-10/")
if (! dir.exists("./03_DEGs")){
  dir.create("./03_DEGs")
}
setwd("./03_DEGs")
library(lance)
library(tidyverse)
