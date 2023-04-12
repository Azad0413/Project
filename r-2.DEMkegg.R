rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/GY0324-12/")
if (! dir.exists("./02_DEMkegg")){
  dir.create("./02_DEMkegg")
}
setwd("./02_DEMkegg")
