rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-320")
if (! dir.exists("./09_cor")){
  dir.create("./09_cor")
}
setwd("./09_cor")
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-320/06_expression/hubgene.xls')

hubgene<-hubgene[-1,]%>%as.data.frame()
dat = read.delim2("/data/nas1/luchunlin/project/BJTC-320/00_rawdata/dat.final.xls", row.names = 1) %>% lc.tableToNum
library(ggcorrplot)
library(corrplot)
## 提取诊断基因表达矩阵
exp<-dat[hubgene$.,]
hub_corr<-round(cor(t(exp)),3)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
hub_p.mat<-round(cor_pmat(t(exp)),3)
col1 <- colorRampPalette(c("#4169E1","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF",
                           "#FDDBC7","#F4A582","#D6604D", "#B2182B","#CC0000","#990000"))
hub_corr_plot<-corrplot(hub_corr,
                        method = "circle",
                        is.corr = T,
                        type = "lower",
                        p.mat = hub_p.mat,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        col = col1(200))
write.table(hub_corr,file = 'hub_corr.xls',
            sep = '\t',
            row.names = T)
