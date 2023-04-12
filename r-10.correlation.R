rm(list = ls())
# 相关性分析-----------
setwd("/data/nas1/luchunlin/project/BJTC-317")
if (! dir.exists("./10_correlation")){
  dir.create("./10_correlation")
}
setwd("./10_correlation")
library(ggcorrplot)
library(corrplot)
dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-317/00_rawdata/dat.final.xls", row.names = 1)  %>% lc.tableToNum
hubgene <- read.delim2('/data/nas1/luchunlin/project/BJTC-317/06_Logistic/hubgene.xls')

hub_exp<-dat[hubgene$sybmol,]%>%log2()
group = read.delim2("/data/nas1/luchunlin/project/BJTC-317/00_rawdata/group.xls")
hub_corr<-round(cor(t(hub_exp)),3)
write.table(hub_corr,file = 'cor.xls',sep = '\t',row.names = T,quote = F)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
hub_p.mat<-round(cor_pmat(t(hub_exp)),3)
write.table(hub_p.mat,file = 'pvalue.xls',sep = '\t',row.names = T,quote = F)
# 筛选一下
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