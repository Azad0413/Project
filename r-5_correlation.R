rm(list = ls())
# 相关性分析-----------
setwd("/data/nas1/luchunlin/project/BJTC-327")
if (! dir.exists("./05_correlation")){
  dir.create("./05_correlation")
}
setwd("./05_correlation")
library(ggcorrplot)
library(corrplot)
# AIH.PBC------
##将基因的矩阵提取出来
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/02_DEIRG/DEIRG.xls')
dat1<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/00_rawdata/AIH.PBC.fpkm.xls',row.names = 1)%>%lc.tableToNum()
df1<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/01_DEGs/DEG_sig(AIH vs.PBC).xls',row.names = 1)
dat1<-na.omit(dat1)
dat1<-log2(dat1+1)
hub_exp<-dat1[hubgene$.,]
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
png('correlation(HIV vs.PBC).png',width = 600,height = 600)
pdf('correlation(HIV vs.PBC).pdf',width = 9,height = 9)
hub_corr_plot<-corrplot(hub_corr,
                        method = "circle",
                        is.corr = T,
                        type = "lower",
                        p.mat = hub_p.mat,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        col = col1(200),
                        tl.col = 'black',
                        tl.offset = 0.5,
                        number.font = 2,
                        number.cex = 0.8)
dev.off()

# AIH.HBV------
##将基因的矩阵提取出来
dat2<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/00_rawdata/AIH.HBV.fpkm.xls',row.names = 1)%>%lc.tableToNum()
df2<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/01_DEGs/DEG_sig(AIH vs.HBV).xls',row.names = 1)
dat2<-na.omit(dat2)
dat2<-log2(dat2+1)
hub_exp2<-dat2[hubgene$.,]
hub_corr2<-round(cor(t(hub_exp2)),3)
write.table(hub_corr2,file = 'cor(AIH vs.HBV).xls',sep = '\t',row.names = T,quote = F)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
hub_p.mat2<-round(cor_pmat(t(hub_exp2)),3)
write.table(hub_p.mat2,file = 'pvalue(AIH vs.HBV).xls',sep = '\t',row.names = T,quote = F)

png('correlation(HIV vs.HBV).png',width = 600,height = 600)
pdf('correlation(HIV vs.HBV).pdf',width = 9,height = 9)
hub_corr_plot<-corrplot(hub_corr2,
                        method = "circle",
                        is.corr = T,
                        type = "lower",
                        p.mat = hub_p.mat,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        col = col1(200),
                        tl.col = 'black',
                        tl.offset = 0.5,
                        number.font = 2,
                        number.cex = 0.8)
dev.off()

