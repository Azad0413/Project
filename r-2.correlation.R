rm(list = ls())
# 相关性分析-----------
setwd("/data/nas1/luchunlin/project/LZZK-503")
if (! dir.exists("./02_correlation")){
  dir.create("./02_correlation")
}
setwd("./02_correlation")
library(ggcorrplot)
library(corrplot)
df<-read.delim2('/data/nas1/luchunlin/project/LZZK-503/01_DEGs/DEG_sig.xls',row.names = 1)
dat<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames<-data.frame(sample=colnames(dat))
colnames$sample<-gsub('.','-',colnames$sample,fixed = T)
colnames(dat)<-colnames$sample
## 'ZEB2-AS1'本身就在差异基因里
df.exp<-dat[rownames(df),]
hub_corr<-round(cor(t(df.exp)),3)
write.table(hub_corr,file = 'cor.xls',sep = '\t',row.names = T,quote = F)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
hub_p.mat<-round(cor_pmat(t(df.exp)),3)
write.table(hub_p.mat,file = 'pvalue.xls',sep = '\t',row.names = T,quote = F)
## 筛选一下
hub_corr<-as.data.frame(hub_corr)
hub_corr<-hub_corr[order(abs(hub_corr$`ZEB2-AS1`),decreasing = T),]
hub_corr<-t(hub_corr)
hub_corr<-as.data.frame(hub_corr)
hub_corr<-hub_corr[order(abs(hub_corr$`ZEB2-AS1`),decreasing = T),]
hub_corr<-hub_corr[,c(1:7)]%>%as.data.frame()

#colnames(hub_corr)<-'ZEB2-AS1'

hub_corr<-hub_corr[c(1:7),]
hub_p.mat<-hub_p.mat[rownames(hub_corr),]%>%as.data.frame()
hub_p.mat<-hub_p.mat[,colnames(hub_corr)]
m<-t(as.matrix(hub_corr))
p<-t(as.matrix(hub_p.mat))

library(dplyr)
tmp = matrix(case_when(p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
source("modified_pheatmap.R")
pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0,
         fontsize_row = 12,
         fontsize_col = 12,
         fontsize = 10)

col1 <- colorRampPalette(c("#4169E1","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF",
                           "#FDDBC7","#F4A582","#D6604D", "#B2182B","#CC0000","#990000"))
hub_corr_plot<-corrplot(m,
                        method = "circle",
                        is.corr = T,
                        type = "lower",
                        p.mat = p,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        col = col1(200))

## 强相关基因 IFFO1 ZEB2 GAB3 MS4A7 GIMAP1 CD37
hubgene<-data.frame(symbol=c('IFFO1','ZEB2','GAB3','MS4A7','GIMAP1','CD37'))
write.table(hubgene,file = 'hubgene.xls',sep = '\t',quote = F,row.names = F)
