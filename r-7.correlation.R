## 08 基因相关性分析-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/LLZK-505")
if (! dir.exists("./07_correlation")){
  dir.create("./07_correlation")
}
setwd("./07_correlation")
## 得到共同强相关基因
library(tidyverse)
library(lance)
dat<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/06_DEGs/normalized.counts.xls',row.names = 1)%>%lc.tableToNum()
colnames<-data.frame(sample=colnames(dat))
colnames$sample<-gsub('.','-',colnames$sample,fixed = T)
colnames(dat)<-colnames$sample
## 08-1 GRIM-19  NDUFS3------
sig.diff<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/06_DEGs/DEG_sig(GRIM-19&NDUFS3).xls',row.names = 1)
hub_exp<-dat[c('NDUFA13','NDUFS3'),]
sig_exp<-dat[rownames(sig.diff),]
hub_exp<-rbind(hub_exp,sig_exp)
hub_corr<-round(cor(t(hub_exp)),3)
hub_corr<-hub_corr[c(1:2),]
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
library(ggcorrplot)
library(corrplot)
hub_p.mat<-round(cor_pmat(t(hub_exp)),3)
hub_p.mat<-hub_p.mat[c(1:2),]
m<-hub_corr
p<-hub_p.mat
write.table(p,file = 'pvalue(GRIM-19&NDUFS3).xls',sep = '\t',row.names = T,quote = F)
p<-read_xlsx('/data/nas1/luchunlin/project/LLZK-505/07_correlation/p(GRIM-19&NDUFS3).xlsx')
p<-column_to_rownames(p,var = 'X1')
p<-t(p)%>%as.matrix()
m<-m[,colnames(p)]
library(Hmisc)
write.table(m,file = 'correlation(GRIM-19&NDUFS3).xls',sep = '\t',row.names = T)
## 筛选相关性最大的10个基因
gene1<-data.frame(symbol=c('NDUFS3','NDUFA13','LY6H','RP11-536O18.1','PCSK2','LINC00592','CTD-3247F14.2','FLG','RPL36AP29','RMRP'))
write.table(gene1,file = 'gene(GRIM-19&NDUFS3).xls',sep = '\t',quote = F,row.names = F)
m<-m[,c('NDUFS3','NDUFA13','LY6H','RP11-536O18.1','PCSK2','LINC00592','CTD-3247F14.2','FLG','RPL36AP29','RMRP')]
p<-p[,colnames(m)]
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

## 08-2 NDUFA4,LRPPRC-----
sig.diff<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/06_DEGs/DEG_sig(NDUFA4&LRPPRC).xls',row.names = 1)
hub_exp<-dat[c('NDUFA4','LRPPRC'),]
sig_exp<-dat[rownames(sig.diff),]
hub_exp<-rbind(hub_exp,sig_exp)
hub_corr<-round(cor(t(hub_exp)),3)
hub_corr<-hub_corr[c(1:2),]
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
library(ggcorrplot)
library(corrplot)
hub_p.mat<-round(cor_pmat(t(hub_exp)),3)
hub_p.mat<-hub_p.mat[c(1:2),]
m<-hub_corr
p<-hub_p.mat
write.table(p,file = 'pvalue(NDUFA4&LRPPRC).xls',sep = '\t',row.names = T,quote = F)
p<-read_xlsx('/data/nas1/luchunlin/project/LLZK-505/07_correlation/p.xlsx')
p<-column_to_rownames(p,var = 'X1')
p<-t(p)%>%as.matrix()
m<-m[,colnames(p)]
library(Hmisc)
write.table(m,file = 'correlation(NDUFA4&LRPPRC).xls',sep = '\t',row.names = T)
## 筛选相关性最大的10个基因
gene2<-data.frame(sample=c('RP11-270C12.3','RP11-20O24.1','AC007009.1','RP11-543P15.1','U3','AC098614.1','AC013404.1','DSC3','MTND1P23','MGAT5B'))
write.table(gene2,file = 'gene(NDUFA4&LRPPRC).xls',sep = '\t',quote = F,row.names = F)
m<-m[,c('RP11-270C12.3','RP11-20O24.1','AC007009.1','RP11-543P15.1','U3','AC098614.1','AC013404.1','DSC3','MTND1P23','MGAT5B')]
p<-p[,colnames(m)]
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



