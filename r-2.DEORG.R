rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-357/")
if (! dir.exists("./02_DEORG")){
  dir.create("./02_DEORG")
}
setwd("./02_DEORG")

oxidative<-read_xlsx('oxidative+stress.xlsx')
sig_diff<-read.delim2("/data/nas1/luchunlin/project/BJTC-357/01_DEGs/DEG_sig.xls", row.names = 1)
dat = read.delim2("/data/nas1/luchunlin/project/BJTC-357/00_rawdata/dat.geo.xls", row.names = 1) %>% lc.tableToNum
#dat1<-read.delim2("/data/nas1/luchunlin/project/BJTC-357/00_rawdata/dat1.xls", row.names = 1) %>% lc.tableToNum
#dat<-dat[,colnames(dat1)]
dat.os<-dat[rownames(dat)%in%oxidative$Symbol,]
dat.df<-dat[rownames(sig_diff),]
### GSVA计算氧化应激相关评分------
# library(GSVA)
# dat2<-as.matrix(dat)
# oxidative$oxidative<-c(rep('oxidative',436))
# gene_list <- split(as.matrix(oxidative)[,1],
#                    oxidative[,2])
# oxstress_score = gsva(dat2,gene_list,
#                     method = "gsva",
#                     ssgsea.norm = TRUE,
#                     verbose = TRUE)
# write.table(oxstress_score,
#             file = "oxstress_score.xls",
#             sep = "\t",
#             quote = F)


## 氧化应激相关基因进行相关性分析（cor>0.4），
library(Hmisc)
nc<-t(rbind(dat.df,dat.os))
m=rcorr(nc)$r[1:nrow(dat.df),(ncol(nc)-length(rownames(dat.os))+1):ncol(nc)]
#m<-t(m)
p=rcorr(nc)$P[1:nrow(dat.df),(ncol(nc)-length(rownames(dat.os))+1):ncol(nc)]

### 
correlation<-paste(signif(m,3),"\n(",
                   signif(p,3),")",sep = '')
dim(correlation)=dim(m)
colnames(correlation)<-colnames(m)
rownames(correlation)<-rownames(m)
write.table(correlation,'correlation.xls',
            sep = '\t',
            row.names = T)

## 筛选
m[abs(m)<0.7]=NA
m<-as.data.frame(t(m))
m<-subset(m,rowSums(m,na.rm = T)>0.7)
m<-m[,-which(apply(m, 2, function(x)all(is.na(x))))]
### p值
p[p>0.05]=NA
p<-as.data.frame(t(p))
p<-subset(p,rowSums(m,na.rm = T)>=0)
#p<-p[,-which(apply(p,2,function(x)all(is.na(x))))]
p<-p[rownames(m),]
p<-p[,colnames(m)]

m<-as.matrix(m)
p<-as.matrix(p)
textMatrix = paste(signif(m, 3), "\n(",
                   signif(p, 3), ")", sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
write.table(correlation,'correlation.xls',
            sep = '\t',
            row.names = T)

DEORG<-data.frame(symbol=colnames(m))
write.table(DEORG,file = 'DEORG.xls',sep = '\t',row.names = F,quote = F)
m[is.na(m)]=0
p[is.na(p)]=1
library(dplyr)
library(dplyr)
tmp = matrix(case_when(p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
library(pheatmap)
source("modified_pheatmap.R")
pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0,
         fontsize_row = 5,
         fontsize_col = 5,
         fontsize = 10)


# lassogene<-read.csv(file = '../04_Lasso/lasso_genes.csv',header = F)
# lassogene$V1<-gsub('_','-',lassogene$V1,fixed = T)
# genecor<-textMatrix[,lassogene$V1]

