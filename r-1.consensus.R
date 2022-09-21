rm(list = ls())
# 一致性聚类----------
setwd("/data/nas1/luchunlin/project/BJTC-258")
if (! dir.exists("./01_consensus")){
  dir.create("./01_consensus")
}
setwd("./01_consensus")
library(lance)
library(tidyverse)
dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-258/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colname<-data.frame(sample=colnames(dat))
colname$sample<-gsub('.','-',colname$sample,fixed = T)
colnames(dat)<-colname$sample
#age<-read.csv(file = 'genage_human.csv')
##307
library(readxl)
age<-read_xlsx('细胞衰老相关基因集-279.xlsx',col_names = F)
age<-data.frame(symbol=age$...1)
dat.age<-dat[rownames(dat)%in%age$symbol,]%>%as.matrix()
library(ConsensusClusterPlus)
#sweep函数减去中位数进行标准化
df<-sweep(dat.age,1, apply(dat.age,1,median,na.rm=T))
maxK <-  4 #最多分成几组
results <- ConsensusClusterPlus(d = as.matrix(df),
                                maxK = maxK,
                                pItem = 0.8,
                                pFeature = 1,
                                reps = 1000,
                                clusterAlg = "pam",
                                distance = "pearson",
                                seed = 1,
                                innerLinkage = "complete",
                                finalLinkage = "complete",
                                corUse = "pairwise.complete.obs",
                                plot = 'pdf',
                                title = "consensus")


icl = calcICL(results,
              title="consensus",
              plot="pdf")

## 筛选最佳聚类数
### 一致性矩阵热图白色块最干净，尽量不掺杂蓝色
### 累积分布曲线下降的坡度最平缓
### delta area 曲线的肘部点横坐标
### 聚类一致性直方图又高又平均

### 可以使用PAC标准进行筛选
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK

for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i

# The optimal K
optK = Kvec[which.min(PAC)]
optK
## [3]
cluster<-results[[3]]$consensusClass
cluster<-as.data.frame(cluster)
cluster<-rownames_to_column(cluster,var = 'sample')
table(cluster$cluster)
cluster$cluster<-ifelse(cluster$cluster=='1','cluster1',
                        ifelse(cluster$cluster=='2','cluster2','cluster3'))
write.table(cluster,file = 'cluster.xls',sep = '\t',quote = F)
