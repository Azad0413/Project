## 06 一致性聚类-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/LLZK-505")
if (! dir.exists("./05_Consensus")){
  dir.create("./05_Consensus")
}
setwd("./05_Consensus")
dat<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames<-data.frame(sample=colnames(dat))
colnames$sample<-gsub('.','-',colnames$sample,fixed = T)
colnames(dat)<-colnames$sample
## 06-1 GRIM-19  NDUFS3------
##根据异基因的表达情况，将TCGA-OV样本聚类。
## 提取表达矩阵
cluster_exp<-dat[c('NDUFA13','NDUFS3'),]
#cluster_exp<-log2(cluster_exp+1)
cluster_exp<-as.matrix(cluster_exp)
#cluster_exp<-cluster_exp[,colnames(cluster_exp)%in%survival$sample]
library(ConsensusClusterPlus)
#sweep函数减去中位数进行标准化
df<-sweep(cluster_exp,1, apply(cluster_exp,1,median,na.rm=T))
maxK <-  6 #最多分成几组
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 1000,              # 抽样次数(一般1000或更多)
                                 pItem = 0.8,             # 抽样比例
                                 pFeature = 1,
                                 clusterAlg = "pam",      # 聚类方法
                                 seed = 100,
                                 title="consensus(GRIM-19&NDUFS3)",
                                 innerLinkage="complete",
                                 plot="pdf")

icl = calcICL(results,
              title="consensus(GRIM-19&NDUFS3)",
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
## [2]
cluster1<-results[[2]]$consensusClass
cluster1<-as.data.frame(cluster1)
colnames(cluster1)<-'cluster'
write.table(cluster1,file = 'cluster(GRIM-19&NDUFS3).xls',sep = '\t',row.names = T,quote = F)
## 06-2 NDUFA4,LRPPRC-----
## 提取表达矩阵
cluster_exp<-dat[c('NDUFA4','LRPPRC'),]
#cluster_exp<-log2(cluster_exp+1)
cluster_exp<-as.matrix(cluster_exp)
#cluster_exp<-cluster_exp[,colnames(cluster_exp)%in%survival$sample]
library(ConsensusClusterPlus)
#sweep函数减去中位数进行标准化
df<-sweep(cluster_exp,1, apply(cluster_exp,1,median,na.rm=T))
maxK <-  6 #最多分成几组
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 1000,              # 抽样次数(一般1000或更多)
                                 pItem = 0.8,             # 抽样比例
                                 pFeature = 1,
                                 clusterAlg = "pam",      # 聚类方法
                                 seed = 100,
                                 title="consensus(NDUFA4&LRPPRC)",
                                 innerLinkage="complete",
                                 plot="pdf")

icl = calcICL(results,
              title="consensus(NDUFA4&LRPPRC)",
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
## [2]
cluster2<-results[[2]]$consensusClass
cluster2<-as.data.frame(cluster2)
colnames(cluster2)<-'cluster'
write.table(cluster2,file = 'cluster(NDUFA4&LRPPRC).xls',sep = '\t',row.names = T,quote = F)

