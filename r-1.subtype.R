rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/HF-0103-1/")
if (! dir.exists("./01_subtype")){
  dir.create("./01_subtype")
}
setwd("./01_subtype")

library(tidyverse)
library(lance)
dat<-read.delim2("../00_rawdata/dat(GSE10846).xls", row.names = 1)%>% lc.tableToNum
dat <- dat^2
geneset<-read.delim2('gene.txt')
dat.cluster<-dat[rownames(dat)%in%geneset$id,]%>%as.matrix()

#dat.cluster<-log2(dat.cluster+1)
#cluster_exp<-log2(cluster_exp+0.0001)
cluster_exp<-as.matrix(dat.cluster)
library(ConsensusClusterPlus)
#sweep函数减去中位数进行标准化
df<-sweep(cluster_exp,1, apply(cluster_exp,1,median,na.rm=T))
maxK <-  9 #最多分成几组
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 1000,              # 抽样次数(一般1000或更多)
                                 pItem = 0.8,             # 抽样比例
                                 pFeature = 1,
                                 clusterAlg = "pam",      # 聚类方法
                                 seed = 123,
                                 title="consensus",
                                 innerLinkage="complete",
                                 plot="pdf")

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
cluster<-results[[2]]$consensusClass
cluster
group <- data.frame(cluster)%>%rownames_to_column(var = 'sample')
group$cluster <- ifelse(group$cluster=='1','cluster 1','cluster 2')
write.table(group,file = 'cluster.xls',sep = '\t',row.names = F,quote = F)


##survival---------
pca_exp<-t(cluster_exp)

cluster_dat<-as.data.frame(pca_exp)
#cluster_dat<-as.data.frame(t(dat_tpm[rownames(dat_tpm)%in%rownames(DE_emtmod),]))
#cluster_dat<-cluster_dat[rownames(cluster_dat)%in%survival$sample,]
cluster_dat$cluster<-as.vector(cluster,)
cluster_dat$sample<-rownames(cluster_dat)
survival<-read_tsv(file = '../00_rawdata/survival(GSE10846).xls')
cluster_dat<-merge(survival,cluster_dat,by='sample')
rownames(cluster_dat)<-cluster_dat$sample
cluster_dat<-cluster_dat[,-1]
library(survival)
library(survminer)
table(cluster_dat$cluster)
kmfit<-survfit(Surv(OS.time, OS) ~ cluster, data =  cluster_dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("Cluster 1","Cluster 2"),
                                      legend.title="cluster",
                                      title="Subtype KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
