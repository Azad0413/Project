rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/YQ444-8/")
if (! dir.exists("./03_cluster")){
  dir.create("./03_cluster")
}
setwd("./03_cluster")

dat<-read.delim2("/data/nas1/luchunlin/project/YQ444-8/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)                                          
dat<-dat[,-c(1:8)]

TRP<-read.delim2('/data/nas1/luchunlin/project/YQ444-8/02_uncox/univariate_cox_result_0.05.xls')
dat.cluster<-dat[rownames(dat)%in%rownames(TRP),]%>%as.matrix()
dat.cluster<-log2(dat.cluster+1)
library(NMF)
##计算最佳rank
rank <- nmf(dat.cluster,2:5,nrun=10,seed=12345,method = 'brunet') 
plot(rank)
### 2类
#使用最佳rank，或者最符合需求的rank
rank2 <- nmf(dat.cluster,2,nrun=10,seed=1234,method = 'brunet') # run parameters,rank=5
##混合系数矩阵
coefmap(rank2,
        annRow = NA,
        annCol = NA,
        main = "Metagene contributions in each sample",
        info = FALSE)
##基底矩阵
##作用：每一行显示主导的基底组分，即每一行有最高负载的基底组分。
basismap(rank2,
         annRow = NA,
         annCol = NA,
         main = "Metagenes",
         info = FALSE)
##一致性矩阵
##作用：基于指定rank评估聚类稳定性的方法，考虑由多个独立NMF运行结果计算得到的连接矩阵。
consensusmap(rank2,
             annRow = NA,
             annCol = NA,
             main = "Consensus matrix",
             info = FALSE)
##查看分组情况
group <- predict(rank2)
group <- as.data.frame(group)
group$group <- paste0('Cluster',group$group)
group$sample <- rownames(group)
group<- group[order(group$group),]
table(group$group)
head(group)
write.table(group,file = 'cluster.xls',sep = '\t',row.names = F,quote = F)
