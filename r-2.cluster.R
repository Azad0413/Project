rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-406-12/")
if (! dir.exists("./02_cluster")){
  dir.create("./02_cluster")
}
setwd("./02_cluster")
library(tidyverse)
library(lance)
library(readxl)
dat<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)                                          
geneset<-read_xlsx('geneset.xlsx')
dat <- dat[,-c(1:59)]
dat.cluster<-dat[rownames(dat)%in%geneset$symbol,]%>%as.matrix()

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
## [2]

cluster<-results[[2]]$consensusClass
cluster
group <- data.frame(cluster)%>%rownames_to_column(var = 'sample')
group$cluster <- ifelse(group$cluster=='1','cluster 1','cluster 2')
write.table(group,file = 'cluster.xls',sep = '\t',row.names = F,quote = F)
# ##PCA------
# ## 肿瘤样本主成分分析
# pca_exp<-t(cluster_exp)
# pca_exp <- log2(pca_exp+1)
# pca1<-prcomp(pca_exp,center = TRUE,scale. = TRUE)
# df1<-pca1$x  ## 提取PC score
# df1 <- as.data.frame(df1)
# summ1 <- summary(pca1)
# summ1
# # 提取主成分的方差贡献率,生成坐标轴标题
# summ1 <- summary(pca1)
# summ1
# xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
# ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
# # 绘制PCA得分图
# cluster<-as.data.frame(cluster)
# cluster$cluster<-ifelse(cluster$cluster=='1','cluster 1','cluster 2')
# library(ggplot2)
# p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2,color = cluster$cluster))+
#   stat_ellipse(aes(fill = cluster$cluster),
#                type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
#   geom_point(size = 3.5)+
#   labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
#   guides(fill = "none")+
#   theme_bw()+
#   scale_fill_manual(values = c("purple","orange","pink"))+
#   scale_colour_manual(values = c("purple","orange","pink"))+
#   theme(plot.title = element_text(hjust = 0.5,size = 15),
#         axis.text = element_text(size = 11),axis.title = element_text(size = 13),
#         legend.text = element_text(size = 11),legend.title = element_text(size = 13),
#         plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
# p.pca1
# ggsave(p.pca1,filename = '04.pca.pdf',w=6,h=5)
# ggsave(p.pca1,filename = '04.pca.png',w=6,h=5)
##survival---------
cluster_dat<-as.data.frame(log2(t(cluster_exp+1)))
#cluster_dat<-as.data.frame(t(dat_tpm[rownames(dat_tpm)%in%rownames(DE_emtmod),]))
#cluster_dat<-cluster_dat[rownames(cluster_dat)%in%survival$sample,]
cluster_dat$cluster<-as.vector(cluster$cluster)
cluster_dat$sample<-rownames(cluster_dat)
survival<-read_tsv(file = '/data/nas1/luchunlin/TCGA_survival/TCGA-LUAD.survival.tsv')
survival<-survival[,-3]
cluster_dat<-merge(survival,cluster_dat,by='sample')
rownames(cluster_dat)<-cluster_dat$sample
cluster_dat<-cluster_dat[,-1]
library(survival)
library(survminer)
table(cluster$cluster)
kmfit<-survfit(Surv(OS.time, OS) ~ cluster, data =  cluster_dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("cluster 1","cluster 2"),
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

