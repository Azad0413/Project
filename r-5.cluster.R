rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./05_cluster")){
  dir.create("./05_cluster")
}
setwd("./05_cluster")
dat.tcga<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
DEARG <- read.delim2('../04_DEARG/DEARG.xls')
survival<-read_tsv(file = '/data/nas1/luchunlin/TCGA_survival/TARGET-OS.survival.tsv')
survival<-survival[,-3]
##一致性聚类------
## 提取表达矩阵
cluster_exp<-dat.tcga[DEARG$symbol,]
#cluster_exp<-log2(cluster_exp+0.0001)
cluster_exp<-as.matrix(cluster_exp)
cluster_exp<-cluster_exp[,colnames(cluster_exp)%in%survival$sample]
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
##PCA
## OV样本主成分分析
pca_exp<-t(cluster_exp)
pca_exp <- log2(pca_exp+1)
pca1<-prcomp(pca_exp,center = TRUE,scale. = TRUE)
df1<-pca1$x  ## 提取PC score
df1 <- as.data.frame(df1)
summ1 <- summary(pca1)
summ1
# 提取主成分的方差贡献率,生成坐标轴标题
summ1 <- summary(pca1)
summ1
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
# 绘制PCA得分图
cluster<-as.data.frame(cluster)
cluster$cluster<-ifelse(cluster$cluster=='1','cluster 1','cluster 2')
library(ggplot2)
p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2,color = cluster$cluster))+
  stat_ellipse(aes(fill = cluster$cluster),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("purple","orange","pink"))+
  scale_colour_manual(values = c("purple","orange","pink"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.pca1
ggsave(p.pca1,filename = '04.pca.pdf',w=6,h=5)
ggsave(p.pca1,filename = '04.pca.png',w=6,h=5)
##survival---------
cluster_dat<-as.data.frame(pca_exp)
#cluster_dat<-as.data.frame(t(dat_tpm[rownames(dat_tpm)%in%rownames(DE_emtmod),]))
#cluster_dat<-cluster_dat[rownames(cluster_dat)%in%survival$sample,]
cluster_dat$cluster<-as.vector(cluster$cluster,)
cluster_dat$sample<-rownames(cluster_dat)
cluster_dat<-merge(survival,cluster_dat,by='sample')
rownames(cluster_dat)<-cluster_dat$sample
cluster_dat<-cluster_dat[,-1]

kmfit<-survfit(Surv(OS.time, OS) ~ cluster, data =  cluster_dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("cluster 1","cluster 2"),
                                      legend.title="cluster",
                                      title="Train KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
##clinical----------
### 画热图（结合临床信息）
clinical<-read_tsv(file ='/data/nas1/luchunlin/TCGA_phenotype/TARGET-OS.clinical.tsv.gz')
#clinical<-GDCquery_clinic(project = "TARGET-OS",type = "clinical")
phenotype<-data.frame(sample=clinical$sample_id,
                      gender=clinical$Gender,
                      age=clinical$`Age at Diagnosis in Days`/365,
                      race=clinical$Race,
                      metastasis=clinical$`Disease at diagnosis`)
phenotype <- phenotype[phenotype$sample%in%colnames(dat.tcga),]
phenotype$age <- round(phenotype$age,digits = 0)
table(phenotype$race)
#phenotype$race <- gsub('Unknown',NA,phenotype$race)
table(phenotype$metastasis)
phenotype$metastasis <- gsub('Metastatic (confirmed)','Metastatic',phenotype$metastasis,fixed = T)
phenotype$metastasis <- gsub('Non-metastatic (confirmed)','Non-metastatic',phenotype$metastasis,fixed = T)
phenotype$metastasis <- gsub('Non-metastatic (Confirmed)','Non-metastatic',phenotype$metastasis,fixed = T)
table(phenotype$age)
phenotype$age<-cut(phenotype$age,breaks = c(1,5,10,15,20,25,30,35),labels = c('1-5','5-10','10-15','15-20','20-25','25-30','30-35'))

heat.group <- merge(group,phenotype,by='sample')
heat.group <- merge(heat.group,survival,by='sample')
heat.group$OS <- ifelse(heat.group$OS==1,'Dead','Alive')
heat.group <- heat.group[order(heat.group$cluster),]
heat.dat <- dat.tcga[,heat.group$sample]
colnames(heat.group)
rownames(heat.group) <- colnames(heat.dat)
heat.group<-dplyr::select(heat.group,c('cluster','gender','age','race','metastasis','OS'))
heat.group$cluster <- as.factor(heat.group$cluster)
heat.group$age<-as.factor(heat.group$age)
heat.group$gender <- as.factor(heat.group$gender)
heat.group$race <- as.factor(heat.group$race)
heat.group$metastasis <- as.factor(heat.group$metastasis)
heat.group$OS <- as.factor(heat.group$OS)
rt_dat<-log2(heat.dat+1)
rt_dat <- rt_dat[DEARG$symbol,]
ann_colors<-list(
  cluster=c('cluster 1'='#FFB7DD','cluster 2'='#77DDFF'),
  gender=c('Male'='#FF8C00','Female'='#20B2AA'),
  age=c('1-5'='#33CCCC','5-10'='#CCCCFF','10-15'='#FF9966','15-20'='#FFCCFF','20-25'='#DDA0DD','25-30'='#FF69B4','30-35'='#40E0D0'),
  race=c('Asian'='#9ACD32','Black or African American'='#FF8C00','White'='blue','Unknown'='green'),
  metastasis=c('Metastatic'='pink','Non-metastatic'='purple'),
  OS=c('Dead'='red','Alive'='green')
)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
pheatmap(rt_dat,
         color = bluered(100),
         border_color = NA,
         annotation_col = heat.group,
         annotation_colors = ann_colors,
         labels_row = NULL,
         clustering_method = 'ward.D2',
         show_rownames = T,
         show_colnames = F,
         fontsize_col = 5,
         cluster_cols = F,
         cluster_rows = T)

## 卡方检验
mydata<-heat.group
##age
age.ka<-xtabs(~mydata$age+mydata$cluster,data = mydata)
chisq.test(age.ka)
# p-value = 0.8884
##gender
gender.ka<-xtabs(~mydata$gender+mydata$cluster,data = mydata)
chisq.test(gender.ka)
# p-value = 1
##race
race.ka<-xtabs(~mydata$race+mydata$cluster,data = mydata)
chisq.test(race.ka)
# p-value = 0.04959
##metastasis
metastasis.ka<-xtabs(~mydata$metastasis+mydata$cluster,data = mydata)
chisq.test(metastasis.ka)
# p-value = 0.9505
##OS
OS.ka<-xtabs(~mydata$OS+mydata$cluster,data = mydata)
chisq.test(OS.ka)
# p-value = 0.941

## 桑葚图

##GSVA
dat<-dat.tcga
library(GSVA)
library(GSEABase)
library(limma)
group2 <- group
cluster1.sample<-group$sample[which(group$cluster=='cluster 1')]
colnames(group2)
colnames(group2)<-c('id','lable')
group2$lable <- ifelse(group2$lable=='cluster 1','cluster1','cluster2')
gsva_exp<-log2(dat[,group2$id]+1)
all(colnames(gsva_exp) == group2$id)
dim(gsva_exp)
## 19712  85
# 分组
group_score <- group2$lable %>% as.factor()
design_score <- model.matrix(~0 + group_score)
rownames(design_score) <- colnames(gsva_exp)
colnames(design_score) <- levels(group_score)
compare_score <- makeContrasts("cluster1-cluster2", levels = design_score)
KEGG_ref <- getGmt("/data/nas1/luchunlin/pipeline/GSVA/c2.cp.kegg.v7.4.symbols.gmt")
es_KEGG <- gsva(as.matrix(gsva_exp), KEGG_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_KEGG <- as.data.frame(es_KEGG)
fit <- lmFit(es_KEGG, design_score)
fit2 <- contrasts.fit(fit ,compare_score)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)

logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > 0)
##27
write.table(allGeneSets,
            file = "GSVA.xls",
            quote = F,
            sep = "\t",
            row.names = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$P.Value),]
dim(DEGeneSets)
write.table(DEGeneSets,file = 'DEGSVA.xls',sep = '\t',row.names = T,quote = F)

## 热图
heat.gsva <- es_KEGG[rownames(DEGeneSets),]
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-group2
group_rt<-group_rt[order(group_rt$lable),]
rt<-heat.gsva[,group_rt$id]
group_rt<-data.frame(group_rt$lable)
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
#x<-log2(heat+1)
x<-t(scale(t(rt)))
ann_colors<-list(
  group = c(cluster1="#00CED1",cluster2="#F08080"))
pdf('07.heatmap.pdf',w=11,h=6)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = T,
         annotation_names_row = F)
dev.off()
png('07.heatmap.png',w=900,h=500)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = T,
         annotation_names_row = F)
dev.off()
