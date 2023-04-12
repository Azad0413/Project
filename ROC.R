# hubgenes <- rownames(hubgenes)
hubgenes <- c("CRY1","ARNTL","PER2","BHLHE41","CRY2","PER3","NPAS2","PER1")

# 09-基于训练集表达矩阵的ROC曲线验证 ----------------------------------------------------
roc_dat <- train
roc_dat <- t(roc_dat)
roc_dat <- as.data.frame(roc_dat)
roc_dat <- roc_dat[,hubgenes]
roc_dat$Sample <- rownames(roc_dat)
conNum <- pheno1$geo_accession[which(pheno1$group == "Control")]
roc_dat$Group <- ifelse(roc_dat$Sample %in% conNum,'Control','HF')
roc_dat <- roc_dat[,c(9:10,1:8)]
rownames(roc_dat) <- NULL
roc_dat$Group <- as.factor(roc_dat$Group)

## 绘制ROC曲线
# install.packages("ROCit")
# devtools::install_github('yikeshu0611/geomROC')
library(plotROC)
library(geomROC)
library(ROCit)
library(pROC)
## *1.CRY1
CRY1 <- rocit(score = roc_dat$CRY1,class = roc_dat$Group,method = 'binormal')
roc.CRY1 <- ggplot()+geom_roc(CRY1,color='red')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.695',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='CRY1',alpha=1,size=6)
roc.CRY1
roc_dat$outome[roc_dat$Group=='Control'] <- 1
roc_dat$outome[roc_dat$Group=='HF'] <- 0
CRY1 <- roc(roc_dat$Group,roc_dat$CRY1,levels=c('Control','HF'))
library(dplyr)
ci.auc(CRY1)
##95% CI: 0.6375-0.7556 (DeLong)
library(verification)
roc.area(as.numeric(as.vector(roc_dat$outome)),CRY1$predictor)
# $p.value
# 1.257963e-09
## *2.ARNTL
ARNTL <- rocit(score = roc_dat$ARNTL,class = roc_dat$Group,method = 'binormal')
roc.ARNTL <- ggplot()+geom_roc(ARNTL,color='red')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.716',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='ARNTL',alpha=1,size=6)
roc.ARNTL
ARNTL <- roc(roc_dat$Group,roc_dat$ARNTL,levels=c('Control','HF'))
ci.auc(ARNTL)
##95% CI: 0.645-0.7603 (DeLong)
roc.area(as.numeric(as.vector(roc_dat$outome)),ARNTL$predictor)
# [1] 3.982132e-10

## *3.PER2
PER2 <- rocit(score = roc_dat$PER2,class = roc_dat$Group,method = 'binormal')
roc.PER2 <- ggplot()+geom_roc(PER2,color='red')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.683',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='PER2',alpha=1,size=6)
roc.PER2
PER2 <- roc(roc_dat$Group,roc_dat$PER2,levels=c('Control','HF'))
ci.auc(PER2)
roc_dat$outome[roc_dat$Group=='Control'] <- 0
roc_dat$outome[roc_dat$Group=='HF'] <- 1
##95% CI:  0.6316-0.7503 (DeLong)
roc.area(as.numeric(as.vector(roc_dat$outome)),PER2$predictor)
# [1] 3.512646e-09


## *4.BHLHE41
BHLHE41 <- rocit(score = roc_dat$BHLHE41,class = roc_dat$Group,method = 'binormal')
roc.BHLHE41 <- ggplot()+geom_roc(BHLHE41,color='red')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.737',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='BHLHE41',alpha=1,size=6)
roc.BHLHE41
BHLHE41 <- roc(roc_dat$Group,roc_dat$BHLHE41,levels=c('Control','HF'))
ci.auc(BHLHE41)

roc_dat$outome[roc_dat$Group=='Control'] <- 0
roc_dat$outome[roc_dat$Group=='HF'] <- 1
##95% CI: 0.6722-0.7822 (DeLong)
roc.area(as.numeric(as.vector(roc_dat$outome)),BHLHE41$predictor)
# 2.783806e-12
## *5.CRY2
CRY2 <- rocit(score = roc_dat$CRY2,class = roc_dat$Group,method = 'binormal')
roc.CRY2 <- ggplot()+geom_roc(CRY2,color='red')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.779',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='CRY2',alpha=1,size=6)
roc.CRY2
CRY2 <- roc(roc_dat$Group,roc_dat$CRY2,levels=c('Control','HF'))
ci.auc(CRY2)
roc_dat$outome[roc_dat$Group=='Control'] <- 0
roc_dat$outome[roc_dat$Group=='HF'] <- 1
##95% CI:0.7321-0.8349 (DeLong)
roc.area(as.numeric(as.vector(roc_dat$outome)),CRY2$predictor)
#  4.070253e-18
## *6.PER3
PER3 <- rocit(score = roc_dat$PER3,class = roc_dat$Group,method = 'binormal')
roc.PER3 <- ggplot()+geom_roc(PER3,color='red')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.739',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='PER3',alpha=1,size=6)
roc.PER3
PER3 <- roc(roc_dat$Group,roc_dat$PER3,levels=c('Control','HF'))
ci.auc(PER3)
roc_dat$outome[roc_dat$Group=='Control'] <- 0
roc_dat$outome[roc_dat$Group=='HF'] <- 1
##95% CI:0.6785-0.7878 (DeLong)
roc.area(as.numeric(as.vector(roc_dat$outome)),PER3$predictor)
#  7.720133e-13

## *7.NPAS2
NPAS2 <- rocit(score = roc_dat$NPAS2,class = roc_dat$Group,method = 'binormal')
roc.NPAS2 <- ggplot()+geom_roc(NPAS2,color='red')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.782',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='NPAS2',alpha=1,size=6)
roc.NPAS2
NPAS2 <- roc(roc_dat$Group,roc_dat$NPAS2,levels=c('Control','HF'))
ci.auc(NPAS2)
roc_dat$outome[roc_dat$Group=='Control'] <- 1
roc_dat$outome[roc_dat$Group=='HF'] <- 0
##95% CI:0.7293-0.8316 (DeLong)
roc.area(as.numeric(as.vector(roc_dat$outome)),NPAS2$predictor)
#  9.130502e-18
## *8.PER1
PER1 <- rocit(score = roc_dat$PER1,class = roc_dat$Group,method = 'binormal')
roc.PER1 <- ggplot()+geom_roc(PER1,color='red')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.643',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='PER1',alpha=1,size=6)
roc.PER1
PER1 <- roc(roc_dat$Group,roc_dat$PER1,levels=c('Control','HF'))
ci.auc(PER1)
roc_dat$outome[roc_dat$Group=='Control'] <- 1
roc_dat$outome[roc_dat$Group=='HF'] <- 0
##95% CI:0.5723-0.6943
roc.area(as.numeric(as.vector(roc_dat$outome)),PER1$predictor)
#  2.641935e-05
# 10-基于验证集的ROC曲线验证 --------------------------------------------------------
# 10-1 数据下载和预处理 -------------------------------------------------------------
## 获取方式:本地下载、上传,文件名:'GSE42955_series_matrix.csv'
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
exp <- read.csv("/data/nas1/liujy/07-BJTC-195/GSE5406_series_matrix.csv")
# 01-样本信息 -----------------------------------------------------------------
library(GEOquery)
eSet <- getGEO('GSE5406', 
               destdir = '.', 
               getGPL = F)
pheno <- pData(eSet[[1]]) ## 获得临床信息,共195个样本
pheno <- pheno[,c(2,10)]
pheno$characteristics_ch1 <- ifelse(pheno$characteristics_ch1 == "normally functioning myocardium from unused donor heart","Control","HF")

# 03-芯片注释信息 ---------------------------------------------------------------
## 芯片注释文件的获取
# library(devtools)
# install_github("jmzeng1314/idmap1",force = T)
# install_github("jmzeng1314/idmap2")
# install_github("jmzeng1314/idmap3")
library(idmap1)
library(idmap2) 
library(idmap3)
## 下载芯片注释信息
ids <- getIDs("GPL96")
colnames(ids)[1] <- "ID_REF"
ids <- ids[,-3]

# 04-训练集 ------------------------------------------------------------------
validationSet <- merge(ids,exp,by = "ID_REF")
validationSet <- validationSet[,-1]
# 19812:211

## 删除Symbol为空的行
validationSet <- validationSet[validationSet$symbol!='', ]
## 查看是否存在表达量都是0的行
View(validationSet[which(rowSums(validationSet[,2:211])==0),])
## 查看symbol列是否存在重复
duplic_list <- validationSet$symbol
aa <- duplic_list[duplicated(validationSet$symbol)]
dupli_expr <- validationSet[which(validationSet$symbol %in% aa),]
## 去重、取均值
validationSet <- aggregate(validationSet[,2:211],
                           by = list(symbol = validationSet$symbol),
                           mean) 
## 18837:314
rownames(validationSet) <- validationSet[,1]
validationSet <- validationSet[,-1]


# 10-5 ROC曲线 --------------------------------------------------------------
# diagnosisGenes <- c("ARNTL","NPAS2","CRY2","BHLHE41","EPHX2")
diagnosisGenes <- c("NPAS2","PER3","CRY2","BHLHE41","ARNTL")
roc_dat <- validationSet[diagnosisGenes,]
roc_dat <- as.data.frame(t(roc_dat))
roc_dat$Sample <- rownames(roc_dat)
conNum <- pheno$geo_accession[which(pheno$characteristics_ch1 == "Control")]
roc_dat$Group <- ifelse(roc_dat$Sample %in% conNum,'Control','HF')
roc_dat <- roc_dat[,c(6:7,1:5)]
rownames(roc_dat) <- NULL
roc_dat$Group <- as.factor(roc_dat$Group)

# 绘制ROC曲线
library(plotROC)
library(geomROC)
library(ROCit)
## 1.ARNTL
ARNTL <- rocit(score = roc_dat$ARNTL,class = roc_dat$Group,method = 'binormal')
roc.ARNTL <- ggplot()+geom_roc(ARNTL,color='red')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.764',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='ARNTL',alpha=1,size=6)
roc.ARNTL
ARNTL <- roc(roc_dat$Group,roc_dat$ARNTL,levels=c('Control','HF'))
ci.auc(ARNTL)
roc_dat$outome[roc_dat$Group=='Control'] <- 1
roc_dat$outome[roc_dat$Group=='HF'] <- 0
##95% CI:0.6646-0.887 (DeLong)
roc.area(as.numeric(as.vector(roc_dat$outome)),ARNTL$predictor)
#   0.0001251673
## 2.PER3
PER3 <- rocit(score = roc_dat$PER3,class = roc_dat$Group,method = 'binormal')
roc.PER3 <- ggplot()+geom_roc(PER3,color='red')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.740',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='PER3',alpha=1,size=6)
roc.PER3
PER3 <- roc(roc_dat$Group,roc_dat$PER3,levels=c('Control','HF'))
ci.auc(PER3)
roc_dat$outome[roc_dat$Group=='Control'] <- 0
roc_dat$outome[roc_dat$Group=='HF'] <- 1
##95% CI:0.6202-0.8976 (DeLong)
roc.area(as.numeric(as.vector(roc_dat$outome)),PER3$predictor)
#  0.0002938796
## 3.CRY2
CRY2 <- rocit(score = roc_dat$CRY2,class = roc_dat$Group,method = 'binormal')
roc.CRY2 <- ggplot()+geom_roc(CRY2,color='red')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.850',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='CRY2',alpha=1,size=6)
roc.CRY2
CRY2 <- roc(roc_dat$Group,roc_dat$CRY2,levels=c('Control','HF'))
ci.auc(CRY2)
roc_dat$outome[roc_dat$Group=='Control'] <- 0
roc_dat$outome[roc_dat$Group=='HF'] <- 1
##95% CI:0.7517-0.959 (DeLong)
roc.area(as.numeric(as.vector(roc_dat$outome)),CRY2$predictor)
#1.183886e-06

## 4.BHLHE41
BHLHE41 <- rocit(score = roc_dat$BHLHE41,class = roc_dat$Group,method = 'binormal')
roc.BHLHE41 <- ggplot()+geom_roc(BHLHE41,color='#840000')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.546',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='BHLHE41',alpha=1,size=6)
roc.BHLHE41
BHLHE41 <- roc(roc_dat$Group,roc_dat$BHLHE41,levels=c('Control','HF'))
ci.auc(BHLHE41)
roc_dat$outome[roc_dat$Group=='Control'] <- 0
roc_dat$outome[roc_dat$Group=='HF'] <- 1
##95% CI:0.1223-0.9944 (DeLong)
roc.area(as.numeric(as.vector(roc_dat$outome)),BHLHE41$predictor)
# 0.002195358
## 5.NPAS2
NPAS2 <- rocit(score = roc_dat$NPAS2,class = roc_dat$Group,method = 'binormal')
roc.NPAS2 <- ggplot()+geom_roc(NPAS2,color='red')+
  roc_diagonal()+
  roc_theme()+
  annotate('text',x=0.75,y=0.25,label='AUC = 0.842',alpha=1,size=6)+
  annotate('text',x=0.15,y=1,label='NPAS2',alpha=1,size=6)
roc.NPAS2
NPAS2 <- roc(roc_dat$Group,roc_dat$NPAS2,levels=c('Control','HF'))
ci.auc(NPAS2)
roc_dat$outome[roc_dat$Group=='Control'] <- 1
roc_dat$outome[roc_dat$Group=='HF'] <- 0
##95% CI:0.7723-0.97 (DeLong)
roc.area(as.numeric(as.vector(roc_dat$outome)),NPAS2$predictor)
# 4.1341e-07