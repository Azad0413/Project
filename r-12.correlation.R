rm(list = ls())
#12 correlation------------
setwd("/data/nas1/luchunlin/project/YCZK-127/")
if (! dir.exists("./12_cor")){
  dir.create("./12_cor")
}
setwd("./12_cor")
## 模型基因与差异免疫浸润细胞的相关性
## 基因表达数据（hub_exp） 和 免疫细胞浸润矩阵 (tiics_result)行名均为样本名称
## Spearman相关分析
tiics_result <- read.delim2('/data/nas1/luchunlin/project/YCZK-127/11_CIBERSORT/tiics_result.xls',row.names = 1)%>%lc.tableToNum()
diff.tiics<-read.delim2('/data/nas1/luchunlin/project/YCZK-127/11_CIBERSORT/diff_cibersort_Table.xls')
tiics_result<-tiics_result[rownames(diff.tiics),]
dat.tcga<-read.delim2("/data/nas1/luchunlin/project/YCZK-127/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
lasso.gene<-data.frame(symbol=c('KRT19','PSME2','HMGB3','MRPL13','SHCBP1'))
hub.exp<-dat.tcga[lasso.gene$symbol,]
hub.exp<-hub.exp[,colnames(tiics_result)]
colname<-data.frame(sample=colnames(hub.exp))
colname$sample<-gsub('.','-',colname$sample,fixed = T)
colnames(hub.exp)<-colname$sample
colnames(tiics_result)<-colname$sample
tiics_exp<-t(tiics_result)
hub.exp<-t(hub.exp)

### KRT19-----
KRT19<-as.numeric(hub.exp[,c('KRT19')])
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation_KRT19 <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),KRT19,method="spearman")
  ## 3.填充
  correlation_KRT19[i,1] = colnames(tiics_exp)[i]
  correlation_KRT19[i,2] = dd$estimate
  correlation_KRT19[i,3] = dd$p.value
}
### 修改名称
colnames(correlation_KRT19) <- c("cell","cor","p.value")
#correlation_KRT19<-correlation_KRT19[order(correlation_KRT19$cor,decreasing = T),]
correlation_KRT19$correlation<-cut(abs(correlation_KRT19$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
correlation_KRT19<-correlation_KRT19[order(correlation_KRT19$cor,decreasing = T),]
### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation_KRT19,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(5,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'KRT19')
## PSME2-----
PSME2<-as.numeric(hub.exp[,c('PSME2')])
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation_PSME2 <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),PSME2,method="spearman")
  ## 3.填充
  correlation_PSME2[i,1] = colnames(tiics_exp)[i]
  correlation_PSME2[i,2] = dd$estimate
  correlation_PSME2[i,3] = dd$p.value
}
### 修改名称
colnames(correlation_PSME2) <- c("cell","cor","p.value")
#correlation_PSME2<-correlation_PSME2[order(correlation_PSME2$cor,decreasing = T),]
correlation_PSME2$correlation<-cut(abs(correlation_PSME2$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
correlation_PSME2<-correlation_PSME2[order(correlation_PSME2$cor,decreasing = T),]
### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation_PSME2,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(5,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'PSME2')
## MRPL13-----
MRPL13<-as.numeric(hub.exp[,c('MRPL13')])
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation_MRPL13 <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),MRPL13,method="spearman")
  ## 3.填充
  correlation_MRPL13[i,1] = colnames(tiics_exp)[i]
  correlation_MRPL13[i,2] = dd$estimate
  correlation_MRPL13[i,3] = dd$p.value
}
### 修改名称
colnames(correlation_MRPL13) <- c("cell","cor","p.value")
#correlation_MRPL13<-correlation_MRPL13[order(correlation_MRPL13$cor,decreasing = T),]
correlation_MRPL13$correlation<-cut(abs(correlation_MRPL13$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
correlation_MRPL13<-correlation_MRPL13[order(correlation_MRPL13$cor,decreasing = T),]
### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation_MRPL13,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(5,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'MRPL13')
## HMGB3-----
HMGB3<-as.numeric(hub.exp[,c('HMGB3')])
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation_HMGB3 <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),HMGB3,method="spearman")
  ## 3.填充
  correlation_HMGB3[i,1] = colnames(tiics_exp)[i]
  correlation_HMGB3[i,2] = dd$estimate
  correlation_HMGB3[i,3] = dd$p.value
}
### 修改名称
colnames(correlation_HMGB3) <- c("cell","cor","p.value")
#correlation_HMGB3<-correlation_HMGB3[order(correlation_HMGB3$cor,decreasing = T),]
correlation_HMGB3$correlation<-cut(abs(correlation_HMGB3$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
correlation_HMGB3<-correlation_HMGB3[order(correlation_HMGB3$cor,decreasing = T),]
### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation_HMGB3,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(5,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'HMGB3')
## SHCBP1-------
SHCBP1<-as.numeric(hub.exp[,c('SHCBP1')])
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation_SHCBP1 <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),SHCBP1,method="spearman")
  ## 3.填充
  correlation_SHCBP1[i,1] = colnames(tiics_exp)[i]
  correlation_SHCBP1[i,2] = dd$estimate
  correlation_SHCBP1[i,3] = dd$p.value
}
### 修改名称
colnames(correlation_SHCBP1) <- c("cell","cor","p.value")
#correlation_SHCBP1<-correlation_SHCBP1[order(correlation_SHCBP1$cor,decreasing = T),]
correlation_SHCBP1$correlation<-cut(abs(correlation_SHCBP1$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))
correlation_SHCBP1<-correlation_SHCBP1[order(correlation_SHCBP1$cor,decreasing = T),]
### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation_SHCBP1,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(5,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'SHCBP1')

