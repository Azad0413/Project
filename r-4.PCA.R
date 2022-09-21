rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-258/")
if (! dir.exists("./04_PCA")){
  dir.create("./04_PCA")
}
setwd("./04_PCA")
library(tidyverse)
library(lance)
dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-258/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colname<-data.frame(sample=colnames(dat))
colname$sample<-gsub('.','-',colname$sample,fixed = T)
colnames(dat)<-colname$sample
DEARG<-read.delim2('/data/nas1/luchunlin/project/BJTC-258/02_DEGs/sig.all.xls')                                              
group<-read.delim2('/data/nas1/luchunlin/project/BJTC-258/01_consensus/cluster.xls')
table(group$cluster)
## 加一个单因素cox-------
library(readxl)                                               
library(readr)
survival<-read_tsv(file = 'TCGA-OV.survival.tsv')
survival<-survival[,-3]
#dat.tcga<-log2(dat.tcga+1)
survival<-survival[survival$sample%in%colnames(dat),]
write.table(survival,file = 'survival.xls',sep = '\t',quote = F,row.names = F)
train_dat<-t(dat)
train_dat<-log2(train_dat+1)
train_dat<-as.data.frame(train_dat)
train_dat<-train_dat[,colnames(train_dat)%in%DEARG$symbol]
train_dat$sample<-rownames(train_dat)
train_dat<-merge(survival,train_dat,by='sample')
train_dat<-column_to_rownames(train_dat,var = 'sample')
train_data<-train_dat
### 单因素cox
library(survival)
library(survminer)
colnames_sum <- colnames(train_data)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(train_data) <- colnames_sum

covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]
#Surv()函数产生一个生存对象  生存时间对生存的影响 对每一个变量构建生存分析公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))  #as.formula(). 将字符串转换成公式。构建formula对象
# coxph函数用于计算cox模型 循环对每一个特征做cox回归分析
univ_models <- lapply(univ_formulas,
                      function(x) {coxph(x, data = train_data)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         #获取HR
                         HR <-signif(x$coef[2], digits=4);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],4)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],4)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", 
                                      HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

## coef是公式中的回归系数b（有时候也叫beta值）。exp(coef)是cox模型中的风险比（HR）
## z代表wald统计量，是coef除以其标准误se(coef)。ower .95 upper .95则是exp(coef)的95%置信区间，可信区间越窄，可信度越高，你的实验越精确，越是真理。
res_mod <- t(as.data.frame(univ_results, check.names = FALSE))
res_mod <- as.data.frame(res_mod)
res_results_0.05 <- res_mod[which(as.numeric(res_mod$p.value) < 0.05),]
res_results_0.05 <- na.omit(res_results_0.05)
# 70
write.table(res_results_0.05,
            file = "univariate_cox_result_0.05.xls",
            quote = F,
            sep = '\t',
            row.names = T)
dim(res_results_0.05)
## 70 2
library(tidyr)
res_results_0.05_2 <- separate(res_results_0.05, "HR (95% CI for HR)",
                               into = c("HR", "HR.95L", "HR.95H"),
                               sep = " ")
res_results_0.05_2 <- separate(res_results_0.05_2, "HR.95L",
                               into = c("HR.95L", "HR.95H"),
                               sep = "\\-")
res_results_0.05_2$HR.95L <- gsub("\\(", "", res_results_0.05_2$HR.95L)
res_results_0.05_2$HR.95H <- gsub("\\)", "", res_results_0.05_2$HR.95H)

res_results_0.05_2[,1:ncol(res_results_0.05_2)] <- as.numeric(unlist(res_results_0.05_2[,1:ncol(res_results_0.05_2)]))
res_results_0.05_2 <- res_results_0.05_2[order(res_results_0.05_2$p.value),]
res_results_0.05_2<-res_results_0.05_2[c(1:10),]
res_results_0.05_2 <- res_results_0.05_2[order(res_results_0.05_2$HR),]
hz <- paste(round(res_results_0.05_2$HR,4),
            "(",round(res_results_0.05_2$HR.95L,4),
            "-",round(res_results_0.05_2$HR.95H,4),")",sep = "")


tabletext <- cbind(c(NA,"Gene",rownames(res_results_0.05_2)),
                   c(NA,"P value",ifelse(res_results_0.05_2$p.value<0.001,
                                         "< 0.001",
                                         round(res_results_0.05_2$p.value,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)

pdf(file = "01.univariate_cox_forest.pdf", height = 7, width = 9, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE,rep(FALSE, 70)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #95%置信区间下限
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # 垂直于x轴的网格线，对应每个刻度
dev.off()
png(filename = "01.univariate_cox_forest.png", height = 400, width = 600)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE,rep(FALSE, 70)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #95%置信区间下限
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # 垂直于x轴的网格线，对应每个刻度

dev.off()
### 70
# PCA------
uncoxgene<-data.frame(symbol=rownames(res_results_0.05))
dat.pca<-log2(dat+1)
#dat.pca<-dat
dat.pca <- dat.pca[uncoxgene$symbol,]
dat.pca <- na.omit(dat.pca)
#mads <- apply(dat.pca,1,mad)
#dat.pca <- dat.pca[rev(order(mads))[1:600],]
dat.pca <- scale(dat.pca,center = T,scale = T)
dat.pca <- as.data.frame(t(dat.pca))
pca<-prcomp(dat.pca
            #,center = T,scale. = T
            )
screeplot(pca,type="lines")
df<-pca$x  ## 提取PC score
df <- as.data.frame(df)
# 提取主成分的方差贡献率,生成坐标轴标题
summ <- summary(pca)
summ
xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

# 绘制PCA得分图
table(group$cluster)


library(ggplot2)
rownames(group)<-group$sample
group<-group[rownames(dat.pca),]
p.pca <- ggplot(data = df,aes(x = PC1,y = PC2,color = group$cluster))+
  stat_ellipse(aes(fill = group$cluster),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 3.5)+
  labs(x = xlab,y = ylab,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("purple","orange","lightblue"))+
  scale_colour_manual(values = c("purple","orange","lightblue"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.pca
ggsave('02.pca.pdf',p.pca,w=6,h=5)
ggsave('02.pca.png',p.pca,w=6,h=5)

y <- eigen(cor(dat.pca))
# 根据累计贡献率大于60%确定主成分
sum(y$values[1:10])/sum(y$values)
s <- df[,1:10]
scores = 0.0
for (i in 1:10)
  scores=(y$values[i]*s[,i])/(sum(y$values[1:10]))+scores
pcaScore <- cbind(s,scores) # 输出综合得分信息
pcaScore <- data.frame(sample=rownames(pcaScore),scores=pcaScore$scores)
write.table(pcaScore,file = 'pcaScore.xls',sep = '\t',row.names = F,quote = F)

### 预后相关基因在不同衰老簇中的表达情况------
uncoxgene$symbol<-gsub("_","-",uncoxgene$symbol)
heat.dat<-dat[uncoxgene$symbol,]
heat.dat<-log2(heat.dat+0.001)
rt.group<-group
rt.group<-rt.group[order(rt.group$cluster),]
heat.dat<-heat.dat[,rt.group$sample]
rt.group<-data.frame(row.names = rownames(rt.group),group=rt.group$cluster)
ann_colors<-list(
  group=c('cluster1'='#FFB7DD','cluster2'='#77DDFF','cluster3'='orange'))
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
pdf(file = "03.heatmap.pdf", height = 10, width = 7)
pheatmap(heat.dat,
         color = bluered(100),
         border_color = NA,
         annotation_col = rt.group,
         annotation_colors = ann_colors,
         labels_row = NULL,
         clustering_method = 'ward.D2',
         show_rownames = T,
         show_colnames = F,
         fontsize_col = 5,
         cluster_cols = F,
         cluster_rows = T,
         fontsize_row = 8)
dev.off()
png(file = "03.heatmap.png", height = 1000, width = 700)
pheatmap(heat.dat,
         color = bluered(100),
         border_color = NA,
         annotation_col = rt.group,
         annotation_colors = ann_colors,
         labels_row = NULL,
         clustering_method = 'ward.D2',
         show_rownames = T,
         show_colnames = F,
         fontsize_col = 5,
         cluster_cols = F,
         cluster_rows = T,
         fontsize_row = 8)
dev.off()

# ## 自写代码-----
# library(ggplot2) 
# library(ggrepel)
# ## 标准化原始数据
# dat.sacle<-scale(dat.pca)
# ## 计算相关系数矩阵
# cor.mat<-cor(dat.sacle)
# ## 特征分解
# rs.mat<-eigen(cor.mat)
# ## 提取特征值，即个主成分的方差
# val<-rs.mat$values
# ## 换算成标准差
# standard.deviation<-sqrt(val)
# ## 计算方差贡献率
# proportion.of.variance<-val/sum(val)
# ## 计算累积贡献率
# cumulative.proportion<-cumsum(proportion.of.variance)
# ## 提取特征向量，即载荷矩阵（loadings）
# load.mat<-as.matrix(rs.mat$vectors)
# ##计算主成分得分
# PC<-dat.sacle%*%load.mat
# colnames(PC)<-paste('PC',1:ncol(PC),sep = '')
# ## 转换成数据框用于绘图
# df2<-as.data.frame(PC)
# ##提取主成分的方差贡献率,生成坐标轴标题
# xlab2 <- paste0("PC1(",round(proportion.of.variance[1]*100,2),"%)")
# ylab2 <- paste0("PC2(",round(proportion.of.variance[2]*100,2),"%)")
# p.pca2 <- ggplot(data = df2,aes(x = PC1,y = PC2,color = group$cluster))+
#   geom_point(size = 3)+
#   theme_bw()+
#   labs(x = xlab2,y = ylab2,color = "Group",title = "Plot of PCA score")+
#   stat_ellipse(aes(fill = group$cluster),
#                type = "norm",geom = "polygon",alpha = 0.2,color = NA)+
#   guides(fill = "none")+
#   theme(plot.title = element_text(hjust = 0.5,size = 15),
#         axis.text = element_text(size = 11),axis.title = element_text(size = 13),
#         legend.text = element_text(size = 11),legend.title = element_text(size = 13),
#         plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
# p.pca2
