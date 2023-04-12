rm(list = ls())
# 04 风险模型构建---------
## 04-1 单因素cox回归----
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./07_univariate_cox")){
  dir.create("./07_univariate_cox")
}
setwd("./07_univariate_cox")
library(lance)
library(readxl)
library(readr)
library(tidyverse)
survival<-read_tsv(file = '/data/nas1/luchunlin/TCGA_survival/TCGA-COAD.survival.tsv')
survival<-survival[,-3]
dat.tcga<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
group<-read.delim2("../01_CIBERSORT/group.xls")
tumor.sample<-group$sample[which(group$group=='Tumor')]
##448
#dat.tcga<-log2(dat.tcga+1)
DEIRG<-read.delim2('../06_DEIRG/DEIRG.xls')
survival<-survival[survival$sample%in%colnames(dat.tcga),]
write.table(survival,file = 'survival.xls',sep = '\t',quote = F,row.names = F)
train_dat<-dat.tcga[,tumor.sample]
train_dat<-t(train_dat)
train_dat<-log2(train_dat+0.1)
train_dat <- t(scale(t(train_dat)))
train_dat<-as.data.frame(train_dat)
train_dat<-train_dat[,colnames(train_dat)%in%DEIRG$symbol]
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
# 7
write.table(res_results_0.05,
            file = "univariate_cox_result_0.05.xls",
            quote = F,
            sep = '\t',
            row.names = T)
dim(res_results_0.05)
## 7 2
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

pdf(file = "01.univariate_cox_forest.pdf", height = 5, width = 9, onefile = F)
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
           xticks = c(0,1,2,3,4,5), #横坐标刻度
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
png(filename = "01.univariate_cox_forest.png", height = 350, width = 600)
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
           xticks = c(0,1,2,3,4,5), #横坐标刻度
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

## 04-2 Lasso回归-----
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./08_Lasso")){
  dir.create("./08_Lasso")
}
setwd("./08_Lasso")
library(glmnet)
x_all <- subset(train_data, select = -c(OS, OS.time))
x_all <- x_all[,rownames(res_results_0.05)]
y_all <- subset(train_data, select = c(OS, OS.time))

# 拟合模型
fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS), 
              family = "cox") 
#dev.new()
png(filename = "01.lasso_model.png", height = 400, width = 500)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
pdf(file = "01.lasso_model.pdf", height = 5)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
# 交叉验证拟合模型
set.seed(8)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),nfold=10,
                  family = "cox") 

png(filename = "02.lasso_verify.png", height = 400, width = 500)
plot(cvfit, las =1)
dev.off()
pdf(file = "02.lasso_verify.pdf", height = 5)
plot(cvfit, las =1)
dev.off()
# 提取指定lambda时特征的系数
coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
#coef.min
#coef.min<-as.data.frame(as.matrix(coef.min))[lasso_geneids$lasso_geneids,]
#rownames(coef.min)<-lasso_geneids$lasso_geneids
#class(coef.min)
cvfit$lambda.min
# [1] 0.00835155
# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)
coef.min
# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids
## 6
## [1] "ZG16"     "MPC1"     "RBM47"    "SMOX"     "CPM"      "DNASE1L3"
write(lasso_geneids, "lasso_genes.csv")
# write.csv(x_all,file = "Lasso_x.csv",quote = F)
# write.csv(y_all,file = "Lasso_y.csv",quote = F)
# 05 风险模型的构建与验证------
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./09_risk")){
  dir.create("./09_risk")
}
setwd("./09_risk")
riskScore=predict(cvfit,newx = as.matrix(x_all),s=cvfit$lambda.min)
coef.min
riskScore<-as.numeric(riskScore)
class(riskScore)
coxGene=lasso_geneids
outCol=c("OS","OS.time",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),0,1))
risk <- as.data.frame(c(cbind(id=rownames(cbind(train_data[,outCol],
                                                riskScore,
                                                risk)),
                              cbind(train_data[,outCol],
                                    riskScore,
                                    risk))))
table(risk$risk)
## 0   1 
## 224 224 
write.table(risk,file = 'risk.xls',quote = F,row.names = F,sep = '\t')
library(ggplot2)
library(ggthemes)
library(Ipaper)
median(riskScore)
risk_dis <- ggplot(risk, aes(x=reorder(id, riskScore), 
                             y=riskScore, 
                             color = factor(risk, 
                                            levels = c(0, 1), 
                                            labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,600,700,800,900,1000,1100)],
                   labels = c(1,100,200,300,400,500,600,700,800,900,1000,1100),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(risk[which(risk$risk==1),]) + 0.5,
             lty = 2) +
  geom_hline(yintercept = median(riskScore),
             lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Train Risk Score Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
risk_dis
ggsave(filename = '01.Train_Risk_Score_Distribution.pdf',risk_dis,w=7,h=5)
ggsave(filename = '01.Train_Risk_Score_Distribution.png',risk_dis,w=7,h=5)
surv_stat <- ggplot(risk, aes(x=reorder(id, riskScore),
                              y=OS.time/365,
                              color = factor(OS,
                                             levels = c(0,1),
                                             labels = c("Alive", "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF","#A73030FF")) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,600,700,800,900,1000,1100)],
                   labels = c(1,100,200,300,400,500,600,700,800,900,1000,1100),
                   expand = c(0.02,0)) +
  ylim(c(0,30))+
  geom_vline(xintercept = nrow(risk[which(risk$risk==1),]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Survival time(Years)",
       title = "Train survival state distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))

surv_stat
ggsave(filename = '02.Train_survival_state_Distribution.pdf',surv_stat,w=7,h=5)
ggsave(filename = '02.Train_survival_state_Distribution.png',surv_stat,w=7,h=5)
## 05-2 KM曲线和ROC曲线------
kmfit<-survfit(Surv(OS.time, OS) ~ risk, data =  risk)
train_survival_median <- ggsurvplot(kmfit,
                                    pval = TRUE, 
                                    
                                    conf.int = F,
                                    legend.labs=c("High risk","Low risk" ),
                                    legend.title="Risk score",
                                    title="Train KM",
                                    font.main = c(15,"bold"),
                                    risk.table = TRUE, 
                                    risk.table.col = "strata", 
                                    linetype = "strata", 
                                    surv.median.line = "hv", 
                                    ggtheme = theme_bw(), 
                                    palette = c("#A73030FF", "#0073C2FF"))
train_survival_median

riskscore <- function(survival_cancer_df,
                      candidate_genes_for_cox,
                      cox_report){
  library("dplyr")
  risk_score_table <- survival_cancer_df[, candidate_genes_for_cox]
  for (each_sig_gene in colnames(risk_score_table)){
    risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*
      (summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table <- cbind(risk_score_table,
                            "total_risk_score"=exp(rowSums(risk_score_table))) %>%
    cbind(survival_cancer_df[,c("OS.time", "OS")])
  risk_score_table <- risk_score_table[,c("OS.time",
                                          "OS",
                                          candidate_genes_for_cox,
                                          "total_risk_score")]
}

multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore,
                           predict.time=single_time,
                           method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
                           risk_score_table = risk)
table(for_multi_ROC$AUC)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)
# 画ROC曲线 
library(scales)
library(geomROC)
library(plotROC)

auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#4682B4", "#FF4040", "#20B2AA")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "Train ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
ggsave('04.Train_ROC.png', ROC,width = 5, height = 4)
ggsave('04.Train_ROC.pdf', ROC,width = 5, height = 4)

### 热图--------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
lasso.dat<-dat.tcga[lasso_geneids,]
group_rt<-data.frame(sample=risk$id,group=risk$riskScore)
group_rt$group<-ifelse(group_rt$group>median(group_rt$group),'High risk','Low risk')
group_rt<-group_rt[order(group_rt$group),]
rt<-lasso.dat[,group_rt$sample]
group_rt<-group_rt$group%>%as.data.frame()
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
x<-log2(rt+1)
#x<-log2(heat+1)
#x<-t(scale(t(heat)))
ann_colors<-list(
  group = c('Low risk'="lightblue",'High risk'="darkorange"))
#annotation_raw<-data.frame(row.names = rownames(heat),
#                           Change=factor(rep(c('Down','Up'),c(10,10))))
library(pheatmap)
pdf('05.heatmap.pdf',w=8,h=6)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         #       annotation_row = annotation_raw,
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 12,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         annotation_names_row = F)
dev.off()
png('05.heatmap.png',w=700,h=600)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         #       annotation_row = annotation_raw,
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 12,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         annotation_names_row = F)
dev.off()

# 外部验证----------
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./15_validation")){
  dir.create("./15_validation")
}
setwd("./15_validation")

library(GEOquery)
library(Biobase)
gset<-getGEO("GSE41258",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL96",destdir = '.')
library(tidyverse)
library(lance)
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat.va<-expr
dat.va$ID<-rownames(dat.va)
dat.va$ID<-as.character(dat.va$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat.va<-dat.va %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
#dat.va <- read_xlsx('GSE138866_130_Omental_Mets_MedGenome_normalized.xlsx')
#dat.va <- dat.va[!duplicated(dat.va$Gene),]%>%column_to_rownames(var = 'Gene')
pd<-pData(a)
table(pd$`tissue:ch1`)
pd <- subset(pd,!`tissue:ch1`=='Normal Colon'|`tissue:ch1`=="Normal Liver"|`tissue:ch1`=='Normal Lung')
#pd <- subset(pd,`tissue:ch1`=='colon tumor')
write.table(dat.va,file = 'dat.va.xls',sep = '\t',quote = F,row.names = T)
survival.va<-data.frame(sample=pd$geo_accession,OS=pd$`fup status (ned,  no evidence of disease, awd, alive with disease, aun, alive unknown, dod, dead of disease, doc, dead of other cause, dun, dead, cause unknown):ch1`,OS.time=pd$`fup interval:ch1`)
#survival.va<-data.frame(sample=pd$geo_accession,OS=pd$`os censor:ch1`,OS.time=pd$`survival (mo):ch1`)
survival.va <- subset(survival.va,OS=='AWD'|OS=='AUN'|OS=='DOC'|OS=='DOD'|OS=='DUN')
#survival.va <- subset(survival.va,OS=='AWD'|OS=='AUN'|OS=='DOD')
table(survival.va$OS)
#survival.va$OS[survival.va$OS=='N/A'] <- NA
#survival.va <- na.omit(survival.va)
#survival.va<-survival.va[which(!is.na(survival.va$OS)),]
survival.va$OS<-ifelse(survival.va$OS=='AUN',0,
                       ifelse(survival.va$OS=='AWD',0,1))
survival.va$OS.time<-as.numeric(survival.va$OS.time)*30
write.table(survival.va,file = 'survival.va.xls',sep = '\t',row.names = F,quote = F)
survival.va$OS<-as.numeric(survival.va$OS)
##GALNTL6 "GALNT15"=====GALNTL17 GALNTL2"
test_data<-t(dat.va[,survival.va$sample])%>%as.data.frame()
#test_data<-log2(test_data+0.001)
test_data<-scale(test_data)%>%as.data.frame()
test_data$sample<-rownames(test_data)
test_data<-merge(survival.va,test_data,by='sample')
test_data<-column_to_rownames(test_data,var='sample')

##多因素-----
# riskScore_out = predict(cox_more_2,type="lp",newdata=test_data)
# risk_out=as.vector(ifelse(riskScore_out>median(riskScore_out),0,1))
# risk_out <- as.data.frame(c(cbind(id=rownames(cbind(test_data[,outCol],
#                                                     riskScore_out,
#                                                     risk_out)),
#                                   cbind(test_data[,outCol],
#                                         riskScore_out,
#                                         risk_out))))
# 
# library(ggplot2)
# library(ggthemes)
# median(riskScore_out)
# table(risk_out$risk_out)

## LASSO------
test_data2<-test_data[,lasso_geneids]
risk_out<-data.frame(test_data2)
risk_out$risk_out<-NA
risk_out$riskScore_out<-NA
cnt<-1

coef.min<-coef.min[lasso_geneids,]

while (cnt < length(rownames(test_data2))+1) {
  risk_out$riskScore_out[cnt]<-sum(coef.min*test_data2[cnt,])
  cnt = cnt + 1
}

dim(test_data2)
cnt<-1
while (cnt < length(rownames(test_data2))+1) {
  risk_out$risk_out[cnt]=as.vector(ifelse(risk_out$riskScore_out[cnt]>median(risk_out$riskScore_out),0,1))
  cnt = cnt + 1
}
riskScore_out<-as.numeric(risk_out$riskScore_out)

risk_out=as.vector(ifelse(riskScore_out>median(riskScore_out),0,1))
risk_out <- as.data.frame(c(cbind(id=rownames(cbind(test_data[,outCol],
                                                    riskScore_out,
                                                    risk_out)),
                                  cbind(test_data[,outCol],
                                        riskScore_out,
                                        risk_out))))
write.table(risk_out,file = 'risk_out.xls',quote = F,row.names = F,sep = '\t')
library(ggplot2)
library(ggthemes)
median(riskScore_out)
table(risk_out$risk_out)
risk_dis_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out), 
                                     y=riskScore_out, 
                                     color = factor(risk_out, 
                                                    levels = c(0, 1), 
                                                    labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,50,100,150,200,250,300,350,400)],
                   labels = c(1,50,100,150,200,250,300,350,400),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
             lty = 2) +
  geom_hline(yintercept = median(riskScore_out),
             lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Validation Risk Score Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
risk_dis_out
ggsave(filename = '01.Validation_Risk_Score_Distribution.pdf',risk_dis_out,w=7,h=5)
ggsave(filename = '01.Validation_Risk_Score_Distribution.png',risk_dis_out,w=7,h=5)
surv_stat_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out),
                                      y=OS.time/365,
                                      color = factor(OS,
                                                     levels = c(0,1),
                                                     labels = c("Alive", "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF", "#A73030FF")) +
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,50,100,150,200,250,300,350,400)],
                   labels = c(1,50,100,150,200,250,300,350,400),
                   expand = c(0.02,0)) +
  ylim(x=c(0,20)) +
  geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Survival time (years)",
       title = "Validation Survival State Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
surv_stat_out
ggsave(filename = '02.Validation_survival_state_Distribution.pdf',surv_stat_out,w=7,h=5)
ggsave(filename = '02.Validation_survival_state_Distribution.png',surv_stat_out,w=7,h=5)
library(survival)
library(survminer)
kmfit_out <- survfit(Surv(OS.time, OS) ~ risk_out, data =  risk_out)
verify_survival_median <- ggsurvplot(kmfit_out,
                                     pval = TRUE, 
                                     conf.int = F,
                                     legend.labs=c("High risk","Low risk" ),
                                     legend.title="Risk score",
                                     title="Validition KM",
                                     font.main = c(15,"bold"),
                                     risk.table = TRUE, 
                                     risk.table.col = "strata", 
                                     linetype = "strata", 
                                     surv.median.line = "hv", 
                                     ggtheme = theme_bw(), 
                                     palette = c("#A73030FF", "#0073C2FF"))

verify_survival_median


library(survivalROC)
multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore,
                           predict.time=single_time,
                           method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
                           risk_score_table = risk_out)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)
#devtools::install_github('yikeshu0611/geomROC')
library(scales)
library(geomROC)
library(plotROC)
library(ggthemes)

auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)


ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365","1095", "1825"),
                     labels = c("1 years", "3 years","5 years"),
                     values = c("#4682B4", "#FF4040", "#20B2AA")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "Test ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
ggsave('04.Test_ROC.png', ROC,width = 5, height = 4)
ggsave('04.Test_ROC.pdf', ROC,width = 5, height = 4)

### 热图--------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
lasso.va<-dat.va[lasso_geneids,]
group_rt<-data.frame(sample=risk_out$id,group=risk_out$riskScore_out)
group_rt$group<-ifelse(group_rt$group>median(group_rt$group),'High risk','Low risk')
group_rt<-group_rt[order(group_rt$group),]
rt<-lasso.va[,group_rt$sample]
group_rt<-group_rt$group%>%as.data.frame()
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
x<-scale(rt)
#x<-log2(heat+1)
#x<-t(scale(t(heat)))
ann_colors<-list(
  group = c('Low risk'="lightblue",'High risk'="darkorange"))
#annotation_raw<-data.frame(row.names = rownames(heat),
#                           Change=factor(rep(c('Down','Up'),c(10,10))))
library(pheatmap)
pdf('05.heatmap.pdf',w=6,h=4)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         #       annotation_row = annotation_raw,
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 12,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         annotation_names_row = F)
dev.off()
png('05.heatmap.png',w=500,h=350)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         #       annotation_row = annotation_raw,
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 12,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         annotation_names_row = F)
dev.off()
