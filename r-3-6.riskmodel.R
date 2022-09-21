rm(list = ls())
# 04 风险模型构建---------
## 04-1 单因素cox回归----

setwd("/data/nas1/luchunlin/project/BJTC-357/")
if (! dir.exists("./03_univariate_cox")){
  dir.create("./03_univariate_cox")
}
setwd("./03_univariate_cox")
library(lance)
library(readxl)
library(readr)
library(tidyverse)
survival<-read.delim2('/data/nas1/luchunlin/project/BJTC-357/00_rawdata/survival.xls')
survival$OS<-as.numeric(survival$OS)
survival$OS.time<-as.numeric(survival$OS.time)

dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-357/00_rawdata/dat.geo(lnc).xls", row.names = 1)%>% lc.tableToNum
dat<-dat[,colnames(dat)%in%survival$sample]
#dat<-t(scale(t(dat)))
#dat<-log2(dat+0.1)
##162
#sig.diff<-read.delim2('/data/nas1/luchunlin/project/BJTC-308/01_DEGs/DEG_sig.xls',row.names = 1)
DEORG<-read.delim2('/data/nas1/luchunlin/project/BJTC-357/02_DEORG/DEORG.xls')
train_dat<-t(dat)%>%as.data.frame()
train_dat<-train_dat[,colnames(train_dat)%in%DEORG$symbol]
train_dat$sample<-rownames(train_dat)
train_dat<-merge(survival,train_dat,by='sample')
train_dat<-column_to_rownames(train_dat,var = 'sample')
#179
### 按照生存状态（OS）进行三七分组。七分作为训练集，三分作为验证集。
library(caret)

set.seed(133)
## 133 
##118  验证集ROC不好
###25  39
all.expr<-train_dat
all.expr$sample<-rownames(all.expr)
sam<-createDataPartition(all.expr$OS,p= .7,list = F)
train_sample<-all.expr$sample[sam]
train_sample<-as.data.frame(train_sample)
test_sample<-all.expr$sample[-sam]
test_sample<-as.data.frame(test_sample)
test_data<-train_dat[rownames(train_dat)%in%test_sample$test_sample,]
train_data<-train_dat[rownames(train_dat)%in%train_sample$train_sample,]
dim(train_data)
##[1] 126 
write.table(train_data,file = 'train.data.xls',sep = '\t',row.names = T,quote = F)
dim(test_data)
##[1] 53
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
                         HR <-signif(x$coef[2], digits=3);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
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
hz <- paste(round(res_results_0.05_2$HR,3),
            "(",round(res_results_0.05_2$HR.95L,3),
            "-",round(res_results_0.05_2$HR.95H,3),")",sep = "")


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
png(filename = "01.univariate_cox_forest.png", height = 400, width = 800)
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
           xticks = c(0,01,2,3,4,5), #横坐标刻度
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

setwd("/data/nas1/luchunlin/project/BJTC-357")
if (! dir.exists("./04_Lasso")){
  dir.create("./04_Lasso")
}
setwd("./04_Lasso")
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
set.seed(3)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),nfold=3,
                  family = "cox") 

png(filename = "02.lasso_verify.png", height = 400, width = 500)
plot(cvfit, las =1)
dev.off()
pdf(file = "02.lasso_verify.pdf", height = 5)
plot(cvfit, las =1)
dev.off()
# 提取指定lambda时特征的系数
coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min
#[1] 0.04719635
# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)
coef.min
# LINC01672    0.850776742
# SNHG1        0.057262976
# KCNMA1_AS1   0.172214585
# CCR5AS      -0.049993331
# LINC01749   -0.000917922
# TMEM220_AS1 -0.068543123
# PCDH9_AS1   -0.023456133
# astive.coefficients <- as.numeric(coef.min)
# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids
## 7
## "LINC01672"   "SNHG1"       "KCNMA1_AS1"  "CCR5AS"      "LINC01749"   "TMEM220_AS1" "PCDH9_AS1"  

write(lasso_geneids, "lasso_genes.csv")
write.csv(x_all,file = "Lasso_x.csv",quote = F)
write.csv(y_all,file = "Lasso_y.csv",quote = F)

# 05 风险模型的构建与验证------
setwd("/data/nas1/luchunlin/project/BJTC-357/")
if (! dir.exists("./05_risk")){
  dir.create("./05_risk")
}
setwd("./05_risk")
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

write.table(risk,file = 'risk.xls',quote = F,row.names = F,sep = '\t')
#0  1 
#30 30  
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
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,120,40,60,80,100,120,140)],
                   labels = c(1,120,40,60,80,100,120,140),
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
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,120,40,60,80,100,120,140)],
                   labels = c(1,120,40,60,80,100,120,140),
                   expand = c(0.02,0)) +
  ylim(c(0,10))+
  geom_vline(xintercept = nrow(risk[which(risk$risk==1),]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Survival time (years)",
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
for_multi_ROC$Time <- factor(for_multi_ROC$Time)
table(for_multi_ROC$AUC)
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
                     values = c("#00468b", "#A73030FF", "#42b540")) +
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
ggsave('04.Train ROC.png', ROC,width = 5, height = 4)
ggsave('04.Train ROC.pdf', ROC,width = 5, height = 4)

## 05-3 外部数据库验证（KM，年ROC）---------
setwd("/data/nas1/luchunlin/project/BJTC-357/")
if (! dir.exists("./06_External_va")){
  dir.create("./06_External_va")
}
setwd("./06_External_va")

### 三七分--------
colnames_sum <- colnames(test_data)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(test_data) <- colnames_sum
test_data2<-test_data[,lasso_geneids]
risk_out<-data.frame(test_data2)
risk_out$risk_out<-NA
risk_out$riskScore_out<-NA
cnt<-1

coef.min<-coef.min[lasso_geneids,]

while (cnt < 54) {
  risk_out$riskScore_out[cnt]<-sum(coef.min*test_data2[cnt,])
  cnt = cnt + 1
}

dim(test_data2)
cnt<-1
while (cnt < 54) {
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

library(ggplot2)
library(ggthemes)
median(riskScore_out)
table(risk_out$risk_out)
## 26 27
risk_dis_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out), 
                                     y=riskScore_out, 
                                     color = factor(risk_out, 
                                                    levels = c(0, 1), 
                                                    labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,10,20,30,40,50,60)],
                   labels = c(1,10,20,30,40,50,60),
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
ggsave(filename = '01.Validation_Score_Distribution.png',risk_dis_out,w=7,h=5)
surv_stat_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out),
                                      y=OS.time/365,
                                      color = factor(OS,
                                                     levels = c(0,1),
                                                     labels = c("Alive", "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF", "#A73030FF")) +
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,10,20,30,40,50,60)],
                   labels = c(1,10,20,30,40,50,60),
                   expand = c(0.02,0)) +
  ylim(x=c(0,10)) +
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
table(for_multi_ROC$AUC)
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
                     values = c("#00468b", "#A73030FF", "#42b540")) +
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
ggsave('04.Test_ROC.png', ROC,width = 5, height = 4)
ggsave('04.Test_ROC.pdf', ROC,width = 5, height = 4)

# ## TCGA--------
# survival.va<-read_tsv(file = 'TCGA-ESCA.survival.tsv')
# survival.va<-survival.va[,-3]
# dat.escc<-read.delim2("/data/nas1/luchunlin/project/BJTC-357/00_rawdata/dat.all(fpkm).xls", row.names = 1)%>% lc.tableToNum
# dat.escc<-dat.escc[,-c(1:11)]
# #dat.tcga<-read.delim2("/data/nas1/luchunlin/project/BJTC-308/00_rawdata/dat.tcga.xls", row.names = 1)%>% lc.tableToNum
# #dat.tcga<-dat.tcga[,colnames(dat.escc)]
# dat.tcga<-log2(dat.escc+1)
# colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
# test_data2<-t(dat.tcga)%>%as.data.frame()
# test_data2<-test_data2[rownames(test_data2)%in%survival.va$sample,]
# test_data2$sample<-rownames(test_data2)
# test_data2<-merge(survival.va,test_data2,by='sample')%>%as.data.frame()
# test_data2<-column_to_rownames(test_data2,var = 'sample')
# colnames(test_data2)<-gsub('-','_',colnames(test_data2),fixed = T)
# colnames(test_data2)<-gsub(' ','_',colnames(test_data2),fixed = T)
# 
# test_data2<-test_data2[,colnames(test_data2)%in%colnames(train_data)]
# # 开始验证
# lasso_geneids<-as.data.frame(lasso_geneids)
# # # 开始验证
# ##手动算
# testdata3<-test_data2[,lasso_geneids$lasso_geneids]
# testdata3<-test_data2[,colnames(test_data2)%in%lasso_geneids$lasso_geneids]
# lasso_geneids
# #ENSG00000246089.5  ENSG00000267302.5
# # "MCPH1_AS1  RNFT1_DT"
# risk_out2<-data.frame(testdata3)
# risk_out2$risk_out<-NA
# risk_out2$riskScore_out<-NA
# cnt<-1
# coef.min<-coef.min[lasso_geneids$lasso_geneids,]
# 
# while (cnt < 81) {
#   risk_out2$riskScore_out[cnt]<-sum(coef.min*testdata3[cnt,])
#   cnt = cnt + 1
# }
# 
# dim(testdata3)
# cnt<-1
# while (cnt < 81) {
#   risk_out2$risk_out[cnt]=as.vector(ifelse(risk_out2$riskScore_out[cnt]>median(risk_out2$riskScore_out),0,1))
#   cnt = cnt + 1
# }
# riskScore_out2<-as.numeric(risk_out2$riskScore_out)
# risk_out2=as.vector(ifelse(riskScore_out2>median(riskScore_out2),0,1))
# risk_out2 <- as.data.frame(c(cbind(id=rownames(cbind(test_data2[,outCol],
#                                                      riskScore_out2,
#                                                      risk_out2)),
#                                    cbind(test_data2[,outCol],
#                                          riskScore_out2,
#                                          risk_out2))))
# table(risk_out2$risk_out2)
# library(ggplot2)
# library(ggthemes)
# median(riskScore_out2)
# 
# risk_dis_out <- ggplot(risk_out2, aes(x=reorder(id, riskScore_out2),
#                                       y=riskScore_out2,
#                                       color = factor(risk_out2,
#                                                      levels = c(0, 1),
#                                                      labels = c("High Risk", "Low Risk")))) +
#   geom_point() +
#   scale_color_manual(values = c("#A73030FF", "#0073C2FF")) +
#   scale_x_discrete(breaks = risk_out2[order(risk_out2$riskScore_out2),]$id[c(1,10,20,30,40,50,60,70,80,90)],
#                    labels = c(1,10,20,30,40,50,60,70,80,90),
#                    expand = c(0.02,0)) +
#   geom_vline(xintercept = nrow(risk_out2[which(risk_out2$risk_out2==1),]) + 0.5,
#              lty = 2) +
#   geom_hline(yintercept = median(riskScore_out2),
#              lty =2) +
#   labs(x = "Patients(increasing risk score)",
#        y = "Risk Score",
#        title = "Validation Risk Score Distribution") +
#   theme_base() +
#   theme(legend.title = element_blank(),
#         legend.key = element_blank(),
#         legend.justification = c(0,1),
#         legend.position = c(0,1),
#         #        legend.margin = margin(c(-5,4,4,3)),
#         legend.background = element_rect(color = "black", size = .3),
#         plot.title = element_text(size = 15, hjust = 0.5))
# risk_dis_out
# risk_out2$OS.time<-as.numeric(risk_out2$OS.time)
# surv_stat_out <- ggplot(risk_out2, aes(x=reorder(id, riskScore_out2),
#                                        y=OS.time/365,
#                                        color = factor(OS,
#                                                       levels = c(0,1),
#                                                       labels = c("Alive", "Dead")))) +
#   geom_point() +
#   scale_color_manual(values = c("#0073C2FF", "#A73030FF")) +
#   scale_x_discrete(breaks = risk_out2[order(risk_out2$riskScore_out2),]$id[c(1,10,20,30,40,50,60,70,80,90)],
#                    labels = c(1,10,20,30,40,50,60,70,80,90),
#                    expand = c(0.02,0)) +
#   ylim(x=c(0,10)) +
#   geom_vline(xintercept = nrow(risk_out2[which(risk_out2$risk_out2==1),]) + 0.5,
#              lty = 2) +
#   labs(x = "Patients(increasing risk score)",
#        y = "Survival time (years)",
#        title = "Validation Survival State Distribution") +
#   theme_base() +
#   theme(legend.title = element_blank(),
#         legend.key = element_blank(),
#         legend.justification = c(0,1),
#         legend.position = c(0,1),
#         #        legend.margin = margin(c(-5,4,4,3)),
#         legend.background = element_rect(color = "black", size = .3),
#         plot.title = element_text(size = 15, hjust = 0.5))
# surv_stat_out
# kmfit_out <- survfit(Surv(OS.time, OS) ~ risk_out2, data =  risk_out2)
# verify_survival_median <- ggsurvplot(kmfit_out,
#                                      pval = TRUE,
#                                      conf.int = F,
#                                      legend.labs=c("High risk","Low risk" ),
#                                      legend.title="Risk score",
#                                      title="Validition KM",
#                                      font.main = c(15,"bold"),
#                                      risk.table = TRUE,
#                                      risk.table.col = "strata",
#                                      linetype = "strata",
#                                      surv.median.line = "hv",
#                                      ggtheme = theme_bw(),
#                                      palette = c("#A73030FF", "#0073C2FF"))
# verify_survival_median
# 
# multi_ROC <- function(time_vector, risk_score_table){
#   library(survivalROC)
#   single_ROC <- function(single_time){
#     for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
#                            status=risk_score_table$OS,
#                            marker=risk_score_table$riskScore_out,
#                            predict.time=single_time,method = 'KM')
#     data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
#                'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
#                'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
#   }
#   multi_ROC_list <- lapply(time_vector, single_ROC)
#   do.call(rbind, multi_ROC_list)
# }
# 
# for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,3,1)),
#                            risk_score_table = risk_out2)
# for_multi_ROC$Time <- factor(for_multi_ROC$Time)
# 
# # 画ROC曲线
# auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
# auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
# auc_y2 <- round(for_multi_ROC[which(for_multi_ROC$Time==730),5][1],2)
# 
# 
# ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
#                                  y=True_positive,
#                                  label=Cut_values,
#                                  color=Time)) +
#   scale_color_manual(breaks = c("365", "730", "1095"),
#                      labels = c("1 years", "2 years", "3 years"),
#                      values = c("#00468b", "#A73030FF", "#42b540")) +
#   geom_roc(labels = F, stat = 'identity') +
#   style_roc() +
#   geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
#   theme_bw() +
#   labs(title = "Test ROC") +
#   theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
#         panel.grid = element_blank(),
#         plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
#   # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) +
#   # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) +
#   # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
#   annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
#            label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
#                      paste('AUC of 2 years =', format(auc_y2,nsmall=2)),
#                      paste('AUC of 3 years =', format(auc_y3,nsmall=2))))
# ROC
# ggsave('Test ROC.png', ROC,width = 5, height = 4)
# ggsave('Test ROC.pdf', ROC,width = 5, height = 4)
# # dat.va<-read.delim2("/data/nas1/luchunlin/project/BJTC-308/00_rawdata/dat1.xls", row.names = 1)%>% lc.tableToNum
# # dat.va<-dat.va[,colnames(dat.va)%in%survival$sample]
# # test_data<-t(dat.va)%>%as.data.frame()
# # colnames_sum <- colnames(test_data)
# # colnames_sum <- gsub("-","_",colnames_sum)
# # colnames_sum <- gsub(" ","_",colnames_sum)
# # colnames(test_data) <- colnames_sum
# # test_data$sample<-rownames(test_data)
# # test_data<-merge(survival,test_data,by='sample')
# # #161
# # test_data<-as.data.frame(test_data)
# # test_data<-column_to_rownames(test_data,var = 'sample')
# # test_data<-test_data[,colnames(test_data)%in%colnames(train_data)]
# # # 开始验证
# # lasso_geneids<-as.data.frame(lasso_geneids)
# # x_all_out <- subset(test_data, select = -c(OS, OS.time))
# # x_all_out <- x_all_out[,rownames(res_results_0.05)]
# # riskScore_out=predict(cvfit,newx = as.matrix(x_all_out),s=cvfit$lambda.min)
# # riskScore_out<-as.numeric(riskScore_out)
# # 
# # risk_out=as.vector(ifelse(riskScore_out>median(riskScore_out),0,1))
# # risk_out <- as.data.frame(c(cbind(id=rownames(cbind(test_data[,outCol],
# #                                                     riskScore_out,
# #                                                     risk_out)),
# #                                   cbind(test_data[,outCol],
# #                                         riskScore_out,
# #                                         risk_out))))
# # 
# # library(ggplot2)
# # library(ggthemes)
# # median(riskScore_out)
# # 
# # risk_dis_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out), 
# #                                      y=riskScore_out, 
# #                                      color = factor(risk_out, 
# #                                                     levels = c(0, 1), 
# #                                                     labels = c("High Risk", "Low Risk")))) +
# #   geom_point() +
# #   scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
# #   scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,25,50,75,100,125,150,175,200)],
# #                    labels = c(1,25,50,75,100,125,150,175,200),
# #                    expand = c(0.02,0)) +
# #   geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
# #              lty = 2) +
# #   geom_hline(yintercept = median(riskScore_out),
# #              lty =2) +
# #   labs(x = "Patients(increasing risk score)",
# #        y = "Risk Score",
# #        title = "Validation Risk Score Distribution") + 
# #   theme_base() +
# #   theme(legend.title = element_blank(),
# #         legend.key = element_blank(),
# #         legend.justification = c(0,1),
# #         legend.position = c(0,1),
# #         #        legend.margin = margin(c(-5,4,4,3)),
# #         legend.background = element_rect(color = "black", size = .3),
# #         plot.title = element_text(size = 15, hjust = 0.5))
# # risk_dis_out
# # surv_stat_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out),
# #                                       y=OS.time/365,
# #                                       color = factor(OS,
# #                                                      levels = c(0,1),
# #                                                      labels = c("Alive", "Dead")))) +
# #   geom_point() +
# #   scale_color_manual(values = c("#0073C2FF", "#A73030FF")) +
# #   scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,25,50,75,100,125,150,175,200)],
# #                    labels = c(1,25,50,75,100,125,150,175,200),
# #                    expand = c(0.02,0)) +
# #   ylim(x=c(0,15)) +
# #   geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
# #              lty = 2) +
# #   labs(x = "Patients(increasing risk score)",
# #        y = "Survival time (years)",
# #        title = "Validation Survival State Distribution") + 
# #   theme_base() +
# #   theme(legend.title = element_blank(),
# #         legend.key = element_blank(),
# #         legend.justification = c(0,1),
# #         legend.position = c(0,1),
# #         #        legend.margin = margin(c(-5,4,4,3)),
# #         legend.background = element_rect(color = "black", size = .3),
# #         plot.title = element_text(size = 15, hjust = 0.5))
# # surv_stat_out
# # kmfit_out <- survfit(Surv(OS.time, OS) ~ risk_out, data =  risk_out)
# # verify_survival_median <- ggsurvplot(kmfit_out,
# #                                      pval = TRUE, 
# #                                      conf.int = F,
# #                                      legend.labs=c("High risk","Low risk" ),
# #                                      legend.title="Risk score",
# #                                      title="Validition KM",
# #                                      font.main = c(15,"bold"),
# #                                      risk.table = TRUE, 
# #                                      risk.table.col = "strata", 
# #                                      linetype = "strata", 
# #                                      surv.median.line = "hv", 
# #                                      ggtheme = theme_bw(), 
# #                                      palette = c("#A73030FF", "#0073C2FF"))
# # verify_survival_median
# # 
# # multi_ROC <- function(time_vector, risk_score_table){
# #   library(survivalROC)   
# #   single_ROC <- function(single_time){
# #     for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
# #                            status=risk_score_table$OS,
# #                            marker=risk_score_table$riskScore_out,
# #                            predict.time=single_time,method = 'KM')
# #     data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
# #                'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
# #                'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
# #   }
# #   multi_ROC_list <- lapply(time_vector, single_ROC)
# #   do.call(rbind, multi_ROC_list)
# # }
# # 
# # for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
# #                            risk_score_table = risk_out)
# # for_multi_ROC$Time <- factor(for_multi_ROC$Time)
# # 
# # # 画ROC曲线 
# # library(scales)
# # library(geomROC)
# # library(plotROC)
# # auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
# # auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
# # auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)
# # 
# # 
# # ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
# #                                  y=True_positive, 
# #                                  label=Cut_values, 
# #                                  color=Time)) + 
# #   scale_color_manual(breaks = c("365", "1095", "1825"),
# #                      labels = c("1 years", "3 years", "5 years"),
# #                      values = c("#4682B4", "#FF4040", "#20B2AA")) +
# #   geom_roc(labels = F, stat = 'identity') + 
# #   style_roc() + 
# #   geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
# #   theme_bw() +
# #   labs(title = "Train ROC") +
# #   theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
# #         panel.grid = element_blank(),
# #         plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
# #   # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
# #   # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) + 
# #   # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
# #   annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
# #            label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
# #                      paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
# #                      paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
# # ROC
# 
# 
# 
# ##GEO----
# survival<-read.delim2('/data/nas1/luchunlin/project/BJTC-357/00_rawdata/survival.xls')
# survival$OS<-as.numeric(survival$OS)
# survival$OS.time<-as.numeric(survival$OS.time)
# 
# dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-357/00_rawdata/dat2.xls", row.names = 1)%>% lc.tableToNum
# dat<-dat[,colnames(dat)%in%survival$sample]
# #dat<-t(scale(t(dat)))
# #dat<-log2(dat+0.01)
# ##162
# #sig.diff<-read.delim2('/data/nas1/luchunlin/project/BJTC-308/01_DEGs/DEG_sig.xls',row.names = 1)
# test_data<-t(dat)%>%as.data.frame()
# colnames(test_data)<-gsub('-','_',colnames(test_data),fixed = T)
# colnames(test_data)<-gsub(' ','_',colnames(test_data),fixed = T)
# test_data<-test_data[,colnames(test_data)%in%colnames(train_data)]
# test_data$sample<-rownames(test_data)
# test_data<-merge(survival,test_data,by='sample')
# test_data<-column_to_rownames(test_data,var = 'sample')
# # test_data2<-test_data[,colnames(test_data)%in%lasso_geneids]
# test_data2<-test_data[,lasso_geneids$lasso_geneids]
# risk_out<-data.frame(test_data2)
# risk_out$risk_out<-NA
# risk_out$riskScore_out<-NA
# cnt<-1
# 
# coef.min<-coef.min[lasso_geneids,]
# 
# while (cnt < 61) {
#   risk_out$riskScore_out[cnt]<-sum(coef.min*test_data2[cnt,])
#   cnt = cnt + 1
# }
# 
# dim(test_data2)
# cnt<-1
# while (cnt < 61) {
#   risk_out$risk_out[cnt]=as.vector(ifelse(risk_out$riskScore_out[cnt]>median(risk_out$riskScore_out),0,1))
#   cnt = cnt + 1
# }
# riskScore_out<-as.numeric(risk_out$riskScore_out)
# 
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
# ## 52
# risk_dis_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out), 
#                                      y=riskScore_out, 
#                                      color = factor(risk_out, 
#                                                     levels = c(0, 1), 
#                                                     labels = c("High Risk", "Low Risk")))) +
#   geom_point() +
#   scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
#   scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,50,100,150,200,250,300,350,400)],
#                    labels = c(1,50,100,150,200,250,300,350,400),
#                    expand = c(0.02,0)) +
#   geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
#              lty = 2) +
#   geom_hline(yintercept = median(riskScore_out),
#              lty =2) +
#   labs(x = "Patients(increasing risk score)",
#        y = "Risk Score",
#        title = "Validation Risk Score Distribution") + 
#   theme_base() +
#   theme(legend.title = element_blank(),
#         legend.key = element_blank(),
#         legend.justification = c(0,1),
#         legend.position = c(0,1),
#         legend.margin = margin(c(-5,4,4,3)),
#         legend.background = element_rect(color = "black", size = .3),
#         plot.title = element_text(size = 15, hjust = 0.5))
# risk_dis_out
# surv_stat_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out),
#                                       y=OS.time/365,
#                                       color = factor(OS,
#                                                      levels = c(0,1),
#                                                      labels = c("Alive", "Dead")))) +
#   geom_point() +
#   scale_color_manual(values = c("#0073C2FF", "#A73030FF")) +
#   scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,50,100,150,200,250,300,350,400)],
#                    labels = c(1,50,100,150,200,250,300,350,400),
#                    expand = c(0.02,0)) +
#   ylim(x=c(0,25)) +
#   geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
#              lty = 2) +
#   labs(x = "Patients(increasing risk score)",
#        y = "Survival time (years)",
#        title = "Validation Survival State Distribution") + 
#   theme_base() +
#   theme(legend.title = element_blank(),
#         legend.key = element_blank(),
#         legend.justification = c(0,1),
#         legend.position = c(0,1),
#         legend.margin = margin(c(-5,4,4,3)),
#         legend.background = element_rect(color = "black", size = .3),
#         plot.title = element_text(size = 15, hjust = 0.5))
# surv_stat_out
# kmfit_out <- survfit(Surv(OS.time, OS) ~ risk_out, data =  risk_out)
# verify_survival_median <- ggsurvplot(kmfit_out,
#                                      pval = TRUE, 
#                                      conf.int = F,
#                                      legend.labs=c("High risk","Low risk" ),
#                                      legend.title="Risk score",
#                                      title="Validition KM",
#                                      font.main = c(15,"bold"),
#                                      risk.table = TRUE, 
#                                      risk.table.col = "strata", 
#                                      linetype = "strata", 
#                                      surv.median.line = "hv", 
#                                      ggtheme = theme_bw(), 
#                                      palette = c("#A73030FF", "#0073C2FF"))
# 
# verify_survival_median
# 
# 
# library(survivalROC)
# multi_ROC <- function(time_vector, risk_score_table){
#   library(survivalROC)
#   single_ROC <- function(single_time){
#     for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
#                            status=risk_score_table$OS,
#                            marker=risk_score_table$riskScore,
#                            predict.time=single_time,
#                            method = 'KM')
#     data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
#                'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
#                'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
#   }
#   multi_ROC_list <- lapply(time_vector, single_ROC)
#   do.call(rbind, multi_ROC_list)
# }
# for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
#                            risk_score_table = risk_out)
# table(for_multi_ROC$AUC)
# for_multi_ROC$Time <- factor(for_multi_ROC$Time)
# #devtools::install_github('yikeshu0611/geomROC')
# library(scales)
# library(geomROC)
# library(plotROC)
# library(ggthemes)
# 
# auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
# auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
# auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)
# 
# 
# ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
#                                  y=True_positive, 
#                                  label=Cut_values, 
#                                  color=Time)) + 
#   scale_color_manual(breaks = c("365", "1095", "1825"),
#                      labels = c("1 years", "3 years", "5 years"),
#                      values = c("#4682B4", "#FF4040", "#20B2AA")) +
#   geom_roc(labels = F, stat = 'identity') + 
#   style_roc() + 
#   geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
#   theme_bw() +
#   labs(title = "Test ROC") +
#   theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
#         panel.grid = element_blank(),
#         plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
#   # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
#   # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) + 
#   # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
#   annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
#            label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
#                      paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
#                      paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
# ROC
# ggsave('Test ROC.png', ROC,width = 5, height = 4)
# ggsave('Test ROC.pdf', ROC,width = 5, height = 4)
