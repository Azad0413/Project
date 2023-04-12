rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./test(TIICs)")){
  dir.create("./test(TIICs)")
}
setwd("./test(TIICs)")

if (! dir.exists("./test7")){
  dir.create("./test7")
}
setwd("./test7")
###ssGSEA + uncox+lasso
### lasso-----
library(tidyverse)
library(lance)
#dat<-read.delim2('/data/nas1/luchunlin/project/JNZK-214-8/01_ssGSEA/ssgsea_result.xls',row.names = 1)%>%lc.tableToNum()
dat<-read.delim2('/data/nas1/luchunlin/project/JNZK-214-8/01_CIBERSORT/cibersort_result.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
# train_data<-t(dat)%>%as.data.frame()
# survival<-read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-COAD.survival.tsv')
# train_data$sample<-rownames(train_data)
# train_data<-merge(survival,train_data,by="sample")
# train_data<-column_to_rownames(train_data,var = 'sample')%>%select(-'X_PATIENT')
# ## 单因素cox--------
library(tidyverse)
library(lance)
#lasso.gene<-read.csv('/data/nas1/luchunlin/project/JNZK-214-8/02_Lasso(TIIC)/lasso_genes.csv',header = F)
#dat<-read.delim2('../01_CIBERSORT/cibersort_result.xls',row.names = 1)%>%lc.tableToNum()
# dat<-read.delim2('/data/nas1/luchunlin/project/JNZK-214-8/01_ssGSEA/ssgsea_result.xls',row.names = 1)%>%lc.tableToNum()
# colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group<-read.delim2('/data/nas1/luchunlin/project/JNZK-214-8/01_ssGSEA/group.xls')
tumor.sample<-group$sample[which(group$group=='Tumor')]
dat<-dat[,tumor.sample]
train_data<-t(dat)%>%as.data.frame()
# #train_data<-train_data[,lasso.gene$V1]
survival<-read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-COAD.survival.tsv')
train_data$sample<-rownames(train_data)
train_data<-merge(survival,train_data,by="sample")
train_data<-column_to_rownames(train_data,var = 'sample')%>%select(-'X_PATIENT')

library(survival)
library(survminer)
colnames_sum <- colnames(train_data)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames_sum<-gsub('(','',colnames_sum,fixed = T)
colnames_sum<-gsub(')','',colnames_sum,fixed = T)
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
res_results_0.2 <- res_mod[which(as.numeric(res_mod$p.value) < 0.2),]
res_results_0.2 <- na.omit(res_results_0.2)
# 3
write.table(res_results_0.2,
            file = "univariate_cox_result_0.2.xls",
            quote = F,
            sep = '\t',
            row.names = T)
dim(res_results_0.2)
## 3 2
library(tidyr)
res_results_0.2_2 <- separate(res_results_0.2, "HR (95% CI for HR)",
                              into = c("HR", "HR.95L", "HR.95H"),
                              sep = " ")
res_results_0.2_2 <- separate(res_results_0.2_2, "HR.95L",
                              into = c("HR.95L", "HR.95H"),
                              sep = "\\-")
res_results_0.2_2$HR.95L <- gsub("\\(", "", res_results_0.2_2$HR.95L)
res_results_0.2_2$HR.95H <- gsub("\\)", "", res_results_0.2_2$HR.95H)

res_results_0.2_2[,1:ncol(res_results_0.2_2)] <- as.numeric(unlist(res_results_0.2_2[,1:ncol(res_results_0.2_2)]))
res_results_0.2_2 <- res_results_0.2_2[order(res_results_0.2_2$HR),]
hz <- paste(round(res_results_0.2_2$HR,4),
            "(",round(res_results_0.2_2$HR.95L,4),
            "-",round(res_results_0.2_2$HR.95H,4),")",sep = "")


tabletext <- cbind(c(NA,"Gene",rownames(res_results_0.2_2)),
                   c(NA,"P value",ifelse(res_results_0.2_2$p.value<0.001,
                                         "< 0.001",
                                         round(res_results_0.2_2$p.value,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)

pdf(file = "univariate_cox_forest.pdf", height = 9, width = 9, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE,rep(FALSE, 70)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.2_2$HR),
           lower=c(NA,NA,res_results_0.2_2$HR.95L), #95%置信区间下限
           upper=c(NA,NA,res_results_0.2_2$HR.95H), #95%置信区间上限
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
png(filename = "univariate_cox_forest.png", height = 600, width = 600)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE,rep(FALSE, 70)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.2_2$HR),
           lower=c(NA,NA,res_results_0.2_2$HR.95L), #95%置信区间下限
           upper=c(NA,NA,res_results_0.2_2$HR.95H), #95%置信区间上限
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
           # hrzl_lines=list(S"2" = gpar(lwd=2, col="black"),
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
##lasso-----
library(glmnet)
x_all <- subset(train_data, select = -c(OS, OS.time))
x_all <- x_all[,rownames(res_results_0.2)]
y_all <- subset(train_data, select = c(OS, OS.time))
# 拟合模型
library(survival)
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
set.seed(2)
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
cvfit$lambda.min
active.min = which(coef.min@i != 0)
coef.min
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids
write(lasso_geneids, "lasso_genes.csv")


# 05 风险模型的构建与验证------
# 
# cox_more_2$coefficients
# riskScore=predict(cox_more_2,type="lp",newdata=train_data)
# coxGene=rownames(as.data.frame(cox_more_2$coefficients))
# outCol=c("OS","OS.time",coxGene)
# risk=as.vector(ifelse(riskScore>median(riskScore),0,1))
# risk <- as.data.frame(c(cbind(id=rownames(cbind(train_data[,outCol],
#                                                 riskScore,
#                                                 risk)),
#                               cbind(train_data[,outCol],
#                                     riskScore,
#                                     risk))))
# table(risk$risk)
# 0   1 
# 224 224  
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
##最佳截断值
# res.cut<-surv_cutpoint(risk,time = 'OS.time',event = 'OS',variables = 'riskScore')
# summary(res.cut)
# cutpoint<-res.cut$cutpoint$cutpoint
# risk$risk<-ifelse(risk$riskScore>cutpoint,0,1)

# write.table(risk,file = 'risk(TIIC).xls',quote = F,row.names = F,sep = '\t')

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
  ylim(c(0,20))+
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
  scale_color_manual(breaks = c("365","1095", "1825"),
                     labels = c("1 years", "3 years","5 years"),
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
ggsave('Train ROC.png', ROC,width = 5, height = 4)
ggsave('Train ROC.pdf', ROC,width = 5, height = 4)

