
# 04 风险模型构建---------
## 04-1 单因素cox回归----
setwd("/data/nas1/luchunlin/project/JNZK-207/")
if (! dir.exists("./10_supplement")){
  dir.create("./10_supplement")
}
setwd("./10_supplement")
library(GEOquery)
library(Biobase)
# gset_va<-getGEO("GSE136337",
#                 destdir = '.',
#                 GSEMatrix = T,
#                 getGPL = F)
# expr_va<-as.data.frame(exprs(gset_va[[1]]))

# gpl1<-getGEO("GPL570",destdir = '.')
# a1=gset_va[[1]]
# gpl1<-Table(gpl1)
# colnames(gpl1)
# probe2symobl1<-gpl1 %>%
#   select('ID','Gene Symbol')%>%
#   filter('Gene Symbol'!='')%>%
#   separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
#   select(-drop)
# probe2symobl1=probe2symobl1[probe2symobl1$symbol!='',]
# 
# a_va=gset_va[[1]]
# pd_va<-pData(a_va)
# dat_va<-expr_va
# dat_va$ID<-rownames(dat_va)
# dat_va$ID<-as.character(dat_va$ID)
# probe2symobl1$ID<-as.character(probe2symobl1$ID)
# dat_va<-dat_va %>%
#   inner_join(probe2symobl1,by='ID')%>%
#   select(-ID)%>%     ## 去除多余信息
#   select(symbol,everything())%>%     ## 重新排列
#   mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
#   arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
#   distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
#   select(-rowMean)%>%     ## 反向选择去除rowMean这一列
#   tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
# dat_final_va<-dat_va
# survival_va<-data.frame(sample=pd_va$geo_accession,
#                         OS=pd_va$`censos.ar:ch1`,
#                         OS.time=pd_va$`monthsos.ar:ch1`)
# table(survival_va$OS)
# survival_va$OS <- ifelse(survival_va$OS=='FALSE',0,1)
# table(survival_va$OS.time)
# survival_va$OS <- as.numeric(survival_va$OS)
# survival_va$OS.time <- as.numeric(survival_va$OS.time)*30
# 
# write.table(dat_va,file = 'dat(GSE136337).xls',sep = '\t',row.names = T,quote = F)
# write.table(survival_va,file = 'survival(GSE136337).xls',sep = '\t',row.names = F,quote = F)
# 

##GSE2658
library(GEOquery)
library(Biobase)
gset_va<-getGEO("GSE2658",
                destdir = '.',
                GSEMatrix = T,
                getGPL = F)
expr_va<-as.data.frame(exprs(gset_va[[1]]))

gpl1<-getGEO("GPL570",destdir = '.')
a1=gset_va[[1]]
gpl1<-Table(gpl1)
colnames(gpl1)
probe2symobl1<-gpl1 %>%
  select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  select(-drop)
probe2symobl1=probe2symobl1[probe2symobl1$symbol!='',]

a_va=gset_va[[1]]
pd_va<-pData(a_va)
dat_va<-expr_va
dat_va$ID<-rownames(dat_va)
dat_va$ID<-as.character(dat_va$ID)
probe2symobl1$ID<-as.character(probe2symobl1$ID)
dat_va<-dat_va %>%
  inner_join(probe2symobl1,by='ID')%>%
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

survival_va<-data.frame(sample=pd_va$geo_accession,
                        OS=pd_va$characteristics_ch1,
                        OS.time=pd_va$characteristics_ch1.2)
table(survival_va$OS)
survival_va$OS <- ifelse(survival_va$OS=='[SURIND=0 (Indicator of disease-related death; integer, 0=alive or death by other cause, 1=disease related death, na=death cause undetermined)]',0,1)
table(survival_va$OS.time)
survival_va$OS.time <- gsub('[SURTIM=','',survival_va$OS.time,fixed = T)
survival_va$OS.time <- gsub(' (Follow-up time in months from Pre-Treatment baseline; integer)]','',survival_va$OS.time,fixed = T)
survival_va$OS <- as.numeric(survival_va$OS)
survival_va$OS.time <- as.numeric(survival_va$OS.time)*30

# write.table(dat_va,file = 'dat(GSE2658).xls',sep = '\t',row.names = T,quote = F)
write.table(survival_va,file = 'survival(GSE2658).xls',sep = '\t',row.names = F,quote = F)

test_data<-t(dat_va[,survival.va$sample])%>%as.data.frame()
test_data<-log2(test_data+1)
#test_data<-scale(test_data)%>%as.data.frame()
test_data$sample<-rownames(test_data)
test_data<-merge(survival.va,test_data,by='sample')
test_data<-column_to_rownames(test_data,var='sample')

###多因素-----
# riskScore_out = predict(cox_more_2,type="lp",newdata=test_data)
# risk_out=as.vector(ifelse(riskScore_out>median(riskScore_out),0,1))
# risk_out <- as.data.frame(c(cbind(id=rownames(cbind(test_data[,outCol],
#                                                     riskScore_out,
#                                                     risk_out)),
#                                   cbind(test_data[,outCol],
#                                         riskScore_out,
#                                         risk_out))))

## LASSO------
test_data2<-test_data[,lasso_geneids$lasso_geneids]
write.table(test_data2,file = 'dat(GSE2658).xls',sep = '\t',row.names = T,quote = F)
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
# 11.88082
table(risk_out$risk_out)
library(ggplot2)
library(ggthemes)
median(riskScore_out2)
risk_dis_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out), 
                                     y=riskScore_out, 
                                     color = factor(risk_out, 
                                                    levels = c(0, 1), 
                                                    labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,50,100,150,200,250,300,350,400,450,500,550,600)],
                   labels = c(1,50,100,150,200,250,300,350,400,450,500,550,600),
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
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,50,100,150,200,250,300,350,400,450,500,550,600)],
                   labels = c(1,50,100,150,200,250,300,350,400,450,500,550,600),
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
surv_pvalue(kmfit_out,method = 'Gehan-Breslow')
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
#auc_y2 <- round(for_multi_ROC[which(for_multi_ROC$Time==730),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365","1095","1825"),
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



###2--------
gset_va<-getGEO("GSE57317",
                destdir = '.',
                GSEMatrix = T,
                getGPL = F)
expr_va<-as.data.frame(exprs(gset_va[[1]]))
a_va=gset_va[[1]]
pd_va<-pData(a_va)
dat_va<-expr_va
dat_va$ID<-rownames(dat_va)
dat_va$ID<-as.character(dat_va$ID)
probe2symobl1$ID<-as.character(probe2symobl1$ID)
dat_va<-dat_va %>%
  inner_join(probe2symobl1,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dat_final_va<-dat_va
survival_va<-data.frame(sample=pd_va$geo_accession,
                        OS=as.numeric(pd_va$`os censored:ch1`),
                        OS.time=as.numeric(pd_va$`OS time:ch1`))


write.table(dat_va,file = 'dat(GSE57317).xls',sep = '\t',row.names = T,quote = F)
write.table(survival_va,file = 'survival(GSE57317).xls',sep = '\t',row.names = F,quote = F)

test_data<-t(dat_va[,survival_va$sample])%>%as.data.frame()
test_data<-log2(test_data+1)
#test_data<-scale(test_data)%>%as.data.frame()
test_data$sample<-rownames(test_data)
test_data<-merge(survival_va,test_data,by='sample')
test_data<-column_to_rownames(test_data,var='sample')

###多因素-----
# riskScore_out = predict(cox_more_2,type="lp",newdata=test_data)
# risk_out=as.vector(ifelse(riskScore_out>median(riskScore_out),0,1))
# risk_out <- as.data.frame(c(cbind(id=rownames(cbind(test_data[,outCol],
#                                                     riskScore_out,
#                                                     risk_out)),
#                                   cbind(test_data[,outCol],
#                                         riskScore_out,
#                                         risk_out))))

## LASSO------
test_data2<-test_data[,lasso_geneids$lasso_geneids]
write.table(test_data2,file = 'dat(GSE57317).xls',sep = '\t',row.names = T,quote = F)
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
write.table(risk_out,file = 'risk_out(GSE57317).xls',quote = F,row.names = F,sep = '\t')


library(ggplot2)
library(ggthemes)
median(riskScore_out)
# 12.54615
table(risk_out$risk_out)
library(ggplot2)
library(ggthemes)
median(riskScore_out2)
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
ggsave(filename = '01.Validation_Risk_Score_Distribution(GSE57317).pdf',risk_dis_out,w=7,h=5)
ggsave(filename = '01.Validation_Risk_Score_Distribution(GSE57317)..png',risk_dis_out,w=7,h=5)
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
ggsave(filename = '02.Validation_survival_state_Distribution(GSE57317).pdf',surv_stat_out,w=7,h=5)
ggsave(filename = '02.Validation_survival_state_Distribution(GSE57317).png',surv_stat_out,w=7,h=5)
kmfit_out <- survfit(Surv(OS.time, OS) ~ risk_out, data =  risk_out)
surv_pvalue(kmfit_out,method = 'Gehan-Breslow')
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
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,3,2)), 
                           risk_score_table = risk_out)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)
#devtools::install_github('yikeshu0611/geomROC')
library(scales)
library(geomROC)
library(plotROC)
library(ggthemes)

auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
#auc_y2 <- round(for_multi_ROC[which(for_multi_ROC$Time==730),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365","1095"),
                     labels = c("1 years", "3 years"),
                     values = c("#4682B4", "#FF4040")) +
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
  annotate("text", x=0.75, y=c(0.25, 0.15),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2))))
ROC
ggsave('04.Test_ROC(GSE57317).png', ROC,width = 5, height = 4)
ggsave('04.Test_ROC(GSE57317).pdf', ROC,width = 5, height = 4)
