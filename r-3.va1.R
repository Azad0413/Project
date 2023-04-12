rm(list = ls())
# 外部验证----------
setwd("/data/nas1/luchunlin/project/JNZK-204(modify)/")
if (! dir.exists("./03_validation")){
  dir.create("./03_validation")
}
setwd("./03_validation")

library(GEOquery)
library(Biobase)
gset<-getGEO("GSE15654",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL8432",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Symbol')%>%
  filter('Symbol'!='')%>%
  separate('Symbol',c('symbol','drop'),sep = '///')%>%
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
unmap <- lasso_geneids[!lasso_geneids%in%rownames(dat.va)]
# table(pd$`tissue:ch1`)
# pd <- subset(pd,!`tissue:ch1`=='Normal Colon'|`tissue:ch1`=="Normal Liver"|`tissue:ch1`=='Normal Lung')
# #pd <- subset(pd,`tissue:ch1`=='colon tumor')
write.table(dat.va,file = 'dat.va.xls',sep = '\t',quote = F,row.names = T)
survival.va<-data.frame(sample=pd$geo_accession,OS=pd$characteristics_ch1.8,OS.time=pd$`days to death:ch1`)
#survival.va<-data.frame(sample=pd$geo_accession,OS=pd$`os censor:ch1`,OS.time=pd$`survival (mo):ch1`)
#survival.va <- subset(survival.va,OS=='AWD'|OS=='AUN'|OS=='DOC'|OS=='DOD'|OS=='DUN')
#survival.va <- subset(survival.va,OS=='AWD'|OS=='AUN'|OS=='DOD')
table(survival.va$OS)
#survival.va$OS[survival.va$OS=='N/A'] <- NA
survival.va <- na.omit(survival.va)
#survival.va<-survival.va[which(!is.na(survival.va$OS)),]
survival.va$OS<-ifelse(survival.va$OS=='death: 0',0,1)
survival.va$OS.time<-as.numeric(survival.va$OS.time)
survival.va <- na.omit(survival.va)
write.table(survival.va,file = 'survival.va.xls',sep = '\t',row.names = F,quote = F)
survival.va$OS<-as.numeric(survival.va$OS)
##GALNTL6 "GALNT15"=====GALNTL17 GALNTL2"
test_data<-t(dat.va[,survival.va$sample])%>%as.data.frame()
#test_data <- na.omit(test_data)
#test_data<-log2(test_data+0.001)
test_data<-scale(test_data)%>%as.data.frame()
test_data$sample<-rownames(test_data)
test_data<-merge(survival.va,test_data,by='sample')
test_data<-column_to_rownames(test_data,var='sample')

df.coef = read.delim2("../02.RiskScore/05.Coefficients.xls")
lasso_geneids <- df.coef$gene
## LASSO------
#unmap <- lasso_geneids[!lasso_geneids%in%colnames(test_data)]
##ADH1A
test_data2<-test_data[,lasso_geneids]%>%lc.tableToNum()
risk_out<-data.frame(test_data2)
risk_out$risk_out<-NA
risk_out$riskScore_out<-NA
cnt<-1

coef.min<-as.numeric(df.coef$coefficient)

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
coxGene=lasso_geneids
outCol=c("OS","OS.time",coxGene)
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
