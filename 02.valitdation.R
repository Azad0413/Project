rm(list=ls())

##验证---------
# library(TCGAbiolinks)
# # 获取TCGA中最新的不同癌种的项目号
# getGDCprojects()$project_id
# query <- GDCquery(project = "CGCI-HTMCP-CC",
#                   data.category = "Transcriptome Profiling",
#                   data.type = "Gene Expression Quantification",
#                   workflow.type = "STAR - Counts")
# 
# GDCdownload(query)
# GDCprepare(query,save=T,save.filename = "HTMCP-CC_mRNA.Rdata")
# load("HTMCP-CC_mRNA.Rdata")
# e_mrna <- data[rowData(data)$gene_type == "protein_coding",]


library(GEOquery)
## GSE24673----
gset<-getGEO("GSE44001",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a<-gset[[1]]
gpl<-getGEO("GPL14951",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','Symbol')%>%
  filter('Symbol'!='')%>%
  separate('Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symbol=probe2symbol[probe2symbol$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
## 敏感和耐药
write.table(dat,file = 'dat.va.xls',sep = '\t',quote = F,row.names = T)
survival.va<-data.frame(sample=pd$geo_accession,DFS=pd$`status_of_dfs:ch1`,DFS.time=pd$`disease_free_survival_(dfs)_(months):ch1`)
#survival.va<-data.frame(sample=pd$geo_accession,OS=pd$`os censor:ch1`,OS.time=pd$`survival (mo):ch1`)
table(survival.va$DFS)
#survival.va<-survival.va[which(!is.na(survival.va$DFS)),]
#survival.va$DFS<-ifelse(survival.va$DFS=='D',1,0)
survival.va$DFS.time<-as.numeric(survival.va$DFS.time)*30
write.table(survival.va,file = 'survival.va.xls',sep = '\t',row.names = F,quote = F)
survival.va$DFS<-as.numeric(survival.va$DFS)
riskmodel<-read.table('/data/nas1/tancj/Project/BJTC-202_repair/multiCoxtest.xls',header=T,row.names=1)

dat <- dat[rownames(riskmodel),]

test_data<-t(dat[,survival.va$sample])%>%as.data.frame()
#test_data<-log2(test_data+0.01)
#test_data<-scale(test_data)%>%as.data.frame()
test_data$sample<-rownames(test_data)
test_data<-merge(survival.va,test_data,by='sample')
test_data<-column_to_rownames(test_data,var='sample')
#input<-read.table('/data/nas1/tancj/Project/BJTC-202_repair/input.txt',header=T,row.names=1)
# colnames(input


astive.coefficients<-riskmodel$coef
lasso_geneids<-rownames(riskmodel)
rt2 <- test_data[,c(3:6)]
rt1<-test_data[,c(1:2)]
risk<-data.frame(rt1)
risk$riskScore<-NA
risk$risk<-NA
cnt<-1
while (cnt <301) {
  risk$riskScore[cnt]<-sum(astive.coefficients*subset(rt2,select=lasso_geneids)[cnt,])
  cnt = cnt + 1
}
cnt<-1
while (cnt <301) {
  risk$risk[cnt]=as.vector(ifelse( risk$riskScore[cnt]>median(risk$riskScore),0,1))
  cnt = cnt + 1
}
risk$Risk<-ifelse(risk$risk==1,'Low','High')
riskScore <- risk$riskScore
table(risk$risk)
# 0   1 
# 150 150 
write.table(risk,file = 'risk.xls',quote = F,row.names = F,sep = '\t')
risk$riskScore <- as.vector(risk$riskScore)
risk <- risk%>%rownames_to_column(var = 'id')

library(ggplot2)
library(ggthemes)
library(Ipaper)
median(riskScore)
# [1] -0.02316101
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
ggsave(filename = '01.Validation_Risk_Score_Distribution.pdf',risk_dis,w=7,h=5)
ggsave(filename = '01.Validation_Risk_Score_Distribution.png',risk_dis,w=7,h=5)
surv_stat <- ggplot(risk, aes(x=reorder(id, riskScore),
                              y=DFS.time/365,
                              color = factor(DFS,
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
       y = "DFS time(Years)",
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

ggsave(filename = '02.Validation_survival_state_Distribution.pdf',surv_stat,w=7,h=5)
ggsave(filename = '02.Validation_survival_state_Distribution.png',surv_stat,w=7,h=5)
## 05-2 KM曲线和ROC曲线------
library(survival)
library(survminer)
kmfit<-survfit(Surv(DFS.time, DFS) ~ risk, data =  risk)
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

riskScore <- function(survival_cancer_df,
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
    cbind(survival_cancer_df[,c("DFS.time", "DFS")])
  risk_score_table <- risk_score_table[,c("DFS.time",
                                          "DFS",
                                          candidate_genes_for_cox,
                                          "total_risk_score")]
}

multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$DFS.time,
                           status=risk_score_table$DFS,
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
ggsave('04.Test ROC.png', ROC,width = 5, height = 4)
ggsave('04.Test ROC.pdf', ROC,width = 5, height = 4)

write.table(pd,file = 'all.clinical.xls',sep = '\t',row.names = F,quote = F)
