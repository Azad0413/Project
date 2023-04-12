rm(list = ls())
# 外部验证----------
setwd("/data/nas1/luchunlin/project/JNZK-204(modify)/")
if (! dir.exists("./03_validation")){
  dir.create("./03_validation")
}
setwd("./03_validation")

library(GEOquery)
library(Biobase)
gset<-getGEO("GSE54236",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL6480",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOL'!='')%>%
  separate('GENE_SYMBOL',c('symbol','drop'),sep = '///')%>%
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
#unmap <- lasso_geneids[!lasso_geneids%in%rownames(dat.va)]
# table(pd$`tissue:ch1`)
# pd <- subset(pd,!`tissue:ch1`=='Normal Colon'|`tissue:ch1`=="Normal Liver"|`tissue:ch1`=='Normal Lung')
# #pd <- subset(pd,`tissue:ch1`=='colon tumor')
write.table(dat.va,file = 'dat.va.xls',sep = '\t',quote = F,row.names = T)
survival.va<-data.frame(sample=pd$geo_accession,OS=pd$status,OS.time=pd$`survival time(months):ch1`)
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
#test_data<-scale(test_data)%>%as.data.frame()
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
  scale_color_manual(values = c("orange","darkgreen"), labels = c("High risk","Low risk")) + 
  geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
             lty = 2) +
  geom_hline(yintercept = median(riskScore_out),
             lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score") +
  theme_classic2() + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), legend.position = c(0.1,0.85),
        axis.line = element_line(size = 0),
        panel.border = element_rect(color = "black", fill = "transparent", size = 1.5))
risk_dis_out
ggsave("01.RiskDistribution.png", width = 10, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("01.RiskDistribution.pdf", width = 10, height = 4, units = "in", dpi = 300, bg = "white")

surv_stat_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out),
                                      y=OS.time,
                                      color = factor(OS,
                                                     levels = c(0,1),
                                                     labels = c("Alive", "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("orange","darkgreen"), labels = c("Dead","Alive")) +
  #ylim(x=c(0,20)) +
  geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Survival time (Days)") +
  theme_classic2() + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), legend.position = c(0.1,0.85),
        axis.line = element_line(size = 0), legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(color = "black", fill = "transparent", size = 1.5))
surv_stat_out 
ggsave("02.SurvDistribution.png", width = 10, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("02.SurvDistribution.pdf", width = 10, height = 4, units = "in", dpi = 300, bg = "white")
library(survival)
kmfit_out <- survfit(Surv(OS.time, OS) ~ risk_out, data =  risk_out)
p.surv = ggsurvplot(kmfit_out, conf.int = F, risk.table = TRUE, break.x.by=20, pval = T, 
                    legend.labs = c("High Risk", "Low Risk")) + xlab("Time (Month)")
p.surv
ggsave("01.KM.png", print(p.surv), width = 10, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("01.KM.pdf", print(p.surv), width = 10, height = 7, units = "in", dpi = 300, bg = "white")


library(survivalROC)
library(tidyverse)
rt = subset(risk_out, select = c(OS, OS.time, riskScore_out))
rt$OS.time <- rt$OS.time / 12
survivalROC_helper <- function(t) {
  
  survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore_out, 
              
              predict.time =t, method="KM")
}
survivalROC_data <- data_frame(t = c(2,4,2)) %>%
  
  mutate(survivalROC = map(t, survivalROC_helper),
         
         ## Extract scalar AUC
         
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         
         ## Put cut off dependent values in a data_frame
         
         df_survivalROC = map(survivalROC, function(obj) {
           
           as_data_frame(obj[c("cut.values","TP","FP")])
           
         })) %>%
  
  dplyr::select(-survivalROC) %>%
  
  unnest() %>%
  
  arrange(t, FP, TP)
survivalROC_data1 <- survivalROC_data %>% 
  
  mutate(auc =sprintf("%.2f",auc))%>% 
  
  unite(year, t,auc,sep = " year \n AUC= ")


year =factor(survivalROC_data1$year)
survivalROC_data1 %>%
  
  ggplot(mapping = aes(x = FP, y = TP)) +
  scale_color_manual(
    values = c("black", "black")) +
  
  geom_path(aes(color= year))+
  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  
  facet_wrap( ~ year) +
  
  theme_bw() +
  labs(x = "FP", y = "TP", title = "Test ROC") +
  theme(axis.text.x = element_text(vjust = 0.5),
        
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 13),
        legend.position = "none",
        panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        text = element_text(face = "bold"))

