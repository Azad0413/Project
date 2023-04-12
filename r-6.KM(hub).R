## 生存分析-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/LZZK-504")
if (! dir.exists("./07_survival")){
  dir.create("./07_survival")
}
setwd("./07_survival")
## 匹配生存数据
dat<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames<-data.frame(sample=colnames(dat))
colnames$sample<-gsub('.','-',colnames$sample,fixed = T)
colnames(dat)<-colnames$sample
survival<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/03_KM/TCGA-OV.survival.tsv') 
survival<-survival[survival$sample%in%colnames(dat),]
hubgene<-c("HIST1H3I","HIST1H1B","HIST1H2AH","HIST1H4F","HIST1H1A","HIST1H2BB","HIST1H3C","HIST1H4L","HIST1H2BM","HIST2H2AB","HIST1H2BI","HIST1H2AJ")
write.table(hubgene,file = 'hubgene.xls',sep = '\t',quote = F,row.names = F)
dat<-t(dat)%>%as.data.frame()
## KM（HIST1H3I）-----
group1<-data.frame(sample=rownames(dat),group=ifelse(dat$HIST1H3I>median(dat$HIST1H3I),'High HIST1H3I','Low HIST1H3I'))
write.table(group1,file = 'group(HIST1H3I).xls',row.names = F,sep = '\t',quote = F)
group1$group<-as.vector(group1$group)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group1$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median1 <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High HIST1H3I","Low HIST1H3I" ),
                                      legend.title="group",
                                      title="KM(HIST1H3I)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median1

## HIST1H1B------
group2<-data.frame(sample=rownames(dat),group=ifelse(dat$HIST1H1B>median(dat$HIST1H1B),'High HIST1H1B','Low HIST1H1B'))
group2$group<-as.vector(group2$group)
write.table(group2,file = 'group(HIST1H1B).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group2$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median2 <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High HIST1H1B","Low HIST1H1B" ),
                                      legend.title="group",
                                      title="KM(HIST1H1B)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median2
## HIST1H2AH------
group3<-data.frame(sample=rownames(dat),group=ifelse(dat$HIST1H2AH>median(dat$HIST1H2AH),'High HIST1H2AH','Low HIST1H2AH'))
group3$group<-as.vector(group3$group)
write.table(group3,file = 'group(HIST1H2AH).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group3$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median3 <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High HIST1H2AH","Low HIST1H2AH" ),
                                      legend.title="group",
                                      title="KM(HIST1H2AH)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median3
## HIST1H4F-----
group4<-data.frame(sample=rownames(dat),group=ifelse(dat$HIST1H4F>median(dat$HIST1H4F),'High HIST1H4F','Low HIST1H4F'))
group4$group<-as.vector(group4$group)
write.table(group4,file = 'group(HIST1H4F).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group4$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median4 <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High HIST1H4F","Low HIST1H4F" ),
                                      legend.title="group",
                                      title="KM(HIST1H4F)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median4
## HIST1H1A-----
group5<-data.frame(sample=rownames(dat),group=ifelse(dat$HIST1H1A>median(dat$HIST1H1A),'High HIST1H1A','Low HIST1H1A'))
group5$group<-as.vector(group5$group)
write.table(group5,file = 'group(HIST1H1A).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group5$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median5 <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High HIST1H1A","Low HIST1H1A" ),
                                      legend.title="group",
                                      title="KM(HIST1H1A)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median5
## HIST1H2BB-------
group6<-data.frame(sample=rownames(dat),group=ifelse(dat$HIST1H2BB>median(dat$HIST1H2BB),'High HIST1H2BB','Low HIST1H2BB'))
group6$group<-as.vector(group6$group)
write.table(group6,file = 'group(HIST1H2BB).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group6$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median6 <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High HIST1H2BB","Low HIST1H2BB" ),
                                      legend.title="group",
                                      title="KM(HIST1H2BB)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median6
## 7HIST1H2BB-------
group7<-data.frame(sample=rownames(dat),group=ifelse(dat$HIST1H3C>median(dat$HIST1H3C),'High HIST1H3C','Low HIST1H3C'))
group7$group<-as.vector(group7$group)
write.table(group7,file = 'group(HIST1H3C).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group7$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median7 <- ggsurvplot(kmfit,
                                       pval = TRUE, 
                                       conf.int = F,
                                       legend.labs=c("High HIST1H3C","Low HIST1H3C" ),
                                       legend.title="group",
                                       title="KM(HIST1H3C)",
                                       font.main = c(15,"bold"),
                                       risk.table = TRUE, 
                                       risk.table.col = "strata", 
                                       linetype = "strata", 
                                       surv.median.line = "hv", 
                                       ggtheme = theme_bw(), 
                                       palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median7
## 8HIST1H2BB-------
group8<-data.frame(sample=rownames(dat),group=ifelse(dat$HIST1H4L>median(dat$HIST1H4L),'High HIST1H4L','Low HIST1H4L'))
group8$group<-as.vector(group8$group)
write.table(group8,file = 'group(HIST1H4L).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group8$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median8 <- ggsurvplot(kmfit,
                                       pval = TRUE, 
                                       conf.int = F,
                                       legend.labs=c("High HIST1H4L","Low HIST1H4L" ),
                                       legend.title="group",
                                       title="KM(HIST1H4L)",
                                       font.main = c(15,"bold"),
                                       risk.table = TRUE, 
                                       risk.table.col = "strata", 
                                       linetype = "strata", 
                                       surv.median.line = "hv", 
                                       ggtheme = theme_bw(), 
                                       palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median8

## 9 HIST1H2BB-------
group9<-data.frame(sample=rownames(dat),group=ifelse(dat$HIST1H2BM>median(dat$HIST1H2BM),'High HIST1H2BM','Low HIST1H2BM'))
group9$group<-as.vector(group9$group)
write.table(group9,file = 'group(HIST1H2BM).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group9$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median9 <- ggsurvplot(kmfit,
                                       pval = TRUE, 
                                       conf.int = F,
                                       legend.labs=c("High HIST1H2BM","Low HIST1H2BM" ),
                                       legend.title="group",
                                       title="KM(HIST1H2BM)",
                                       font.main = c(15,"bold"),
                                       risk.table = TRUE, 
                                       risk.table.col = "strata", 
                                       linetype = "strata", 
                                       surv.median.line = "hv", 
                                       ggtheme = theme_bw(), 
                                       palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median9

## 10 HIST1H2BB-------
group10<-data.frame(sample=rownames(dat),group=ifelse(dat$HIST2H2AB>median(dat$HIST2H2AB),'High HIST2H2AB','Low HIST2H2AB'))
group10$group<-as.vector(group10$group)
write.table(group10,file = 'group(HIST2H2AB).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group10$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median10 <- ggsurvplot(kmfit,
                                       pval = TRUE, 
                                       conf.int = F,
                                       legend.labs=c("High HIST2H2AB","Low HIST2H2AB" ),
                                       legend.title="group",
                                       title="KM(HIST2H2AB)",
                                       font.main = c(15,"bold"),
                                       risk.table = TRUE, 
                                       risk.table.col = "strata", 
                                       linetype = "strata", 
                                       surv.median.line = "hv", 
                                       ggtheme = theme_bw(), 
                                       palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median10

## 11 HIST1H2BB-------
group11<-data.frame(sample=rownames(dat),group=ifelse(dat$HIST1H2BI>median(dat$HIST1H2BI),'High HIST1H2BI','Low HIST1H2BI'))
group11$group<-as.vector(group11$group)
write.table(group11,file = 'group(HIST1H2BI).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group11$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median11 <- ggsurvplot(kmfit,
                                        pval = TRUE, 
                                        conf.int = F,
                                        legend.labs=c("High HIST1H2BI","Low HIST1H2BI" ),
                                        legend.title="group",
                                        title="KM(HIST1H2BI)",
                                        font.main = c(15,"bold"),
                                        risk.table = TRUE, 
                                        risk.table.col = "strata", 
                                        linetype = "strata", 
                                        surv.median.line = "hv", 
                                        ggtheme = theme_bw(), 
                                        palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median11

## 12 HIST1H2BB-------
group12<-data.frame(sample=rownames(dat),group=ifelse(dat$HIST1H2AJ>median(dat$HIST1H2AJ),'High HIST1H2AJ','Low HIST1H2AJ'))
group12$group<-as.vector(group12$group)
write.table(group12,file = 'group(HIST1H2AJ).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group12$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median12 <- ggsurvplot(kmfit,
                                        pval = TRUE, 
                                        conf.int = F,
                                        legend.labs=c("High HIST1H2AJ","Low HIST1H2AJ" ),
                                        legend.title="group",
                                        title="KM(HIST1H2AJ)",
                                        font.main = c(15,"bold"),
                                        risk.table = TRUE, 
                                        risk.table.col = "strata", 
                                        linetype = "strata", 
                                        surv.median.line = "hv", 
                                        ggtheme = theme_bw(), 
                                        palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median12

