## 甲状腺癌------
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-320")
if (! dir.exists("./11_cancer")){
  dir.create("./11_cancer")
}
setwd("./11_cancer")
library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)
pheno<-read_tsv('TCGA-THCA.GDC_phenotype.tsv.gz')
## 读取从xena下载的数据
tcga.expr<-read_tsv(file = 'TCGA-THCA.htseq_fpkm.tsv.gz')
tcga.expr<-as.data.frame(tcga.expr)
rownames(tcga.expr)<-tcga.expr[,1]
tcga.expr<-tcga.expr[,-1]
## xena下载的数据经过了log2+1转化，需要将其还原
tcga.expr<-2^tcga.expr-1
## 对数据进行id转化
genecode<-read.table(file = 'gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat.tcga<-tcga.expr
dat.tcga$ID <- rownames(dat.tcga)
dat.tcga$ID<-as.character(dat.tcga$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat.tcga<-dat.tcga %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat.tcga)
## 筛选癌症组织，去掉癌旁组织。01-09为肿瘤，10-19为正常对照
mete=data.frame(colnames(dat.tcga))  # 取第一行样本id
for (i in 1:length(mete[,1])) {
  num=as.numeric(as.character(substring(mete[i,1],14,15)))
  if(num %in% seq(1,9)){mete[i,2]="T"}
  if(num %in% seq(10,29)){mete[i,2]="N"}
}
names(mete)=c("id","group")
mete$group=as.factor(mete$group)
mete=subset(mete,mete$group=="T")
exp_tumor<-dat.tcga[,which(colnames(dat.tcga)%in%mete$id)]
exp_tumor<-as.data.frame(exp_tumor)
# 510
exp_control<-dat.tcga[,which(!colnames(dat.tcga)%in%mete$id)]
exp_control<-as.data.frame(exp_control)
control.sample<-colnames(exp_control)
tumor.sample<-colnames(exp_tumor)
# 58
dat.final<-cbind(exp_control,exp_tumor)
pheno<-read_tsv('TCGA-THCA.GDC_phenotype.tsv.gz')
pheno<-pheno[pheno$submitter_id.samples%in%colnames(dat.final),]
table(pheno$primary_diagnosis.diagnoses)
pheno<-subset(pheno,primary_diagnosis.diagnoses=='Papillary adenocarcinoma, NOS'|primary_diagnosis.diagnoses=='Papillary carcinoma, columnar cell'|primary_diagnosis.diagnoses=='Papillary carcinoma, follicular variant')
dat.final<-dat.final[,colnames(dat.final)%in%pheno$submitter_id.samples]

write.table(dat.final,file = 'dat.tcga.xls',sep = '\t',quote = F,row.names = T)
clinical<-GDCquery_clinic(project = "TCGA-THCA",type = "clinical")
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-320/06_expression/hubgene.xls')
hubgene<-hubgene[-1,]%>%as.data.frame()
hub_exp<-dat.final[hubgene$.,]
hub_exp<-log2(hub_exp+1)
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','PTC')

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
exp_plot <- ggplot(hub_exp2,aes(x = Symbol, y = expr, fill = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 't.test') +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+facet_wrap(~Symbol,scales = "free",nrow = 3) 
exp_plot
ggsave(filename = 'expression.pdf',w=8,h=8)
ggsave(filename = 'expression.png',w=8,h=8)
### PFI（无进展生存期）-----
survival <- read_xlsx("/data/nas1/luchunlin/project/BJTC-320/11_cancer/THCA_survival_dat.xlsx") %>% 
  subset(select = c("sample", "PFI", "PFI.time"))
##TNFAIP3
exp_tumor<-exp_tumor[,colnames(exp_tumor)%in%pheno$submitter_id.samples]
km.dat<-t(exp_tumor)%>%as.data.frame()
km.dat<-log2(km.dat+1)
group<-data.frame(sample=rownames(km.dat),group=ifelse(km.dat$TNFAIP3>median(km.dat$TNFAIP3),'High TNFAIP3','Low TNFAIP3'))
group$group<-as.vector(group$group)
km.dat$group<-as.vector(group$group)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
km.dat<-column_to_rownames(km.dat,var = 'sample')
table(km.dat$group)
library(survival)
library(survminer)
kmfit<-survfit(Surv(PFI.time, PFI) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High TNFAIP3","Low TNFAIP3"),
                                      legend.title="group",
                                      title="TNFAIP3 KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF","#0073C2FF"))

cluster_survival_median
## PTK2B
km.dat<-t(exp_tumor)%>%as.data.frame()
km.dat<-log2(km.dat+1)
group<-data.frame(sample=rownames(km.dat),group=ifelse(km.dat$PTK2B>median(km.dat$PTK2B),'High PTK2B','Low PTK2B'))
group$group<-as.vector(group$group)

km.dat$group<-as.vector(group$group)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
km.dat<-column_to_rownames(km.dat,var = 'sample')
table(km.dat$group)
kmfit<-survfit(Surv(PFI.time, PFI) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High PTK2B","Low PTK2B"),
                                      legend.title="group",
                                      title="PTK2B KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF","#0073C2FF"))
cluster_survival_median
## STAT1
km.dat<-t(exp_tumor)%>%as.data.frame()
km.dat<-log2(km.dat+1)
group<-data.frame(sample=rownames(km.dat),group=ifelse(km.dat$STAT1>median(km.dat$STAT1),'High STAT1','Low STAT1'))
group$group<-as.vector(group$group)

km.dat$group<-as.vector(group$group)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
km.dat<-column_to_rownames(km.dat,var = 'sample')
table(km.dat$group)
kmfit<-survfit(Surv(PFI.time, PFI) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High STAT1","Low STAT1"),
                                      legend.title="group",
                                      title="STAT1 KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF","#0073C2FF"))
cluster_survival_median
## MMP9
km.dat<-t(exp_tumor)%>%as.data.frame()
km.dat<-log2(km.dat+1)
group<-data.frame(sample=rownames(km.dat),group=ifelse(km.dat$MMP9>median(km.dat$MMP9),'High MMP9','Low MMP9'))
group$group<-as.vector(group$group)

km.dat$group<-as.vector(group$group)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
km.dat<-column_to_rownames(km.dat,var = 'sample')
table(km.dat$group)
kmfit<-survfit(Surv(PFI.time, PFI) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High MMP9","Low MMP9"),
                                      legend.title="group",
                                      title="MMP9 KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF","#0073C2FF"))
cluster_survival_median
