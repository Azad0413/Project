rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/HF-0103-1/")
if (! dir.exists("./18_IC50")){
  dir.create("./18_IC50")
}
setwd("./18_IC50")
library(lance)
library(tidyverse)
dat<-read.delim2('../00_rawdata/dat(GSE10846).xls',row.names = 1)%>%lc.tableToNum()
# colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
risk<-read.delim2('../07_risk/risk.xls')
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
dat<-dat[,risk$id]
library(pRRophetic)
library(ggplot2)
set.seed(12345)
model_geneids<-read.delim2('../06_Lasso/lasso_genes.csv',header = F)
model_expr<-dat[model_geneids$V1, risk$id]
#model_expr <- log2(model_expr+1)
riskscore<-data.frame(risk$id,risk$risk)
colnames(riskscore)<-c('sample','risk')
riskscore$risk[which(riskscore$risk==1)] <-'Low risk'
riskscore$risk[which(riskscore$risk==0)] <-'High risk'
head(riskscore)
library(oncoPredict)

library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(reshape2)
library(ggpubr)


th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='./DataFiles/Training Data/'
GDSC1_Expr = readRDS(file=file.path(dir,'GDSC1_Expr (RMA Normalized and Log Transformed).rds'))
GDSC1_Res = readRDS(file = file.path(dir,"GDSC1_Res.rds"))
GDSC1_Res <- exp(GDSC1_Res) 

testExpr<- as.matrix(model_expr)
dim(testExpr)  

calcPhenotype(trainingExprData = GDSC1_Expr,
              trainingPtype = GDSC1_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 5, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
result1 <- read.csv('calcPhenotype_Output/DrugPredictions.csv',row.names = 1)

drug <- data.frame(drug=colnames(result1))
drug_plot1 <- result1[c('LFM.A13_192','QL.XII.47_235','QL.XII.61_1203','QL.X.138_331','Lenalidomide_1020')]
drug_plot1 <- log2(drug_plot1+1)
colnames(drug_plot1) <- c('LFM-A13','QL-XII-47','QL-XII-61','QL-X-138','Lenalidomide')
#LFM.A13_192	QL.XII.47_235 QL.XII.61_1203  QL.X.138_331		Lenalidomide_1020



###GDSC2 Ibrutinib_1799------
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='./DataFiles/Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

testExpr<- as.matrix(model_expr)
dim(testExpr)  

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 5, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
result2 <- read.csv('calcPhenotype_Output/DrugPredictions.csv',row.names = 1)

drug <- data.frame(drug=colnames(result2))
drug_plot2 <- result2[c('Ibrutinib_1799')]
drug_plot2 <- log2(drug_plot2+1)
colnames(drug_plot2) <- c('Ibrutinib')

### 差异分析--------
drug_plot <- cbind(drug_plot1,drug_plot2)%>%as.data.frame()%>%
  rownames_to_column(var = 'sample')

drug_plot <- merge(riskscore,drug_plot,by='sample')
drug_plot2 <- tidyr::gather(drug_plot,Drug,IC50,-c('sample','risk'))

library(rstatix)
library(ggplot2)
library(ggpubr)
res.drug <- drug_plot2%>%
  group_by(Drug)%>%
  wilcox_test(IC50 ~ risk) %>%
  adjust_pvalue(method = 'BH')%>%
  add_significance('p')

write.table(res.drug,file = 'res.drug.xls',sep = '\t',row.names = F,quote = F)
DE.drug<-res.drug[which(res.drug$p<0.05),]
# write.table(DE.drug,file = 'DE.drug.xls',sep = '\t',row.names = F,quote = F)
colnames(drug_plot2)
violin<-drug_plot2[drug_plot2$Drug%in%res.drug$Drug[which(res.drug$p<0.05)],]
drug_plot <- ggplot(violin, aes(x=risk,
                                y=IC50,
                                fill=risk)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="", x="", y = "IC50",size=20) +
  stat_compare_means(data = violin,
                     mapping = aes(group = risk),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.4) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=0,hjust=,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+facet_wrap(~Drug,scales = "free",nrow = 2) +
  guides(fill='none')
drug_plot
ggsave(filename = '01.drug.plot.pdf',drug_plot,w=6,h=5)
ggsave(filename = '01.drug.plot.png',drug_plot,w=6,h=5)

