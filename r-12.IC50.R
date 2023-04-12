rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-441-3/")
if (! dir.exists("./12_IC50")){
  dir.create("./12_IC50")
}
setwd("./12_IC50")

library(tidyverse)

dat<-read.csv('../00_rawdata/dat.fpkm.xls',sep = '\t',row.names = 1)
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group<-read.delim2('../05_survival/group(UGCG).xls')
high.sample<-group$sample[which(group$group=='High UGCG')]

model_expr <- log2(dat+1)
model_expr<-model_expr[,group$sample]

###GDSC2 venetoclax-----
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='/data/nas1/luchunlin/GDSC/DataFiles/Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

testExpr<- as.matrix(model_expr)
dim(testExpr)  
library(oncoPredict)
class(GDSC2_Res)
GDSC2_Res <- GDSC2_Res[,c('Venetoclax_1909',"Camptothecin_1003")]
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
result2 <- log2(result2+1)
drug_plot2 <- result2[c('Venetoclax_1909')]%>%rownames_to_column(var='sample')
# drug_plot2 <- log2(drug_plot2+1)
colnames(drug_plot2)[2] <- c('Venetoclax')

### 差异分析--------
drug_plot <- merge(group,drug_plot2,by='sample')
drug_plot2 <- tidyr::gather(drug_plot,Drug,IC50,-c('sample','group'))

library(rstatix)
library(ggplot2)
library(ggpubr)
res.drug <- drug_plot2%>%
  wilcox_test(IC50 ~ group) %>%
  adjust_pvalue(method = 'BH')%>%
  add_significance('p')

write.table(res.drug,file = 'res.drug.xls',sep = '\t',row.names = F,quote = F)
# DE.drug<-res.drug[which(res.drug$p<0.05),]
# write.table(DE.drug,file = 'DE.drug.xls',sep = '\t',row.names = F,quote = F)
colnames(drug_plot2)
# violin<-drug_plot2[drug_plot2$Drug%in%res.drug$Drug[which(res.drug$p<0.05)],]
violin<-drug_plot2
drug_plot <- ggplot(violin, aes(x=group,
                                y=IC50,
                                fill=group)) +
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
                     mapping = aes(group = group),
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
ggsave(filename = '01.drug.plot.pdf',drug_plot,w=5,h=5)
ggsave(filename = '01.drug.plot.png',drug_plot,w=5,h=5)
