rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-399-11/")
if (! dir.exists("./18_IC50")){
  dir.create("./18_IC50")
}
setwd("./18_IC50")
library(tidyverse)
library(lance)
dat<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
risk<-read.delim2('../08_risk/risk.xls')
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
dat<-log2(dat[,risk$id]+1)
library(pRRophetic)
library(ggplot2)
set.seed(12345)
lasso_geneids<-read.delim2('../07_Lasso/lasso_genes.csv',header = F)
model_expr<-dat[lasso_geneids$V1, risk$id]
riskscore<-data.frame(risk$id,risk$risk)
colnames(riskscore)<-c('sample','risk')
riskscore$risk[which(riskscore$risk==1)] <-'Low risk'
riskscore$risk[which(riskscore$risk==0)] <-'High risk'
head(riskscore)

# 
# install.packages("oncoPredict")
library(oncoPredict)

library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(reshape2)
library(ggpubr)
#Temozolomide_1375(替莫唑胺) Vincristine_1818(长春新碱)Carmustine_1807（卡莫司汀）

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
result <- read.csv('calcPhenotype_Output/DrugPredictions.csv',row.names = 1)

drug_plot <- result[c('Temozolomide_1375','Vincristine_1818','Carmustine_1807')]
drug_plot <- log2(drug_plot+1)
colnames(drug_plot) <- c('Temozolomide','Vincristine','Carmustine')
# Etoposide(依托泊苷) Cisplatin(顺铂)
#drug<-read.table(file = '/data/nas1/luchunlin/pipeline/Medicinal_Sensity/drugs.txt',sep='\t',header=F)
drug <- data.frame(V1=c('Cisplatin','Etoposide'))
ic50<-data.frame(riskscore$sample)
a<-data.frame(row.names=riskscore$sample,riskscore$risk)
colnames(a)<-'risk'
cnt<-1
while (cnt < 3) {
  predictedPtype <- pRRopheticPredict(as.matrix(model_expr), drug[cnt,],selection=1)
  Tipifarnib<-data.frame(predictedPtype)
  colnames(Tipifarnib)<-drug[cnt,]
  a<-cbind(a,Tipifarnib)
  cnt = cnt + 1
}
#write.table(a,'IC50.xls',sep='\t',quote=F)
b<-a
b <- cbind(b,drug_plot)
# b[b<0]<-NA
# 先写成函数的形式，方便调用
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
c<-removeColsAllNa(b)
na_flag <- apply(is.na(c), 2, sum)
x <- c[, which(na_flag == 0)]
dim(x)
medicinal_result <- t(subset(x, select = -risk)) 
high_group <- risk$id[which(risk$risk==0)]
low_group <- risk$id[which(risk$risk==1)]
pvalue = padj = log2FoldChange <- matrix(0, nrow(medicinal_result), 1)
for (i in 1:nrow(medicinal_result)){
  pvalue[i, 1] = p.value = wilcox.test(medicinal_result[i, high_group],
                                       medicinal_result[i, low_group])$p.value
  log2FoldChange[i, 1] = mean(medicinal_result[i, high_group]) - 
    mean(medicinal_result[i, low_group])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(medicinal_result))
high_group_res <- signif(apply(medicinal_result[rownames(rTable), high_group], 
                               1,
                               median), 4)
low_group_res <- signif(apply(medicinal_result[rownames(rTable), low_group], 
                              1, 
                              median), 4)
rTable <- data.frame(high_group_res, 
                     low_group_res,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$drugs <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$drugs, "****",  sep = ""),
                                          paste(rTable$drugs, "***", sep = "")),
                                   paste(rTable$drugs, "**", sep = "")),
                            paste(rTable$drugs, "*",  sep = "")), 
                     rTable$drugs)

write.table(rTable,
            file = "drugs_wilcox_test.xls",
            quote = F,
            row.names = F)
DE.drug<-rTable[which(rTable$pvalue<0.05),]


library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
all(rownames(rTable) == rownames(medicinal_result))
drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
drugs_res <- drugs_res[which(rownames(drugs_res) %in% rownames(DE.drug)),]
drugs_res <- drugs_res[-c(1,2),]
violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                            "High risk", "Low risk") 
violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High risk", "Low risk"))
head(violin_dat)
drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                                color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                add = "jitter",
                                short.panel.labs = T,
                                ggtheme = theme_bw()) +
  stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                           facet.by = "drugs",
                           short.panel.labs = T,
                           panel.labs.background = list(fill = "white"),
                           ncol = 3,
                           scales = "free_y") + xlab("") + ylab("IC(50)") +
  # geom_text(data=data_text,
  #           mapping=aes(x=x,y=y,label=label),nudge_x=0.1,nudge_y=0.1)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 13, face = "bold"))
drugs_hub_boxplot
ggsave(filename = "01.drugs.plot.pdf", height = 3, width = 7)
ggsave(filename = "01.drugs.plot.png", height = 3, width = 7)

