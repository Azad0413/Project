rm(list = ls())
#09 checkpoint-------------
setwd("/data/nas1/luchunlin/project/BJTC-308/")
if (! dir.exists("./14_IC50")){
  dir.create("./14_IC50")
}
setwd("./14_IC50")

dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-308/00_rawdata/dat2.xls", row.names = 1)%>% lc.tableToNum
risk<-read.delim2('/data/nas1/luchunlin/project/BJTC-308/05_risk/risk.xls')
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
dat<-dat[,risk$id]
library(pRRophetic)
library(ggplot2)
set.seed(12345)
lasso_geneids<-data.frame(lasso_geneids=c('ACAA1','ACOT11','DDHD1','B4GALNT1'))
model_expr<-dat[lasso_geneids$lasso_geneids, risk$id]
riskscore<-data.frame(risk$id,risk$risk)
colnames(riskscore)<-c('sample','risk')
riskscore$risk[which(riskscore$risk==1)] <-'Low risk'
riskscore$risk[which(riskscore$risk==0)] <-'High risk'
head(riskscore)
drug<-read.table(file = 'drugs.txt',sep='\t',header=F)
ic50<-data.frame(riskscore$sample)

a<-data.frame(row.names=riskscore$sample,riskscore$risk)

colnames(a)<-'risk'

cnt<-1

while (cnt < 139) {
  
  predictedPtype <- pRRopheticPredict(as.matrix(model_expr), drug[cnt,],selection=1)
  
  Tipifarnib<-data.frame(predictedPtype)
  
  colnames(Tipifarnib)<-drug[cnt,]
  
  a<-cbind(a,Tipifarnib)
  
  cnt = cnt + 1
}

write.table(a,'IC50.xls',sep='\t',quote=F)

b<-a
b[b<0]<-NA
# 先写成函数的形式，方便调用
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
c<-removeColsAllNa(b)
na_flag <- apply(is.na(c), 2, sum)
x <- c[, which(na_flag == 0)]
View(x)
dim(x)
# [1] 60 108
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
rTable$sig <- ifelse(rTable$padj < 0.05,
                     ifelse(rTable$padj < 0.01, 
                            ifelse(rTable$padj < 0.001,
                                   ifelse(rTable$padj < 0.0001,
                                          paste(rTable$drugs, "****",  sep = ""),
                                          paste(rTable$drugs, "***", sep = "")),
                                   paste(rTable$drugs, "**", sep = "")),
                            paste(rTable$drugs, "*",  sep = "")), 
                     rTable$drugs)

write.table(rTable,
            file = "drugs_wilcox_test.xls",
            quote = F,
            row.names = F)

### 发散条形图绘制

#install.packages('ggprism')
library(ggprism)
## 横坐标药物，纵坐标：IC50(H)/IC50(L)-1
dat_plot<-data.frame(drug=rownames(rTable),
                     'IC50(H)/IC50(L)-1'=(rTable$high_group_res/rTable$low_group_res-1),
                     pvalue=rTable$pvalue,
                     padj=rTable$padj)
dat_plot$threshold=factor(ifelse(dat_plot$padj<0.05&dat_plot$pvalue<0.05,'P<0.05 & FDR<0.05',ifelse(dat_plot$pvalue<0.05&dat_plot$padj>0.05,'P<0.05 & FDR>0.05','P>0.05 & FDR>0.05')))
dat_plot<-dat_plot%>%arrange(desc(dat_plot$IC50.H..IC50.L..1))

dat_plot$drug<-factor(dat_plot$drug,levels=dat_plot$drug)

p <- ggplot(data = dat_plot,aes(x = drug,y = IC50.H..IC50.L..1,fill = threshold)) +
  geom_col()+
  scale_fill_manual(values = c('P<0.05 & FDR<0.05'= '#5F9EA0','P>0.05 & FDR>0.05'='#cccccc','P<0.05 & FDR>0.05'='#FFD700')) +
  xlab('') + 
  ylab('Exp(Median IC50(H))/Exp(Median IC50(L))-1') + 
  theme_prism(border = T) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 65,size = 7,
                               hjust = 1,vjust = 1),
    axis.text.y = element_text(size = 13),
    legend.position = c(0.85,0.85),
    legend.text = element_text(size = 8,face = 'bold')
  )
p
