rm(list = ls())
#08 风险评分与临床指标相关性分析-------------
setwd("/data/nas1/luchunlin/project/BJTC-258/")
if (! dir.exists("./13_Medicinal_Sensity")){
  dir.create("./13_Medicinal_Sensity")
}
setwd("./13_Medicinal_Sensity")
dat.tcga<-read.delim2("/data/nas1/luchunlin/project/BJTC-258/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colname<-data.frame(sample=colnames(dat.tcga))
colname$sample<-gsub('.','-',colname$sample,fixed = T)
colnames(dat.tcga)<-colname$sample
group<-read.delim2('/data/nas1/luchunlin/project/BJTC-258/05_survival/clinical.xls')%>%
  dplyr::select(c('sample','group'))
uncoxgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-258/04_PCA/univariate_cox_result_0.05.xls')

library(pRRophetic)
library(ggplot2)
set.seed(12345)
model_expr<-dat.tcga[rownames(uncoxgene),group$sample]
agingscore<-group
colnames(agingscore)<-c('sample','group')
head(agingscore)
drug<-read.table(file = 'drugs.txt',sep='\t',header=F)
ic50<-data.frame(agingscore$sample)
a<-data.frame(row.names=agingscore$sample,agingscore$group)
colnames(a)<-'score'
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
# b[b<0]<-NA
# 先写成函数的形式，方便调用
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
c<-removeColsAllNa(b)
na_flag <- apply(is.na(c), 2, sum)
x <- c[, which(na_flag == 0)]
                                                                                                                                                                                                                                                                                                                                                                                                       dim(x)
medicinal_result <- t(subset(x, select = -score)) 
high_group <- group$sample[which(group$group=='High')]
low_group <- group$sample[which(group$group=='Low')]
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
ggsave(filename = '01.all_medicinal.pdf',w=15,h=8)
ggsave(filename = '01.all_medicinal.png',w=15,h=8)
rTable2 <- rTable[which(rTable$high_group_res > rTable$low_group_res &
                          rTable$padj < 0.05),]
dim(rTable2)
# [1]  51  7
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
all(rownames(rTable) == rownames(medicinal_result))
drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, padj=rTable$padj)
# drugs_res <- drugs_res[rownames(rTable2),]
drugs_res <- drugs_res[which(rownames(drugs_res) %in% 
                               c("Paclitaxel", "Docetaxel",
                                 "Cisplatin","Gemcitabine","Mitomycin.C","Doxorubicin")),]
violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","padj"))
violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                            "High score", "Low score") 
violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High score", "Low score"))
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
                           scales = "free_y") + xlab("") + ylab("Estimated ln(IC50)") +
  # geom_text(data=data_text,
  #           mapping=aes(x=x,y=y,label=label),nudge_x=0.1,nudge_y=0.1)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 13, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 13, face = "bold"))
drugs_hub_boxplot
ggsave(filename = "02.drugs.plot.pdf", height = 6, width = 10)
ggsave(filename = "02.drugs.plot.png", height = 6, width = 10)
