rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-327/")
if (! dir.exists("./07_CIBERSORT")){
  dir.create("./07_CIBERSORT")
}
setwd("./07_CIBERSORT")
# AIH.PBC------
library(tidyverse)
library(lance)
dat1<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/00_rawdata/AIH.PBC.fpkm.xls',row.names = 1)%>%lc.tableToNum()
dat1<-na.omit(dat1)
group1<-read.delim2("/data/nas1/luchunlin/project/BJTC-327/01_DEGs/group1.xls")
table(group1$group)
AIH.sample<-group1$sample[which(group1$group=='AIH')]
PBC.sample<-group1$sample[which(group1$group=='PBC')]
library(readxl)
LM22gene<-read_xlsx('/data/nas1/luchunlin/pipeline/CIBERSORT/LM22gene.xlsx')%>%as.data.frame()
## 首先，所有的数据要求非负，无缺失值，未经过log转化
## 1、Affymetrix microarrays芯片数据要求 MAS5 or RMA标准化
## 2、Illumina Beadchip 需要经limma package处理
## 3、RNA-seq数据可以使用FPKM和TPM

CIBERSORT_exp<-dat1[rownames(dat1)%in%LM22gene$Gene,]

setwd('/data/nas1/luchunlin/pipeline/CIBERSORT')

write.table(CIBERSORT_exp,
            file = "CIBERSORT_exp.txt",
            quote = F,
            sep = "\t",
            row.names = T)
{ 
  source("Cibersort.R")
  
  result <- CIBERSORT('/data/nas1/luchunlin/pipeline/CIBERSORT/LM22.txt',
                      'CIBERSORT_exp.txt', 
                      perm = 1000, ##Permutations for significance analysis是用来计算单个样本估算免疫浸润的p值，大多数文章会采用1000次。数值越大，运行时间越久，
                      QN = F)
  cibersort_raw <- read.table("CIBERSORT-Results.txt",
                              header = T,
                              sep = "\t",
                              row.names = 1,
                              check.names = F)
  cibersort_result <- t(cibersort_raw[,-c(23,24,25)])
}
setwd("/data/nas1/luchunlin/project/BJTC-327/07_CIBERSORT/")
{
  tiics_result <- cibersort_result
  pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
  for (i in 1:nrow(tiics_result)){
    pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, PBC.sample],
                                         tiics_result[i, AIH.sample])$p.value
    log2FoldChange[i, 1] = mean(tiics_result[i, PBC.sample]) - 
      mean(tiics_result[i, AIH.sample])
  }
  padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
  rTable <- data.frame(log2FoldChange, 
                       pvalue, 
                       padj,
                       row.names = rownames(tiics_result))
  PBC <- signif(apply(tiics_result[rownames(rTable), PBC.sample], 
                       1,
                       mean), 4)
  AIH <- signif(apply(tiics_result[rownames(rTable), AIH.sample], 
                      1, 
                      mean), 4)
  rTable <- data.frame(PBC, 
                       AIH,
                       rTable[, c("padj", "pvalue", "log2FoldChange")])
  rTable$immune_cell <- rownames(rTable)
  rTable$sig <- ifelse(rTable$padj < 0.05,
                       ifelse(rTable$padj < 0.01, 
                              ifelse(rTable$padj < 0.001,
                                     ifelse(rTable$padj < 0.0001,
                                            paste(rTable$immune_cell, "****",  sep = ""),
                                            paste(rTable$immune_cell, "***", sep = "")),
                                     paste(rTable$immune_cell, "**", sep = "")),
                              paste(rTable$immune_cell, "*",  sep = "")), 
                       rTable$immune_cell)
  
  write.table(rTable,
              file = "cibersort_tiics_wilcox_test(PBC).xls",
              quote = F,
              row.names = F,
              sep = '\t')
}
write.table(tiics_result,file = 'tiics_result(PBC).xls',sep = '\t',row.names = T,quote = F)
#diff_cibersort_Table<-rTable[which(rTable$pvalue<0.05),]
##13
# write.table(diff_cibersort_Table,
#             file = 'diff_cibersort_Table.xls',
#             quote = F,
#             sep = '\t')
## 13-1 柱状堆叠图--------
cibersort_result<-as.data.frame(cibersort_result)
cibersort_result$cell<-rownames(cibersort_result)
box_dat <- gather(cibersort_result, key=sample, value='score', -c("cell"))
box_dat$group<-ifelse(box_dat$sample%in%PBC.sample,'PBC','AIH')

library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(7,"Paired"))

box_plot <- ggplot(box_dat, aes(x=sample, 
                                y=100*score,
                                fill=cell)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='Relative Percent',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(22))+
  facet_grid(~box_dat$group,scales= "free",space= "free")
box_plot
ggsave('cibersort(AIH vs.PBC).pdf',box_plot,width = 10,height = 7)
ggsave('cibersort(AIH vs.PBC).png',box_plot,width = 10,height = 7)
violin_dat<-cibersort_result%>%as.data.frame()
violin_dat$cell<-rownames(violin_dat)
violin_dat <- gather(violin_dat, key=sample, value=score, -c("cell"))
violin_dat$group <- ifelse(violin_dat$sample%in%PBC.sample,'PBC','AIH')

head(violin_dat)

violin_plot <- ggplot(violin_dat, aes(x=cell,
                                      y=score,
                                      fill=group)) +
 # geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
              width=0.1,
              position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
              position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="Immune Cell", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())#+facet_wrap(~cell,scales = "free",nrow = 3) 
violin_plot

ggsave("boxplot(AIH vs.PBC).pdf",width = 10,height = 6)
ggsave("boxplot(AIH vs.PBC).png",width = 10,height = 6)


## AIH HBV-------
dat2<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/00_rawdata/AIH.HBV.fpkm.xls',row.names = 1)%>%lc.tableToNum()
dat2<-na.omit(dat2)
group2<-read.delim2("/data/nas1/luchunlin/project/BJTC-327/01_DEGs/group2.xls")
table(group2$group)
AIH.sample<-group2$sample[which(group2$group=='AIH')]
HBV.sample<-group2$sample[which(group2$group=='HBV')]
library(readxl)
LM22gene<-read_xlsx('/data/nas1/luchunlin/pipeline/CIBERSORT/LM22gene.xlsx')%>%as.data.frame()
## 首先，所有的数据要求非负，无缺失值，未经过log转化
## 1、Affymetrix microarrays芯片数据要求 MAS5 or RMA标准化
## 2、Illumina Beadchip 需要经limma package处理
## 3、RNA-seq数据可以使用FPKM和TPM
CIBERSORT_exp<-dat2[rownames(dat2)%in%LM22gene$Gene,]
setwd('/data/nas1/luchunlin/pipeline/CIBERSORT')

write.table(CIBERSORT_exp,
            file = "CIBERSORT_exp.txt",
            quote = F,
            sep = "\t",
            row.names = T)
{ 
  source("Cibersort.R")
  
  result <- CIBERSORT('/data/nas1/luchunlin/pipeline/CIBERSORT/LM22.txt',
                      'CIBERSORT_exp.txt', 
                      perm = 1000, ##Permutations for significance analysis是用来计算单个样本估算免疫浸润的p值，大多数文章会采用1000次。数值越大，运行时间越久，
                      QN = F)
  cibersort_raw <- read.table("CIBERSORT-Results.txt",
                              header = T,
                              sep = "\t",
                              row.names = 1,
                              check.names = F)
  cibersort_result <- t(cibersort_raw[,-c(23,24,25)])
}
setwd("/data/nas1/luchunlin/project/BJTC-327/07_CIBERSORT/")
{
  tiics_result <- cibersort_result
  pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
  for (i in 1:nrow(tiics_result)){
    pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, HBV.sample],
                                         tiics_result[i, AIH.sample])$p.value
    log2FoldChange[i, 1] = mean(tiics_result[i, HBV.sample]) - 
      mean(tiics_result[i, AIH.sample])
  }
  padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
  rTable <- data.frame(log2FoldChange, 
                       pvalue, 
                       padj,
                       row.names = rownames(tiics_result))
  HBV <- signif(apply(tiics_result[rownames(rTable), HBV.sample], 
                      1,
                      mean), 4)
  AIH <- signif(apply(tiics_result[rownames(rTable), AIH.sample], 
                      1, 
                      mean), 4)
  rTable <- data.frame(HBV, 
                       AIH,
                       rTable[, c("padj", "pvalue", "log2FoldChange")])
  rTable$immune_cell <- rownames(rTable)
  rTable$sig <- ifelse(rTable$padj < 0.05,
                       ifelse(rTable$padj < 0.01, 
                              ifelse(rTable$padj < 0.001,
                                     ifelse(rTable$padj < 0.0001,
                                            paste(rTable$immune_cell, "****",  sep = ""),
                                            paste(rTable$immune_cell, "***", sep = "")),
                                     paste(rTable$immune_cell, "**", sep = "")),
                              paste(rTable$immune_cell, "*",  sep = "")), 
                       rTable$immune_cell)
  
  write.table(rTable,
              file = "cibersort_tiics_wilcox_test(HBV).xls",
              quote = F,
              row.names = F,
              sep = '\t')
}
write.table(tiics_result,file = 'tiics_result(HBV).xls',sep = '\t',row.names = T,quote = F)
#diff_cibersort_Table<-rTable[which(rTable$pvalue<0.05),]
##13
# write.table(diff_cibersort_Table,
#             file = 'diff_cibersort_Table.xls',
#             quote = F,
#             sep = '\t')
## 13-1 柱状堆叠图--------
cibersort_result<-as.data.frame(cibersort_result)
cibersort_result$cell<-rownames(cibersort_result)
box_dat <- gather(cibersort_result, key=sample, value='score', -c("cell"))
box_dat$group<-ifelse(box_dat$sample%in%HBV.sample,'HBV','AIH')

library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(7,"Paired"))

box_plot <- ggplot(box_dat, aes(x=sample, 
                                y=100*score,
                                fill=cell)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='Relative Percent',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(22))+
  facet_grid(~box_dat$group,scales= "free",space= "free")
box_plot
ggsave('cibersort(AIH vs.HBV).pdf',box_plot,width = 10,height = 7)
ggsave('cibersort(AIH vs.HBV).png',box_plot,width = 10,height = 7)

violin_dat<-cibersort_result%>%as.data.frame()
violin_dat$cell<-rownames(violin_dat)
violin_dat <- gather(violin_dat, key=sample, value=score, -c("cell"))
violin_dat$group <- ifelse(violin_dat$sample%in%HBV.sample,'HBV','AIH')

head(violin_dat)

violin_plot <- ggplot(violin_dat, aes(x=cell,
                                      y=score,
                                      fill=group)) +
  # geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="Immune Cell", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())#+facet_wrap(~cell,scales = "free",nrow = 3) 
violin_plot

ggsave("boxplot(AIH vs.HBV).pdf",width = 10,height = 6)
ggsave("boxplot(AIH vs.HBV).png",width = 10,height = 6)
