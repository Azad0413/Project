rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/BJTC-356/")
if (! dir.exists("./05_DEMetas")){
  dir.create("./05_DEMetas")
}
setwd("./05_DEMetas")

## 正离子和负离子2类
pos.metabo <- read.csv('pos-metaboAnalystInput.csv',row.names = 1)
neg.metabo <- read.csv('neg-metaboAnalystInput.csv',row.names = 1)
pos.metabo <- pos.metabo[-1,]
neg.metabo <- neg.metabo[-1,]
colnames(pos.metabo)
pos.metabo <- pos.metabo[,-c(23:30)]
neg.metabo <- neg.metabo[,-c(23:30)]
diff.metabo.pos <- read.delim2('diff_pos_quant_identification.txt')
diff.metabo.neg <- read.delim2('diff_neg_quant_identification.txt')

diff.metabo.neg <- diff.metabo.neg[!duplicated(diff.metabo.neg$ID),]
diff.metabo.pos <- diff.metabo.pos[!duplicated(diff.metabo.pos$ID),]

##本项目使用的是差异倍数大于等于 2 或小于等于 0.5，p-value 值小于 0.01作为筛选条件。
diff.metabo.neg$t.test_p.value_BHcorrect <- as.numeric(diff.metabo.neg$t.test_p.value_BHcorrect)
diff.metabo.pos$t.test_p.value_BHcorrect <- as.numeric(diff.metabo.pos$t.test_p.value_BHcorrect)

length(unique(diff.metabo.neg$ID))
length(unique(diff.metabo.pos$ID))
# sig.diff.pos <- subset(diff.metabo.pos,diff.metabo.pos$ratio>=2 & diff.metabo.pos$t.test_p.value_BHcorrect<0.01 |
#                          diff.metabo.pos$ratio<=0.5 & diff.metabo.pos$t.test_p.value_BHcorrect<0.01)
# 
# sig.diff.neg <- subset(diff.metabo.neg,diff.metabo.neg$ratio>=2 & diff.metabo.neg$t.test_p.value_BHcorrect<0.01 |
#                          diff.metabo.neg$ratio<=0.5 & diff.metabo.neg$t.test_p.value_BHcorrect<0.01)
# 
# length(unique(sig.diff.neg$ID))
# length(unique(sig.diff.pos$ID))

# sig.diff.pos <- read.delim2('diff_pos_quant_remove_no_adduct.txt')
# sig.diff.neg <- read.delim2('diff_neg_quant_remove_no_adduct.txt')


diff.neg.dat <- neg.metabo[diff.metabo.neg$ID,]
diff.pos.dat <- pos.metabo[diff.metabo.pos$ID,]

### 热图--------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)

dat_rep<-rbind(head(diff.metabo.neg[order(diff.metabo.neg$ratio,decreasing = T),],267),
                               head(diff.metabo.neg[order(diff.metabo.neg$ratio,decreasing = F),],307))
df.group <- data.frame(sample=colnames(neg.metabo),group=c(rep('Persistent',22),rep('Paroxysmal',22)))
group_rt<-df.group%>%as.data.frame()
group_rt<-group_rt[order(group_rt$group),]
rt<-diff.neg.dat[,group_rt$sample]
group_rt<-data.frame(group_rt$group)
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
heat<-rt[dat_rep$ID,]%>%lc.tableToNum()

x<-heat
#x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(Paroxysmal="#00CED1",Persistent="#F08080"))
pdf('01.neg.heatmap.pdf',w=6,h=6)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = F,
         annotation_names_row = F)
dev.off()
png('01.neg.heatmap.png',w=500,h=500)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = F,
         annotation_names_row = F)
dev.off()


### 热图pos--------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)

dat_rep<-rbind(head(diff.metabo.pos[order(diff.metabo.pos$ratio,decreasing = T),],477),
               head(diff.metabo.pos[order(diff.metabo.pos$ratio,decreasing = F),],555))
df.group <- data.frame(sample=colnames(pos.metabo),group=c(rep('Persistent',22),rep('Paroxysmal',22)))
group_rt<-df.group%>%as.data.frame()
group_rt<-group_rt[order(group_rt$group),]
rt<-diff.pos.dat[,group_rt$sample]
group_rt<-data.frame(group_rt$group)
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
heat<-rt[dat_rep$ID,]%>%lc.tableToNum()

x<-heat
#x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(Paroxysmal="#00CED1",Persistent="#F08080"))
pdf('02.pos.heatmap.pdf',w=6,h=6)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = F,
         annotation_names_row = F)
dev.off()
png('02.pos.heatmap.png',w=500,h=500)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = F,
         annotation_names_row = F)
dev.off()
