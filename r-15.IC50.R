rm(list = ls())
#09 checkpoint-------------
setwd("/data/nas1/luchunlin/project/NN-0118-2/")
if (! dir.exists("./15_IC50")){
  dir.create("./15_IC50")
}
setwd("./15_IC50")

library(lance)
library(tidyverse)
dat<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group <- read.delim2('../10_KM/group(CCNB1).xls')
table(group$group)
high.sample<-group$sample[which(group$group=='High CCNB1')]
low.sample<-group$sample[which(group$group=='Low CCNB1')]
dat<-dat[,group$sample]
library(pRRophetic)
library(ggplot2)
set.seed(12345)
hubgene <- read.delim2('../06_PPI/hubgene.xls')
model_expr<-dat[hubgene$symbol,]
drug<-read.table(file = '/data/nas1/luchunlin/pipeline/Medicinal_Sensity/drugs.txt',sep='\t',header=F)
ic50<-data.frame(group$sample)
a<-data.frame(row.names=group$sample,group=group$group)
colnames(a)<-'group'

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
medicinal_result <- t(subset(x, select = -group)) 
high_group <- high.sample
low_group <- low.sample
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
ggsave(p,filename = '01.allplot.pdf',w=10,h=6)
ggsave(p,filename = '01.allplot.png',w=10,h=6)


### 差异最显著的5个-----
top5 <- c('Shikonin','MG.132','RDEA119','AZ628','FH535')
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
all(rownames(rTable) == rownames(medicinal_result))
drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
drugs_res <- drugs_res[top5,]
violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
table(group$group)
violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                            "High CCNB1", "Low CCNB1") 
violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High CCNB1", "Low CCNB1"))
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
ggsave(filename = "02.drugs.plot.pdf", height = 6, width = 7)
ggsave(filename = "02.drugs.plot.png", height = 6, width = 7)

## 基因与药物的相关性--------
hubgene <- 'CCNB1'
hub_exp <- dat[hubgene,]%>%lc.tableToNum()
DEdrug <- drugs_res[,-c(1,1008)]
DEdrug <- t(DEdrug)%>%as.data.frame()
rownames(DEdrug) <- gsub('.','-',rownames(DEdrug),fixed = T)
cnt<-1
while(cnt<2){
  library(psych)
  x <-as.numeric(hub_exp[cnt,])
  y <-DEdrug
  library(psych)
  d <- corr.test(y,x,use="complete",method = 'spearman')
  r <- data.frame(d$r)
  p <- data.frame(d$p)
  write.table(p,file=paste(rownames(hub_exp)[cnt], "p.xls",sep='_'),sep="\t",quote=F)
  write.table(r,file=paste(rownames(hub_exp)[cnt], "r.xls",sep='_'),sep="\t",quote=F)
  correlation<-data.frame(rownames(p),r$d.r,p$d.p)
  colnames(correlation) <- c("cell","Correlation","p.value")
  correlation<-correlation[correlation$p.value<0.05,]
  correlation$sig[correlation$p.value>0.05] <- "ns"   
  correlation$sig[correlation$p.value<0.05&correlation$p.value>0.01] <- "*"   
  correlation$sig[correlation$p.value<0.01&correlation$p.value>0.001] <- "**"  
  correlation$sig[correlation$p.value<0.001&correlation$p.value>0.0001] <- "***"   
  correlation$sig[correlation$p.value<0.0001] <- "****"
  correlation$'correlation'<-abs(correlation$Correlation)
  #correlation_sig<-correlation[c(3,9,10,11,23,29,7,12,18,31,1,14,15,25,28),]
  library("viridis")
  #trace(ggdotchart,edit=T)
  p0<-ggdotchart(correlation, x = "cell", y = "Correlation",
                 
                 dot.size ='correlation',
                 # filled.contour(p.value,col=Lab.palette(20)),
                 color ='p.value',
                 # label='sig',
                 
                 #  font.label = list( size = 9),
                 
                 #ylab="mgp = c(3, 2, 0)",
                 # font.label = list(size=10,face="bold",color='black',position="dodge",
                 #                   vjust=0.5),
                 #  y.label=0.7,
                 # palette = palette,
                 # 按照cyl填充颜色
                 # palette = palette, # 修改颜色
                 sorting = "descending",
                 add = "segments",                             # 添加棒子
                 rotate = TRUE,
                 ggtheme = theme_pubr(),                        # 改变主题
                 ylab="Correlation Coefficient (r)",
                 xlab='',
                 # dot.size = 6 ,
                 title=rownames(hub_exp)[cnt]
                 
  )+scale_colour_gradient( high = "#4682B4",low = "#66CDAA")
  
  p10<-p0+theme(legend.position = "right",
                panel.background = element_blank())+geom_hline(aes(yintercept=0),linetype="dashed",lwd = 0.2)+
    theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
          axis.text.x =element_text(size=12,family = "Times", face = "bold"),
          axis.title.y =element_text(size=20,family = "Times", face = "bold"),
          axis.text.y=element_text(size=12,family = "Times", face = "bold"),
          plot.title=element_text(size=15,family = "Times", face = "bold",hjust=0.5),
          legend.text = element_text(size = 13, family = "Times"),
          legend.title = element_text(size = 16, family = "Times",face = "bold"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  p10
  ggsave(paste(cnt+2,rownames(hub_exp)[cnt],'pdf',sep='.'),w=6,h=4)
  ggsave(paste(cnt+2,rownames(hub_exp)[cnt],'png',sep='.'),w=6,h=4)
  
  cnt<-cnt+1
}

