rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-441-3/")
if (! dir.exists("./13_GeneMANIA")){
  dir.create("./13_GeneMANIA")
}
setwd("./13_GeneMANIA")


gene <- read.delim2('genemania-genes.txt')
diff <- read.csv('../01_DEGs/DEG_sig.xls',sep = '\t',row.names = 1)
intersect <- data.frame(symbol=intersect(rownames(diff),gene$Symbol)%>%as.data.frame()) 
                        
##4
write.table(intersect,file = 'intersect.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('DEGs'=rownames(diff),'GeneMANIA'=gene$Symbol)
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEGs','GeneMANIA'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png('01.venn.png',w=400,h=400)
ggvenn(mydata,c('DEGs','GeneMANIA'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()


###相关性
library(ggcorrplot)
dat<-read.csv('../00_rawdata/dat(GSE9476).xls',sep = '\t',row.names = 1)
all.gene <- data.frame(symbol=c(intersect$.,'UGCG'))
hub.exp<-dat[all.gene$symbol,]
exp_corr<-round(cor(t(hub.exp),method = 'spearman'),3)
write.table(exp_corr,file = 'gene.correlation.xls',
            sep = '\t',
            row.names = T)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著

hub_p.mat<-round(cor_pmat(t(hub.exp)),3)

pdf('02.corrheatmap.pdf',w=6,h=6)
col1 <- colorRampPalette(c("#4169E1","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF",
                                    "#FDDBC7","#F4A582","#D6604D", "#B2182B","#CC0000","#990000"))
hub_corr_plot<-corrplot(exp_corr,
                       method = "circle",
                       is.corr = T,
                       type = "lower",
                       p.mat = hub_p.mat,
                       insig = "blank",
                       outline = "white",
                       addCoef.col ="black",
                       col = col1(200))
dev.off()
png('02.corrheatmap.png',w=500,h=500)
col1 <- colorRampPalette(c("#4169E1","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF",
                                    "#FDDBC7","#F4A582","#D6604D", "#B2182B","#CC0000","#990000"))
hub_corr_plot<-corrplot(exp_corr,
                       method = "circle",
                       is.corr = T,
                       type = "lower",
                       p.mat = hub_p.mat,
                       insig = "blank",
                       outline = "white",
                       addCoef.col ="black",
                       col = col1(200))
dev.off()

##expression---------
library(tidyverse)

group<-read.delim2('../00_rawdata/group(GSE9476).xls')
control.sample<-group$sample[which(group$group=='control')]
aml.sample<-group$sample[which(group$group=='AML')]
## UGCG-------

hub_exp<-dat[intersect$.,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','AML')

hub_exp2$Group <- factor(hub_exp2$Group,levels = c('control','AML'))
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
  labs(title="", x="", y = "expression level",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test') +
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
        panel.grid.minor = element_blank())+facet_wrap(~Symbol,scales = "free",nrow = 2) 
exp_plot
ggsave(filename = '03.expression.pdf',exp_plot,w=7,h=7)
ggsave(filename = '03.expression.png',exp_plot,w=7,h=7)


##KM
dat<-read.csv('../00_rawdata/dat.fpkm.xls',sep = '\t',row.names = 1)
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
hub_exp<-dat[intersect$.,]

## KM curve
survival<-read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-LAML.survival.tsv')%>%
  column_to_rownames(var = 'sample')

dat.km <- dat[intersect$.,]
dat.km <- log2(dat.km+1)
survival <- survival[colnames(dat.km),]
dat.km <- t(dat.km)%>%as.data.frame()
dat.km <- cbind(dat.km,survival)
dat.km$OS <- as.numeric(dat.km$OS)
dat.km$OS.time <- as.numeric(dat.km$OS.time)
library(survival)
library(survminer)
for (i in c(1:4)) {
  library(survival)
  dat.km$group <- factor(ifelse(dat.km[,i]>median(dat.km[,i]),'High','Low'),levels = c('High','Low'))
  kmfit <- survfit(Surv(OS.time,OS) ~ group, data = dat.km)
  KM <- ggsurvplot(kmfit,
                   pval = TRUE, 
                   conf.int = F,
                   legend.labs=c("High","Low" ),
                   legend.title="group",
                   title=colnames(dat.km[i]),
                   font.main = c(15,"bold"),
                   risk.table = TRUE, 
                   risk.table.col = "strata", 
                   linetype = "strata", 
                   surv.median.line = "hv", 
                   ggtheme = theme_bw(), 
                   palette = c("#A73030FF", "#0073C2FF"))
  print(KM)
  fn1 = paste0(i+3,'.',colnames(dat.km)[i],'.pdf')
  fn2 = paste0(i+3,'.',colnames(dat.km)[i],'.png')
  ggsave(fn1,KM$plot,width = 4,h=3,units = 'in',limitsize = 300)
  ggsave(fn2,KM$plot,width = 4,h=3,units = 'in',limitsize = 300)
  i <- i+1
}
