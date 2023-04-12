rm(list = ls())
# 02 差异分析---------
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./04_DEARG")){
  dir.create("./04_DEARG")
}
setwd("./04_DEARG")
##差异基因------
diff1 <- read.delim2('../01_DEGs(GSE16088)/DEG_sig(GSE16088).xls')
up1 <- data.frame(symbol=rownames(diff1[which(diff1$change=='UP'),]))
down1 <- data.frame(symbol=rownames(diff1[which(diff1$change=='DOWN'),]))

diff2 <- read.delim2('../02_DEGs(GSE99671)/DEG_sig(GSE99671).xls',row.names = 1)
up2 <- data.frame(symbol=rownames(diff2[which(diff2$change=='UP'),]))
down2 <- data.frame(symbol=rownames(diff2[which(diff2$change=='DOWN'),]))

diff3 <- read.delim2('../03_DEGs(GSE19276)/DEG_sig(GSE19276).xls')
up3 <- data.frame(symbol=rownames(diff3[which(diff3$change=='UP'),]))
down3 <- data.frame(symbol=rownames(diff3[which(diff3$change=='DOWN'),]))

up <- data.frame(symbol=intersect(up1$symbol,up2$symbol))
up <- data.frame(symobl=intersect(up$symbol,up3$symbol))
##11
down <- data.frame(symbol=intersect(down1$symbol,down2$symbol))
down <- data.frame(symobl=intersect(down$symbol,down3$symbol))
##57
diff <- rbind(up,down)
##68
write.table(diff,file = 'DEGs.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)

mydata1<-list('UP(GSE16088)'=up1$symbol,'UP(GSE99671)'=up2$symbol,'UP(GSE19276)'=up3$symbol)
pdf(file = '01.up.pdf',w=6,h=6)
ggvenn(mydata1,c('UP(GSE16088)','UP(GSE99671)','UP(GSE19276)'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
png(file = '01.up.png',w=400,h=400)
ggvenn(mydata1,c('UP(GSE16088)','UP(GSE99671)','UP(GSE19276)'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()

mydata2<-list('DOWN(GSE16088)'=down1$symbol,'DOWN(GSE99671)'=down2$symbol,'DOWN(GSE19276)'=down3$symbol)
pdf(file = '02.down.pdf',w=6,h=6)
ggvenn(mydata2,c('DOWN(GSE16088)','DOWN(GSE99671)','DOWN(GSE19276)'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
png(file = '02.down.png',w=400,h=400)
ggvenn(mydata2,c('DOWN(GSE16088)','DOWN(GSE99671)','DOWN(GSE19276)'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()


## anoikis
geneset <- read.csv('Anoikis.csv')
colnames(geneset)
geneset <- data.frame(symbol=geneset[geneset$Relevance.score>0.4,]$Gene.Symbol)
##501

DEARG <- data.frame(symbol=intersect(geneset$symbol,diff$symobl))
##7
write.table(DEARG,file = 'DEARG.xls',sep = '\t',row.names = F,quote = F)

mydata3<-list('DEGs'=diff$symobl,'Anoikis'=geneset$symbol)
pdf(file = '03.DEARG.pdf',w=6,h=6)
ggvenn(mydata3,c('DEGs','Anoikis'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png(file = '03.DEARG.png',w=400,h=400)
ggvenn(mydata3,c('DEGs','Anoikis'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()

## expression
### 1------
dat<-read.delim2("../00_rawdata/dat(GSE16088).xls", row.names = 1)%>% lc.tableToNum
group = read.delim2("../00_rawdata/group(GSE16088).xls")
table(group$group)
control.sample<-group$sample[which(group$group=='control')]
df<-diff1
## 
hub_exp<-dat[DEARG$symbol,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','OS')

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  t_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
df<-df[stat.test$Symbol,]
stat.test$p.adj<-df$adj.P.Val%>%as.numeric()%>%round(digits = 3)
stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
                        ifelse(stat.test$p.adj<0.05,"**",
                               ifelse(stat.test$p.adj<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700"), name = "Group")+
  labs(title="Expression(GSE16088)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(12,12,9,11,8,8,12.5),
                     size = 3.2,
                     family = "Times",
                     label = "p.adj",
                     #parse = T,
                     face = "bold")+
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
        panel.grid.minor = element_blank())+
  guides(fill='none')+
  facet_wrap(~Symbol,scales = "free",nrow = 1) 
exp_plot
ggsave('04.expression(GSE16088).pdf',exp_plot,width = 12,height = 3.5)
ggsave('04.expression(GSE16088).png',exp_plot,width = 12,height = 3.5)

### 2------
dat<-read.delim2("../00_rawdata/dat(GSE99671).xls", row.names = 1)%>% lc.tableToNum
library(edgeR)
dat <- cpm(dat)
dat <- log2(dat+1)%>%as.data.frame()
group = read.delim2("../00_rawdata/group(GSE99671).xls")
table(group$group)
control.sample<-group$sample[which(group$group=='control')]
df<-diff2
## 
hub_exp<-dat[DEARG$symbol,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','OS')

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  t_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
df<-df[stat.test$Symbol,]
stat.test$p.adj<-df$padj%>%as.numeric()%>%round(digits = 3)
stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
                        ifelse(stat.test$p.adj<0.05,"**",
                               ifelse(stat.test$p.adj<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700"), name = "Group")+
  labs(title="Expression(GSE99671)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(7,7,9,13,13,3,10),
                     size = 3.2,
                     family = "Times",
                     label = "p.adj",
                     #parse = T,
                     face = "bold")+
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
        panel.grid.minor = element_blank())+
  guides(fill='none')+
  facet_wrap(~Symbol,scales = "free",nrow = 1) 
exp_plot
ggsave('05.expression(GSE99671).pdf',exp_plot,width = 12,height = 3.5)
ggsave('05.expression(GSE99671).png',exp_plot,width = 12,height = 3.5)


### 3------
dat<-read.delim2("../00_rawdata/dat(GSE19276).xls", row.names = 1)%>% lc.tableToNum
group = read.delim2("../00_rawdata/group(GSE19276).xls")
table(group$group)
control.sample<-group$sample[which(group$group=='control')]
df<-diff3
## 
hub_exp<-dat[DEARG$symbol,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','OS')

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  t_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
df<-df[stat.test$Symbol,]
stat.test$p.adj<-df$adj.P.Val%>%as.numeric()%>%round(digits = 3)
stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
                        ifelse(stat.test$p.adj<0.05,"**",
                               ifelse(stat.test$p.adj<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700"), name = "Group")+
  labs(title="Expression(GSE19276)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(4,5,3,8,8.5,4,5),
                     size = 3.2,
                     family = "Times",
                     label = "p.adj",
                     #parse = T,
                     face = "bold")+
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
        panel.grid.minor = element_blank())+
  guides(fill='none')+
  facet_wrap(~Symbol,scales = "free",nrow = 1) 
exp_plot
ggsave('06.expression(GSE19276).pdf',exp_plot,width = 12,height = 3.5)
ggsave('06.expression(GSE19276).png',exp_plot,width = 12,height = 3.5)
