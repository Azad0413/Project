rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./01_CIBERSORT")){
  dir.create("./01_CIBERSORT")
}
setwd("./01_CIBERSORT")
library(tidyverse)
library(lance)
dat<-read.delim2("/data/nas1/luchunlin/project/JNZK-214-8/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group<-data.frame(sample=colnames(dat),group=c(rep('control',41),rep('Tumor',448)))
write.table(group,file = 'group.xls',sep = '\t',row.names = F,quote = F)
tumor.sample<-group$sample[which(group$group=='Tumor')]
control.sample<-group$sample[which(group$group=='control')]
library(immunedeconv)
set_cibersort_binary("CIBERSORT.R")
set_cibersort_mat("LM22.txt")
res.cibersort<-deconvolute(as.matrix(dat),method = 'cibersort')
write.table(res.cibersort,'cibersort.txt',sep = '\t',col.names = T,row.names = F,quote = F)
save(file = 'cibersort.Rdata',res.cibersort)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(7,"Paired"))
##画图
pdf('01.cibersort.box.pdf',w=10,h=8)
res.cibersort %>%
  gather(sample, fraction, -cell_type) %>%
  merge(group,by='sample')%>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
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
  facet_grid(~group,scales= "free",space= "free")
dev.off()
png('01.cibersort.box.png',w=700,h=500)
res.cibersort %>%
  gather(sample, fraction, -cell_type) %>%
  merge(group,by='sample')%>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
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
  facet_grid(~group,scales= "free",space= "free")
dev.off()

## 差异-------
dat.cibersort <- res.cibersort %>% 
  tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
dat.cibersort <- merge(group, dat.cibersort, by = "sample")
dat.cibersort<-dat.cibersort[,-23]
dat.cibersort2 <- tidyr::gather(dat.cibersort, ImmuneCell, Score, -c("sample", "group"))
library(rstatix)
library(ggplot2)
library(ggpubr)
colnames(dat.cibersort2)
stat_cibersort <- dat.cibersort2 %>% 
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
write.table(stat_cibersort,file = 'stat.cibersort.xls',sep = '\t',row.names = F,quote = F)
DE.cibersort<-stat_cibersort[which(stat_cibersort$p<0.05),]
write.table(DE.cibersort,file = 'DE.cibersort.xls',sep = '\t',row.names = F,quote = F)
colnames(dat.cibersort2)
violin.cibersort<-dat.cibersort2[dat.cibersort2$ImmuneCell%in%stat_cibersort$ImmuneCell[which(stat_cibersort$p<0.05)],]
cibersort_plot <- ggplot(violin.cibersort, aes(x=ImmuneCell,
                                               y=Score,
                                               fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  #stat_boxplot(geom="errorbar", 
  #             width=0.1,
  #             position = position_dodge(0.9)) +
  #geom_boxplot(width=0.7,
  #             position=position_dodge(0.9),
  #             outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="Immune Cell", x="", y = "Score",size=20) +
  stat_compare_means(data = violin.cibersort,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
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
        panel.grid.minor = element_blank())+facet_wrap(~ImmuneCell,scales = "free",nrow = 3) 
cibersort_plot
ggsave(filename = '02.cibersort.plot.pdf',cibersort_plot,w=10,h=8)
ggsave(filename = '02.cibersort.plot.png',cibersort_plot,w=11,h=8)
