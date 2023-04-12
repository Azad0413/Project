rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-302")
if (! dir.exists("./15_immune")){
  dir.create("./15_immune")
}
setwd("./15_immune")

### GSE182616 烧伤过程免疫水平变化
library(GEOquery)
library(Biobase)
library(tidyverse)
library(dplyr)
gset<-getGEO("GSE182616",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
pd<-pData(a)
gpl<-getGEO("GPL17077",destdir = '.')
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOLl'!='')%>%
  separate('GENE_SYMBOL',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symbol=probe2symbol[probe2symbol$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
### 提取0、2、4、8、12、24h的表达数据；
colnames(pd)
clinical<-data.frame(Sample=pd$geo_accession,
                     time=pd$characteristics_ch1#,
                     #  statu=pd$`mortality:ch1`,
                    # TBSA=pd$`tbsa:ch1`
                    )
table(clinical$time)
clinical$time<-ifelse(clinical$time=='time point: hr0','0',
                      ifelse(clinical$time=='time point: hr2','2',
                             ifelse(clinical$time=='time point: hr4','4',
                                    ifelse(clinical$time=='time point: hr8','8',
                                           ifelse(clinical$time=='time point: hr12','12',
                                                  ifelse(clinical$time=='time point: hr24','24',NA))))))
table(clinical$time)
clinical$time<-as.numeric(clinical$time)
clinical<-clinical[order(clinical$time),]
clinical<-clinical[c(1:305),]
## 把表达矩阵提取出来
data.final<-dat[,clinical$Sample]
write.table(clinical,file = 'group.time.xls',sep = '\t',quote = F,row.names = F)
write.table(data.final,file = 'data.time.xls',sep = '\t',quote = F,row.names = T)

### ssGSEA---------
library(GSVA)
gene_set <- read.table("/data/nas1/luchunlin/pipeline/ssGSEA/mmc3.txt",
                       header = T,
                       sep ="\t")
dat.final2 <- as.matrix(data.final)
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])

ssgsea_score = gsva(dat.final2, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.table(ssgsea_score,
            file = "ssgsea_result.xls",
            sep = "\t",
            quote = F)

group<-clinical
## 28个免疫浸润细胞画折线图
violin_dat<-ssgsea_score%>%as.data.frame()
violin_dat$cell<-rownames(violin_dat)
violin_dat <- gather(violin_dat, key=sample, value=score, -c("cell"))
table(group$time)
violin_dat$group <- c(rep('0h',1680),rep('2h',1456),rep('4h','1428'),rep('8h',1372),rep('12h',1344),rep('24h',1260))
#library(rstatix)
#stat.test<-violin_dat%>%
#  group_by(cell)%>%
#  wilcox_test(score~group)%>%
#  adjust_pvalue(method = 'fdr')%>%
#  add_significance("p.adj")
violin_dat$group<-factor(violin_dat$group,levels = c('0h','2h','4h','8h','12h','24h'))
head(violin_dat)

violin_plot <- ggplot(violin_dat, aes(x=group,
                                      y=score,
                                      fill=group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
              width=0.1,
              position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
              position=position_dodge(0.9),
              outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄

  scale_fill_manual(values= c("#FF6A6A", "#20B2AA",'#FFCCCC','#77FFCC','#D28EFF','#FA8072'), name = "Group")+
  labs(title="Changes of immune infiltration level at different timepoint", x="", y = "Score",size=20) +
  # stat_compare_means(data = violin_dat,
  #                    mapping = aes(group = group),
  #                    label ="p.signif",
  #                    method = 'wilcox.test',
  #                    paired = F,
  #                    position = 'identity') +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=10), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 8),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill='none')+
  facet_wrap(~cell,scales = "free",nrow = 6) 
  
violin_plot
ggsave('boxplot.pdf',violin_plot,w=11,h=11)
ggsave('boxplot.png',violin_plot,w=12,h=11)
