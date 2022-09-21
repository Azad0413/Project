rm(list = ls())
#08 HLA_exp-------------
setwd("/data/nas1/luchunlin/project/YCZK-127/")
if (! dir.exists("./10_HLA")){
  dir.create("./10_HLA")
}
setwd("./10_HLA")
## 分析不同risk score组中24个HLA家族基因的表达
dat.tcga<-read.delim2("/data/nas1/luchunlin/project/YCZK-127/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colname<-data.frame(sample=colnames(dat.tcga))
colname$sample<-gsub('.','-',colname$sample,fixed = T)
colnames(dat.tcga)<-colname$sample
risk<-read.delim2('/data/nas1/luchunlin/project/YCZK-127/05_risk/risk.xls')
dat.tcga<-dat.tcga[,risk$id]
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
risk2 <- risk

risk2$risk_label <- ifelse(risk$risk == 0, "High", "Low")

dim(dat.tcga)
## 把HLA基因表达矩阵提取出来
HLA<-read_xlsx(file = 'HLA_family.xlsx')
HLA.dat<-dat.tcga[HLA$Symbol,]
HLA.dat<-log2(HLA.dat+1)
HLA.dat$gene<-rownames(HLA.dat)
violin_dat <- gather(HLA.dat, key=sample, value='expr', -c("gene"))
head(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% high.sample,
                           "High", "Low") 
head(violin_dat)
colnames(violin_dat)
library(rstatix)
stat.test<-violin_dat%>%
  group_by(gene)%>%
  wilcox_test(expr ~ group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'HLA_exp_result.xls',
            sep = '\t',
            row.names = F)
violin_plot <- ggplot(violin_dat, aes(x=gene, 
                                      y=expr,
                                      fill=group)) +
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#A73030FF", "#0073C2FF"), name = "Group")+
  labs(title="Expression of HLA family genes", x="", y = "log2(expr+1)",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  #  geom_signif(comparisons = my_comparisons,
  #              test = t.test,
  #              map_signif_level = T)+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
violin_plot
