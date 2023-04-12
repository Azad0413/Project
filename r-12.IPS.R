rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-255-2/")
if (! dir.exists("./12_ips")){
  dir.create("./12_ips")
}
setwd("./12_ips")

risk <- read.delim2('../06_risk/risk.xls')%>%dplyr::select(c('id','risk'))
risk$barcode <- risk$id
risk$barcode <- substr(risk$barcode,1,12)
ips<-read.delim2('TCIA-ClinicalData (1).tsv')
ips<-ips[ips$barcode%in%risk$barcode,]
colnames(ips)
ips <- ips[,c(1,23:26)]

ips <- ips[,c(1,2,4,3,5)]
colnames(ips) <- c('barcode','IPS','IPS-CTLA4 blocker','IPS-PD1/PDL1/PDL2 blocker','IPS-CTLA4-and PD1/PDL1/PDL2 blocker')

IPS<-merge(ips,risk,by='barcode')%>%column_to_rownames(var = 'id')

IPS <- na.omit(IPS)
IPS$risk <- ifelse(IPS$risk==0,'High risk','Low risk')
High.sample <- rownames(IPS[which(IPS$risk=='High risk'),])

IPS <- IPS[,-c(1,6)]
ips.dat <- t(IPS)%>%as.data.frame()
ips.dat$ips <- rownames(ips.dat)

violin_dat <- gather(ips.dat, key=sample, value='score', -c("ips"))
head(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% High.sample,
                           "High risk", "Low risk") 
head(violin_dat)
colnames(violin_dat)
library(rstatix)
stat.test<-violin_dat%>%
  group_by(ips)%>%
  wilcox_test(score ~ group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox.result.xls',sep = '\t',row.names = F,quote = F)
violin_plot <- ggplot(violin_dat, aes(x=group, 
                                      y=score,
                                      fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.2,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#DB3587","#35B3B2"), name = "Group")+
  labs(title="IPS results", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.45) +
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
        panel.grid.minor = element_blank())+facet_wrap(~ips,scales = "free",nrow = 3) +
  guides(fill='none')
violin_plot
ggsave('01.IPS.pdf',violin_plot,w=7,h=6)
ggsave('01.IPS.png',violin_plot,w=7,h=6)
