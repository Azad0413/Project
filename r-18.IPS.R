rm(list = ls())
setwd("/data/nas1/luchunlin/project/HZ0301-3/")
if (! dir.exists("./18_ips")){
  dir.create("./18_ips")
}
setwd("./18_ips")

ips<-read_delim('TCIA-ClinicalData.tsv')%>%data.frame(.)
ips<-ips[,c('barcode','ips_ctla4_neg_pd1_neg')]
rownames(ips)<-ips$barcode
risk <- read.delim2('../07_risk/risk.xls')

risk$sample<-substr(risk$id,1,12)
risk<-risk[!duplicated(risk$sample),]
risk$risk <- ifelse(risk$risk==0,'High risk','Low risk')
rownames(risk)<-risk$sample
data<-merge(risk,ips,by = 'row.names')
rownames(data)<-data$Row.names
write.csv(data,'ips_data.csv',quote=F) 

library(ggplot2)
library(ggpubr)
table(data$risk)
my_comparisons = list( c('High risk','Low risk'))
ggplot(data, aes(risk, ips_ctla4_neg_pd1_neg)) + 
  ggdist::stat_halfeye(aes(color=risk,fill=risk),adjust = .5, width = .7, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(aes(color=risk),width = .1, outlier.shape = NA) +
  # gghalves::geom_half_point(aes(color=Species),side = "l", range_scale = .4, alpha = .5) +
  #ggsci::scale_color_nejm()+
  #ggsci::scale_fill_nejm() +
  geom_jitter(aes(color=risk),width = .05, alpha = .3) +
  scale_fill_manual(values=c("#DB3587","#35B3B2")) +
  scale_color_manual(values=c("#DB3587","#35B3B2")) + 
  coord_flip()+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",label.y = c(10.5),
                     color='black',family="Times",face = "bold")+
  labs(x = "Risk", y = "IPS", title = "") + 
  theme(axis.title.x = element_text(size = 20, face = "bold", family = "Times",color='black'),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
        axis.text.y = element_text(size = 16, face = "bold", family = "Times",color='black'),
        legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
        legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
        legend.position = "bottom"
  )+
  theme(panel.grid =element_blank()) +   
  theme(axis.text = element_blank()) +  
  theme(axis.ticks = element_blank())+
  theme(panel.background = element_blank())
ggsave("01.ips_boxplot.pdf",w=7,h=6)
ggsave("01.ips_boxplot.png",w=7,h=6)

