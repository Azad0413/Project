rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-399-11/")
if (! dir.exists("./16_IPS")){
  dir.create("./16_IPS")
}
setwd("./16_IPS")

risk <- read.delim2('../08_risk/risk.xls')%>%dplyr::select(c('id','riskScore'))
risk$barcode <- risk$id
risk$barcode <- substr(risk$barcode,1,12)
ips<-read.delim2('IPS.tsv')
ips<-ips[ips$barcode%in%risk$barcode,]
colnames(ips)
ips <- ips[,c(1,16:19)]

ips$IPS <- rowMeans(ips[,c(2:5)])

IPS<-merge(ips,risk,by='barcode')

IPS<-IPS[,c(7,6,8)]%>%column_to_rownames(var = 'id')

colnames(IPS)<-c('IPS','Risk')

IPS$riskScore <- as.numeric(IPS$riskScore)

IPS <- na.omit(IPS)
IPS$Risk <- ifelse(IPS$Risk>median(IPS$Risk),'High risk','Low risk')
# Pairwise t-test between groups
stat.test <- IPS %>%
#  group_by(Risk) %>%  ##????
  wilcox_test(IPS ~ Risk) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance("p")


palette <- c("#6495ed","#dc143c")
ggviolin(IPS, x = "Risk", y = "IPS",add = c("jitter", "mean_sd"),
         fill = "Risk",  palette = palette,
         size=0.5,xlab="", ylab="IPS",ggtheme = theme_bw(), legend = "top",font.axis = 2)+ 
  stat_pvalue_manual(stat.test, label = "p.signif",y.position = 13)+
  theme(axis.title.x =element_text(size=14,colour="black"),
        axis.text.x =element_text(size=13,colour="black"),
        axis.title.y =element_text(size=14,colour="black"),
        axis.text.y=element_text(size=13,colour="black"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
ggsave(width=5,height=5,'01.IPS.pdf')
ggsave(width=5,height=5,'01.IPS.png')
