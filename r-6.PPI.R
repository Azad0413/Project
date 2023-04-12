rm(list = ls())
# 02 差异分析---------
setwd("/data/nas1/luchunlin/project/LZZK-519-10/")
if (! dir.exists("./06_PPI")){
  dir.create("./06_PPI")
}
setwd("./06_PPI")

# library(readxl)
# ## top degree 直方图
# top10degree<-read_xlsx('topdegree.xlsx')
# #colnames(top10degree)<-c('Symbol','Degree')
# top10degree<-as.data.frame(top10degree)
# top10degree<-top10degree[order(top10degree$Degree),]
# # 
# class(top10degree$Degree)
# top10degree$Degree <- factor(top10degree$Degree)
# top10degree$Symbol
# top10degree$Symbol <- factor(top10degree$Symbol,levels = c('CD4','GAD2','GOT1','ATP5A1','RAC2','ATP5B','SDHB','BTK','FYB'))
# library(RColorBrewer)
# # top10degree$group<-c(rep('low',7),rep('mid',1),rep('high',1))
# degree<-ggplot(data = top10degree,aes(x=Symbol,y=Degree,fill=Degree))+
#   geom_bar(stat = 'identity')+
#   scale_fill_manual(values= brewer.pal(10,'Set3'))+
#   coord_flip()+
#   ggtitle('Top10 Degree')+
#   guides(fill='none')+
#   theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
#         axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
#   labs(x='.')
# 
# degree
# ggsave(filename = '02.TOPdegree.pdf',w=5,h=5)
# ggsave(filename = '02.TOPdegree.png',w=5,h=5)
hubgene <- data.frame(Symbol=c('SDHB','SUCLA2','SLC25A3','ATP5B','ATP5A1'))
write.table(hubgene,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)
