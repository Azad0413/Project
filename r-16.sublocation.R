rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./16_sublocation")){
  dir.create("./16_sublocation")
}
setwd("./16_sublocation")

# lncLocator
hubgene <- read.delim2('../07_hubgene/hubgene.xls')
hubgene

CTD_2626G11.2.loc <- data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                               score=as.numeric(c('0.818','0.129','0.010','0.035','0.008')))
### 直方图
CTD_2626G11.2<-ggplot(data = CTD_2626G11.2.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'darkblue','lightgreen','purple'))+
  # coord_flip()+
  ggtitle('CTD-2626G11.2 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
CTD_2626G11.2

AP000696.2.loc <- data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                                score=as.numeric(c('0.791','0.105','0.035','0.053','0.016')))
### 直方图
AP000696.2<-ggplot(data = AP000696.2.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'darkblue','lightgreen','purple'))+
  # coord_flip()+
  ggtitle('AP000696.2 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
AP000696.2

RP11_528A4.2.loc <- data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                             score=as.numeric(c('0.665','0.276','0.009','0.026','0.024')))
### 直方图
 RP11_528A4.2<-ggplot(data =  RP11_528A4.2.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'darkblue','lightgreen','purple'))+
  # coord_flip()+
  ggtitle(' RP11-528A4.2 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
RP11_528A4.2


LINC00645.loc <- data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                               score=as.numeric(c('0.874','0.021','0.021','0.082','0.002')))
### 直方图
LINC00645<-ggplot(data =  LINC00645.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'darkblue','lightgreen','purple'))+
  # coord_flip()+
  ggtitle(' LINC00645 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
LINC00645

RP4_655J12.4.loc <- data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                               score=as.numeric(c('0.752','0.157','0.012','0.059','0.019')))
### 直方图
RP4_655J12.4<-ggplot(data = RP4_655J12.4.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'darkblue','lightgreen','purple'))+
  # coord_flip()+
  ggtitle('RP4-655J12.4 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
RP4_655J12.4
RP11_321G12.1.loc <- data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                                score=as.numeric(c('0.781','0.054','0.065','0.060','0.040')))
### 直方图
RP11_321G12.1<-ggplot(data = RP11_321G12.1.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'darkblue','lightgreen','purple'))+
  # coord_flip()+
  ggtitle('RP11-321G12.1 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
RP11_321G12.1


RP11_195B3.1.loc <- data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                            score=as.numeric(c('0.244','0.180','0.124','0.420','0.032')))
### 直方图
RP11_195B3.1<-ggplot(data =  RP11_195B3.1.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'darkblue','lightgreen','purple'))+
  # coord_flip()+
  ggtitle(' RP11-195B3.1 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
RP11_195B3.1

PTCSC3.loc <- data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                            score=as.numeric(c('0.758','0.052','0.025','0.051','0.114')))
### 直方图
PTCSC3<-ggplot(data =  PTCSC3.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'darkblue','lightgreen','purple'))+
  # coord_flip()+
  ggtitle(' PTCSC3 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
PTCSC3

library(patchwork)
all <- CTD_2626G11.2+AP000696.2+RP11_528A4.2+LINC00645+RP4_655J12.4+RP11_321G12.1+RP11_195B3.1+PTCSC3+
  plot_layout(ncol = 3)&
  theme_bw()&
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold",size = 8),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=7), )
all
ggsave('01.location.pdf',all,w=8,h=7)
ggsave('01.location.png',all,w=8,h=7)
