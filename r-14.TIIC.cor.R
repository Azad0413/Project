rm(list = ls())
#part3 与差异免疫细胞的相关性--------
## spearman 生物标志物与差异免疫浸润细胞的相关性。
setwd("/data/nas1/luchunlin/project/HF-0106-2/")
if (! dir.exists("./14_tiic_cor")){
  dir.create("./14_tiic_cor")
}
setwd("./14_tiic_cor")
res.cibersort <- read.delim2('../13_CIBERSORT/cibersort.txt',row.names = 1)%>%lc.tableToNum()
dat <- read.delim2('../00_rawdata/dat(GSE169568).xls',row.names = 1)%>%lc.tableToNum()
hubgene <- read.delim2('../07_features/features.xls')
hub_exp <- dat[hubgene$symbol,]

DEtiic <- res.cibersort

DEtiic <- t(DEtiic)%>%as.data.frame()
cnt<-1
while(cnt<5){
  library(psych)
  x <-as.numeric(hub_exp[cnt,])
  y <-DEtiic
  library(psych)
  d <- corr.test(y,x,use="complete",method = 'pearson')
  r <- data.frame(d$r)
  p <- data.frame(d$p)
  write.table(p,file=paste(rownames(hub_exp)[cnt], "p.xls",sep='_'),sep="\t",quote=F)
  write.table(r,file=paste(rownames(hub_exp)[cnt], "r.xls",sep='_'),sep="\t",quote=F)
  correlation<-data.frame(rownames(p),r$d.r,p$d.p)
  colnames(correlation) <- c("cell","Correlation","p.value")
  correlation<-correlation[correlation$p.value<0.05,]
  correlation$sig[correlation$p.value>0.05] <- "ns"   
  correlation$sig[correlation$p.value<0.05&correlation$p.value>0.01] <- "*"   
  correlation$sig[correlation$p.value<0.01&correlation$p.value>0.001] <- "**"  
  correlation$sig[correlation$p.value<0.001&correlation$p.value>0.0001] <- "***"   
  correlation$sig[correlation$p.value<0.0001] <- "****"
  correlation$'correlation'<-abs(correlation$Correlation)
  #correlation_sig<-correlation[c(3,9,10,11,23,29,7,12,18,31,1,14,15,25,28),]
  library("viridis")
  #trace(ggdotchart,edit=T)
  p0<-ggdotchart(correlation, x = "cell", y = "Correlation",
                 
                 dot.size ='correlation',
                 # filled.contour(p.value,col=Lab.palette(20)),
                 color ='p.value',
                 # label='sig',
                 
                 #  font.label = list( size = 9),
                 
                 #ylab="mgp = c(3, 2, 0)",
                 # font.label = list(size=10,face="bold",color='black',position="dodge",
                 #                   vjust=0.5),
                 #  y.label=0.7,
                 # palette = palette,
                 # 按照cyl填充颜色
                 # palette = palette, # 修改颜色
                 sorting = "descending",
                 add = "segments",                             # 添加棒子
                 rotate = TRUE,
                 ggtheme = theme_pubr(),                        # 改变主题
                 ylab="Correlation Coefficient (r)",
                 xlab='',
                 # dot.size = 6 ,
                 title=rownames(hub_exp)[cnt]
                 
  )+scale_colour_gradient( high = "#4682B4",low = "#66CDAA")
  
  p10<-p0+theme(legend.position = "right",
                panel.background = element_blank())+geom_hline(aes(yintercept=0),linetype="dashed",lwd = 0.2)+
    theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
          axis.text.x =element_text(size=12,family = "Times", face = "bold"),
          axis.title.y =element_text(size=20,family = "Times", face = "bold"),
          axis.text.y=element_text(size=12,family = "Times", face = "bold"),
          plot.title=element_text(size=20,family = "Times", face = "bold",hjust=0.5),
          legend.text = element_text(size = 16, family = "Times"),
          legend.title = element_text(size = 18, family = "Times",face = "bold"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  p10
  ggsave(paste(cnt,rownames(hub_exp)[cnt],'pdf',sep='.'),w=8,h=6)
  ggsave(paste(cnt,rownames(hub_exp)[cnt],'png',sep='.'),w=8,h=6)
  
  cnt<-cnt+1
}


