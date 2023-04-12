rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/BJTC-356/")
if (! dir.exists("./10_ssGSEA")){
  dir.create("./10_ssGSEA")
}
setwd("./10_ssGSEA")

dat <- read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
group <- read.delim2("../01_GSVA/group.xls")
table(group$group)

## 03-1 ssGSEA--------
library(GSVA)
gene_set <- read.table("/data/nas1/luchunlin/pipeline/ssGSEA/mmc3.txt",
                       header = T,
                       sep ="\t")
dat.final2 <- as.matrix(dat)
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

## 03-2 差异免疫细胞鉴定-------
## 差异-------
dat.ssgsea <- ssgsea_score %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
dat.ssgsea <- merge(group, dat.ssgsea, by = "sample")
dat.ssgsea2 <- tidyr::gather(dat.ssgsea, ImmuneCell, Score, -c("sample", "group"))
library(rstatix)
library(ggplot2)
library(ggpubr)
colnames(dat.ssgsea2)
stat_ssgsea <- dat.ssgsea2 %>% 
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
write.table(stat_ssgsea,file = 'stat.ssgsea.xls',sep = '\t',row.names = F,quote = F)
DE.ssgsea<-stat_ssgsea[which(stat_ssgsea$p<0.05),]
write.table(DE.ssgsea,file = 'DE.ssgsea.xls',sep = '\t',row.names = F,quote = F)
colnames(dat.ssgsea2)
violin.ssgsea<-dat.ssgsea2[dat.ssgsea2$ImmuneCell%in%stat_ssgsea$ImmuneCell[which(stat_ssgsea$p<0.05)],]
ssgsea_plot <- ggplot(violin.ssgsea, aes(x=group,
                                               y=Score,
                                               fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="Immune Cell", x="", y = "Score",size=20) +
  stat_compare_means(data = violin.ssgsea,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.45) +
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
        panel.grid.minor = element_blank())+facet_wrap(~ImmuneCell,scales = "free",nrow = 1) 
ssgsea_plot
ggsave(filename = '01.ssgsea.plot.pdf',ssgsea_plot,w=6,h=4)
ggsave(filename = '01.ssgsea.plot.png',ssgsea_plot,w=6,h=4)


## spearman 生物标志物与差异免疫浸润细胞的相关性。
hubgene <- read.delim2('../06_quadrant/final.gene.xls')
hub_exp <- dat[hubgene$symbol,]%>%lc.tableToNum()
DEtiic <- t(ssgsea_score)%>%as.data.frame()
#DEtiic <- DEtiic[,stat_ssgsea$ImmuneCell[which(stat_ssgsea$p<0.05)]]

cnt<-1
while(cnt<10){
  library(psych)
  x <-as.numeric(hub_exp[cnt,])
  y <-DEtiic
  library(psych)
  d <- corr.test(y,x,use="complete",method = 'spearman')
  r <- data.frame(d$r)
  p <- data.frame(d$p)
  write.table(p,file=paste(cnt+1,rownames(hub_exp)[cnt], "p.xls",sep='_'),sep="\t",quote=F)
  write.table(r,file=paste(cnt+1,rownames(hub_exp)[cnt], "r.xls",sep='_'),sep="\t",quote=F)
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
  ggsave(paste(cnt+1,rownames(hub_exp)[cnt],'pdf',sep='.'),w=6,h=4)
  ggsave(paste(cnt+1,rownames(hub_exp)[cnt],'png',sep='.'),w=6,h=4)
  
  cnt<-cnt+1
}
