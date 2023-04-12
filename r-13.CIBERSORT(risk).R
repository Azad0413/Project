rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/HF-0103-1/")
if (! dir.exists("./13_CIBERSORT(risk)")){
  dir.create("./13_CIBERSORT(risk)")
}
setwd("./13_CIBERSORT(risk)")

dat<-read.delim2("../00_rawdata/dat(GSE10846).xls", row.names = 1)%>% lc.tableToNum
# colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group<-read.delim2('../07_risk/risk.xls')%>%dplyr::select(c('id','risk'))
colnames(group)<-c('sample','group')
group$group <- ifelse(group$group=='1','Low risk','High risk')
Low.sample<-group$sample[which(group$group=='Low risk')]
High.sample<-group$sample[which(group$group=='High risk')]
dat <- dat[,group$sample]
dat <- dat^2
library(immunedeconv)
set_cibersort_binary("CIBERSORT.R")
set_cibersort_mat("LM22.txt")
res.cibersort<-deconvolute(as.matrix(dat),method = 'cibersort')

write.table(res.cibersort,'cibersort.txt',sep = '\t',col.names = T,row.names = F,quote = F)
save(file = 'cibersort.Rdata',res.cibersort)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(7,"Paired"))
##画图
# pdf('01.cibersort.box.pdf',w=10,h=8)
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
# dev.off()
# png('01.cibersort.box.png',w=700,h=500)
# res.cibersort %>%
#   gather(sample, fraction, -cell_type) %>%
#   merge(group,by='sample')%>%
#   # 绘制堆积条形图
#   ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
#   geom_bar(position = 'stack',stat = 'identity')+
#   scale_y_continuous(expand = c(0,0))+
#   theme_bw()+
#   labs(x='',
#        y='Relative Percent',
#        fill='')+
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = 'top') +
#   scale_fill_manual(values = mypalette(22))+
#   facet_grid(~group,scales= "free",space= "free")
# dev.off()

## 差异-------
##去掉在75%的样本中结果为0的细胞
res.cibersort2 <- res.cibersort%>%column_to_rownames(var = 'cell_type')
keep <- rowSums(res.cibersort2>0)>=floor(0.75*ncol(res.cibersort2))
res.cibersort2 <- res.cibersort2[keep,]

dat.cibersort <- res.cibersort2 %>% 
  # tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
dat.cibersort <- merge(group, dat.cibersort, by = "sample")
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
cibersort_plot <- ggplot(violin.cibersort, aes(x=group,
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
  stat_compare_means(data = violin.cibersort,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.4) +
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
        panel.grid.minor = element_blank())+facet_wrap(~ImmuneCell,scales = "free",nrow = 2) 
cibersort_plot
ggsave(filename = '01.cibersort.plot.pdf',cibersort_plot,w=9,h=6)
ggsave(filename = '01.cibersort.plot.png',cibersort_plot,w=9.5,h=6)
### 风险分数与差异免疫细胞相关性----------
riskscore <- read.delim2('../07_risk/risk.xls')%>%select(c('id','riskScore'))%>%column_to_rownames(var = 'id')%>%t()

DEtiic <- res.cibersort%>%column_to_rownames(var = 'cell_type')
DEtiic <- DEtiic[DE.cibersort$ImmuneCell,]
DEtiic <- t(DEtiic)%>%as.data.frame()
cnt<-1
while(cnt<2){
  library(psych)
  x <-as.numeric(riskscore[cnt,])
  y <-DEtiic
  library(psych)
  d <- corr.test(y,x,use="complete",method = 'pearson')
  r <- data.frame(d$r)
  p <- data.frame(d$p)
  write.table(p,file=paste(rownames(riskscore)[cnt], "p.xls",sep='_'),sep="\t",quote=F)
  write.table(r,file=paste(rownames(riskscore)[cnt], "r.xls",sep='_'),sep="\t",quote=F)
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
                  title='Correlation between riskScore and DEcells'
                 
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
  ggsave(paste(cnt+1,rownames(riskscore)[cnt],'correlation','pdf',sep='.'),w=7,h=5)
  ggsave(paste(cnt+1,rownames(riskscore)[cnt],'correlation','png',sep='.'),w=7,h=5)
  
  cnt<-cnt+1
}

#### 基因与免疫浸润细胞相关性------
hubgene <- read.delim2('../06_Lasso/lasso_genes.csv',header = F)
diff_score <- t(DEtiic)

library(Hmisc)
hub_exp <- dat[hubgene$V1,]
nc<-t(rbind(diff_score,hub_exp))
m=rcorr(nc,type = 'spearman')$r[1:nrow(diff_score),(ncol(nc)-length(hubgene$V1)+1):ncol(nc)]
m<-t(m)
p=rcorr(nc,type = 'spearman')$P[1:nrow(diff_score),(ncol(nc)-length(hubgene$V1)+1):ncol(nc)]
p<-t(p)
library(dplyr)
library(dplyr)

tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- round(cor,digits = 3)
cor[abs(cor)<0.1] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
write.table(textMatrix,file = 'correlation.xls',sep = '\t',row.names = T,quote = F)
library(WGCNA)
pdf(file = '03.correlation.pdf',w=8,h=6)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("Biomarker-TIICs correlation"))
dev.off()
png(file = '03.correlation.png',w=600,h=400)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("Biomarker-TIICs correlation"))
dev.off()



