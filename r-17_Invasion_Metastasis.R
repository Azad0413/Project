rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-399-11/")
if (! dir.exists("./17_Invasion_Metastasis")){
  dir.create("./17_Invasion_Metastasis")
}
setwd("./17_Invasion_Metastasis")
#PARG 1 相关性----------
#风险评分与侵袭相关基因、EMT相关基因、血管生成相关基因的表达相关性。
risk <- read.delim2('../08_risk/risk.xls')
dat <- read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat) <- gsub('.','-',colnames(dat),fixed = T)
dat <- dat[,risk$id]
riskScore <- risk%>%select(c('id','riskScore'))%>%column_to_rownames(var = 'id')
riskScore$riskScore <- as.numeric(riskScore$riskScore)
riskScore <- t(riskScore)%>%as.data.frame()
###侵袭------
invasion <- read_xlsx('invasion.xlsx')
invasion_exp <- dat[invasion$GeneName,]
invasion_exp <- na.omit(invasion_exp)
library(Hmisc)
gene.exp <- invasion_exp

nc<-t(rbind(riskScore,gene.exp))
m=rcorr(nc)$r[1:nrow(riskScore),(ncol(nc)-length(rownames(gene.exp))+1):ncol(nc)]
#m<-t(m)

p=rcorr(nc)$P[1:nrow(riskScore),(ncol(nc)-length(rownames(gene.exp))+1):ncol(nc)]
#p<-t(p)

library(dplyr)
library(dplyr)
cor.all <- t(rbind(m,p))
colnames(cor.all) <- c('correlation','pvalue')
write.table(cor.all,file = 'correlation(invasion).xls',sep = '\t',row.names = F,quote = F)
##top20
cor.final <- data.frame(cor.all)
cor.final <- cor.final[order(abs(cor.final$correlation),decreasing = T),]
cor.final <- cor.final[c(1:20),]

m <- t(m[rownames(cor.final)])
p <- t(p[rownames(cor.final)])
rownames(m) <- 'riskScore'
rownames(p) <- 'riskScore'

tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- signif(cor,3)
#cor[abs(cor)<0.1] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
library(WGCNA)
pdf(file = '01.correlation(Invasion).pdf',w=10,h=3.5)
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
               main = paste("riskScore-Invasion correlation"))
dev.off()
png(file = '01.correlation(Invasion).png',w=800,h=250)
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
               main = paste("riskScore-Invasion correlation"))
dev.off()

###EMT---------
EMT <- read_xlsx('EMT.xlsx')
EMT_exp <- dat[EMT$GeneSymbol,]
EMT_exp <- na.omit(EMT_exp)
library(Hmisc)
gene.exp <- EMT_exp
nc<-t(rbind(riskScore,gene.exp))
m=rcorr(nc)$r[1:nrow(riskScore),(ncol(nc)-length(rownames(gene.exp))+1):ncol(nc)]
#m<-t(m)

p=rcorr(nc)$P[1:nrow(riskScore),(ncol(nc)-length(rownames(gene.exp))+1):ncol(nc)]
#p<-t(p)

library(dplyr)
library(dplyr)
cor.all <- t(rbind(m,p))
colnames(cor.all) <- c('correlation','pvalue')
write.table(cor.all,file = 'correlation(EMT).xls',sep = '\t',row.names = F,quote = F)
##top20
cor.final <- data.frame(cor.all)
cor.final <- cor.final[order(abs(cor.final$correlation),decreasing = T),]
cor.final <- cor.final[c(1:20),]

m <- t(m[rownames(cor.final)])
p <- t(p[rownames(cor.final)])
rownames(m) <- 'riskScore'
rownames(p) <- 'riskScore'

tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- signif(cor,3)
#cor[abs(cor)<0.1] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
library(WGCNA)
pdf(file = '02.correlation(EMT).pdf',w=10,h=3.5)
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
               main = paste("riskScore-EMT correlation"))
dev.off()
png(file = '02.correlation(EMT).png',w=800,h=250)
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
               main = paste("riskScore-EMT correlation"))
dev.off()

###血管生成---------
library(readxl)
Angiogenesis <- read_xlsx('angiogenesis.xlsx')
Angiogenesis_exp <- dat[Angiogenesis$Symbol,]
Angiogenesis_exp <- na.omit(Angiogenesis_exp)
library(Hmisc)
gene.exp <- Angiogenesis_exp
nc<-t(rbind(riskScore,gene.exp))
m=rcorr(nc)$r[1:nrow(riskScore),(ncol(nc)-length(rownames(gene.exp))+1):ncol(nc)]
#m<-t(m)

p=rcorr(nc)$P[1:nrow(riskScore),(ncol(nc)-length(rownames(gene.exp))+1):ncol(nc)]
#p<-t(p)

library(dplyr)
library(dplyr)
cor.all <- t(rbind(m,p))
colnames(cor.all) <- c('correlation','pvalue')
write.table(cor.all,file = 'correlation(Angiogenesis).xls',sep = '\t',row.names = F,quote = F)
##top20
cor.final <- data.frame(cor.all)
cor.final <- cor.final[which(cor.final$pvalue<0.05),]
cor.final <- cor.final[order(abs(cor.final$correlation),decreasing = T),]
#cor.final <- cor.final[c(1:20),]

m <- t(m[rownames(cor.final)])
p <- t(p[rownames(cor.final)])
rownames(m) <- 'riskScore'
rownames(p) <- 'riskScore'

tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- signif(cor,3)
#cor[abs(cor)<0.1] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
library(WGCNA)
pdf(file = '03.correlation(Angiogenesis).pdf',w=10,h=3.5)
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
               main = paste("riskScore-Angiogenesis correlation"))
dev.off()
png(file = '03.correlation(Angiogenesis).png',w=800,h=250)
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
               main = paste("riskScore-Angiogenesis correlation"))
dev.off()

#PART2 评分---------
###侵袭------
library(GSVA)
geneset <- invasion%>%select('GeneName')
geneset$invasion<- c(rep('invasion',97))

group<-risk%>%select(c('id','riskScore'))
group$riskScore <- ifelse(group$riskScore>median(group$riskScore),'High risk','Low risk')
colnames(group) <- c('sample','group')
table(group$group)

gene_list <- split(as.matrix(geneset)[,1],
                   geneset[,2])
dat2 <- as.matrix(dat)
score = gsva(dat2, gene_list,
             method = "ssgsea",
             ssgsea.norm = TRUE,
             verbose = TRUE)
write.table(score,
            file = "Invasion_score.xls",
            sep = "\t",
            quote = F)

Invasion_score<-t(score)%>%as.data.frame()
violin_dat <- Invasion_score
High_sample <- group$sample[which(group$group=='High risk')]
violin_dat$group <- ifelse(rownames(violin_dat)%in%High_sample,'High risk','Low risk')
colnames(violin_dat)
library(rstatix)
stat.test<-violin_dat%>%
 # group_by(gene)%>%
  wilcox_test(invasion ~ group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox_result(invasion).xls',
            sep = '\t',
            row.names = F)

violin_plot <- ggplot(violin_dat, aes(x=group, 
                                      y=invasion,
                                      fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#F08080", "#48D1CC"), name = "Group")+
  labs(title="Invasion", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.5) +
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
        panel.grid.minor = element_blank())#+facet_wrap(~gene,scales = "free",nrow = 3) 
violin_plot
ggsave('04.Invasion.pdf',violin_plot,w=4,h=4.5)
ggsave('04.Invasion.png',violin_plot,w=4,h=4.5)

###EMT------
library(GSVA)
geneset <- EMT%>%select('GeneSymbol')
geneset$EMT<- c(rep('EMT',1184))
gene_list <- split(as.matrix(geneset)[,1],
                   geneset[,2])
score = gsva(dat2, gene_list,
             method = "ssgsea",
             ssgsea.norm = TRUE,
             verbose = TRUE)
write.table(score,
            file = "EMT_score.xls",
            sep = "\t",
            quote = F)

EMT_score<-t(score)%>%as.data.frame()
violin_dat <- EMT_score
High_sample <- group$sample[which(group$group=='High risk')]
violin_dat$group <- ifelse(rownames(violin_dat)%in%High_sample,'High risk','Low risk')
colnames(violin_dat)
library(rstatix)
stat.test<-violin_dat%>%
  # group_by(gene)%>%
  wilcox_test(EMT ~ group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox_result(EMT).xls',
            sep = '\t',
            row.names = F)

violin_plot <- ggplot(violin_dat, aes(x=group, 
                                      y=EMT,
                                      fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#F08080", "#48D1CC"), name = "Group")+
  labs(title="EMT", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.5) +
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
        panel.grid.minor = element_blank())#+facet_wrap(~gene,scales = "free",nrow = 3) 
violin_plot
ggsave('05.EMT.pdf',violin_plot,w=4,h=4.5)
ggsave('05.EMT.png',violin_plot,w=4,h=4.5)


###Angiogenesis------
library(GSVA)
geneset <- Angiogenesis%>%select('Symbol')
geneset$Angiogenesis<- c(rep('Angiogenesis',36))
gene_list <- split(as.matrix(geneset)[,1],
                   geneset[,2])
score = gsva(dat2, gene_list,
             method = "ssgsea",
             ssgsea.norm = TRUE,
             verbose = TRUE)
write.table(score,
            file = "Angiogenesis_score.xls",
            sep = "\t",
            quote = F)

Angiogenesis_score<-t(score)%>%as.data.frame()
violin_dat <- Angiogenesis_score
violin_dat$group <- ifelse(rownames(violin_dat)%in%High_sample,'High risk','Low risk')
colnames(violin_dat)
library(rstatix)
stat.test<-violin_dat%>%
  # group_by(gene)%>%
  wilcox_test(Angiogenesis ~ group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox_result(Angiogenesis).xls',
            sep = '\t',
            row.names = F)

violin_plot <- ggplot(violin_dat, aes(x=group, 
                                      y=Angiogenesis,
                                      fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#F08080", "#48D1CC"), name = "Group")+
  labs(title="Angiogenesis", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.5) +
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
        panel.grid.minor = element_blank())#+facet_wrap(~gene,scales = "free",nrow = 3) 
violin_plot
ggsave('06.Angiogenesis.pdf',violin_plot,w=4,h=4.5)
ggsave('06.Angiogenesis.png',violin_plot,w=4,h=4.5)
