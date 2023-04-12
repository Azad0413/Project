rm(list = ls())
setwd("/data/nas1/luchunlin/project/SJZZK-428-10/")
if (! dir.exists("./11_PCA")){
  dir.create("./11_PCA")
}
setwd("./11_PCA")

hubgene<-read.delim2('../08_Lasso/lasso_genes.csv',header = F)
dat = read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1) %>% lc.tableToNum
colnames(dat) <- gsub('.','-',colnames(dat),fixed = T)

group<-read.delim2("../09_risk/risk.xls")%>%select(c('id','riskScore'))
colnames(group) <- c('sample','group')
group$group <- as.numeric(group$group)
group$group <- ifelse(group$group>median(group$group),'High risk','Low risk')

dat <- dat[,group$sample]
dat <- log2(dat+1)
hub.exp<-dat[hubgene$V1,]
hub.exp<-hub.exp[order(rownames(hub.exp),decreasing = F),]  
cluster_exp<-t(hub.exp)
pca<-prcomp(cluster_exp,center = TRUE,scale. = TRUE)
df<-pca$x  ## 提取PC score
df <- as.data.frame(df)
summ <- summary(pca)
summ
# 提取主成分的方差贡献率,生成坐标轴标题
summ <- summary(pca)
summ
xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
# 绘制PCA得分图
table(group$group)
library(ggplot2)
p.pca <- ggplot(data = df,aes(x = PC1,y = PC2,color = group$group))+
  stat_ellipse(aes(fill = group$group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 3.5)+
  labs(x = xlab,y = ylab,color = "Group",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("purple","orange","pink"))+
  scale_colour_manual(values = c("purple","orange","pink"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.pca
ggsave('01.pca.pdf',p.pca,w=6,h=5)
ggsave('01.pca.png',p.pca,w=6,h=5)


## tSNE------
library(Rtsne)
set.seed(123123)
tsne_out <- Rtsne(cluster_exp,pca=FALSE,perplexity=10,theta=0.0)
# 获取tSNE的坐标值
str(tsne_out)
# 其中在Y中存储了画图坐标
tsnes=tsne_out$Y
colnames(tsnes) <- c("tSNE1", "tSNE2") #为坐标添加列名
# 在此基础上添加颜色分组信息，首先还是将tsnes这个矩阵变成数据框，然后增加一列group信息，最后映射在geom_point中
tsnes=as.data.frame(tsnes)
rownames(tsnes)
#group=c(rep('group1',cell_num),rep('group2',cell_num))
tsnes$group=group$group
p.tsne <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=group))+
  labs(color = "Group",title = "tSNE Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("purple","orange","pink"))+
  scale_colour_manual(values = c("purple","orange","pink"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.tsne
ggsave('02.tsne.pdf',p.tsne,w=6,h=5)
ggsave('02.tsne.png',p.tsne,w=6,h=5)


