rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-321")
if (! dir.exists("./18_PCA")){
  dir.create("./18_PCA")
}
setwd("./18_PCA")
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/08_va/hub_final.xls')
dat = read.delim2("/data/nas1/luchunlin/project/BJTC-321/00_rawdata/dat.final.xls", row.names = 1) %>% lc.tableToNum
hub.exp<-dat[hubgene$hubgene,]
group<-read.delim2("/data/nas1/luchunlin/project/BJTC-321/00_rawdata/group.xls")
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
  labs(x = xlab,y = ylab,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("purple","orange","pink"))+
  scale_colour_manual(values = c("purple","orange","pink"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.pca
ggsave('pca.pdf',p.pca,w=6,h=5)
ggsave('pca.png',p.pca,w=6,h=5)
