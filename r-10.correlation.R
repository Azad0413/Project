rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-321")
if (! dir.exists("./10_cor")){
  dir.create("./10_cor")
}
setwd("./10_cor")
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/08_va/hub_final.xls')
dat = read.delim2("/data/nas1/luchunlin/project/BJTC-321/00_rawdata/dat.final.xls", row.names = 1) %>% lc.tableToNum
library(ggcorrplot)
library(corrplot)
## 提取诊断基因表达矩阵
exp<-dat[hubgene$hubgene,]
hub_corr<-round(cor(t(exp)),3)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
hub_p.mat<-round(cor_pmat(t(exp)),3)
col1 <- colorRampPalette(c("#4169E1","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF",
                           "#FDDBC7","#F4A582","#D6604D", "#B2182B","#CC0000","#990000"))
pdf('correlatin.pdf',w=6,h=6)
hub_corr_plot<-corrplot(hub_corr,
                        method = "circle",
                        is.corr = T,
                        type = "lower",
                        p.mat = hub_p.mat,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        tl.col = 'black',
                        col = col1(200))
dev.off()
png('correlatin.png',w=500,h=500)
hub_corr_plot<-corrplot(hub_corr,
                        method = "circle",
                        is.corr = T,
                        type = "lower",
                        p.mat = hub_p.mat,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        tl.col = 'black',
                        col = col1(200))
dev.off()
write.table(hub_corr,file = 'hub_corr.xls',
            sep = '\t',
            row.names = T)
## 挑取两对最正相关的基因 共表达关系图
##  CST1  CST4
library(ggpubr)
cor_1<-data.frame(t(exp[c('CST1','CST4'),]))
cor1_plot<-ggplot(cor_1,aes(x=CST1,y=CST4 ))+
  geom_point(color ="#B51F2E")+
  geom_smooth(method = 'lm', formula = y ~ x, se = F,color='red',size=0.3,)+
  stat_cor(data=cor_1, method = "pearson")+
  scale_fill_manual(values = c("#2774C4", "#BEBEBE")) +theme_bw()+
  theme(axis.text.x=element_text(hjust = 1,colour="black",size=16), 
        axis.text.y=element_text(size=16,face="plain"), 
        axis.title.y=element_text(size = 15,face="plain"),
        axis.title.x=element_text(size = 15,face="plain"))+
  labs(x="CST1",y="CST4")+
  theme(panel.grid =element_blank()) +   ## 删去网格线
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank())+
  theme(panel.background = element_blank(),axis.text.x= element_blank(),axis.text.y= element_blank())+#去除背景
  guides(fill=F)  ##去除图例
cor1_plot
ggsave('cor_plot.pdf',cor1_plot,w=4,h=4)
ggsave('cor_plot.png',cor1_plot,w=4,h=4)
# cor_2<-data.frame(t(hub_exp[c('FCER1G','C1QB'),]))
# cor2_plot<-ggplot(cor_2,aes(x=FCER1G,y=C1QB ))+
#   geom_point(color ="#B51F2E")+
#   geom_smooth(method = 'lm', formula = y ~ x, se = F,color='red',size=0.3,)+
#   stat_cor(data=cor_2, method = "spearman")+
#   scale_fill_manual(values = c("#2774C4", "#BEBEBE")) +theme_bw()+
#   theme(axis.text.x=element_text(hjust = 1,colour="black",size=16), 
#         axis.text.y=element_text(size=16,face="plain"), 
#         axis.title.y=element_text(size = 15,face="plain"),
#         axis.title.x=element_text(size = 15,face="plain"))+
#   labs(x="FCER1G",y="C1QB")+
#   theme(panel.grid =element_blank()) +   ## 删去网格线
#   theme(axis.text = element_blank()) +   ## 删去所有刻度标签
#   theme(axis.ticks = element_blank())+
#   theme(panel.background = element_blank(),axis.text.x= element_blank(),axis.text.y= element_blank())+#去除背景
#   guides(fill=F)  ##去除图例
# cor2_plot