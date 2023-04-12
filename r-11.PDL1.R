rm(list = ls())
#09 checkpoint-------------
setwd("/data/nas1/luchunlin/project/BJTC-308/")
if (! dir.exists("./11_PDL1")){
  dir.create("./11_PDL1")
}
setwd("./11_PDL1")

# IMvigor210 --------------------------------------------------------------
library(IMvigor210CoreBiologies)
data(cds)   ##clinic 
expMatrix <- counts(cds)  ##count
eff_length2 <- fData(cds)[,c("entrez_id","length","symbol")]
rownames(eff_length2) <- eff_length2$entrez_id
head(eff_length2)
feature_ids <- rownames(expMatrix)
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2),]  
mm <- match(rownames(expMatrix),rownames(eff_length2))
eff_length2 <- eff_length2[mm,]
x <- expMatrix/eff_length2$length
eset <- t(t(x)/colSums(x))*1e6
summary(duplicated(rownames(eset)))   ##fpkm
library(IOBR)
eset <- IOBR::anno_eset(eset = eset,
                        annotation = eff_length2,
                        symbol = "symbol",
                        probe = "entrez_id",
                        method = "mean")   ##symbol 

# fpkm转TPM
# FPKM2TPM <- function(fpkm){
#   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
# }
# eset[which(eset<0)] <- 0
#eset <- apply(eset,2,FPKM2TPM)   ##tpm
pdata <- pData(cds)
colnames(pdata) <- gsub(colnames(pdata),pattern = " ",replacement = "_")
write.csv(eset,'eset.csv',quote=F)  ##TPM表达矩阵
write.csv(pdata,'pdata.csv',quote=F)  ##样本临床信息
colnames(pdata)
pdata <- rownames_to_column(pdata[,c("binaryResponse",
                                     "Immune_phenotype",
                                     "FMOne_mutation_burden_per_MB",
                                     "Neoantigen_burden_per_MB",
                                     "censOS","os")],var = "ID")
colnames(pdata)<-c("ID","BOR_binary","Immune_phenotype","TMB","TNB","fustat","futime")
table(pdata$BOR_binary)
pdata$response<-ifelse(pdata$BOR_binary=='NA','FALSE','TURE')
pdata$response[is.na(pdata$response)]<-'FALSE'
table(pdata$response)
#pdata<-pdata[!is.na(pdata$BOR_binary),]
# pdata$BOR_binary<-ifelse(pdata$BOR_binary=="CR/PR","R","NR")

save.image("03.IMvigor210.Rdata")
clinic<-data.frame(pdata)
#clinic<-clinic[!is.na(clinic$BOR_binary),]
rownames(clinic)<-clinic[,1]
clinic$futime<-round(clinic$futime*30)
load('/data/nas1/luchunlin/project/BJTC-308/model.RData')
exp<-log2(eset+1)
lasso_geneids<-as.character(lasso_geneids$lasso_geneids)
rt2<-exp[lasso_geneids,]
rt2<-subset(rt2,select=rownames(clinic))
rt2<-t(rt2)
rt1<-clinic
rt1$riskscore<-NA
rt1$risk<-NA
cnt<-1
while (cnt < 349) {
  clinic$riskscore[cnt]<-sum(astive.coefficients*subset(rt2,select=lasso_geneids)[cnt,])
  cnt = cnt + 1
}
cnt<-1
while (cnt < 349) {
  clinic$risk[cnt]=as.vector(ifelse(clinic$riskscore[cnt]>median(clinic$riskscore),0,1))
  cnt = cnt + 1
}
clinic$Risk<-ifelse(clinic$risk==1,'Low','High')
write.table(clinic,'clinic.txt',sep='\t',quote=F)


# CR/PR  SD/PD  小提琴图 ------------------------------------------------------
stat.test <- clinic %>%
  #group_by(dose) %>%  ##分组
  wilcox_test(riskscore~ BOR_binary) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")%>%
  mutate(y.position =0)
violin_dat<-clinic[!is.na(clinic$BOR_binary),]
palette <- c("#E89591","#75B0D6")
ggviolin(violin_dat, x = "BOR_binary", y = "riskscore",add = c("mean_sd"),
         fill = "BOR_binary",  palette = palette,
         size=0.5,xlab="", ylab="Risk score",ggtheme = theme_bw(), legend = "top",font.axis = 2)+ 
  stat_pvalue_manual(stat.test, label = "p.adj",family='Times',size=6)+
  theme(axis.title.x =element_text(size=20,family = "Times", face = "bold"),
        axis.text.x =element_text(size=16,family = "Times", face = "bold"),
        axis.title.y =element_text(size=20,family = "Times", face = "bold"),
        axis.text.y=element_text(size=16,family = "Times", face = "bold"),
        legend.text=element_text(size=14,family = "Times", face = "bold"),
        legend.title=element_text(size=15,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())#+
#  facet_wrap(~Immune_phenotype,scales = "free",nrow = 1) 
ggsave(width=7,height=8,'05.violin.pdf')
ggsave(width=7,height=8,'05.violin.png')

# KM ----------------------------------------------------------------------
kmfit<-survfit(Surv(futime,fustat) ~ risk, data =  clinic)
pdf('06.KM.pdf',w=7,h=7,onefile=F,family='Times')
s_surv<-ggsurvplot(kmfit,
                   pval = TRUE, conf.int = F,legend.labs=c("High","Low" ),legend.title="Risk", title="IMVigor210",
                   risk.table = TRUE,
                   risk.table.col = "strata",
                   #linetype = "strata",
                   # surv.median.line = "hv",
                   ggtheme = theme_bw(),
                   palette=c("#dc143c","#6495ed"))

s_surv$table <- s_surv$table + 
  labs(x = "Time (days)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 16, face = "bold", family = "Times"),
        plot.title = element_text(size = 20, face = "bold", family = "Times"),
        text = element_text(size = 18, face = "bold", family = "Times"))
s_surv$plot <- s_surv$plot + 
  labs(x = "Time (days)") +
  labs(y = "Survival probability") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
        axis.text = element_text(size = 14, face = "bold", family = "Times"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Times"),
        legend.text = element_text(size = 16, face = "bold", family = "Times"),
        text = element_text(size = 20, face = "bold", family = "Times"))
print(s_surv)
dev.off()

png('06.KM.png',w=7,h=7,family='Times',units='in',res=600 )
s_surv<-ggsurvplot(kmfit,
                   pval = TRUE, conf.int = F,legend.labs=c("High","Low" ),legend.title="Risk", title="IMVigor210",
                   risk.table = TRUE,
                   risk.table.col = "strata",
                   #linetype = "strata",
                   # surv.median.line = "hv",
                   ggtheme = theme_bw(),
                   palette=c("#dc143c","#6495ed"))

s_surv$table <- s_surv$table + 
  labs(x = "Time (days)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 16, face = "bold", family = "Times"),
        plot.title = element_text(size = 20, face = "bold", family = "Times"),
        text = element_text(size = 18, face = "bold", family = "Times"))
s_surv$plot <- s_surv$plot + 
  labs(x = "Time (days)") +
  labs(y = "Survival probability") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
        axis.text = element_text(size = 14, face = "bold", family = "Times"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Times"),
        legend.text = element_text(size = 16, face = "bold", family = "Times"),
        text = element_text(size = 20, face = "bold", family = "Times"))
print(s_surv)
dev.off()

# ROC (1,3,5)---------------------------------------------------------------------
library(survivalROC)
ROC <- clinic
cutoff_1<- 365*1
cutoff_2 <- 365*2
cutoff_3 <- 365*3
#cutoff_7 <- 365*7
year_1= survivalROC(Stime=ROC$futime,
                    status=ROC$fustat,
                    marker = ROC$riskscore,
                    predict.time = cutoff_1,
                    method = 'KM')
year_2= survivalROC(Stime=ROC$futime,
                    status=ROC$fustat,
                    marker = ROC$riskscore, 
                    predict.time = cutoff_2,
                    method = 'KM')

year_3= survivalROC(Stime=ROC$futime,
                    status=ROC$fustat,
                    marker = ROC$riskscore, 
                    predict.time = cutoff_3,
                    method = 'KM')

pdf('01.ROC.pdf',w=9,h=8.5)
par(pin = c(4,4), mar = c(5,6,6,1),family="serif")      ##par(family="serif")  :all
# pin: 以英寸表示的图形尺寸（宽和高）
# mai: 以数值向量表示的边界大小，顺序为“下、左、上、右”，单位为英寸
# mar: 以数值向量表示的边界大小，顺序为“下、左、上、右”，单位为英分。默认（5,4,4,2）+ 0.1
plot(year_1$FP, year_1$TP,
     type="l",col="red",xlim=c(0,1), ylim=c(0,1),
     xlab="FP",
     ylab="TP",
     lwd=2,  #lwd曲线宽度
     #main="TCGA-GBM, Method=KM\n Year = 1, 3, 5",
     main="ESCC, Method=KM\n Year = 1, 2, 3",
     cex.axis=1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=2.0,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=2.0,   ##标题的缩放倍数。
     #family = "serif", 
     font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)  

abline(0,1,col="gray",lty=2)
lines(year_2$FP, year_2$TP, type="l",col="green",xlim=c(0,1), ylim=c(0,1))
lines(year_3$FP, year_3$TP, type="l",col="blue",xlim=c(0,1), ylim=c(0,1), lwd=2)
#lines(year_7$FP, year_5$TP, type="l",col="blue",xlim=c(0,1), ylim=c(0,1), lwd=2)

legend(0.55,0.2,c(
  paste("AUC of 1 year =",round(year_1$AUC,3)),
  paste("AUC of 2 year =",round(year_2$AUC,3)),
  paste("AUC of 3 year =",round(year_3$AUC,3))
  #paste("AUC of 7 year =",round(year_7$AUC,3))
),
x.intersp=1, y.intersp=0.8,
lty= 1 ,lwd= 2,col=c("red","green",'blue'),
bty = "n",
seg.len=1,cex=1.5,
text.font = 2)
dev.off()

png('01.ROC.png',w=9,h=8.5,units='in',res=600)
par(pin = c(4,4), mar = c(5,6,6,1),family="serif")      ##par(family="serif")  :all
# pin: 以英寸表示的图形尺寸（宽和高）
# mai: 以数值向量表示的边界大小，顺序为“下、左、上、右”，单位为英寸
# mar: 以数值向量表示的边界大小，顺序为“下、左、上、右”，单位为英分。默认（5,4,4,2）+ 0.1
plot(year_1$FP, year_1$TP,
     type="l",col="red",xlim=c(0,1), ylim=c(0,1),
     xlab="FP",
     ylab="TP",
     lwd=2,  #lwd曲线宽度
     #main="TCGA-GBM, Method=KM\n Year = 1, 3, 5",
     main="TCGA-GBM, Method=KM\n Year = 1, 3, 5",
     cex.axis=1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=2.0,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=2.0,   ##标题的缩放倍数。
     #family = "serif", 
     font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)  

abline(0,1,col="gray",lty=2)
lines(year_2$FP, year_2$TP, type="l",col="green",xlim=c(0,1), ylim=c(0,1))
lines(year_3$FP, year_3$TP, type="l",col="blue",xlim=c(0,1), ylim=c(0,1), lwd=2)
#lines(year_7$FP, year_5$TP, type="l",col="blue",xlim=c(0,1), ylim=c(0,1), lwd=2)

legend(0.55,0.2,c(
  paste("AUC of 1 year =",round(year_1$AUC,3)),
  paste("AUC of 2 year =",round(year_2$AUC,3)),
  paste("AUC of 3 year =",round(year_3$AUC,3))
  #paste("AUC of 7 year =",round(year_7$AUC,3))
),
x.intersp=1, y.intersp=0.8,
lty= 1 ,lwd= 2,col=c("red","green",'blue'),
bty = "n",
seg.len=1,cex=1.5,
text.font = 2)
dev.off()


library(RColorBrewer)
library(ROCR)
library(pROC)
Affair<-clinic
Affair$outcome<-ifelse(Affair$fustat==1,'Dead','Alive')
outcome <-Affair$outcome
roc <- roc(outcome ~ riskscore, Affair, levels=c('Alive', 'Dead'))
pdf('07.ROC.pdf',w=8,h=8,family='Times')
plot(roc, 
     print.auc=TRUE, #设置是否添加AUC值标签
     auc.polygon=TRUE, #设置是否添加AUC值面积多边形
     # max.auc.polygon=TRUE, #设置是否添加最大AUC值面积多边形
     auc.polygon.col="skyblue", #设置AUC值面积多边形的填充色
     # grid=c(0.1, 0.2), #添加网格线
     # grid.col=c("green", "red"), #设置网格线颜色
     print.thres=TRUE,
     cex.axis=1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
dev.off()

png('07.ROC.png',w=8,h=8,family='Times', units='in',res=600)
plot(roc, 
     print.auc=TRUE, #设置是否添加AUC值标签
     auc.polygon=TRUE, #设置是否添加AUC值面积多边形
     # max.auc.polygon=TRUE, #设置是否添加最大AUC值面积多边形
     auc.polygon.col="skyblue", #设置AUC值面积多边形的填充色
     # grid=c(0.1, 0.2), #添加网格线
     # grid.col=c("green", "red"), #设置网格线颜色
     print.thres=TRUE,
     cex.axis=1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.lab=1.8,   ##坐标轴刻度文字的缩放倍数。类似cex。
     cex.main=1.8,   ##标题的缩放倍数。
     font.axis = 2, 
     font.lab = 2, 
     font.main = 2, 
     font.sub =2)
dev.off()
## 柱状堆叠图-----
library(RColorBrewer)
colnames(clinic)
box_dat<-clinic[,c('response','Risk','riskscore')]
box_dat$percentage<-1/174

box_dat$percentage<-as.numeric(box_dat$percentage)
box_plot <- ggplot(box_dat, aes(x=Risk, 
                                y=percentage,
                                fill=response)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='Percentage',
       fill='')+
  theme(
    # axis.text.x = element_blank(),
    #     axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = c("#E89591","#75B0D6"))#+
#  facet_grid(~box_dat$Risk,scales= "free",space= "free")
box_plot
## 小提琴图---------
colnames(box_dat)
violin_plot <- ggplot(box_dat, aes(x=response, 
                                      y=riskscore,
                                      fill=response)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#E89591","#75B0D6"), name = "Group")+
  labs(title="", x="", y = "riskScore",size=20) +
  stat_compare_means(data = box_dat,
                     mapping = aes(group = response),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
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
        panel.grid.minor = element_blank())
violin_plot
