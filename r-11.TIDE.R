rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-255-2/")
if (! dir.exists("./11_TIDE")){
  dir.create("./11_TIDE")
}
setwd("./11_TIDE")

dat.tcga<-read.csv("../00_rawdata/dat.fpkm.xls", row.names = 1,sep = '\t')
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
risk<-read.delim2('../06_risk/risk.xls')
dat.tcga<-dat.tcga[,risk$id]
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
# tide_dat[which(tide_dat<0)] <- 0
tide_dat <- apply(dat.tcga,2,FPKM2TPM)
tide_dat <- log2(tide_dat + 1)
rownmean <- apply(tide_dat,1,mean)
tide_dat2 <- sweep(tide_dat, 1, rownmean)
dim(tide_dat2)
write.table(tide_dat2,
            file ="tide_dat.txt",
            sep = "\t",
            quote = F,
            row.names = T)


tide_result <- read.csv("02.tide_score.csv",header = T)
colnames(tide_result)
#View(tide_result)
tide_result2 <- subset(tide_result, select = c("Patient", "TIDE","Dysfunction","Exclusion"))
tide_plot_dat <- data.frame(Patient=risk$id,
                            riskScore=risk$riskScore,
                            risk_group=risk$risk)
tide_plot_dat <- merge(tide_plot_dat, tide_result2, by = "Patient")
tide_plot_dat$risk_group <- factor(tide_plot_dat$risk_group,
                                   level = c(0, 1),
                                   labels = c("High Risk", "Low Risk"))

dim(tide_plot_dat)
colnames(tide_plot_dat)
class(tide_plot_dat$riskScore)
tide_plot_dat$riskScore <- as.numeric(tide_plot_dat$riskScore)
p1<-ggplot(tide_plot_dat,aes(x=riskScore,y=Exclusion))+geom_point(aes(x=riskScore,y=Exclusion,color =risk_group), size=2,alpha=0.9)+geom_smooth(method = 'lm', formula = y ~ x, se = T,color='black')+ stat_cor(data=tide_plot_dat, method = "spearman")+
  scale_fill_manual(values =c("black")) +
  theme_bw()+
  scale_color_manual(values = c("#DB3587","#35B3B2"), name = "Risk") +
  theme(axis.title = element_text(size = 20, face = "bold", family = "Times",color='black'),
        axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
        axis.text.y = element_text(size = 16, face = "bold", family = "Times",color='black'),
        legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
        legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
        legend.position = "bottom"
  )+
  labs(x="Risk score",y="Exclusion score")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p1
library(ggsci)
table(tide_plot_dat$risk_group)
my_comparisons = list( c("High Risk", "Low Risk"))
p2<-ggplot(tide_plot_dat, aes(x=risk_group, y=Exclusion,fill=risk_group)) + 
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",
                     color='black',family="Times",face = "bold")+
  scale_fill_npg()+
  scale_fill_manual(values=c("#DB3587","#35B3B2")) + 
  labs(x = "", y = "", title = "") + 
  #theme_bw() + 
  # geom_text(aes(label = B, vjust = 1.1, hjust = -0.5, angle = 45), show_guide = FALSE) + 
  theme(panel.grid =element_blank()) +   
  theme(axis.text = element_blank()) +  
  theme(axis.ticks = element_blank())+
  theme(panel.background = element_blank())+
  guides(fill=F) 
p2
library(patchwork)
p<-p1+p2+plot_layout(widths = c(2, 0.3))
p
pdf('01.Exclusion_boxplot.pdf',width=7,height=6,family='Times')
print(p)
dev.off()

png('01.Exclusion_boxplot.png',width=7,height=6,family='Times',units='in',res=600)
print(p)
dev.off()


colnames(tide_plot_dat)
p1<-ggplot(tide_plot_dat,aes(x=riskScore,y=Dysfunction))+geom_point(aes(x=riskScore,y=Dysfunction,color =risk_group), size=2,alpha=0.9)+geom_smooth(method = 'lm', formula = y ~ x, se = T,color='black')+ stat_cor(data=tide_plot_dat, method = "spearman")+
  scale_fill_manual(values =c("black")) +
  theme_bw()+
  scale_color_manual(values = c("#DB3587","#35B3B2"), name = "Risk") +
  theme(axis.title = element_text(size = 20, face = "bold", family = "Times",color='black'),
        axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
        axis.text.y = element_text(size = 16, face = "bold", family = "Times",color='black'),
        legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
        legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
        legend.position = "bottom"
  )+
  labs(x="Risk score",y="Dysfunction score")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p1
library(ggsci)
table(tide_plot_dat$risk_group)
my_comparisons = list( c("High Risk", "Low Risk"))
p2<-ggplot(tide_plot_dat, aes(x=risk_group, y=Dysfunction,fill=risk_group)) + 
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",
                     color='black',family="Times",face = "bold")+
  scale_fill_npg()+
  scale_fill_manual(values=c("#DB3587","#35B3B2")) + 
  labs(x = "", y = "", title = "") + 
  #theme_bw() + 
  # geom_text(aes(label = B, vjust = 1.1, hjust = -0.5, angle = 45), show_guide = FALSE) + 
  theme(panel.grid =element_blank()) +   
  theme(axis.text = element_blank()) +  
  theme(axis.ticks = element_blank())+
  theme(panel.background = element_blank())+
  guides(fill=F) 
p2
library(patchwork)
p<-p1+p2+plot_layout(widths = c(2, 0.3))
p
pdf('02.Dysfunction_boxplot.pdf',width=7,height=6,family='Times')
print(p)
dev.off()

png('02.Dysfunction_boxplot.png',width=7,height=6,family='Times',units='in',res=600)
print(p)
dev.off()

colnames(tide_plot_dat)
p1<-ggplot(tide_plot_dat,aes(x=riskScore,y=TIDE))+geom_point(aes(x=riskScore,y=TIDE,color =risk_group), size=2,alpha=0.9)+geom_smooth(method = 'lm', formula = y ~ x, se = T,color='black')+ stat_cor(data=tide_plot_dat, method = "spearman")+
  scale_fill_manual(values =c("black")) +
  theme_bw()+
  scale_color_manual(values = c("#DB3587","#35B3B2"), name = "Risk") +
  theme(axis.title = element_text(size = 20, face = "bold", family = "Times",color='black'),
        axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
        axis.text.y = element_text(size = 16, face = "bold", family = "Times",color='black'),
        legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
        legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
        legend.position = "bottom"
  )+
  labs(x="Risk score",y="TIDE score")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p1
library(ggsci)
table(tide_plot_dat$risk_group)
my_comparisons = list( c("High Risk", "Low Risk"))
p2<-ggplot(tide_plot_dat, aes(x=risk_group, y=TIDE,fill=risk_group)) + 
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",
                     color='black',family="Times",face = "bold")+
  scale_fill_npg()+
  scale_fill_manual(values=c("#DB3587","#35B3B2")) + 
  labs(x = "", y = "", title = "") + 
  #theme_bw() + 
  # geom_text(aes(label = B, vjust = 1.1, hjust = -0.5, angle = 45), show_guide = FALSE) + 
  theme(panel.grid =element_blank()) +   
  theme(axis.text = element_blank()) +  
  theme(axis.ticks = element_blank())+
  theme(panel.background = element_blank())+
  guides(fill=F) 
p2
library(patchwork)
p<-p1+p2+plot_layout(widths = c(2, 0.3))
p
pdf('03.TIDE_boxplot.pdf',width=7,height=6,family='Times')
print(p)
dev.off()

png('03.TIDE_boxplot.png',width=7,height=6,family='Times',units='in',res=600)
print(p)
dev.off()
