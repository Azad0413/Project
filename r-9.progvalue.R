rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./09_progvalue")){
  dir.create("./09_progvalue")
}
setwd("./09_progvalue")
dat<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
## 保留tumor
dat <- dat[,-c(1:72)]
dat <- log2(dat+1)
hubgene <- read.delim2('../07_hubgene/hubgene.xls')
dat <- t(dat[hubgene$symbol,])%>%as.data.frame()%>%rownames_to_column(var = 'sample')
survival<-read.delim2('../08_KM/survival.xls')%>%rownames_to_column(var = 'sample')
phenotype<-read.delim2('/data/nas1/luchunlin/TCGA_phenotype/TCGA-KIRC.GDC_phenotype.tsv.gz')
train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                            age=phenotype$age_at_initial_pathologic_diagnosis,
                            gender=phenotype$gender.demographic,
                            grade=phenotype$neoplasm_histologic_grade,
                            stage=phenotype$tumor_stage.diagnoses,
                            T.stage=phenotype$pathologic_T,
                            N.stage=phenotype$pathologic_N,
                            M.stage=phenotype$pathologic_M)
train_phenotype<-merge(train_phenotype,dat,by='sample')
train_phenotype<-merge(train_phenotype,survival,by='sample')
write.table(train_phenotype,file = 'phenotype.xls',sep = '\t',row.names = F,quote = F)
train_phenotype$OS<-as.numeric(train_phenotype$OS)
train_phenotype$OS.time<-as.numeric(train_phenotype$OS.time)
train_phenotype2<-train_phenotype
table(train_phenotype2$stage)
train_phenotype2$stage<-gsub('not reported',NA,train_phenotype2$stage)

train_phenotype2$stage<-gsub('stage iv','4',train_phenotype2$stage)
train_phenotype2$stage<-gsub('stage iii','3',train_phenotype2$stage)
train_phenotype2$stage<-gsub('stage ii','2',train_phenotype2$stage)
train_phenotype2$stage<-gsub('stage i','1',train_phenotype2$stage)

table(train_phenotype2$age)
train_phenotype2$age<-ifelse(train_phenotype2$age>60,'1','0')

table(train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('a','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('b','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('c','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T4','3',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T3','3',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T2','2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T1','1',train_phenotype2$T.stage)

table(train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('NX',NA,train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N0','0',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N1','1',train_phenotype2$N.stage)

table(train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('MX',NA,train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('M0','0',train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('M1','1',train_phenotype2$M.stage)
train_phenotype2$M.stage[train_phenotype2$M.stage==''] <- NA

table(train_phenotype2$grade)
train_phenotype2$grade<-gsub('GX',NA,train_phenotype2$grade)
train_phenotype2$grade<-gsub('G4','4',train_phenotype2$grade)
train_phenotype2$grade<-gsub('G3','3',train_phenotype2$grade)
train_phenotype2$grade<-gsub('G2','1/2',train_phenotype2$grade)
train_phenotype2$grade<-gsub('G1','1/2',train_phenotype2$grade)
train_phenotype2$grade[train_phenotype2$grade==''] <- NA

table(train_phenotype2$gender)
train_phenotype2$gender <- ifelse(train_phenotype2$gender=='male',1,0)
colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','Age','Gender','Grade','Stage','T.stage','N.stage','M.stage','CTD-2626G11.2','AP000696.2','RP11-528A4.2','LINC00645','RP4-655J12.4','RP11-321G12.1','RP11-195B3.1','PTCSC3','OS','OS.time')

train_risk_clinical <- train_phenotype2%>%column_to_rownames(var = 'id')

dim(train_risk_clinical)
colnames_train <- colnames(train_risk_clinical)
covariates_train <- colnames_train[-which(colnames_train %in% c("OS", "OS.time"))]

train_risk_clinical$Stage<-factor(train_risk_clinical$Stage)
train_risk_clinical$Grade <- factor(train_risk_clinical$Grade)
train_risk_clinical$T.stage <- factor(train_risk_clinical$T.stage)
train_risk_clinical$N.stage <- factor(train_risk_clinical$N.stage)
train_risk_clinical$M.stage <- factor(train_risk_clinical$M.stage)
library(survival)
res.age = coxph(Surv(time = OS.time, event = OS) ~ Age, data = train_risk_clinical) %>% summary
res.age = c(res.age$conf.int[-2], res.age$coefficients[5])
res.gender = coxph(Surv(time = OS.time, event = OS) ~ Gender, data = train_risk_clinical) %>% summary
res.gender = c(res.gender$conf.int[-2], res.gender$coefficients[5])
res.stage = coxph(Surv(time = OS.time, event = OS) ~ Stage, data = train_risk_clinical) %>% summary
res.stage = cbind(res.stage$conf.int[,-2], res.stage$coefficients[,5])
res.tstage = coxph(Surv(time = OS.time, event = OS) ~ T.stage, data = train_risk_clinical) %>% summary
res.tstage = cbind(res.tstage$conf.int[,-2], res.tstage$coefficients[,5])
res.nstage = coxph(Surv(time = OS.time, event = OS) ~ N.stage, data = train_risk_clinical) %>% summary
res.nstage = c(res.nstage$conf.int[-2], res.nstage$coefficients[5])
res.mstage = coxph(Surv(time = OS.time, event = OS) ~ M.stage, data = train_risk_clinical) %>% summary
res.mstage = c(res.mstage$conf.int[-2], res.mstage$coefficients[5])
res.grade = coxph(Surv(time = OS.time, event = OS) ~ Grade, data = train_risk_clinical) %>% summary
res.grade = cbind(res.grade$conf.int[,-2], res.grade$coefficients[,5])
colnames(train_risk_clinical) <- gsub('-','_',colnames(train_risk_clinical))

res.CTD_2626G11.2 = coxph(Surv(time = OS.time, event = OS) ~ CTD_2626G11.2, data = train_risk_clinical) %>% summary
res.CTD_2626G11.2 = c(res.CTD_2626G11.2$conf.int[-2], res.CTD_2626G11.2$coefficients[5])
res.AP000696.2 = coxph(Surv(time = OS.time, event = OS) ~ AP000696.2, data = train_risk_clinical) %>% summary
res.AP000696.2 = c(res.AP000696.2$conf.int[-2], res.AP000696.2$coefficients[5])
res.RP11_528A4.2 = coxph(Surv(time = OS.time, event = OS) ~ RP11_528A4.2, data = train_risk_clinical) %>% summary
res.RP11_528A4.2 = c(res.RP11_528A4.2$conf.int[-2], res.RP11_528A4.2$coefficients[5])
res.LINC00645 = coxph(Surv(time = OS.time, event = OS) ~ LINC00645, data = train_risk_clinical) %>% summary
res.LINC00645 = c(res.LINC00645$conf.int[-2], res.LINC00645$coefficients[5])
res.RP4_655J12.4 = coxph(Surv(time = OS.time, event = OS) ~ RP4_655J12.4, data = train_risk_clinical) %>% summary
res.RP4_655J12.4 = c(res.RP4_655J12.4$conf.int[-2], res.RP4_655J12.4$coefficients[5])
res.RP11_321G12.1 = coxph(Surv(time = OS.time, event = OS) ~ RP11_321G12.1, data = train_risk_clinical) %>% summary
res.RP11_321G12.1 = c(res.RP11_321G12.1$conf.int[-2], res.RP11_321G12.1$coefficients[5])
res.RP11_195B3.1 = coxph(Surv(time = OS.time, event = OS) ~ RP11_195B3.1, data = train_risk_clinical) %>% summary
res.RP11_195B3.1 = c(res.RP11_195B3.1$conf.int[-2], res.RP11_195B3.1$coefficients[5])
res.PTCSC3 = coxph(Surv(time = OS.time, event = OS) ~ PTCSC3, data = train_risk_clinical) %>% summary
res.PTCSC3 = c(res.PTCSC3$conf.int[-2], res.PTCSC3$coefficients[5])

colnames(train_risk_clinical)
res.ref = c(1,1,1,NA)
res = rbind(res.age,res.gender,res.ref,res.grade,res.ref,res.stage,res.ref,res.tstage,res.nstage,res.mstage,res.CTD_2626G11.2,res.AP000696.2,
            res.RP11_528A4.2,res.LINC00645,res.RP4_655J12.4,res.RP11_321G12.1,res.RP11_195B3.1,res.PTCSC3) %>% as.data.frame()
rownames(res)
res$Indicators = c("Age","Gender","Grade1/2(Reference)","Grade3","Grade4","Stage1(Reference)",'Stage2','Stage3','Stage4',
                   "T.stage1(Reference)","T.stage2","T.stage3","N1 vs.N0","M1 vs.M0","CTD_2626G11.2","AP000696.2","RP11_528A4.2","LINC00645",
                   "RP4_655J12.4","RP11_321G12.1","RP11_195B3.1","PTCSC3")
colnames(res) = c("hr","low","up","pv","Indicator")
res$p = signif(res$pv, 2) %>% paste0("p = ", .)
res$p[is.na(res$pv)] = NA
res$Indicator = factor(res$Indicator, levels = rev(res$Indicator))
rownames(res) <- res$Indicator
res2 <- data.frame(p.value=res$pv,
                   HR=res$hr,
                   HR.95L=res$low,
                   HR.95H=res$up,
                   Indicator=res$Indicator)
rownames(res2) <- res2$Indicator
write.table(res2, file = "univariate_cox_prog_forest.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))
library(tidyr)
hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
hz[c(3,6,10)] <- ""

tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.0001,
                                      "< 0.0001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
pdf(file = '01.uncoxforest.pdf',height = 7, width = 12, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,rep(FALSE, 11)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,1,2,3,4,5,6,7,8,9,10,11), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate") # 垂直于x轴的网格线，对应每个刻度
dev.off()

png(file = '01.uncoxforest.png',height = 600, width = 1000)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,rep(FALSE, 11)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,1,2,3,4,5,6,7,8,9,10,11), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate") # 垂直于x轴的网格线，对应每个刻度
dev.off()
## 07-2 多因素Cox----------
res.mul = coxph(Surv(time = OS.time, event = OS) ~ Age+Grade+Stage+T.stage+N.stage+M.stage+AP000696.2+PTCSC3, data = train_risk_clinical)%>% summary
res.mul = cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
rownames(res.mul)
res.mul = rbind(res.mul[c(1),], res.ref, res.mul[c(2:3),],res.ref,res.mul[c(4:6),],res.ref,res.mul[c(7:12),])

res.mul$Indicators = c("Age","Grade1/2(Reference)","Grade3","Grade4","Stage1(Reference)",'Stage2','Stage3','Stage4',
                       "T.stage1(Reference)","T.stage2","T.stage3","N1 vs.N0","M1 vs.M0","AP000696.2","PTCSC3")
colnames(res.mul) = c("hr","low","up","pv","Indicator")
res.mul$p = signif(res.mul$pv, 2) %>% paste0("p = ", .)
res.mul$p[is.na(res.mul$pv)] = NA
res.mul$Indicator = factor(res.mul$Indicator, levels = rev(res.mul$Indicator))
rownames(res.mul) <- res.mul$Indicator

multi_res <- data.frame(p.value=res.mul$pv,
                        HR=res.mul$hr,
                        HR.95L=res.mul$low,
                        HR.95H=res.mul$up,
                        Indicator=res.mul$Indicator)
rownames(multi_res) <- multi_res$Indicator
multi_res
write.csv(multi_res,
          file = "multivariate_cox_prog_result.csv",
          quote = F,
          row.names = T)
multi_res <- subset(multi_res, select = -c(Indicator))

library(tidyr)
hz <- paste(round(multi_res$HR,3),
            "(",round(multi_res$HR.95L,3),
            "-",round(multi_res$HR.95H,3),")",sep = "")
hz

hz[c(2,5,9)] <- ""
tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.0001,
                                      "< 0.0001",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
pdf(file = '02.mulcoxforest.pdf',height = 7, width = 12, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #95%置信区间下限
           upper=c(NA,multi_res$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,50,100,150,200,250,300), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate") # 垂直于x轴的网格线，对应每个刻度
dev.off()

png(file = '02.mulcoxforest.png',height = 600, width = 1000)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #95%置信区间下限
           upper=c(NA,multi_res$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,50,100,150,200,250,300), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate") # 垂直于x轴的网格线，对应每个刻度
dev.off()
##07-3 构建COX模型，绘制列线图---------
multi_cov<-c("Age","Grade","Stage","T.stage","N.stage","M.stage","AP000696.2","PTCSC3")
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~',
                                   paste(multi_cov,
                                         sep = '',
                                         collapse = '+')))


cox_more_prog <- coxph(cox_data_prog,
                       data = as.data.frame(train_risk_clinical))

# Nomogram
library(rms)
ddist <- datadist(train_risk_clinical)
options(datadist='ddist')

# 构建COX模型，绘制列线图
train_risk_clinical$Age <- ifelse(train_risk_clinical$Age==1,'>60','<=60')
res.cox <- psm(cox_data_prog,
               data = train_risk_clinical, dist = 'lognormal')
surv <- Survival(res.cox) # 构建生存概率函数
function(x) surv(365, x) # 1年事件发生概率
function(x) surv(1095, x) # 3年事件发生概率
function(x) surv(1825, x) # 5年事件发生概率


nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(365, x),
                               function(x) surv(1095, x),
                               function(x) surv(1825, x)),
                    funlabel=c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"),
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99),
                    lp=F)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
png(filename = "03.nomogram_line_points.png", height = 700, width = 1300)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
pdf(file = "03.nomogram_line_points.pdf", height = 9, width = 17)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
##07-4 构建校准曲线---------

coxm_1 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 1*365)
cal_1 <-calibrate(coxm_1,u=1*365,cmethod='KM',m=50,B=100)
coxm_3 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 3*365)
cal_3 <-calibrate(coxm_3,u=3*365,cmethod='KM',m=100,B=100)

coxm_5 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,
              y=T,
              time.inc = 365*5)
cal_5<-calibrate(coxm_5,u=365*5,cmethod='KM',m=100,B=100)

##绘制3年生存期校曲线
##time.in 和 u 要是一样的，都是要评价的时间节点

pdf(file = '04.calibrate.pdf',w=9,h=8)
par(mar=c(7,4,4,3),cex=1.5)

plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability)',#便签
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#00468b",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#便签
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#ed0000",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#便签
     ylab='Actual 1-5-year Survival Probability',#标签
     col="#42b540",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围

#加上图例
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()

png(file = '04.calibrate.png',w=700,h=600)
par(mar=c(7,4,4,3),cex=1.5)

plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability)',#便签
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#00468b",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#便签
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#ed0000",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#便签
     ylab='Actual 1-5-year Survival Probability',#标签
     col="#42b540",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围

#加上图例
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()


## 热图------

## DCA-----
library(rmda)
dca.dat <- train_risk_clinical
dca.dat$Age <- ifelse(dca.dat$Age=='>60',1,0)
complex<-decision_curve(OS ~Age+Grade+Stage+T.stage+N.stage+M.stage+AP000696.2+PTCSC3,data =dca.dat,family = binomial(link ='logit'),
                        thresholds = seq(0,1, by = 0.01),
                        confidence.intervals= 0.95,
                        study.design = 'case-control',
                        population.prevalence= 0.3
)

age.dca <- decision_curve(OS ~Age,data =dca.dat,family = binomial(link ='logit'),
                          thresholds = seq(0,1, by = 0.01),
                          confidence.intervals= 0.95,
                          study.design = 'case-control',
                          population.prevalence= 0.3
)
grade.dca <- decision_curve(OS ~Grade,data =dca.dat,family = binomial(link ='logit'),
                            thresholds = seq(0,1, by = 0.01),
                            confidence.intervals= 0.95,
                            study.design = 'case-control',
                            population.prevalence= 0.3
)

stage.dca <- decision_curve(OS ~Stage,data =dca.dat,family = binomial(link ='logit'),
                            thresholds = seq(0,1, by = 0.01),
                            confidence.intervals= 0.95,
                            study.design = 'case-control',
                            population.prevalence= 0.3
)

tstage.dca <- decision_curve(OS ~T.stage,data =dca.dat,family = binomial(link ='logit'),
                            thresholds = seq(0,1, by = 0.01),
                            confidence.intervals= 0.95,
                            study.design = 'case-control',
                            population.prevalence= 0.3
)

nstage.dca <- decision_curve(OS ~N.stage,data =dca.dat,family = binomial(link ='logit'),
                             thresholds = seq(0,1, by = 0.01),
                             confidence.intervals= 0.95,
                             study.design = 'case-control',
                             population.prevalence= 0.3
)
mstage.dca <- decision_curve(OS ~M.stage,data =dca.dat,family = binomial(link ='logit'),
                             thresholds = seq(0,1, by = 0.01),
                             confidence.intervals= 0.95,
                             study.design = 'case-control',
                             population.prevalence= 0.3
)
colnames(dca.dat)
AP000696.2.dca <- decision_curve(OS ~AP000696.2,data =dca.dat,family = binomial(link ='logit'),
                             thresholds = seq(0,1, by = 0.01),
                             confidence.intervals= 0.95,
                             study.design = 'case-control',
                             population.prevalence= 0.3
)
PTCSC3.dca <- decision_curve(OS ~PTCSC3,data =dca.dat,family = binomial(link ='logit'),
                                 thresholds = seq(0,1, by = 0.01),
                                 confidence.intervals= 0.95,
                                 study.design = 'case-control',
                                 population.prevalence= 0.3
)


dca.list <- list(age.dca,grade.dca,stage.dca,tstage.dca,nstage.dca,mstage.dca,AP000696.2.dca,PTCSC3.dca,complex)

pdf("05.DCA.pdf", width = 6, height = 6)
plot_decision_curve(dca.list,
                    curve.names=c('Age','Grade','Stage','T.stage','N.stage','M.stage','AP000696.2','PTCSC3','Nomogram'),
                    cost.benefit.axis =FALSE,col= c("#FDC086","#7FC97F","#BEAED4","#74CCBE","#ED7474","red",'pink','orange','blue'),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()
png("05.DCA.png", width = 500, height = 500)
plot_decision_curve(dca.list,
                    curve.names=c('Age','Grade','Stage','T.stage','N.stage','M.stage','AP000696.2','PTCSC3','Nomogram'),
                    cost.benefit.axis =FALSE,col= c("#FDC086","#7FC97F","#BEAED4","#74CCBE","#ED7474","red",'pink','orange','blue'),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()

##临床影响曲线CIC-------
pdf(file = '06.CIC.pdf', width = 6, height = 5)
plot_clinical_impact(complex,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T,ylim=c(0,1000),
                     legend.position= "topright")
dev.off()
#红色曲线（Numberhigh risk）表示，在各个阈概率下，被simple或complex模型划分为阳性（高风险）的人数；
#蓝色曲线（Number high risk with outcome）为各个阈概率下真阳性的人数。意义一目了然。
png(file = '06.CIC.png', width = 400, height = 350)
plot_clinical_impact(complex,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T,ylim=c(0,1000),
                     legend.position= "topright")
dev.off()
