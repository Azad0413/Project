rm(list = ls())
library(glmnet)
library(survival)
library(survminer)
library(tidyverse)
library(magrittr)
library(lance)

df.tmp = read.delim2("clean.data/fpkm.tsv", row.names = 1) %>% lc.tableToNum %>% add(1) %>% log2 %>% as.data.frame
df.surv = read.delim2("clean.data/survival.tsv")
gene.inter = read.delim2("01.DEG/05.Invasion.DEG.xls")

smp = intersect(df.surv$sample, colnames(df.tmp))
smp = smp[!str_ends(smp, "1..")]
df.tmp = df.tmp[smp]
df.surv = subset(df.surv, sample %in% smp)
df.tmp = df.tmp[as.character(gene.inter$GeneID),]
rownames(df.tmp) = gene.inter$Gene

res.univar = apply(df.tmp, 1, function(x){
  res.sub = coxph(Surv(time = df.surv$OS.time, event = df.surv$OS) ~ x) %>% summary
  p.value = signif(res.sub$wald["pvalue"], digits=2)
  hr = signif(res.sub$coef[2], digits=3)
  hr.low = signif(res.sub$conf.int[,"lower .95"],digits = 3)
  hr.high = signif(res.sub$conf.int[,"upper .95"],digits = 3)
  return(c(hr, p.value, hr.low, hr.high))
}) %>% t %>% as.data.frame

colnames(res.univar) = c("HR","p.value","CI.low","CI.high")
res.univar$EntrezID = gene.inter$GeneID
res.univar$Gene = rownames(res.univar)
res.univar = res.univar[c(5,6,1:4)]
write.table(res.univar, "02.RiskScore/01.UnivarCox.xls", sep = "\t", col.names = T, row.names = F, quote = F)

res.univar.plot = subset(res.univar, p.value < 0.01)
res.univar.plot = res.univar.plot[order(res.univar.plot$HR, decreasing = T),]
res.univar.plot$Gene = factor(res.univar.plot$Gene, levels = res.univar.plot$Gene)
res.univar.plot$plabel = res.univar.plot$p.value %>% signif(3)
res.univar.plot$hlabel = paste0(res.univar.plot$HR,"(", res.univar.plot$CI.low, "-", res.univar.plot$CI.high, ")")

library(aplot)
p = ggplot(data = res.univar.plot) + 
  geom_vline(xintercept = 1, color = "grey60", linetype = 1) + 
  geom_segment(aes(x = CI.low, xend = CI.high, y = Gene, yend = Gene), color = "darkblue", size = 1) +
  geom_point(aes(y = Gene, x = HR), size = 3, shape = 15, color = "red") +
  theme_minimal() + 
  ylab(NULL) + xlab("Harzard Ratio") + 
  scale_x_continuous(limits = c(0.1,2)) +
  scale_color_manual(values = c("green","red"), guide = guide_none()) +
  theme(axis.line.x = element_line(size = 1), axis.title = element_text(size = 14), 
        axis.text.x = element_text(size = 12), axis.text.y = element_blank(),
        panel.grid = element_blank(), plot.background = element_rect(fill = "white", colour = "transparent"))
dt = ggplot(data = res.univar.plot) + 
  geom_text(aes(y = Gene, x = 1, label = plabel), hjust = 0.5, fontface = "bold") +
  theme_minimal() + ggtitle("p.value") + xlab(NULL) + ylab(NULL) +
  theme(axis.text = element_blank(), plot.background = element_rect(fill = "white", colour = "transparent"),
        plot.title = element_text(hjust = 0.5), line = element_blank())
dt2 = ggplot(data = res.univar.plot) + 
  geom_text(aes(y = Gene, x = 1, label = hlabel), hjust = 0.5, fontface = "bold") +
  theme_minimal() + ggtitle("Harzad Ratio") + xlab(NULL) + ylab(NULL) +
  theme(axis.text = element_blank(), plot.background = element_rect(fill = "white", colour = "transparent"),
        plot.title = element_text(hjust = 0.5), line = element_blank())
dt3 = ggplot(data = res.univar.plot) + 
  geom_text(aes(y = Gene, x = 1, label = Gene), hjust = 0, fontface = "bold") +
  theme_minimal() + xlab(NULL) + ylab(NULL) +
  theme(axis.text= element_blank(), plot.background = element_rect(fill = "white", colour = "transparent"),
        plot.title = element_text(hjust = 0.5), line = element_blank())

p.final = insert_left(p, dt2, width = 0.35) %>% insert_left(dt, width = 0.2) %>% insert_left(dt3, width = 0.4)
p.final
ggsave("02.RiskScore/02.UnivarCox.png", p.final, width = 10, height = 8, units = "in", dpi = 300, bg = "white")
ggsave("02.RiskScore/02.UnivarCox.pdf", p.final, width = 10, height = 8, units = "in", dpi = 300, bg = "white")

gs = res.univar.plot$Gene %>% levels
df.tmp = df.tmp[gs,]
df.fpkm.sub = df.tmp
identical(colnames(df.fpkm.sub),df.surv$sample)

set.seed(1234567)
res.lasso = cv.glmnet(t(df.fpkm.sub), Surv(time = df.surv$OS.time, event = df.surv$OS), family = "cox",
                      nfolds = 10, alpha = 1, type.measure = "C")
plot(res.lasso)
coef(res.lasso, s = res.lasso$lambda.min)
ggsave("02.RiskScore/03.Lasso.CV.png", plot(res.lasso), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("02.RiskScore/03.Lasso.CV.pdf", plot(res.lasso), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("02.RiskScore/04.Lasso.Coef.png", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 8, height = 7, dpi = 300, units = "in", bg = "white")
ggsave("02.RiskScore/04.Lasso.Coef.pdf", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 8, height = 7, dpi = 300, units = "in", bg = "white")

df.coef = coef(res.lasso, s = res.lasso$lambda.min)
df.coef = cbind(gene = rownames(df.coef), coefficient = df.coef[,1])
df.coef = lc.tableToNum(df.coef) %>% as.data.frame
df.coef = subset(df.coef, coefficient != 0) %>% as.data.frame %>% lc.tableToNum()
write.table(df.coef, "02.RiskScore/05.Coefficients.xls", sep = "\t", quote = F, col.names = T, row.names = F)

identical(df.surv$sample, colnames(df.fpkm.sub))
df.surv$risk = predict(res.lasso, s = res.lasso$lambda.min, newx = t(df.fpkm.sub), type = "link")[,1]
df.surv$risk.group = ifelse(df.surv$risk > median(df.surv$risk), "High", "Low")
df.surv = df.surv[order(df.surv$risk, decreasing = F),]
df.surv$rn = 1:nrow(df.surv)

ggplot(df.surv, aes(x = rn, y = risk, color = risk.group)) + 
  geom_hline(yintercept = median(df.surv$risk), linetype = 2) +
  geom_vline(xintercept = sum(df.surv$risk.group=="Low"), linetype = 2) +
  geom_point() +
  scale_color_manual(values = c("orange","lightgreen"), labels = c("High risk","Low risk")) +
  xlab("Patients (increasing riskscore)") +
  ylab("Risk Score") +
  theme_classic2() + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), legend.position = c(0.1,0.85),
        axis.line = element_line(size = 0),
        panel.border = element_rect(color = "black", fill = "transparent", size = 1.5))
ggsave("02.RiskScore/06.RiskDistribution.png", width = 10, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("02.RiskScore/06.RiskDistribution.pdf", width = 10, height = 4, units = "in", dpi = 300, bg = "white")

df.surv$surv = ifelse(df.surv$OS == 0, "Alive", "Dead")
ggplot(df.surv, aes(x = rn, y = OS.time, color = surv)) + 
  geom_vline(xintercept = sum(df.surv$risk.group=="Low"), linetype = 2) +
  geom_point() +
  scale_color_manual(values = c("lightgreen","orange")) +
  xlab("Patients (increasing riskscore)") +
  ylab("Survival Time (Days)") +
  theme_classic2() + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), legend.position = c(0.1,0.85),
        axis.line = element_line(size = 0), legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(color = "black", fill = "transparent", size = 1.5))
ggsave("02.RiskScore/07.SurvDistribution.png", width = 10, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("02.RiskScore/07.SurvDistribution.pdf", width = 10, height = 4, units = "in", dpi = 300, bg = "white")

surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk.group, data = df.surv)
pv = survminer::surv_pvalue(surv.fit, data = df.surv)[2] %>% signif(3)
df.surv$risk.group = factor(df.surv$risk.group, levels = c("Low", "High"))
surv.diff = survdiff(Surv(time = OS.time,event = OS) ~ risk.group, data = df.surv)
HR = ((surv.diff$obs[2]/surv.diff$exp[2])/(surv.diff$obs[1]/surv.diff$exp[1])) %>% signif(3)
up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/surv.diff$exp[2]+1/surv.diff$exp[1])) %>% signif(3)
low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/surv.diff$exp[2]+1/surv.diff$exp[1])) %>% signif(3)
output.text = paste0("HR ", HR, "(95% CI ", low95,"-", up95, "); logrank p=", pv)

p.surv = ggsurvplot(surv.fit, data = df.surv, conf.int = F, risk.table = TRUE, break.x.by=1200,
                    pval = T, legend.labs = c("High Risk", "Low Risk")) + xlab("Time (Day)")
p.surv
ggsave("02.RiskScore/08.KM.png", print(p.surv), width = 10, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("02.RiskScore/08.KM.pdf", print(p.surv), width = 10, height = 7, units = "in", dpi = 300, bg = "white")

library(pROC)
df.train.1year = subset(df.surv, OS == 1 | OS.time >= 365*1)
df.train.2year = subset(df.surv, OS == 1 | OS.time >= 365*2)
df.train.3year = subset(df.surv, OS == 1 | OS.time >= 365*3)
df.train.5year = subset(df.surv, OS == 1 | OS.time >= 365*5)

df.train.1year$event = ifelse(df.train.1year$OS.time >= 365*1, 0, 1)
df.train.2year$event = ifelse(df.train.2year$OS.time >= 365*2, 0, 1)
df.train.3year$event = ifelse(df.train.3year$OS.time >= 365*3, 0, 1)
df.train.5year$event = ifelse(df.train.5year$OS.time >= 365*5, 0, 1)

roc.1year = roc(df.train.1year$event, df.train.1year$risk)
roc.2year = roc(df.train.2year$event, df.train.2year$risk)
roc.3year = roc(df.train.3year$event, df.train.3year$risk)
roc.5year = roc(df.train.5year$event, df.train.5year$risk)

auc.1year = auc(roc.1year) %>% round(2) %>% paste0("1 year: ", .)
auc.2year = auc(roc.2year) %>% round(2) %>% paste0("2 year: ", .)
auc.3year = auc(roc.3year) %>% round(2) %>% paste0("3 year: ", .)
auc.5year = auc(roc.5year) %>% round(2) %>% paste0("5 year: ", .)
auc.annot = paste("Area Under Curve", auc.1year, auc.3year, auc.5year, sep = "\n")

ggroc(list("1 year" = roc.1year, "3 year" = roc.3year,"5 year" = roc.5year), legacy.axes = T, size = 0.8) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey30") +
  geom_text(aes(x = 0.8, y = 0.25, label = auc.annot), color = "grey20", size = 5) +
  theme_pubclean(base_size = 16) + 
  guides(color = guide_legend("Survival Time")) +
  theme(panel.grid.major.x = element_line(color = "grey", linetype = 3), 
        panel.border = element_rect(color = "black", fill = "transparent"), 
        aspect.ratio = 1)
ggsave("02.RiskScore/09.ROC.png", width = 7, height = 8, units = "in", dpi = 600)
ggsave("02.RiskScore/09.ROC.pdf", width = 7, height = 8, units = "in", dpi = 600)

write.table(df.surv, "full.table.tsv", quote = F, row.names = F, col.names = T, sep = "\t")





