rm(list = ls())
library(glmnet)
library(survival)
library(survminer)
library(tidyverse)
library(magrittr)
library(lance)
library(caret)

df.tmp = read.delim2("/data/nas1/yangly/Project/JNZK204/clean.data/fpkm.tsv", row.names = 1) %>% lc.tableToNum %>% add(1) %>% log2 %>% as.data.frame
df.surv = read.delim2("/data/nas1/yangly/Project/JNZK204/clean.data/survival.tsv")
gene.inter = read.delim2("../01.DEG/05.Selected.DEG.xls")

df.tmp = df.tmp[as.character(gene.inter$EntrezID),]
rownames(df.tmp) = gene.inter$Gene

smp = intersect(df.surv$sample, colnames(df.tmp))
smp = smp[!str_ends(smp, "1..")]
df.surv = subset(df.surv, sample %in% smp)
set.seed(123); train.index = createDataPartition(df.surv$OS, 1, 0.7)[[1]]
smp.train = smp[train.index]
smp.test = setdiff(smp, smp.train)
df.surv.train = subset(df.surv, sample %in% smp.train)
df.surv.test = subset(df.surv, sample %in% smp.test)
df.tmp.train = df.tmp[df.surv.train$sample]
df.tmp.test = df.tmp[df.surv.test$sample]

res.univar = apply(df.tmp.train, 1, function(x){
  res.sub = coxph(Surv(time = df.surv.train$OS.time, event = df.surv.train$OS) ~ x) %>% summary
  p.value = signif(res.sub$wald["pvalue"], digits=2)
  hr = signif(res.sub$coef[2], digits=3)
  hr.low = signif(res.sub$conf.int[,"lower .95"],digits = 3)
  hr.high = signif(res.sub$conf.int[,"upper .95"],digits = 3)
  return(c(hr, p.value, hr.low, hr.high))
}) %>% t %>% as.data.frame

colnames(res.univar) = c("HR","p.value","CI.low","CI.high")
res.univar$EntrezID = gene.inter$EntrezID
res.univar$Gene = rownames(res.univar)
res.univar = res.univar[c(5,6,1:4)]
write.table(res.univar, "../02.RiskScore/01.UnivarCox.xls", sep = "\t", col.names = T, row.names = F, quote = F)

res.univar.plot = subset(res.univar, p.value < 0.05)
res.univar.plot = res.univar.plot[order(res.univar.plot$HR, decreasing = T),]
res.univar.plot$Gene = factor(res.univar.plot$Gene, levels = res.univar.plot$Gene)
res.univar.plot$plabel = res.univar.plot$p.value %>% signif(2)
res.univar.plot$hlabel = paste0(res.univar.plot$HR,"(", res.univar.plot$CI.low, "-", res.univar.plot$CI.high, ")")

library(aplot)
p = ggplot(data = res.univar.plot) + 
  geom_vline(xintercept = 1, color = "grey60", linetype = 1) + 
  geom_segment(aes(x = CI.low, xend = CI.high, y = Gene, yend = Gene), color = "darkblue", size = 1) +
  geom_point(aes(y = Gene, x = HR), size = 3, shape = 15, color = "red") +
  theme_minimal() + 
  ylab(NULL) + xlab("Harzard Ratio") + 
  #scale_x_continuous(limits = c(0,8)) +
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
ggsave("02.RiskScore/02.UnivarCox.png", p.final, width = 10, height = 10, units = "in", dpi = 300, bg = "white")
ggsave("02.RiskScore/02.UnivarCox.pdf", p.final, width = 10, height = 10, units = "in", dpi = 300, bg = "white")

gs = res.univar.plot$Gene %>% levels()
df.fpkm.sub = df.tmp.train[gs,]
identical(colnames(df.fpkm.sub),df.surv.train$sample)

#set.seed(123456)
set.seed(77)
res.lasso = cv.glmnet(t(df.fpkm.sub), Surv(time = df.surv.train$OS.time, event = df.surv.train$OS), family = "cox",
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
write.table(df.coef, "../02.RiskScore/05.Coefficients.xls", sep = "\t", quote = F, col.names = T, row.names = F)

identical(df.surv.train$sample, colnames(df.fpkm.sub))
df.surv.train$risk = predict(res.lasso, s = res.lasso$lambda.min, newx = t(df.fpkm.sub), type = "link")[,1]
df.surv.train$risk.group = ifelse(df.surv.train$risk > median(df.surv.train$risk), "High", "Low")
df.surv.train = df.surv.train[order(df.surv.train$risk, decreasing = F),]
df.surv.train$rn = 1:nrow(df.surv.train)
write.table(df.surv.train,file = 'risk.xls',sep = '\t',row.names = F,quote = F)

ggplot(df.surv.train, aes(x = rn, y = risk, color = risk.group)) + 
  geom_hline(yintercept = median(df.surv.train$risk), linetype = 2) +
  geom_vline(xintercept = sum(df.surv.train$risk.group=="Low"), linetype = 2) +
  geom_point() +
  scale_color_manual(values = c("orange","darkgreen"), labels = c("High risk","Low risk")) +
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

df.surv.train$surv = ifelse(df.surv.train$OS == 0, "Alive", "Dead")
ggplot(df.surv.train, aes(x = rn, y = OS.time, color = surv)) + 
  geom_vline(xintercept = sum(df.surv.train$risk.group=="Low"), linetype = 2) +
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

surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk.group, data = df.surv.train)
p.surv = ggsurvplot(surv.fit, data = df.surv.train, conf.int = F, risk.table = TRUE, break.x.by=800,
                    pval = T, legend.labs = c("High Risk", "Low Risk")) + xlab("Time (Day)")
p.surv
ggsave("02.RiskScore/08.KM.png", print(p.surv), width = 10, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("02.RiskScore/08.KM.pdf", print(p.surv), width = 10, height = 7, units = "in", dpi = 300, bg = "white")

library(pROC)
df.train.7year = subset(df.surv.train, OS == 1 | OS.time >= 365*7)
df.train.3year = subset(df.surv.train, OS == 1 | OS.time >= 365*3)
df.train.5year = subset(df.surv.train, OS == 1 | OS.time >= 365*5)
df.train.1year = subset(df.surv.train, OS == 1 | OS.time >= 365*1)

df.train.7year$event = ifelse(df.train.7year$OS.time >= 365*7, 0, 1)
df.train.3year$event = ifelse(df.train.3year$OS.time >= 365*3, 0, 1)
df.train.5year$event = ifelse(df.train.5year$OS.time >= 365*5, 0, 1)
df.train.1year$event = ifelse(df.train.1year$OS.time >= 365*1, 0, 1)

roc.7year = roc(df.train.7year$event, df.train.7year$risk)
roc.3year = roc(df.train.3year$event, df.train.3year$risk)
roc.5year = roc(df.train.5year$event, df.train.5year$risk)
roc.1year = roc(df.train.1year$event, df.train.1year$risk)

auc.7year = auc(roc.7year) %>% round(3) %>% paste0("7 year: ", .)
auc.3year = auc(roc.3year) %>% round(3) %>% paste0("3 year: ", .)
auc.5year = auc(roc.5year) %>% round(3) %>% paste0("5 year: ", .)
auc.1year = auc(roc.1year) %>% round(3) %>% paste0("1 year: ", .)
auc.annot = paste("Area Under Curve", auc.1year, auc.3year, auc.5year, auc.7year, sep = "\n")

ggroc(list("1 year" = roc.1year, "3 year" = roc.3year, "5 year" = roc.5year, "7 year" = roc.7year), legacy.axes = T, size = 0.8) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey30") +
  geom_text(aes(x = 0.8, y = 0.25, label = auc.annot), color = "grey20", size = 5) +
  theme_pubclean(base_size = 16) + 
  guides(color = guide_legend("Survival Time")) +
  theme(panel.grid.major.x = element_line(color = "grey", linetype = 3), 
        panel.border = element_rect(color = "black", fill = "transparent"), 
        aspect.ratio = 1)
ggsave("02.RiskScore/09.ROC.png", width = 7, height = 8, units = "in", dpi = 600)
ggsave("02.RiskScore/09.ROC.pdf", width = 7, height = 8, units = "in", dpi = 600)

identical(df.surv.test$sample, colnames(df.tmp.test))
df.tmp.test = df.tmp.test[gs,]
df.surv.test$risk = predict(res.lasso, s = res.lasso$lambda.min, newx = t(df.tmp.test), type = "link")[,1]
df.surv.test$risk.group = ifelse(df.surv.test$risk > median(df.surv.test$risk), "High", "Low")
df.surv.test = df.surv.test[order(df.surv.test$risk, decreasing = F),]
df.surv.test$rn = 1:nrow(df.surv.test)

ggplot(df.surv.test, aes(x = rn, y = risk, color = risk.group)) + 
  geom_hline(yintercept = median(df.surv.test$risk), linetype = 2) +
  geom_vline(xintercept = sum(df.surv.test$risk.group=="Low"), linetype = 2) +
  geom_point() +
  scale_color_manual(values = c("orange","darkgreen"), labels = c("High risk","Low risk")) +
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
ggsave("03.Validation//01.RiskDistribution.png", width = 10, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Validation//01.RiskDistribution.pdf", width = 10, height = 4, units = "in", dpi = 300, bg = "white")

df.surv.test$surv = ifelse(df.surv.test$OS == 0, "Alive", "Dead")
ggplot(df.surv.test, aes(x = rn, y = OS.time, color = surv)) + 
  geom_vline(xintercept = sum(df.surv.test$risk.group=="Low"), linetype = 2) +
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
ggsave("03.Validation//02.SurvDistribution.png", width = 10, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Validation//02.SurvDistribution.pdf", width = 10, height = 4, units = "in", dpi = 300, bg = "white")

surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk.group, data = df.surv.test)
p.surv = ggsurvplot(surv.fit, data = df.surv.test, conf.int = F, risk.table = TRUE, break.x.by=800,
                    pval = T, legend.labs = c("High Risk", "Low Risk")) + xlab("Time (Day)")
p.surv
ggsave("03.Validation//04.KM.png", print(p.surv), width = 10, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("03.Validation//04.KM.pdf", print(p.surv), width = 10, height = 7, units = "in", dpi = 300, bg = "white")

library(pROC)
df.test.7year = subset(df.surv.test, OS == 1 | OS.time >= 365*7)
df.test.3year = subset(df.surv.test, OS == 1 | OS.time >= 365*3)
df.test.5year = subset(df.surv.test, OS == 1 | OS.time >= 365*5)
df.test.1year = subset(df.surv.test, OS == 1 | OS.time >= 365*1)

df.test.7year$event = ifelse(df.test.7year$OS.time >= 365*7, 0, 1)
df.test.3year$event = ifelse(df.test.3year$OS.time >= 365*3, 0, 1)
df.test.5year$event = ifelse(df.test.5year$OS.time >= 365*5, 0, 1)
df.test.1year$event = ifelse(df.test.1year$OS.time >= 365*1, 0, 1)

roc.7year = roc(df.test.7year$event, df.test.7year$risk)
roc.3year = roc(df.test.3year$event, df.test.3year$risk)
roc.5year = roc(df.test.5year$event, df.test.5year$risk)
roc.1year = roc(df.test.1year$event, df.test.1year$risk)

auc.1year = auc(roc.1year) %>% round(3) %>% paste0("1 year: ", .)
auc.3year = auc(roc.3year) %>% round(3) %>% paste0("3 year: ", .)
auc.5year = auc(roc.5year) %>% round(3) %>% paste0("5 year: ", .)
auc.7year = auc(roc.7year) %>% round(3) %>% paste0("7 year: ", .)
auc.annot = paste("Area Under Curve", auc.1year, auc.3year, auc.5year, auc.7year, sep = "\n")

ggroc(list("1 year" = roc.1year, "3 year" = roc.3year,"5 year" = roc.5year, "7 year" = roc.7year), legacy.axes = T, size = 0.8) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey30") +
  geom_text(aes(x = 0.8, y = 0.25, label = auc.annot), color = "grey20", size = 5) +
  theme_pubclean(base_size = 16) + 
  guides(color = guide_legend("Survival Time")) +
  theme(panel.grid.major.x = element_line(color = "grey", linetype = 3), 
        panel.border = element_rect(color = "black", fill = "transparent"), 
        aspect.ratio = 1)
ggsave("03.Validation//04.ROC.png", width = 7, height = 8, units = "in", dpi = 600)
ggsave("03.Validation//04.ROC.pdf", width = 7, height = 8, units = "in", dpi = 600)

df.surv = rbind(df.surv.train, df.surv.test)
df.surv$risk.group = ifelse(df.surv$risk > median(df.surv$risk), "High", "Low")
write.table(df.surv, "full.table.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

rm(list = ls())
df.coef = read.delim2("02.RiskScore/05.Coefficients.xls")
df.coef$id = AnnotationDbi::select(org.Hs.eg.db, df.coef$gene, "ENTREZID", "SYMBOL")$ENTREZID
df.exp = read.delim2("/data/nas1/yangly/Project/JNZK204/clean.data/fpkm.tsv", row.names = 1) %>% lc.tableToNum()
df.exp = df.exp[as.character(df.coef$id),]
rownames(df.exp) = df.coef$gene
df.deg = read.delim2("01.DEG/03.DEG.xls", row.names = 1)
df.deg = df.deg[as.character(df.coef$id),]
rownames(df.deg) = df.deg$Gene

grps = ifelse(colnames(df.exp) %>% str_ends(".1.."),"normal","tumor")
df.meta = data.frame(sample = colnames(df.exp), Group = grps)
rownames(df.meta) = df.meta$sample
df.meta$Group = factor(df.meta$Group, levels = c("normal","tumor"), labels = c("Normal", "Tumor"))
df.meta = df.meta[order(df.meta$Group, decreasing = F),]
df.exp = df.exp[df.meta$sample]

p = pheatmap::pheatmap(log2(df.exp+1), annotation_col = df.meta[-1], cluster_rows = T, cluster_cols = F, scale = "row", 
                       color = colorRampPalette(c("green","green","green","black","red","red","red"))(100),
                       annotation_row = df.deg['Direction'], annotation_names_col = F, annotation_names_row = F,
                       annotation_colors = list(Group = c(Normal = "darkgreen", Tumor = "darkorange"),
                                                Direction = c(Up = "red", Down = "blue")),
                       show_colnames = F)
ggsave("02.RiskScore/10.Model.Gene.png", p, width = 12, height = 4.5, units = "in", dpi = 300, bg = "white")
ggsave("02.RiskScore/10.Model.Gene.pdf", p, width = 12, height = 4.5, units = "in", dpi = 300, bg = "white")

df.exp = log2(df.exp + 1) %>% t %>% as.data.frame()

library(Rtsne)
df.meta.sub = subset(df.meta, Group == "Tumor")
df.exp.sub = df.exp[df.meta.sub$sample,]
set.seed(1234567);res.tsne = Rtsne(df.exp.sub, pca = TRUE, pca_scale = TRUE, perplexity = 45)
ggplot(as.data.frame(res.tsne$Y), aes(x = V1, y = V2)) + geom_point()














