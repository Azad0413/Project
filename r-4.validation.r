rm(list = ls())
library(GEOquery)
library(org.Hs.eg.db)
library(survival)
library(survminer)

options(timeout = 10^7)
Sys.setenv("VROOM_CONNECTION_SIZE" = 10^7)
gs = getGEO("GSE19234", destdir = ".", getGPL = F)
df.exp = gs[[1]] %>% exprs %>% as.data.frame
df.pheno = gs[[1]] %>% pData %>% as.data.frame
boxplot(df.exp[1:5])
df.exp = log2(df.exp+1)
df.exp = limma::normalizeQuantiles(df.exp)
boxplot(df.exp[1:5])

df.annot = getGEO("GPL570", AnnotGPL = T, destdir = ".") %>% Table
df.annot = subset(df.annot, !str_detect(`Gene ID`, "/"))
df.exp$gene = df.annot$`Gene ID`[match(rownames(df.exp), as.character(df.annot$ID))]
df.exp = subset(df.exp, gene != "")
df.exp = na.omit(df.exp)
df.exp = aggregate(.~gene, FUN = max, data = df.exp)
rownames(df.exp) = as.character(df.exp$gene)
write.table(df.exp, "clean.data/GSE19234.expression.tsv", sep = "\t", row.names = F, quote = F)
write.table(df.pheno, "clean.data/GSE19234.pheno.tsv", sep = "\t", row.names = F, quote = F)

df.exp = df.exp[-1]

df.coef = read.delim2("02.RiskScore/05.Coefficients.xls")
coef = as.numeric(df.coef$coefficient)
genes = AnnotationDbi::select(org.Hs.eg.db, df.coef$gene, "ENTREZID", "SYMBOL")
df.exp = df.exp[genes$ENTREZID,]
df.exp = df.exp[df.pheno$geo_accession]

df.pheno$risk = apply(df.exp, 2, function(x){x %*% coef})
df.pheno$OS.time = df.pheno$`days since initial diagnosis:ch1` %>% as.numeric()
df.pheno$OS = (df.pheno$`staus dead or alive:ch1` == "1") %>% as.numeric()
df.pheno = df.pheno[order(df.pheno$risk, decreasing = F),]
df.pheno = subset(df.pheno, !is.na(OS.time))
df.pheno$risk.group = ifelse(df.pheno$risk >= median(df.pheno$risk), "High", "Low")
df.pheno$rn = 1:nrow(df.pheno)

surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk.group, data = df.pheno)
p.surv = ggsurvplot(surv.fit, data = df.pheno, conf.int = F, risk.table = TRUE, break.x.by=800, pval = T, 
                    legend.labs = c("High Risk", "Low Risk")) + xlab("Time (Month)")
p.surv

ggplot(df.pheno, aes(x = rn, y = risk, color = risk.group)) + 
  geom_hline(yintercept = median(df.pheno$risk), linetype = 2) +
  geom_vline(xintercept = sum(df.pheno$risk.group=="Low"), linetype = 2) +
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
ggsave("03.Validation/01.RiskDistribution.png", width = 10, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Validation/01.RiskDistribution.pdf", width = 10, height = 4, units = "in", dpi = 300, bg = "white")

df.pheno$surv = ifelse(df.pheno$OS == 0, "Alive", "Dead")
ggplot(df.pheno, aes(x = rn, y = OS.time, color = surv)) + 
  geom_vline(xintercept = sum(df.pheno$risk.group=="Low"), linetype = 2) +
  geom_point() +
  scale_color_manual(values = c("lightgreen","orange")) +
  xlab("Patients (increasing riskscore)") +
  ylab("Survival Time (Months)") +
  theme_classic2() + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), legend.position = c(0.1,0.85),
        axis.line = element_line(size = 0), legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(color = "black", fill = "transparent", size = 1.5))
ggsave("03.Validation/03.SurvDistribution.png", width = 10, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Validation/03.SurvDistribution.pdf", width = 10, height = 4, units = "in", dpi = 300, bg = "white")

surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk.group, data = df.pheno)
p.surv = ggsurvplot(surv.fit, data = df.pheno, conf.int = F, risk.table = TRUE, break.x.by=800, pval = T, 
                    legend.labs = c("High Risk", "Low Risk")) + xlab("Time (Day)")
p.surv
ggsave("03.Validation/03.KM.png", print(p.surv), width = 10, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("03.Validation/03.KM.pdf", print(p.surv), width = 10, height = 7, units = "in", dpi = 300, bg = "white")

library(pROC)
df.surv = df.pheno
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
auc.annot
ggroc(list("1 year" = roc.1year,"3 year" = roc.3year, "5 year" = roc.5year), legacy.axes = T, size = 0.8) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey30") +
  geom_text(aes(x = 0.8, y = 0.25, label = auc.annot), color = "grey20", size = 5) +
  theme_pubclean(base_size = 16) + 
  guides(color = guide_legend("Survival Time")) +
  theme(panel.grid.major.x = element_line(color = "grey", linetype = 3), 
        panel.border = element_rect(color = "black", fill = "transparent"), 
        aspect.ratio = 1)
ggsave("03.Validation/04.ROC.png", width = 7, height = 8, units = "in", dpi = 600)
ggsave("03.Validation/04.ROC.pdf", width = 7, height = 8, units = "in", dpi = 600)

rm(list = ls())
library(org.Hs.eg.db)
library(survival)
library(survminer)

df.exp = read.delim2("uvm/fpkm.tsv", row.names = 1) %>% lc.tableToNum %>% add(1) %>% log2 %>% as.data.frame
df.coef = read.delim2("02.RiskScore/05.Coefficients.xls")
coef = as.numeric(df.coef$coefficient)
genes = AnnotationDbi::select(org.Hs.eg.db, df.coef$gene, "ENTREZID", "SYMBOL")
df.exp = df.exp[genes$ENTREZID,]
df.pheno = read.delim2("uvm/survival.tsv")
df.exp = df.exp[df.pheno$sample]

df.pheno$risk = apply(df.exp, 2, function(x){x %*% coef})
df.pheno = df.pheno[order(df.pheno$risk, decreasing = F),]
df.pheno$risk.group = ifelse(df.pheno$risk >= median(df.pheno$risk), "High", "Low")
df.pheno$rn = 1:nrow(df.pheno)

surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk.group, data = df.pheno)
p.surv = ggsurvplot(surv.fit, data = df.pheno, conf.int = F, risk.table = TRUE, break.x.by=800, pval = T, 
                    legend.labs = c("High Risk", "Low Risk")) + xlab("Time (Month)")
p.surv

ggplot(df.pheno, aes(x = rn, y = risk, color = risk.group)) + 
  geom_hline(yintercept = median(df.pheno$risk), linetype = 2) +
  geom_vline(xintercept = sum(df.pheno$risk.group=="Low"), linetype = 2) +
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
ggsave("03.Validation/01.RiskDistribution.png", width = 10, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Validation/01.RiskDistribution.pdf", width = 10, height = 4, units = "in", dpi = 300, bg = "white")

df.pheno$surv = ifelse(df.pheno$OS == 0, "Alive", "Dead")
ggplot(df.pheno, aes(x = rn, y = OS.time, color = surv)) + 
  geom_vline(xintercept = sum(df.pheno$risk.group=="Low"), linetype = 2) +
  geom_point() +
  scale_color_manual(values = c("lightgreen","orange")) +
  xlab("Patients (increasing riskscore)") +
  ylab("Survival Time (Months)") +
  theme_classic2() + 
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), legend.position = c(0.1,0.85),
        axis.line = element_line(size = 0), legend.background = element_rect(fill = "transparent"),
        panel.border = element_rect(color = "black", fill = "transparent", size = 1.5))
ggsave("03.Validation/03.SurvDistribution.png", width = 10, height = 4, units = "in", dpi = 300, bg = "white")
ggsave("03.Validation/03.SurvDistribution.pdf", width = 10, height = 4, units = "in", dpi = 300, bg = "white")

surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk.group, data = df.pheno)
p.surv = ggsurvplot(surv.fit, data = df.pheno, conf.int = F, risk.table = TRUE, break.x.by=800, pval = T, 
                    legend.labs = c("High Risk", "Low Risk")) + xlab("Time (Day)")
p.surv
ggsave("03.Validation/03.KM.png", print(p.surv), width = 10, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("03.Validation/03.KM.pdf", print(p.surv), width = 10, height = 7, units = "in", dpi = 300, bg = "white")

library(pROC)
df.surv = df.pheno
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
auc.annot
ggroc(list("1 year" = roc.1year,"3 year" = roc.3year, "5 year" = roc.5year), legacy.axes = T, size = 0.8) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "grey30") +
  geom_text(aes(x = 0.8, y = 0.25, label = auc.annot), color = "grey20", size = 5) +
  theme_pubclean(base_size = 16) + 
  guides(color = guide_legend("Survival Time")) +
  theme(panel.grid.major.x = element_line(color = "grey", linetype = 3), 
        panel.border = element_rect(color = "black", fill = "transparent"), 
        aspect.ratio = 1)
ggsave("03.Validation/04.ROC.png", width = 7, height = 8, units = "in", dpi = 600)
ggsave("03.Validation/04.ROC.pdf", width = 7, height = 8, units = "in", dpi = 600)