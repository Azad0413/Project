rm(list = ls())
library(lance)
library(magrittr)
library(stringr)
library(ggplot2)
library(org.Hs.eg.db)
library(survival)
library(rms)
library(ggpubr)

df.surv = read.delim2("full.table.tsv")
df.pheno = read.delim2("clean.data/pheno.tsv")
df.pheno$sample = make.names(df.pheno$submitter_id.samples)
df.pheno = df.pheno[c("sample","age_at_initial_pathologic_diagnosis",
                      "pathologic_M","pathologic_N","pathologic_T",
                      "gender.demographic")]
colnames(df.pheno) = c("sample","Age","M_Stage","N_Stage","T_Stage","Gender")
df.pheno = lc.tableToNum(df.pheno)
df.full = merge(df.surv, df.pheno, by = "sample")
df.full = lc.tableToNum(df.full)

res.risk = coxph(Surv(time = OS.time, event = OS) ~ risk, data = df.full) %>% summary
res.risk = c(res.risk$conf.int[-2], res.risk$coefficients[5])

res.age = coxph(Surv(time = OS.time, event = OS) ~ Age, data = df.full) %>% summary
res.age = c(res.age$conf.int[-2], res.age$coefficients[5])

res.gender = coxph(Surv(time = OS.time, event = OS) ~ Gender, data = df.full) %>% summary
res.gender = c(res.gender$conf.int[-2], res.gender$coefficients[5])

df.full$T_Stage = str_match(df.full$T_Stage, "T[1-4]")[,1]
df.full$T_Stage = c("T1"="T1-T2","T2"="T1-T2","T3"="T3-T4","T4"="T3-T4")[df.full$T_Stage]
res.t = coxph(Surv(time = OS.time, event = OS) ~ T_Stage, data = df.full) %>% summary
res.t = c(res.t$conf.int[-2], res.t$coefficients[5])

df.full$N_Stage = str_match(df.full$N_Stage, "N[0-2]")[,1]
df.full$N_Stage = c("N0"="N0","N1"="N1-N2","N2"="N1-N2")[df.full$N_Stage]
res.n = coxph(Surv(time = OS.time, event = OS) ~ N_Stage, data = df.full) %>% summary
res.n = c(res.n$conf.int[-2], res.n$coefficients[5])

df.full$M_Stage = str_match(df.full$M_Stage, "M[0-1]")[,1]
res.m = coxph(Surv(time = OS.time, event = OS) ~ M_Stage, data = df.full) %>% summary
res.m = c(res.m$conf.int[-2], res.m$coefficients[5])

res = rbind(res.risk,res.age,res.gender,res.t,res.n,res.m) %>% as.data.frame()
res$Indicators = c("Risk", "Age", "Gender", "T Stage\n(T3-4 vs T1-2)", "N Stage\n(N1-2 vs N0)", "M Stage\n(M0 vs M1)")
colnames(res) = c("hr","low","up","pv","Indicator")
res$p = signif(res$pv, 2) %>% paste0("p = ", .)
res$p[is.na(res$pv)] = NA
res$Indicator = factor(res$Indicator, levels = rev(res$Indicator))

ggplot(data = res) +
  geom_point(mapping = aes(x = hr, y = Indicator), pch = 15, size = 3) +
  geom_errorbar(mapping = aes(xmin = low, xmax = up, y = Indicator), width = 0.1) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60") +
  geom_text(aes(x = 0, y = Indicator, label = p), hjust = 0) +
  theme_pubclean(base_size = 14) + 
  xlab("Harzard Ratio")
ggsave("05.Prognostic/01.Univar.Cox.png", width = 9, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("05.Prognostic/01.Univar.Cox.pdf", width = 9, height = 5, units = "in", dpi = 300, bg = "white")

res.mul = coxph(Surv(OS.time,event=OS)~risk+Age+T_Stage+N_Stage+M_Stage, data = df.full) %>% summary
res.mul = cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
res.mul$Indicators = c("Risk", "Age", "T Stage\n(T3-4 vs T1-2)", "N Stage\n(N1-2 vs N0)", "M Stage\n(M0 vs M1)")
colnames(res.mul) = c("hr","low","up","pv","Indicator")
res.mul$p = signif(res.mul$pv, 2) %>% paste0("p = ", .)
res.mul$p[is.na(res.mul$pv)] = NA
res.mul$Indicator = factor(res.mul$Indicator, levels = rev(res.mul$Indicator))
ggplot(data = res.mul) +
  geom_point(mapping = aes(x = hr, y = Indicator), pch = 15, size = 3) +
  geom_errorbar(mapping = aes(xmin = low, xmax = up, y = Indicator), width = 0.1) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60") +
  geom_text(aes(x = -1, y = Indicator, label = p), hjust = 0) +
  theme_pubclean(base_size = 14) + 
  scale_x_continuous(breaks = c(1,5,10,15), limits = c(-1,10)) +
  xlab("Harzard Ratio")
ggsave("05.Prognostic/02.Multivar.Cox.png", width = 9, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("05.Prognostic/02.Multivar.Cox.pdf", width = 9, height = 5, units = "in", dpi = 300, bg = "white")

library(rms)
colnames(df.full)[5] = "Risk"
model4nomo = rms::cph(Surv(time = OS.time, event = OS)~Risk+Age+T_Stage+N_Stage+M_Stage, data = df.full, surv = T, x = T, y = T)
data.dist = datadist(df.full)
options(datadist = "data.dist")

surv <- Survival(model4nomo)
surv4 = function(x)surv(365*7,lp=x)
surv1 <- function(x)surv(365*3,lp=x) 
surv2 <- function(x)surv(365*5,lp=x) 
surv3 <- function(x)surv(365*1,lp=x) 

nomoplot = nomogram(model4nomo, fun = list(surv3,surv1,surv2), lp = F, funlabel = c("1-year Survival", "3-year Survival", "5-year Survival"),
                    maxscale = 100, fun.at=c('0.95', '0.9','0.85','0.8','0.7','0.6','0.5','0.3','0.1'))
png("05.Prognostic/03.Nomogram.png", width = 10, height = 8, res = 1200, bg = "white", units = "in")
plot(nomoplot, xfrac=0.2)
dev.off()
pdf("05.Prognostic/03.Nomogram.pdf", width = 10, height = 8)
plot(nomoplot, xfrac=0.2)
dev.off()

model4cal1 = rms::cph(Surv(time = OS.time, event = OS) ~ Risk+Age+T_Stage+N_Stage+M_Stage, 
                      data = df.full, time.inc = 365,
                      surv = T, x = T, y = T)
set.seed(123456)
cal1 = calibrate(model4cal1, cmethod="KM", method="boot", u=365, m=80, B=500)

png("05.Prognostic/04.CaliberationCurve_1year.png", width = 7, height = 7, units = "in", res = 1200)
plot(cal1, subtitles = F, 
     xlab = "Predicted 1 year Survival", ylab = "Actual 1 year Survival", cex.lab = 1.2)
dev.off()
pdf("05.Prognostic/04.CaliberationCurve_1year.pdf", width = 7, height = 7)
plot(cal1, subtitles = F, 
     xlab = "Predicted 1 year Survival", ylab = "Actual 1 year Survival", cex.lab = 1.2)
dev.off()

model4cal2 = rms::cph(Surv(time = OS.time, event = OS) ~ Risk+Age+T_Stage+N_Stage+M_Stage, 
                      data = df.full, time.inc = 365*3,
                      surv = T, x = T, y = T)
set.seed(123456)
cal2 = calibrate(model4cal2, cmethod="KM", method="boot", u=365*3, m=80, B=500)

png("05.Prognostic/05.CaliberationCurve_3year.png", width = 7, height = 7, units = "in", res = 1200)
plot(cal2, subtitles = F, 
     xlab = "Predicted 3 year Survival", ylab = "Actual 3 year Survival", cex.lab = 1.2)
dev.off()
pdf("05.Prognostic/05.CaliberationCurve_3year.pdf", width = 7, height = 7)
plot(cal2, subtitles = F, 
     xlab = "Predicted 3 year Survival", ylab = "Actual 3 year Survival", cex.lab = 1.2)
dev.off()

model4cal3 = rms::cph(Surv(time = OS.time, event = OS) ~ Risk+Age+T_Stage+N_Stage+M_Stage, 
                      data = df.full, time.inc = 365*5,
                      surv = T, x = T, y = T)
set.seed(123456)
cal3 = calibrate(model4cal3, cmethod="KM", method="boot", u=365*5, m=80, B=500)

png("05.Prognostic/06.CaliberationCurve_5year.png", width = 7, height = 7, units = "in", res = 1200)
plot(cal3, subtitles = F, 
     xlab = "Predicted 5 year Survival", ylab = "Actual 5 year Survival", cex.lab = 1.2)
dev.off()
pdf("05.Prognostic/06.CaliberationCurve_5year.pdf", width = 7, height = 7)
plot(cal3, subtitles = F, 
     xlab = "Predicted 5 year Survival", ylab = "Actual 5 year Survival", cex.lab = 1.2)
dev.off()
