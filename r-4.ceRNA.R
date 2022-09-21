rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-336/")
if (! dir.exists("./04_ceRNA")){
  dir.create("./04_ceRNA")
}
setwd("./04_ceRNA")

### 3个差异基因
DEmRNA<-read.delim2('../01_DEmRNAs/DEG_sig(mRNA).xls',row.names = 1)
DElncRNA<-read.delim2('../02_DElncRNAs/DEG_sig(lncRNA).xls',row.names = 1)
DEmiRNA<-read.delim2('../03_DEmiRNAs/DEG_sig(miRNA).xls',row.names = 1)
sig.all<-rbind(DEmRNA,DElncRNA,DEmiRNA)
dim(sig.all)
##7386
write.table(sig.all,file = 'DEG.sig.xls',sep = '\t',row.names = T,quote = F)
