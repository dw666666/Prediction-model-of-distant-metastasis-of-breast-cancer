rm(list = ls()) 
load(file = "step2-2output.Rdata")

library(limma)
design=model.matrix(~Group)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
#为deg数据框添加几列

library(dplyr)
deg <- mutate(deg,probe_id=rownames(deg))
head(deg)
#加上探针注释
table(!duplicated(ids$probe_id))
table(!duplicated(ids$symbol))
#按symbol列去重
ids= ids[!duplicated(ids$symbol),]


deg <- inner_join(deg,ids,by="probe_id")
head(deg)
nrow(deg)

#加change列,标记上下调基因
logFC_t=0.5
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg <- mutate(deg,change)
table(deg$change)
deg42837<-deg
#加ENTREZID列，用于富集分析（symbol转entrezid，然后inner_join）
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(deg$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)#人类

dim(deg)
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
dim(deg)
length(unique(deg$symbol))
save(Group,deg,logFC_t,P.Value_t,gse_number,file = "step3-3DEGoutput.Rdata")
save(deg42837,file = "step3-3_42837DEGoutput.Rdata")
