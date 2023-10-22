rm(list = ls())  
load(file = "step2_2output.Rdata")

table(!duplicated(ids$probe_id))
table(!duplicated(ids$symbol))

#ids = ids[order(deg$P.Value,-abs(deg$logFC)),]
library(dplyr)
exp<-as.data.frame(exp)
exp <- mutate(exp,probe_id=rownames(exp))
ids= ids[!duplicated(ids$symbol),]
ids$probe_id<-as.character(ids$probe_id)
exp <- inner_join(exp,ids,by="probe_id")
nrow(exp)
exp<-exp[-which(exp$symbol==""),]
rownames(exp)<-exp$symbol
exp<-exp[,-c(156,157)]


exp<-as.matrix(exp)
#筛选基因并转置矩阵
exp=t(exp)



exp<-as.data.frame(exp)
for(i in 1:ncol(exp)){
  exp[,i]<-as.numeric(exp[,i])
}


class(exp$JRKL)
dim(exp)
#exp=as.data.frame(lapply(exp,as.numeric))
#exp <- as.numeric(as.character(unlist(lapply(exp, as.numeric))))
save(exp,file = "stepdataForWGCNAoutput.Rdata")
