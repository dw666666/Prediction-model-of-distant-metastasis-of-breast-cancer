rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
gse_number = "GSE43837"
eSet <- getGEO(gse_number, 
               destdir = '.', 
               getGPL = F)

class(eSet)
length(eSet)
eSet = eSet[[1]]
#(1)提取表达矩阵exp
exp<- exprs(eSet)
exp = log2(exp+1)
exp[1:4,1:4]

boxplot(exp)
#(2)提取临床信息
pd <- pData(eSet)
#pd<-pd[which(pd$characteristics_ch1.1%in%c("metastasis (1: yes, 0: no/censored): 0","metastasis (1: yes, 0: no/censored): 1")),]


#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl_number <- eSet@annotation



save(gse_number,pd,exp,gpl_number,file = "step1-1output.Rdata")
