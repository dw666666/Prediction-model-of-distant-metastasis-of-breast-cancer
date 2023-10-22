rm(list = ls())
options(stringsAsFactors = F)
library(GEOquery)
require(stats)
gse_number = "GSE9893"
eSet <- getGEO(gse_number, 
               destdir = '.', 
               getGPL = F)

class(eSet)
length(eSet)
eSet = eSet[[1]]
#(1)提取表达矩阵exp
exp<- exprs(eSet)

boxplot(exp)

#(2)提取临床信息
pd <- pData(eSet)

#(3)调整pd的行名顺序与exp列名完全一致
p = identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
#(4)提取芯片平台编号
gpl_number <- eSet@annotation


save(pd,file = "steppdForWGCNAoutput.Rdata")
save(gse_number,pd,exp,gpl_number,file = "step1__1output.Rdata")
