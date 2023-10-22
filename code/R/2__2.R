# Group(实验分组)和ids(探针注释)
rm(list = ls())  
load(file = "step1__1output.Rdata")
library(stringr)

Group=ifelse(str_detect(pd$characteristics_ch1.15,"no"),"control","metastasis")

#设置参考水平，指定levels，对照组在前，处理组在后
Group = factor(Group,
               levels = c("control","metastasis"))
Group
table(Group)
# 注意levels与因子内容必须对应一致


gpl_number 

if(!require(hgu133a2))BiocManager::install("hgu133a2")
library(hgu133a2.db)
ls("package:hgu133a2.db")
ids <- toTable(hgu133a2SYMBOL)
head(ids)


gpl=getGEO(filename = 'GPL5049.soft.gz')
gpl=gpl@dataTable@table
ids<-gpl[,c("ID","Gene_Symbol")]
colnames(ids)<-c('probe_id','symbol')
save(exp,Group,ids,gse_number,file = "step2_2output.Rdata")
