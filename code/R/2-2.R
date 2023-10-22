# Group(实验分组)和ids(探针注释)
rm(list = ls())  
load(file = "step1-1output.Rdata")
library(stringr)

Group=ifelse(str_detect(pd$`tissue:ch1`,"primary"),"control","metastasis")

#设置参考水平，指定levels，对照组在前，处理组在后
Group = factor(Group,
               levels = c("control","metastasis"))
Group
table(Group)
# 注意levels与因子内容必须对应一致


gpl_number 

if(!require(hgu133plus2))BiocManager::install("hgu133plus2")
library(u133x3p.db)
ls("package: u133x3p.db")
ids <- toTable( u133x3pSYMBOL)
head(ids)
save(exp,Group,ids,gse_number,file = "step2-2output.Rdata")
