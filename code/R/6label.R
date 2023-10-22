rm(list = ls()) 
load(file = "step1__1output.Rdata")
load(file = "step6.xoutput.Rdata")
library(stringr)
library(dplyr)


label<-pd[,c("characteristics_ch1.15","characteristics_ch1")]
label<-select(label,-2)
label$characteristics_ch1.15=ifelse(str_detect(pd$characteristics_ch1.15,"no"),"control","metastasis")
colnames(label)="label"
p = identical(rownames(data),rownames(label));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]

save(label,file = "step6.youtput.Rdata")
