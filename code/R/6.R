rm(list = ls()) 
load(file = "step6.xoutput.Rdata")
load(file = "step6.youtput.Rdata")
library(dplyr)


data<-as.data.frame(data)
data <- mutate(data,label$label)
colnames(data)[22]<-"label"
sort<-colnames(data)
sort<-sort[-22]
#colnames(data)<-c("CP","NQO1","KRT17","ESR1","PNMT","KRT14","label")
write.csv(data,file = "data.csv")
save(sort,file = "step_datasort_output.Rdata")
