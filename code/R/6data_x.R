rm(list = ls())  
load(file = "step2_2output.Rdata")
load(file = "step_lasso_output.Rdata")





id<-deglaaso[,c("probe_id","symbol")]
data<-exp[which(rownames(exp)%in%id$probe_id),]
data<-t(data)
a<-id[match(colnames(data),id$probe_id),"symbol"]
a
colnames(data)<-a
save(data,file = "step6.xoutput.Rdata")

