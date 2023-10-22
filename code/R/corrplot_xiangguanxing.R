rm(list = ls()) 
#install.packages("corrplot")#安装包
library("corrplot")#加载包
load("step_lasso_output.Rdata")
load("step2_2output.Rdata")
library(dplyr)
exp<-as.data.frame(exp)
exp <- mutate(exp,probe_id=rownames(exp))
#加上探针注释
table(!duplicated(ids$probe_id))
table(!duplicated(ids$symbol))
#按symbol列去重
#ids = ids[order(deg$P.Value,-abs(deg$logFC)),]
ids= ids[!duplicated(ids$symbol),]
ids$probe_id<-as.character(ids$probe_id)
exp <- inner_join(exp,ids,by="probe_id")
exp<-exp[which(exp$probe_id%in%deglaaso$probe_id),]
rownames(exp)=exp$symbol
exp<-exp[,-c(156,157)]
#exp<-exp[Group=="metastasis"]



exp <- as.matrix(exp)
exp<-t(exp)
qc <- cor(exp)
qc
col3 <- colorRampPalette(c("blue", "white","red"))
corrplot(qc, order = "hclust", addrect = 2, col = col3(20))
corrplot(qc, type = "lower", tl.pos = "tp",tl.cex=.7,col = col3(20),tl.col = "black")
corrplot(qc, add = TRUE, type = "upper", method = "number",
         col = col3(20), diag = FALSE, tl.pos = "n", cl.pos = "n", number.cex = 0.3)
save(exp,file = "dataForT_Test.Rdata")
