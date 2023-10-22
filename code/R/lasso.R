rm(list = ls()) 
load(file = "step2_2output.Rdata")
load(file = "limma&WGCNA.RData")
library(glmnet)


exp<-exp[rownames(exp)%in%degsLimmaANDWGCNA$probe_id,]
exp<-as.matrix(exp)
x<-t(exp)
y <- Group
alpha1_fit <- glmnet(x,y,alpha=1,family="binomial")
plot(alpha1_fit,xvar="lambda",label=FALSE)
set.seed(101)
alpha1.fit <- cv.glmnet(x,y,type.measure = "class",alpha=1,family="binomial")
plot(alpha1.fit)
print(alpha1.fit)

coef(alpha1_fit,s=alpha1.fit$lambda.min)
#90nothing<-c(1266,1738,2276,2330,5467,5755,7253,8068,8900,9071,10334,10335,10763,12619,13139,15009,19415,19435)
nothing<-c(2483,4060,4911,5546,5610,5768,7575,7891,8292,8325,8750,10629,11839,13031,13637,13713,13829,14321,14762,16408,19413)
deglaaso<-degsLimmaANDWGCNA[(which(degsLimmaANDWGCNA$probe_id%in%nothing)),]




save(deglaaso,file = "step_lasso_output.Rdata")
write.table(deglaaso$symbol,file="TF.txt")
