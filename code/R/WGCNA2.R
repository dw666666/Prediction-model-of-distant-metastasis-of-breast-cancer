#构建表达网络
library(WGCNA)
options(stringsAsFactors = FALSE)
rm(list = ls()) 
lnames=load(file = "stepWGCNA_1output.Rdata")
lnames


#选择软阈值
powers=c(c(1:10),seq(from=12,to=20,by=2))
powers
sft=pickSoftThreshold(exp,powerVector=powers,verbose=5)
sft
#绘图准备
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1=0.9#字符大小

plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",ylab = "Scale Free Topology Model Fit,signed R~2",type = "n",
     main = paste("Scale independemce"));
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers,cex=cex1,col = "red");
     abline(h=0.90,col="red")#查看位于0.9以上的点，官网是0.85
#平均连接度

plot(sft$fitIndices[,1],sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",ylab = "Mean Connectivity",type = "n",
     main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels = powers,cex=cex1,col="red")

sft$powerEstimate

#if(is.na(power)){
#  power=ifelse(nSamples<20,ifelse(type=="unsigned",9,18),
 #              ifelse(nSamples<30,ifelse(type=="unsigned",8,16),
  #                    ifelse(nSamples<40,ifelse(type=="unsigned",7,14),
   #                          ifelse(type=="unsigned",6,12))
    #           )
  #)
#}


nGenes=ncol(exp)
load(file = "WGCNA_net.RData")
cor=WGCNA::cor
net=blockwiseModules(exp,power = 6,maxBlockSize = nGenes,
                     TOMType = "unsigned",minModuleSize = 30,
                     reassignThreshold = 0,mergeCutHeight = 0.5,
                     numericLabels = TRUE,pamRespectsDendro = FALSE,
                     saveTOMs = TRUE,
                     saveTOMFileBase = "femaleMouseTOM",
                     verbose = 3)
cor=stats::cor

##查看划分的模块数和每个模块里面包含的基因个数table(net$colors)
table(net$colors)





#模块标识的层次聚类树状图
sizeGrWindow(12,9)
mergedColors=labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,hang=0.03,
                    addGuide = TRUE,guideHang = 0.05)

#保存分配模块和模块包含的基因信息
moduleLabels=net$colors
modulColors=labels2colors(net$colors)
MEs=net$MEs;
geneTree=net$dendrograms[[1]];
save(MEs,moduleLabels,modulColors,geneTree,
     file = "FemaleLiver-02-networkConstruction-auto.RData")
save(net,file = "WGCNA_net.RData")