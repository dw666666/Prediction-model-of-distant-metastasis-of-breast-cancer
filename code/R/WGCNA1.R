#数据预处理
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("WGCNA")
#nBiocManager::install("Rcpp")
#BiocManager::install("GenomeInfoDbData")
#library(GenomeInfoDbData)
#library('Rcpp')
#BiocManager::install("GO.db")
library(WGCNA)
options(stringsAsFactors = FALSE)
rm(list = ls()) 
load(file = "stepdataForWGCNAoutput.Rdata")
#检查缺失值和识别离群值（异常值）


dim(exp)
gsg=goodSamplesGenes(exp,verbose = 3)
gsg$allOK#如果是false则需要删除缺失值


#聚类所有样本，观察是否有离群值或者异常值
sampleTree=hclust(dist(exp),method = "average")
sizeGrWindow(12,10)
par(cex=0.6)
par(mar=c(0,6,2,0))
plot(sampleTree,main = "Sample clustering to detect outliers",sub = "",xlab = "",
     cex.lab=1.5,cex.axis=1.5,cex.main=2)


abline(h=150,col="red")
clust=cutreeStatic(sampleTree,cutHeight = 150,minSize =10)
table(clust)
keepSamples=(clust==1)
exp=exp[keepSamples,]
nGenes=ncol(exp)
nSamples=nrow(exp)

#载入表型数据
load(file = "steppdForWGCNAoutput.Rdata")
dim(pd)
names(pd)
allPd=pd[,c(24,25)]
dim(allPd)
names(allPd)
femaleSamples=rownames(exp)
traitROWs=match(femaleSamples,rownames(allPd))
dataTraits=allPd[traitROWs,]
table(dataTraits$characteristics_ch1.15)
library(stringr)
dataTraits$characteristics_ch1.14=ifelse(str_detect(dataTraits$characteristics_ch1.14,"no"),0,1)
dataTraits$characteristics_ch1.15=ifelse(str_detect(dataTraits$characteristics_ch1.15,"no"),0,1)
collectGarbage()
#可视化表型数据与基因表达量的联系，重构样本聚类树
sampleTree2=hclust(dist(exp),method = "average")
traitColors=numbers2colors(dataTraits$characteristics_ch1.15,signed = FALSE)
plotDendroAndColors(sampleTree2,traitColors,
                    groupLabels = "distant_metastases",
                    main="Sample dendrogram and trait heatmap")
save(exp,dataTraits,file = "stepWGCNA_1output.Rdata")
