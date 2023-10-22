library(WGCNA)
options(stringsAsFactors = FALSE)
rm(list = ls()) 
lnames=load("stepWGCNA_1output.Rdata")
lnames
lnames=load("FemaleLiver-02-networkConstruction-auto.RData")
lnames

dataTraits$No_metastasis=dataTraits$characteristics_ch1.15
dataTraits=dataTraits[,-1]
colnames(dataTraits)=c("Metastasis","No_metastasis")
library(stringr)
dataTraits$No_metastasis=ifelse(str_detect(dataTraits$Metastasis,"1"),0,1)

nGenes=ncol(exp)
nSamples=nrow(exp)
MEs0=moduleEigengenes(exp,modulColors)$eigengenes
MEs=orderMEs(MEs0)
moduleTraitCor=cor(MEs,dataTraits,use="p");
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples);

sizeGrWindow(10,6)
#展示模块与表型数据的相关系数和P值
textMatrix=paste(signif(moduleTraitCor,2),"\n(",
                 signif(moduleTraitPvalue,1),")",sep="");
dim(textMatrix)=dim(moduleTraitCor)
par(mar=c(6,8.5,3,3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(dataTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim=c(-1,1),
               main=paste("Module-trait relationshaips"))


metastasis=as.data.frame(dataTraits$Metastasis)
names(metastasis)="metastasis"
modNames=substring(names(MEs),3)
geneModuleMembership=as.data.frame(cor(exp,MEs,use="p"));
MMPvalue=as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),nSamples));
names(geneModuleMembership)=paste("MM",modNames,sep="");
names(MMPvalue)=paste("p.MM",modNames,sep="");
geneTraitSignificance=as.data.frame(cor(exp,metastasis,use = "p"));#和转移的关联
GSPvalue=as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples));
names(geneTraitSignificance)=paste("GS.",names(metastasis),sep="");
names(GSPvalue)=paste("p.GS.",names(metastasis),sep="");



module = "midnightblue"
column = match(module,modNames);
moduleGenes = modulColors==module;
sizeGrWindow(7,7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[ moduleGenes,column]),
                   abs(geneTraitSignificance[ moduleGenes,1]),
                   xlab = paste("Module Membership in",module,"module"),
                   ylab = "Gene significance for metastasis",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2,cex.lab = 1.2,cex.axis = 1.2,col = module)
table(moduleGenes)
moduleGenesMidnightblue=moduleGenes



module = "blue"
column = match(module,modNames);
moduleGenes = modulColors==module;
sizeGrWindow(7,7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[ moduleGenes,column]),
                   abs(geneTraitSignificance[ moduleGenes,1]),
                   xlab = paste("Module Membership in",module,"module"),
                   ylab = "Gene significance for metastasis",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2,cex.lab = 1.2,cex.axis = 1.2,col = module)
table(moduleGenes)
moduleGenesblue=moduleGenes


module = "greenyellow"
column = match(module,modNames);
moduleGenes = modulColors==module;
sizeGrWindow(7,7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[ moduleGenes,column]),
                   abs(geneTraitSignificance[ moduleGenes,1]),
                   xlab = paste("Module Membership in",module,"module"),
                   ylab = "Gene significance for metastasis",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2,cex.lab = 1.2,cex.axis = 1.2,col = module)
table(moduleGenes)
moduleGenesgreenyellow=moduleGenes



class(exp)
exp1<-as.matrix(exp)
exp1<-t(exp1)
exp1<-as.data.frame(exp1)
exp1$symbol=rownames(exp1)

#degsMidnightblue<-exp1[moduleGenesMidnightblue,"symbol"]
degsblue<-exp1[moduleGenesblue,"symbol"]
degsgreenyellow<-exp1[moduleGenesgreenyellow,"symbol"]
degs2<-c(degsgreenyellow,degsblue)
save(degs2,file = "degs2.RData")




names(exp)
names(exp)[modulColors=="midnightblue"]
annot=read.csv(file = "")






