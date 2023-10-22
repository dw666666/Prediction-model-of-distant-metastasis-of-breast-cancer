rm(list = ls()) 
load("limma&WGCNA.RData")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("DOSE")
library("ggnewscale")

genes=as.vector(degsLimmaANDWGCNA$symbol)
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs <- as.character(entrezIDs)
degsLimmaANDWGCNA$entrezID=entrezIDs
degsLimmaANDWGCNA=degsLimmaANDWGCNA[degsLimmaANDWGCNA$entrezID!="NA",]
gene=degsLimmaANDWGCNA$entrezID
kk <- enrichGO(gene = gene,
               
               OrgDb = org.Hs.eg.db,
               
               pvalueCutoff =0.15,
               
               qvalueCutoff = 0.1,
               
               ont="all",
               
               readable =T)


dotplot(kk, split="ONTOLOGY",showCategory = 10,label_format=100) + facet_grid(ONTOLOGY~., scale="free")+
  theme(axis.text.y = element_text(size = 8))

#bar<-barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale='free')+ theme_bw(base_size = 10)+
#  theme(plot.title = element_text(hjust = 0.5),  axis.text.y = element_text(size = 5))
#barplot(kk, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
KEGG<-enrichKEGG(gene,#KEGG富集分析
                 organism = "hsa",
                 pvalueCutoff = 0.3,
                 qvalueCutoff = 0.3)
dotplot(KEGG,label_format=100)+
  theme(axis.text.y = element_text(size = 8))
KEGG2 <- pairwise_termsim(KEGG)
enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk")


