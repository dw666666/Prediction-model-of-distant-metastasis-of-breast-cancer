rm(list = ls()) 
lnames=load("limma&WGCNA.RData")

library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
library(GOplot)
library(DOSE)
library(ggnewscale)
#library(ggridges)



library(topGO)
library(circlize)
library(ComplexHeatmap)


GO_database <- 'org.Hs.eg.db' 
KEGG_database <- 'hsa' 

#gene ID转换
gene <- bitr(degsLimmaANDWGCNA$symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
gene<-gene$SYMBOL

GO<-enrichGO(gene,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENTREZID",
              ont = "ALL",
              pvalueCutoff = 0.5,
              qvalueCutoff = 0.5,
              readable = T)
GO<-enrichGO(gene = gene,
         
         OrgDb = org.Hs.eg.db,
         
         pvalueCutoff =0.5,
         
         qvalueCutoff = 0.5,
         
         ont="all",
         
         readable =T)

library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = "hsa",
                 pvalueCutoff = 0.1,
                 qvalueCutoff = 0.1)

#barplot(KEGG,showCategory = 40,title = 'KEGG Pathway')
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
dotplot(KEGG)



enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)
enrichplot::cnetplot(KEGG,circular=FALSE,colorEdge = TRUE)


#富集到的功能集/通路集之间的关联网络图
GO2 <- pairwise_termsim(GO)
KEGG2 <- pairwise_termsim(KEGG)
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")
enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk")



#GSEA
colnames(degsLimmaANDWGCNA)[8] <- 'SYMBOL'
degsLimmaANDWGCNA_merge <- merge(degsLimmaANDWGCNA,gene,by='SYMBOL')
GSEA_input <- degsLimmaANDWGCNA_merge$logFC
names(GSEA_input) = degsLimmaANDWGCNA_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE)
GSEA_GO <- gseGO(GSEA_input, OrgDb ="org.Hs.eg.db", pvalueCutoff = 0.5)
KEGG_kk_entrez <- gseKEGG(geneList     = GSEA_input,
                          organism     = 'hsa', 
                          pvalueCutoff = 3)  

#GSEA富集图
library(ggridges)
ridgeplot(GO_kk_entrez) 
gseaplot2(GO_kk_entrez,1)
gseaplot2(GO_kk_entrez,1:3)



