rm(list = ls()) 
load("step3__3DEGoutput.Rdata")
load("step3-3_42837DEGoutput.Rdata")
load("degs2.RData")
degs<-deg[(deg$change)!="stable",]
deg42837<-deg42837[(deg42837$change)!="stable",]

degs<-degs[which(deg42837$symbol%in%degs$symbol),]
degsLimmaANDWGCNA<-degs[which(degs$symbol%in%degs2),]



#degsLimmaANDWGCNA<-degsLimmaANDWGCNA[1:85,]



#degsLimmaANDWGCNA<-degsLimmaANDWGCNA[which(degsLimmaANDWGCNA$symbol!="Dynamic range 4"),]
save(degsLimmaANDWGCNA,file = "limma&WGCNA.RData")


ppiGene=degsLimmaANDWGCNA$symbol
write.csv(ppiGene,"ppi.csv",row.names = FALSE)






library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

deg1=degs$symbol
deg2=deg42837$symbol
deg3=degs2
x = list(deg1, deg2,deg3)

venn.diagram(
  x,
  category.names = c("DIFF(GSE9893)" , "DIFF(GSE43837)","WGCNA"),
  filename = 'venn_plot2.png',
  output=TRUE,
  
 
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300,
  compression = "lzw",
  
 
  lwd = 2, 
  lty = "blank",  
  fill = myCol,  
  
 
  cex = .5,  ；
  fontface = "bold",  
  fontfamily = "sans",  
  
  
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",  
  cat.pos = c(-27, 27, 180),  
  cat.dist = c(0.055, 0.055, 0.055),  
  cat.fontfamily = "sans",
  rotation = 1  
)

inter <- get.venn.partitions(x)
inter=inter[1,5][["1"]]
inter






