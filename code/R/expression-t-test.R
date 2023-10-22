rm(list = ls()) 
library(stringr)
library(ggplot2)
library(gghalves)
library(tidyverse)
library(ggpubr)
load("dataForT_Test.Rdata")
load("steppdForWGCNAoutput.Rdata")
p = identical(rownames(pd),rownames(exp));p
if(!p) exp = exp[,match(rownames(pd),rownames(exp))]
Group=ifelse(str_detect(pd$characteristics_ch1.15,"no"),"control","metastasis")
table(Group)
exp=as.data.frame(exp)
exp$group=Group



for (i in 1:21) {
  print(i)
  metastasis=exp[exp$group=="metastasis",][,i]
  control=exp[exp$group=="control",][,i]
  x <- c(metastasis,control)
  group <- c(rep("metastasis",48),rep("control",107))
  shapiro.test(metastasis) 
  shapiro.test(control) 
  bartlett.test(x~group)
  print(t.test(metastasis,control,paired = FALSE,var.equal = F))
}
#shapiro.test(metastasis) 

#Shapiro-Wilk normality test



for (i in 1:nrow(exp)) {
  rownames(exp)[i]=paste(rownames(exp)[i], "_",exp$group[i])
}

exp=t(exp)
exp=as.data.frame(exp)

exp=exp[-22,]
exp$gene=rownames(exp)


exp <- exp %>% 
  pivot_longer(cols = !gene, 
               names_to = "Samples", 
               values_to = "Values")

colnames(exp)[1] <- "Genes"



exp$group <- str_split(exp$Samples, "_", simplify = T)[,2]

head(exp)
exp$Values<-as.double(exp$Values)
ggplot()+
  geom_half_violin(
    data = exp %>% filter(group == " control"),
    aes(x = Genes,y = Values),colour="white",fill="#1ba7b3",side = "l"
  )+
  geom_half_violin(
    data = exp %>% filter(group == " metastasis"),
    aes(x = Genes,y = Values),colour="white",fill="#dfb424",side = "r"
  )+
  theme_bw()+
  xlab("")+
  ylab("Expression")+
  geom_point(data = exp, aes(x = Genes,y = Values, fill = group),
             stat = 'summary', fun=mean,
             position = position_dodge(width = 0.2),color="black")+

  stat_summary(data = exp, aes(x = Genes,y = Values, group= group),
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar', color='black',
               width=0.01,size=0.5,
               position = position_dodge(width = 0.2))+
  stat_compare_means(data = exp, aes(x = Genes,y = Values, group = group),
                    
                     method = "t.test",
                     symnum.args=list(cutpoints = c(1e-03, 1e-04, 1e-05, 1e-06,1e-7,1e-8,1e-9,1e-10),
                                     symbols = c("p<1e-03","p<1e-04", "p<1e-05", "p<1e-06","p<1e-07","p<1e-08","p<1e-09")),
                     label = "p.signif",
                     label.y = max(exp$Values),
                     hide.ns = F)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "top",
        panel.grid=element_blank(),
        legend.justification = "right",
        
        )
  


ggsave("violin_plot.pdf", height = 5, width = 10)




GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

ggplot(exp, aes(x = Genes,y = Values, fill = group))+
  geom_split_violin(trim = T,colour="white")+
  geom_point(stat = 'summary',fun=mean,
             position = position_dodge(width = 0.2))+
  scale_fill_manual(values = c("#1ba7b3","#dfb424"))+
  stat_summary(fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar',color='black',
               width=0.01,size=0.5,
               position = position_dodge(width = 0.2))+
  stat_compare_means(data = exp, aes(x = Genes,y = Values),
                     method = "t.test",
                     
                     symnum.args=list(cutpoints = c(1e-03, 1e-04, 1e-05, 1e-06,1e-7,1e-8,1e-9,1e-10),
                                      symbols = c("p<1e-03","p<1e-04", "p<1e-05", "p<1e-06","p<1e-07","p<1e-08","p<1e-09")),
                     label = "p.signif",
                     label.y = max(exp$Values),
                     hide.ns = F)+
  theme_bw()+
  xlab("")+
  ylab("Expression")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "top",
        #legend.key = element_rect(fill = c("#1ba7b3","#dfb424")),
        legend.justification = "right")
