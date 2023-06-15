library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final  
library(Seurat)
table(Idents(sce))

deg = FindMarkers(sce,ident.1 = 'NK',
            ident.2 = 'B')
head(deg[order(deg$p_val),])
table(Idents(sce))

library(EnhancedVolcano)
res=deg
head(res)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'avg_log2FC',
                y = 'p_val_adj')

plot(deg$avg_log2FC,(deg$pct.1 - deg$pct.2))
barplot()


CD14_deg = FindMarkers(sce,ident.1 = 'CD14+ Mono',
                  ident.2 = 'B',
                  logfc.threshold = 0,
                  min.pct = 0
                  )
head(CD14_deg[order(CD14_deg$p_val),])
FCGR3A_deg = FindMarkers(sce,ident.1 = 'FCGR3A+ Mono',
                       ident.2 = 'B',
                       logfc.threshold = 0,
                       min.pct = 0)
head(FCGR3A_deg[order(FCGR3A_deg$p_val),])


ids=c(rownames(CD14_deg[abs(CD14_deg$avg_log2FC)>2,]),
              rownames(FCGR3A_deg[abs(FCGR3A_deg$avg_log2FC)>2,]))
ids
df= data.frame(
  FCGR3A_deg = FCGR3A_deg[ids,'avg_log2FC'],
  CD14_deg = CD14_deg[ids,'avg_log2FC']
)
library(ggpubr)
ggscatter(df, x = "FCGR3A_deg", y = "CD14_deg",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson",  label.sep = "\n")
)

library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final  
library(Seurat)
table(Idents(sce))

sce$celltype = Idents(sce)
sce$group = sample(1:2,ncol(sce),replace = T)
table(sce$celltype,sce$group )

Idents(sce) = paste0('c',sce$group )
table(Idents(sce))
degs = lapply(unique(sce$celltype), function(x){
  FindMarkers(sce[,sce$celltype==x],ident.1 = 'c1',
              ident.2 = 'c2')
})
x=degs[[1]]
do.call(rbind,lapply(degs, function(x){
  table(x$avg_log2FC > 0 )
}))



