library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final  

library(Seurat)
genes_to_check = c("CD3D","CD3E" )
#genes_to_check = c("CD4","CD8A" )
FeaturePlot(sce,genes_to_check)
library(stringr)  
library(ggplot2)  
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  ) +theme(axis.text.x = element_text(angle = 90))

p
mat = sce@assays$RNA@counts[ genes_to_check ,]
table(
  mat[1,]>0 ,mat[2,]>0 
)
sce$ok = mat[1,]>0  | mat[2,]>0 
table(sce$ok )
FeaturePlot(sce,'ok')


library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A',
                   'CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'C1QA',  'C1QB',  # mac
                   'S100A9', 'S100A8', 'MMP19',# monocyte
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'KLRB1','NCR1', # NK 
                   'FGF7','MME', 'ACTA2', ## fibo 
                   'DCN', 'LUM',  'GSN' , ## mouse PDAC fibo 
                   'Amy1' , 'Amy2a2', # Acinar_cells
                   'PECAM1', 'VWF',  ## endo 
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )
genes_to_check=genes_to_check[genes_to_check %in% rownames(sce)]
library(Seurat)
library(ggplot2) 
pl = lapply(genes_to_check, function(cg){  FeaturePlot(sce, cg,) + NoLegend() + NoAxes() })
ps <- cowplot::plot_grid(plotlist = pl)
ps  
ggsave("FeaturePlot_umap.pdf",width = 16,height = 15)

pl = lapply(genes_to_check, function(cg){  FeaturePlot(sce, cg,order = T,raster = T) + NoLegend() + NoAxes() })
ps <- cowplot::plot_grid(plotlist = pl)
ps  
ggsave("FeaturePlot_umap2.pdf",width = 16,height = 15)


