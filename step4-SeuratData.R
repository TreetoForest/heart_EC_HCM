# install.packages('devtools')
# devtools::install_github('satijalab/seurat-data')
library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final  

library(scRNAstat) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)

DimPlot(sce,label = T)

x='check_pbmc3k_by_scRNAstat'
dir.create( x )
sce = basic_qc(sce=sce,org='human',
               dir = x)  
sce = basic_filter(sce)  
sce = basic_workflow(sce,dir = x)   
markers_figures <- basic_markers(sce,
                                 org='human',
                                 group='seurat_clusters',
                                 dir = x)
p_umap = DimPlot(sce,reduction = 'umap',  
                 group.by = 'seurat_clusters',
                 label.box = T,  label = T,repel = T)
p_umap+markers_figures[[1]]

ggsave(paste0('umap_markers_for_',x,'.pdf'),width = 12)



