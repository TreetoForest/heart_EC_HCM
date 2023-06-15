# BiocManager::install('TENxPBMCData') # (286 KB)
library(TENxPBMCData)
args(TENxPBMCData) 

tenx_pbmc4k <- TENxPBMCData(dataset = "pbmc4k")
tenx_pbmc4k
counts(tenx_pbmc4k)


library(scRNAstat) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)

ct = as.data.frame( assay(tenx_pbmc4k))
ct[1:4,1:4]
colnames(ct) = 1:ncol(ct)

library(AnnoProbe) 
ag=annoGene(rownames(ct),
            ID_type = 'ENSEMBL',species = 'human'
)
head(ag)
ag=ag[!duplicated(ag$SYMBOL),]
ag=ag[!duplicated(ag$ENSEMBL),]

pos = match(ag$ENSEMBL,rownames(ct))
ct = ct[pos,]
rownames(ct) = ag$SYMBOL
ct[1:4,1:4]
 
library(scRNAstat) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
sce = CreateSeuratObject(counts =  ct)
x='tenx_pbmc4k'
dir.create(x)
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



