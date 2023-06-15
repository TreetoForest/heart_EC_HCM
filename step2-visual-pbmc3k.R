sce=pbmc
sce$celltype=Idents(sce)

p1=FeaturePlot(sce,'CD4')
p2=DimPlot(sce, reduction = "umap", 
        label = TRUE, repel = T,pt.size = 0.5) + NoLegend()
p3=VlnPlot(sce,'CD4',group.by = 'celltype')
library(patchwork)
p1+p2
p1+p3 

pos=sce@reductions$umap@cell.embeddings
pos=pos[sce@assays$RNA@counts['CD4',]>1,]
head(pos)
library(ggplot2)
p2+geom_point(aes(x=UMAP_1,y=UMAP_2), 
              shape = 21, colour = "black",
              fill = "blue", size = 0.5,  
              data = as.data.frame(pos))

