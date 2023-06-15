library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final   
library(Seurat)
table(Idents(sce))
p1=DimPlot(sce,label = T)
sce$celltype = Idents(sce)

sce.all = sce
colnames(sce.all@meta.data)

library(patchwork)
library(ggplot2)
colnames(sce.all@meta.data)
FeaturePlot(sce.all,'LYZ')+DimPlot(sce.all,
                                    group.by = "seurat_clusters",label = T,repel = T)
ggsave('LYZ-and-umap-0.8.pdf',width = 5,height = 8)

av <-AverageExpression(sce.all,
                      group.by = "seurat_clusters",
                       assays = "RNA") 
av=av[[1]]
head(av)
write.csv(av,file = 'AverageExpression-0.8.csv')
cg=names(tail(sort(apply(av, 1, sd)),1000))
pheatmap::pheatmap(cor(av[cg,]),
                   file = 'AverageExpression-0.8.pdf')
dev.off()

# load(file = 'phe-by-markers.Rdata')
# sce.all@meta.data = phe
# sce.all
# table(Idents(sce.all))  
# Idents(sce.all)=sce.all$celltype

table(Idents(sce.all))  
library(patchwork)
FeaturePlot(sce.all,'LYZ')+DimPlot(sce.all,label = T,repel = T)
ggsave('LYZ-and-umap-celltype.pdf',width = 5,height = 8)

av <-AverageExpression(sce.all,
                       # group.by = "celltype",
                       assays = "RNA") 
av=av[[1]]
head(av)
write.csv(av,file = 'AverageExpression-celltype.csv')
cg=names(tail(sort(apply(av, 1, sd)),1000))
pheatmap::pheatmap(cor(av[cg,]),
                   file = 'AverageExpression-celltype.pdf')
dev.off()

library(patchwork)
VlnPlot(sce.all,'LYZ',pt.size = 0) + 
  VlnPlot(sce.all,'LYZ',    group.by = "seurat_clusters",pt.size = 0) 
ggsave('LYZ-VlnPlot.pdf',width = 10,height = 8)


cg= c('ACADS','CDIPT','CYB5R3','CYB5R4','DCXR','G6PD','IDH3G','MLYCD','PFKL','PIK3CA')
cg = cg[cg %in% rownames(sce.all)]
cg

for (g in cg) {
  #g=cg[1]
  library(patchwork)
  p = VlnPlot(sce.all, g, pt.size = 0) + 
    VlnPlot(sce.all,g,    group.by = "seurat_clusters",pt.size = 0) 
  print(p)
  ggsave(paste0('choose-',g,'-VlnPlot.pdf'),width = 10,height = 8) 
}

genes_to_check = cg
pl = lapply(genes_to_check, function(cg){  FeaturePlot(sce.all, cg,) + NoLegend() + NoAxes() })
ps <- cowplot::plot_grid(plotlist = pl)
ps  
ggsave("FeaturePlot_umap.pdf",width = 16,height = 15)

pl = lapply(genes_to_check, function(cg){  FeaturePlot(sce.all, cg,order = T,raster = T) + NoLegend() + NoAxes() })
ps <- cowplot::plot_grid(plotlist = pl)
ps  
ggsave("FeaturePlot_umap2.pdf",width = 16,height = 15)

pro='pbmc3k'

# 基本上代替热图的小提琴图
table(Idents(sce.all))
p_all_markers=DotPlot(sce.all, 
                      group.by  = "celltype",
                      features = genes_to_check,
                      scale = T,assay='RNA' ) + coord_flip() + 
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers
ggplot2::ggsave(paste0(pro,'_all_genes_to_check-DotPlot.pdf'),height = 4,width =5)


p1 <- VlnPlot(sce.all,features =  genes_to_check,
              group.by  = "celltype",
              flip = T,stack = T )
p1 + NoLegend()
ggplot2::ggsave(paste0(pro,'_genes_to_check-VlnPlot-heatmap.pdf'),height = 4,width =5)

sce.Scale <- ScaleData(subset(sce.all,downsample=100),features =  genes_to_check  )  
DoHeatmap(sce.Scale,
          features =  genes_to_check,
          group.by = "celltype",
          assay = 'RNA', label = T)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))

ggsave(paste0(pro,filename = "_genes_to_check-pheatmap.pdf") )





 


