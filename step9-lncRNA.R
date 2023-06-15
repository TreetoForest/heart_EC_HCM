library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final   
table(Idents(sce))

ct=sce@assays$RNA@counts 
ct[1:4,1:4]
dim(ct)
library(AnnoProbe) 
ids=annoGene( rownames(ct),'SYMBOL','human')
head(ids) 
tail(sort(table(ids$biotypes)))

lncGenes= unique(ids[ids$biotypes=='lncRNA',1])
pbGenes = unique(ids[ids$biotypes=='protein_coding',1])

sce=PercentageFeatureSet(sce, features = lncGenes ,
                         col.name = "percent_lncGenes")
fivenum(sce@meta.data$percent_lncGenes)

C=sce@assays$RNA@counts
dim(C)
C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
# 这里的C 这个矩阵，有一点大，可以考虑随抽样
C=C[,sample(1:ncol(C),1000)]
boxplot(as.matrix(Matrix::t(C[lncGenes, ])),
        cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
        col = (scales::hue_pal())(50)[50:1], 
        horizontal = TRUE)


#可视化细胞的上述比例情况 
feats <- c("nFeature_RNA", "nCount_RNA","percent_lncGenes")
p1=VlnPlot(sce, group.by = "orig.ident", 
           features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
p1
library(ggplot2) 
ggsave(filename="Vlnplot1.pdf",plot=p1) 


sce = CreateSeuratObject(
  counts = ct[pbGenes,] 
)
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce )
sce <- ScaleData(sce ) 
sce <- RunPCA(sce,  verbose = FALSE)
sce <- FindNeighbors(sce, dims = 1:10, verbose = FALSE)
sce <- FindClusters(sce, resolution = 0.5, verbose = FALSE)
sce <- RunUMAP(sce, dims = 1:10, umap.method = "uwot", metric = "cosine")
table(sce$seurat_clusters)
umap_pbGenes= DimPlot(sce,label = T,repel = T) 
sce_pbGenes = sce

sce <- pbmc3k.final    
sce = CreateSeuratObject(
  counts = ct[lncGenes,] 
)
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce )
sce <- ScaleData(sce ) 
sce <- RunPCA(sce, features = VariableFeatures(object = sce), 
              verbose = FALSE)
sce <- FindNeighbors(sce, dims = 1:10, verbose = FALSE)
sce <- FindClusters(sce, resolution = 0.5, verbose = FALSE)
sce <- RunUMAP(sce, dims = 1:10, umap.method = "uwot", metric = "cosine")
table(sce$seurat_clusters)
sce_lncGenes = sce
umap_lncGenes = DimPlot(sce,label = T,repel = T) 

library(patchwork)
umap_lncGenes + umap_pbGenes + 
identical(colnames(sce_lncGenes),colnames(sce_pbGenes))
library(gplots)
balloonplot(table(
  sce_lncGenes$seurat_clusters,
  Idents(pbmc3k.final)
))
balloonplot(table(
  sce_pbGenes$seurat_clusters,
  Idents(pbmc3k.final)
))


rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
# BiocManager::install('scRNAseq')
library(scRNAseq)
fluidigm <- ReprocessedFluidigmData()
fluidigm
ct <- floor(assays(fluidigm)$rsem_counts)
ct[1:4,1:4] 
dim(ct)
ct=ct[rowSums(ct>1)>10,]
dim(ct)
ct[1:4,1:4]
dim(ct)
library(AnnoProbe) 
ids=annoGene( rownames(ct),'SYMBOL','human')
head(ids) 
tail(sort(table(ids$biotypes)))

lncGenes= unique(ids[ids$biotypes=='lncRNA',1])

sce = CreateSeuratObject(
  counts = ct[lncGenes,] 
)
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce )
sce <- ScaleData(sce ) 
sce <- RunPCA(sce, features = VariableFeatures(object = sce), 
              verbose = FALSE)
sce <- FindNeighbors(sce, dims = 1:10, verbose = FALSE)
sce <- FindClusters(sce, resolution = 0.5, verbose = FALSE)
sce <- RunUMAP(sce, dims = 1:10, umap.method = "uwot", metric = "cosine")
table(sce$seurat_clusters)
sce_lncGenes = sce
umap_lncGenes = DimPlot(sce,label = T,repel = T) 
umap_lncGenes

table(
  colData(fluidigm)$Biological_Condition,
  sce$seurat_clusters
)

