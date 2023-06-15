library(Seurat)
# https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
## Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/") 
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", 
                           min.cells = 3, min.features = 200)  
ct=pbmc@assays$RNA@counts
ct
ct[ct>0]=1 
ct

pbmc@assays$RNA@counts=ct
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 10000) 
pbmc <- NormalizeData(pbmc) 
## Identify the 2000 most highly variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
## In addition we scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), 
               verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindClusters(pbmc, resolution = 0.5, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:10, umap.method = "uwot", metric = "cosine")
table(pbmc$seurat_clusters)
# pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25,  logfc.threshold = 0.25, verbose = FALSE)

library(patchwork)
library(ggplot2)
p1=DimPlot(pbmc, reduction = "umap", group.by = 'seurat_clusters',
        label = TRUE, pt.size = 0.5) 
p2=DotPlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", 
                           "CD14", "FCER1A", "FCGR3A", 
                           "LYZ", "PPBP", "CD8A"),
        group.by = 'seurat_clusters')+theme(axis.text.x = element_text(angle = 45, 
                                                                       vjust = 0.5, hjust=0.5))

p1+p2 
phe=pbmc@meta.data
save(phe,file = 'phe-by-0-1-matrix.Rdata')


load(file = 'phe-by-basic-seurat.Rdata')
phe_basic=phe
load(file = 'phe-by-0-1-matrix.Rdata')
phe_0_1=phe
identical(rownames(phe_0_1),rownames(phe_basic))
library(gplots)
balloonplot(table(phe_basic$seurat_clusters,phe_0_1$seurat_clusters))

phe_0_1$type_by_0_1 = ifelse(phe_0_1$seurat_clusters %in% c(0,1,2,5,8),'Tcells',
       ifelse(phe_0_1$seurat_clusters %in% c(3),'Bcells','myeoloid'
       ))
table(phe_0_1$type_by_0_1)
phe_basic$type_by_basic = ifelse(phe_basic$seurat_clusters %in% c(0,2,4,6),'Tcells',
                             ifelse(phe_basic$seurat_clusters %in% c(3),'Bcells','myeoloid'
                             ))
table(phe_basic$type_by_basic)
table(phe_basic$type_by_basic,phe_0_1$type_by_0_1)





