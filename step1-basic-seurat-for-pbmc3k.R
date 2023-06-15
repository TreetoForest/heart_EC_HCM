library(Seurat)
# https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
## Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

## Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", 
                           min.cells = 3, min.features = 200)
# ## Identification of mithocondrial genes
# pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# 
# ## Filtering cells following standard QC criteria.
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
#                  percent.mt < 5)

## Normalizing the data
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
phe=pbmc@meta.data
save(phe,file = 'phe-by-basic-seurat.Rdata')

# pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25,  logfc.threshold = 0.25, verbose = FALSE)
DimPlot(pbmc, reduction = "umap", group.by = 'seurat_clusters',
        label = TRUE, pt.size = 0.5) 
DotPlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", 
                           "CD14", "FCER1A", "FCGR3A", 
                           "LYZ", "PPBP", "CD8A"),
        group.by = 'seurat_clusters')
## Assigning cell type identity to clusters
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T",
                     "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

pbmc$cluster_by_counts=Idents(pbmc)
table(pbmc$cluster_by_counts)





