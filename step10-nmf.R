library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final   

library(Seurat)
table(Idents(sce))
p1=DimPlot(sce,label = T);p1

sub_sce = sce[,grepl('Mono',Idents(sce))]
dim(sub_sce)
table(Idents(sub_sce))
sub_sce$celltype = Idents(sub_sce)
p2=DimPlot(sub_sce,label = T)
library(patchwork)
p1+p2



sub_sce=CreateSeuratObject(
  counts = sub_sce@assays$RNA@counts,
  meta.data = sub_sce@meta.data
)
library(dplyr)
sub_sce = NormalizeData(sub_sce) %>% FindVariableFeatures() %>% ScaleData(do.center = F)

suppressPackageStartupMessages(library(NMF))
vm <- sub_sce@assays$RNA@scale.data
res <- nmf(vm,rank=2,method = "snmf/r") 

# 默认交替最小二乘法(Alternating Least Squares(ALS))——snmf/r  
# 参数rank=2，是期望的细胞亚群数量

plot(t(coef(res) ),col=sub_sce$celltype)
head(basis(res))

# 期望单细胞分成2群，拿出来每个单细胞亚群各自的特征基因
fs <- extractFeatures(res,10L)
fs <- lapply(fs,function(x)rownames(res)[x])
names(fs)=paste0('gp',1:2) 

library(ggplot2) 
th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
DotPlot(sce[,grepl('Mono',Idents(sce))],
        features =   fs)  + coord_flip() +th


## 选择后续分析的因子，使用 NMF 运行的结果进行降维和聚类 
sub_sce
sub_sce <- RunPCA(sub_sce,verbose = F)
sub_sce@reductions$nmf <- sub_sce@reductions$pca
sub_sce@reductions$nmf@cell.embeddings <- t(coef(res) )
sub_sce@reductions$nmf@feature.loadings <- basis(res) 
sub_sce <- RunUMAP(sub_sce,reduction = "nmf",dims = 1:2) 

sub_sce <- FindNeighbors(sub_sce,reduction = "nmf",dims = 1:2) %>% FindClusters(resolution = 0.1)

sub_sce$cluster <- apply(NMF::coefficients(res)[1:2,],2,which.max)

table( sub_sce$celltype ,sub_sce$cluster)
table(Idents(sub_sce) ,sub_sce$cluster)

Idents(sub_sce) <- sub_sce$cluster
table(Idents(sub_sce))
DimPlot(sub_sce,reduction = 'umap',label = T)


