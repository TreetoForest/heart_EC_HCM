#### 此流程说明如何产生 SWNE 的结果图（以经典的pbmc3k数据进行演示） ####

## step 0：清空变量，设置环境，加载R包 ####
rm(list = ls())
getOption('timeout')
options(timeout = 10000)

if(F){
  if(!require(remotes)){ install.packages("remotes") }  # If not already installed;
  remotes::install_github("linxihui/NNLM")
  remotes::install_github("yanwu2014/swne")
  
}

suppressPackageStartupMessages(library(NNLM))
suppressPackageStartupMessages(library(swne))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(SeuratData))

## step 1：加载数据 ####
data("pbmc3k")
obj <- pbmc3k.final
DimPlot(obj,reduction = "umap",label = T,label.size = 3,label.box = T,repel = T)&NoLegend()
ggsave("umap.pdf",width = 5,height = 5)

## 大多scRNA-seq数据都使用高变基因进行分析，这里也取出高变基因
## Pull out overdispersed genes as defined by Seurat
var.genes <- VariableFeatures(obj)
length(var.genes)
# [1] 2000

## 取出Seurat定义的细胞类型
## Pull out overdispersed genes as defined by Seurat
cell.clusters <- Idents(obj)
levels(cell.clusters)
 

## step 2：SWNE降维 ####
## 生成SWNE嵌入最简单的方法就是使用RunSWNE函数
## Run SWNE
genes.embed <- c("MS4A1", "GNLY", "CD3E", "CD14",
                 "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
swne.embedding <- RunSWNE(obj, k = 16, genes.embed = genes.embed, 
                          sample.groups = cell.clusters)
 

## 可视化
## Plot SWNE
PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = cell.clusters,
         do.label = T, label.size = 3.5, pt.size = 1, show.legend = F,
         seed = 42)
ggsave("swne.pdf",width = 5,height = 5)

## RunSWNE函数也可以返回一个Seurat对象，将SWNE设置为一种降维方式（同tSNE和UMAP）
# 这样后续分析以及可视化都是可以直接使用 seurat 的各种函数
## Run SWNE
obj <- RunSWNE(obj, k = 16, genes.embed = genes.embed,
               return.format = "seurat") 
DimPlot(obj, reduction = "umap",label = T,label.size = 3,label.box = T,repel = T)&NoLegend()
DimPlot(obj, reduction = "swne",label = T,label.size = 3,label.box = T,repel = T)&NoLegend()
ggsave("swneplot.pdf",width = 5,height = 5)


## step 3：一步步了解SWNE嵌入过程 ####
## 首先，提取counts、标准化和调整基因方差，同时保持标准化矩阵的非负性（非负矩阵分解）
norm.counts <- ExtractNormCounts(obj, obj.type = "seurat", rescale.method = "log")
# calculating variance fit ... using gam Warning message:
#   In pf(exp(df$res), n.obs, n.obs, lower.tail = F, log.p = F) : 产生了NaNs
dim(norm.counts)
# [1] 13714  2638

## 我们使用FindNumFactors函数来确定要使用的最佳因子数量(跟选择合适的PC降维差不多的概念)。
## 对于大型数据集，这个函数可能会很慢，因为它迭代不同的k值，所以一个简单的“hack”就是让k等于重要主成分(PC)的数量。
k.range <- seq(2,20,2) ## Range of factors to iterate over
k.err <- FindNumFactors(norm.counts[var.genes,], k.range = k.range, n.cores = 8, do.plot = T)
ggsave("factors.pdf",width = 6,height = 5)

## 然后运行NMF分解。我们可以使用独立成分分析(ICA)、非负SVD (nnsvd)或完全随机的初始化来初始化NMF。
## ICA被推荐用于大多数数据集。RunNMF的输出是基因load(W)和NMF嵌入(H)的列表。
k <- 16  # 这里我们选择前16个因子
nmf.res <- RunNMF(norm.counts[var.genes,], k = k)
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|

## 我们可以使用Seurat预先计算的共享最近邻(SNN)矩阵，也可以自己重新计算
# pc.scores <- t(GetCellEmbeddings(se.obj, reduction.type = "pca", dims.use = 1:k))
# snn <- CalcSNN(pc.scores)
obj <- FindNeighbors(obj, k = 10, prune.SNN = 1/15)
# Computing nearest neighbor graph
# Computing SNN

## 我们可以通过使用PAGA图删除集群中没有共享统计上显著数量的边缘的细胞之间的边缘来修剪SNN矩阵
snn <- as(obj@graphs$RNA_snn, "dgCMatrix")
knn <- as(obj@graphs$RNA_nn, "dgCMatrix") ## Extract kNN matrix
snn <- PruneSNN(snn, knn, clusters = cell.clusters, qval.cutoff = 1e-3)

## 运行SWNE嵌入：三个关键参数是alpha,exp snn.Exp和n_pull，它们控制因子和相邻单元格如何影响单元格坐标
alpha.exp <- 1.25 # Increase this > 1.0 to move the cells closer to the factors. Values > 2 start to distort the data.
snn.exp <- 0.25 # Lower this < 1.0 to move similar cells closer to each other
n_pull <- 3 # The number of factors pulling on each cell. Must be at least 3.
swne.embedding <- EmbedSWNE(nmf.res$H, SNN = snn, alpha.exp = alpha.exp, snn.exp = snn.exp,
                            n_pull = n_pull)
# Initial stress        : 0.15146
# stress after  10 iters: 0.04553, magic = 0.338
# stress after  20 iters: 0.04276, magic = 0.500
# stress after  30 iters: 0.04268, magic = 0.500

## 先将因子名称设置为空
swne.embedding$H.coords$name <- ""

## 为了帮助解释这些细胞群，让我们挑选一些关键的PBMC基因嵌入
genes.embed <- c("MS4A1", "GNLY", "CD3E", "CD14",
                 "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")

## 由于我们只在overdispered的基因上运行了NMF，我们需要将其余的基因投射到NMF投影上，以获得所有基因的基因load
nmf.res$W <- ProjectFeatures(norm.counts, nmf.res$H, n.cores = 8)

## 将关键的PBMC基因嵌入到可视化图像中，并重新可视化
swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, genes.embed, n_pull = n_pull)

## 让我们用嵌入的关键基因绘制SWNE图。细胞或细胞簇离基因越近，表达水平越高。
## 我们为可复制的簇颜色设置了一种种子，以便每个地块将使用相同的颜色来标记簇。
color.seed <- 42
PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = cell.clusters, do.label = T,
         label.size = 3.5, pt.size = 1, show.legend = F, seed = color.seed)
ggsave("swne_gene.pdf",width = 5,height = 5)

## 我们可以通过将其中一个关键基因的表达叠加到图上来验证嵌入的基因。
gene.use <- "CD8A"
gene.expr <- norm.counts[gene.use,]
FeaturePlotSWNE(swne.embedding, gene.expr, gene.use, alpha.plot = 0.4, label.size = 3.5, pt.size = 1.25)
ggsave("CD8A_gene.pdf",width = 5.5,height = 5)

## 我们也可以做一个t-SNE图来比较
obj <- RunTSNE(obj)
tsne.scores <- Embeddings(obj, "tsne")
PlotDims(tsne.scores, sample.groups = cell.clusters, pt.size = 1, label.size = 3.5, alpha = 0.4,
         show.legend = F, seed = color.seed, show.axes = F)
ggsave("tsne.pdf",width = 5,height = 5)

## 我们也可以使用基因load矩阵来解释这些因子。在这里，我们按基因load提取了每个因子的前3个基因。
## 由于NMF倾向于创建基于部件的数据表示，这些因子通常对应于解释数据的关键生物过程或基因模块
gene.loadings <- nmf.res$W
top.factor.genes.df <- SummarizeAssocFeatures(gene.loadings, features.return = 3)
head(top.factor.genes.df) 

## 可以做一个热图来显示每个因子的最重要的基因
gene.loadings.heat <- gene.loadings[unique(top.factor.genes.df$feature),]
ggHeat(gene.loadings.heat, clustering = "col")
ggsave("factor_heatmap.pdf",width = 5,height = 5)

## 我们可以提取与其他绘图方法(如Monocle)兼容的聚类颜色
color.mapping <- ExtractSWNEColors(swne.embedding, sample.groups = cell.clusters, seed = color.seed)
color.mapping

