# 1.加载R包 ####
rm(list = ls())
library(Seurat)
library(tidyverse)
library(ggsci)
library(pagoda2) 
library(Matrix) 
library(conos)
library(dplyr)
library(entropy)
library(conosPanel)

library(help="igraph");help("igraph") 

# remotes::install_github('kharchenkolab/conosPanel')
# install.packages('conosPanel', repos='https://kharchenkolab.github.io/drat/', type='source')
# install.packages("conos")

# sudo mkdir -p /Volumes/Builds/packages/high-sierra-x86_64/Rlib/4.1/igraph/libs/
# sudo cp  /Library/Frameworks/R.framework/Versions/4.1/Resources/library/igraph/libs/igraph.so /Volumes/Builds/packages/high-sierra-x86_64/Rlib/4.1/igraph/libs/ 
  
# 2.设置路径 ####
getwd()
dir.create("step7_pagoda2")
setwd("step7_pagoda2/")

# 3.加载数据 ####
## 3.1 测试数据 ----
library(conosPanel)
panel <- conosPanel::panel 
# panel是一个List，包含4个单细胞样本的表达量稀疏矩阵
lapply(panel, dim)

cm <- panel[[1]]
str(cm)
par(mfrow=c(1,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0)
hist(log10(colSums(cm)+1),main='molecules per cell',col='cornsilk',xlab='log10(molecules per cell)')
hist(log10(rowSums(cm)+1),main='molecules per gene',col='cornsilk',xlab='log10(molecules per gene])')
counts <- gene.vs.molecule.cell.filter(cm,min.cell.size=500)
hist(log10(rowSums(counts)+1),main='Molecules per gene',xlab='molecules (log10)',col='cornsilk')
abline(v=1,lty=2,col=2)


## 3.2 真实数据 ----

library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final  
library(Seurat)
table(Idents(sce))

cm=sce@assays$RNA@counts
str(cm) 
par(mfrow=c(1,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0)
hist(log10(colSums(cm)+1),main='molecules per cell',col='cornsilk',xlab='log10(molecules per cell)')
hist(log10(rowSums(cm)+1),main='molecules per gene',col='cornsilk',xlab='log10(molecules per gene])')

hist(log10(rowSums(counts)+1),main='Molecules per gene',xlab='molecules (log10)',col='cornsilk')
abline(v=1,lty=2,col=2)

counts <- counts[rowSums(counts)>=10,]
dim(counts) 
rownames(counts) <- make.unique(rownames(counts))

r <- Pagoda2$new(counts,log.scale=TRUE, n.cores=2) 
r$adjustVariance(plot=T,gam.k=10) 
r$calculatePcaReduction(nPcs=50,n.odgenes=3e3) 
r$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine') 
r$getKnnClusters(method=infomap.community,type='PCA')

M <- 30; 
r$getEmbedding(type='PCA',embeddingType = 'largeVis', 
               M=M,perplexity=30,gamma=1/M,alpha=1)
r$plotEmbedding(type='PCA',show.legend=F,
                shuffle.colors=F,
              alpha=0.1 ) + ggtitle('clusters (largeVis)')


r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=F,n.cores=30)
r$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,
                min.group.size=1,
                shuffle.colors=F,mark.cluster.cex=1,
                alpha=0.1,main='clusters (tSNE)')

gene <-"HBB" 
r$plotEmbedding(type='PCA',embeddingType='tSNE',
                colors=r$counts[,gene],shuffle.colors=F,
              alpha=0.1,main=gene)
 
### 3.可视化uncertainty结果及最终注释结果 ####
# 此函数返回属于每个组的每个单元格的数据集中的概率、不确定性分数和最终标签
con$plotPanel(colors=new.label.info$uncertainty, 
              show.legend=TRUE,
              legend.title="Uncertainty",
              legend.pos=c(1, 0))
ggsave("cellannot_all_sample_uncertainty.pdf")

con$plotPanel(groups=new.label.info$labels, show.legend=FALSE)
ggsave("cellannot_all_sample_labble.pdf") 



r$getDifferentialGenes(type='PCA',verbose=T,clusterType='community')

de <- r$diffgenes$PCA[[1]][['1']]; 
r$plotGeneHeatmap(genes=rownames(de)[1:15],groups=r$clusters$PCA[[1]])

gene <-"IL32" 
r$plotEmbedding(type='PCA',embeddingType='tSNE',
                colors=r$counts[,gene],shuffle.colors=F,
              alpha=0.1,main=gene)

r$getHierarchicalDiffExpressionAspects(type='PCA',
                                       clusterName='community',z.threshold=3)

# app <- p2.make.pagoda1.app(r, inner.clustering=TRUE,
#                            embeddingType='tSNE',
#                            clusterType='community' , 
#                            row.clustering=list(order=rev(1:nrow(r$misc$pathwayOD$xv))) )
# show.app(app,'pbmc',browse=T)



## 3.3 数据检查 ----
str(panel,1) 
head(colnames(panel[[1]])) 
# 快速检查细胞名称是否唯一
any(duplicated(unlist(lapply(panel,colnames)))) 
# [1] FALSE

## 3.4 数据预处理 ----

### 3.4.1.用 pagoda2 进行预处理 ####
set.seed(1)
panel.preprocessed <- lapply(panel, basicP2proc, 
                             n.cores=1, 
                             min.cells.per.gene=0, # 不删除低表达基因
                             n.odgenes=2e3, 
                             get.largevis=FALSE, 
                             make.geneknn=FALSE)

typeof(panel.preprocessed)
# [1] "list"
names(panel.preprocessed)
# [1] "MantonBM1_HiSeq_1" "MantonBM2_HiSeq_1" "MantonCB1_HiSeq_1"
# [4] "MantonCB2_HiSeq_1"

### 3.4.2.用 Seurat 进行预处理 ####
library(Seurat)
panel.preprocessed.seurat <- lapply(panel, basicSeuratProc)

# 这里可以是之前创建好的每个样本构成Seurat对象list。
# 可以通过Seurat包中的函数将矩阵过滤，去除细胞周期的影响。

# 4.构建Conos对象 ####
con <- Conos$new(panel.preprocessed, n.cores=1) 
# 使用conos生成图形可以利用并行处理，根据个人电脑进行配置
con

## 4.1 单样本可视化(Joint Graph) ####
space='PCA' # PCA, CPCA,CCA
con$buildGraph(k=30, k.self=5, 
               space=space,  # PCA, CPCA,CCA
               ncomps=30, 
               n.odgenes=2000, 
               matching.method='mNN', 
               metric='angular', 
               score.component.variance=TRUE, 
               verbose=TRUE)

plotComponentVariance(con, space=space) 
# space 应与buildGraph函数中的space相对应
ggsave(paste0("plotComponentVariance_",space,".pdf"))

# # 计算结果储存在con$pairs$PCA
# con$pairs$PCA <- NULL # 清空缓存

###1. plotPanel ####
con$plotPanel(clustering="multilevel",
              use.local.clusters=TRUE, 
              title.size=6)
ggsave(paste0("plotPanel_",space,".pdf"))


###2. leiden.community法分群 ####
resolution = 1 # 可以适当修改分群
con$findCommunities(method=leiden.community, resolution=resolution) # 相当于Seurat包中的FindClusters函数
con$plotPanel(font.size=4) # 绘图
table(con$clusters$leiden$groups) 
ggsave(paste0("plotPanel_",resolution,".pdf"))


# 多样本展示
table(con$clusters$leiden$groups)
con$embedGraph(method='largeVis')
con$plotGraph(clustering='leiden')
ggsave(paste0("plotGraph_leiden.pdf"),width = 6,height = 6)
con$plotGraph(color.by='sample', mark.groups=FALSE, alpha=0.1, show.legend=TRUE)
ggsave(paste0("plotGraph_leiden_by_sample.pdf"),width = 6,height = 6)

###3. plotClusterBarplots ####
plotClusterBarplots(con, legend.height = 0.1) 
# 可视化每个cluster不同样本的占比
ggsave(paste0("plotClusterBarplots_",resolution,".pdf"))

###4.gene expression ####
# 可视化基因表达情况，相当于FeaturePlot
# gene_to_check <- c("GZMK","MKI67","TOP2A") # 可检查已知的maker gene
# for (gene in gene_to_check) {
#   gene=gene_to_check[1]
#   print(gene)
#   con$plotPanel(gene = gene,title=paste0(gene,' expression'))
#   ggsave(paste0("plotPanel_gene_expression_",gene,".pdf"),
#          width = 6,height = 6)
# }

## 4.2 数据中存在NA ####
# 如果存在NA，可视化将不会出现

## create new variable for panel.preprocessed
preprocessed_panel_example <- panel.preprocessed

## set NAs within all cell sin MantonCB2_HiSeq_1
preprocessed_panel_example$MantonCB2_HiSeq_1$clusters$PCA$multilevel <- NA

## create new Conos object
con_example <- Conos$new(preprocessed_panel_example, n.cores=1)

## construct joint graph
con_example$buildGraph()

## now plot the panel with NA values as none
con_example$plotPanel(clustering="multilevel", 
                      use.local.clusters=TRUE, 
                      title.size=6,plot.na=TRUE)

## 4.3 多样本可视化 ####

### 1.embedGraph by cluster####
con$embedGraph(method='largeVis') # 嵌入（默认使用'largeVis'）
con$plotGraph(alpha=0.1)
ggsave("allsample_PCA.pdf")

### 2.embedGraph by sample ####
con$plotGraph(color.by='sample', 
              mark.groups=FALSE, 
              alpha=0.1, # 调整点的透明度
              show.legend=TRUE)
### 3.gene expression ####
# 可视化基因表达情况，相当于FeaturePlot
gene_to_check <- c("GZMK","MKI67","TOP2A") # 可检查已知的maker gene
for (gene in gene_to_check) {
  print(gene)
  con$plotGraph(gene = gene,title=paste0(gene,' expression'),alpha=0.5)
  ggsave(paste0("plotGraph_gene_expression_",gene,".pdf"),width = 6,height = 6)
}

### 4.walktrap.community 法分群 ####
step = 7 # 可以适当修改分群建议8~10之间，相当于Seurat包中的FindClusters函数
con$findCommunities(method=igraph::walktrap.community,steps=step) # 相当于Seurat包中的FindClusters函数
con$plotPanel(font.size=4) # 绘图
ggsave(paste0("plotPanel_",step,".pdf"))

# 可视化 --- walktrap法
table(con$clusters$walktrap$groups)
con$plotGraph(clustering='walktrap')
ggsave(paste0("plotGraph_walktrap.pdf"),width = 6,height = 6)

# 可视化 --- 单样本展示
table(con$clusters$walktrap$groups)
con$plotPanel(clustering='walktrap')
ggsave(paste0("plotPanel_walktrap.pdf"),width = 6,height = 6)


##4.4 更改embedding参数 ####
###1. largeVis ####
# 设置embedGraph 主要参数
alpha=0.001 # 值的大小影响分群数
sgd_batched=1e8 #  值的大小影响离散程度
# 计算 
con$embedGraph(method = 'largeVis', # largeVis 为默认计算方法
               alpha=alpha, 
               embedding.name="example_embedding", 
               sgd_batched=sgd_batched)   
table(con$clusters$walktrap$groups)
con$plotGraph(clustering='walktrap', size=0.1)
ggsave(paste("embedGraph_largeVis_",alpha,"_",sgd_batched,".pdf"))

### 2. UMAP ####
# 设置主要参数 
spread = 15
min.dist = 0.01
min.prob.lower = 1e-3
# 计算
con$embedGraph(method="UMAP", 
               min.dist=min.dist, 
               spread=spread, 
               min.prob.lower=min.prob.lower) # 默认min.prob.lower = 1e-3

con$plotGraph(clustering='walktrap', size=0.1)
ggsave(paste("plotGraph_UMAP_",spread,"_",min.dist,".pdf"))

# 单样本展示
con$plotPanel(clustering='walktrap', size=0.1, use.common.embedding=TRUE)
ggsave(paste("plotPanel_UMAP_",spread,"_",min.dist,".pdf"))

### 3.greedyModularityCut ---- 
# 类似于Seurat包中的clustree函数
fc <- greedyModularityCut(con$clusters$walktrap$result, 40)
con$plotGraph(groups=fc$groups, size=0.1)
dend <- as.dendrogram(fc$hc)

pdf(file = "greedyModularityCut_40.pdf")
plot(dend)
dev.off()

## 4.5 细胞亚群注释 ####
### 1.构建1个样本的注释文件 ####
cellannot <- read.table(file.path(find.package('conos'), 
                                  'extdata', 'cellannot.txt'),
                        header=FALSE, sep='\t')
head(cellannot)
table(cellannot)
cellannot <- setNames(cellannot[,2], cellannot[,1])
con$plotPanel(groups = cellannot)
ggsave("cellannot_one_sample.pdf")

### 2.将标签从一个带注释的样本传播到其他样本 ####
new.label.info <- con$propagateLabels(labels = cellannot,
                                      verbose=TRUE)

### 3.可视化uncertainty结果及最终注释结果 ####
# 此函数返回属于每个组的每个单元格的数据集中的概率、不确定性分数和最终标签
con$plotPanel(colors=new.label.info$uncertainty, 
              show.legend=TRUE,
              legend.title="Uncertainty",
              legend.pos=c(1, 0))
ggsave("cellannot_all_sample_uncertainty.pdf")

con$plotPanel(groups=new.label.info$labels, show.legend=FALSE)
ggsave("cellannot_all_sample_labble.pdf") 


con$plotPanel(clustering='walktrap', size=0.1, 
              groups=new.label.info$labels
              , show.legend=FALSE, 
              use.common.embedding=TRUE)
ggsave("cellannot_all_sample_labble_umap.pdf")


head(new.label.info$label.distribution)

# 5.细胞类型&细胞亚群差异表达分析 ####

## 1.指定细胞类型信息 ####
new.annot <- new.label.info$labels

## 2.计算 ####
de.info <- con$getDifferentialGenes(groups=new.annot, append.auc=TRUE)

## 3.查看指定细胞亚群的差异基因结果 ####

# 这里以B cells 举例
head(de.info$`B cells`)

#设置主要参数
new.annot= new.label.info$labels
  gene="CD74"
pdf(paste0("expression_gene",gene,".pdf"))
cowplot::plot_grid(con$plotGraph(groups=new.annot), 
                   con$plotGraph(gene=gene) )
dev.off()

gene="CD79A"
pdf(paste0("expression_gene",gene,".pdf"))
cowplot::plot_grid(con$plotGraph(groups=new.annot), 
                   con$plotGraph(gene=gene) )
dev.off()

# 筛选指定细胞亚群中，AUC > 0.75 的基因，按照Precision排序
# 这里以 monocytes 举例
de.info$monocytes %>% filter(AUC > 0.75) %>% arrange(-Precision) %>% head()
head(de.info$monocytes)

## 4.可视化基因表达情况 ####
con$plotGraph(gene="CD14")

## 5.绘制Top gene表达量热图 ####
# 设置主要参数 
n.genes.per.cluster = 5 

# 绘图
pdf(file =paste0("DEG_Top_marker_gene_",n.genes.per.cluster,".pdf") )
plotDEheatmap(con,as.factor(new.annot),
              de.info,
              n.genes.per.cluster = n.genes.per.cluster, 
              column.metadata=list(samples=con$getDatasetPerCell()), 
              row.label.font.size = 10) # 设置基因字体大小
dev.off()

## 6.绘制指定基因热图 ####
gns <- c("GZMB","IL32","CD3E","LYZ","HLA-DRA","IGHD","GNLY","IGHM","GZMK")
n <- length(gns)
pdf(file =paste0("DEG_gene_marker_",n,".pdf") )
plotDEheatmap(con,new.annot,de.info[-c(3,10)], 
              n.genes.per.cluster = 30, 
              column.metadata=list(samples=con$getDatasetPerCell()), 
              row.label.font.size = 7, 
              labeled.gene.subset = gns)
dev.off()


## 7.使用 Leiden 或 walktrap 绘制 DE聚类热图 ####
### 7.1 Leiden法绘制 DE聚类 热图 ####
# leiden.community 聚类分析
con$findCommunities(method = leiden.community, resolution = 1.0)
# leiden.community 找差异基因
leiden.de <- con$getDifferentialGenes(clustering = "leiden", 
                                      append.auc = TRUE, 
                                      groups=con$clusters$leiden$groups)
# 设置主要参数
n.genes.per.cluster=5

# 绘图
pdf(paste0("leiden_DEG_Top_",n.genes.per.cluster,".pdf"))

plotDEheatmap(con = con, 
              groups = as.factor(con$clusters$leiden$groups), 
              de =leiden.de, 
              n.genes.per.cluster = 5, 
              column.metadata=list(samples=con$getDatasetPerCell()), 
              row.label.font.size = 7)
dev.off()


### 7.2 walktrap法绘制 DE聚类 热图 ####
# walktrap 聚类分析
con$findCommunities(method = igraph::walktrap.community, steps = 10)
# walktrap.community 找差异基因
walktrap.de <- con$getDifferentialGenes(clustering = "walktrap", 
                                      append.auc = TRUE, 
                                      groups=con$clusters$leiden$groups)
# 设置主要参数
n.genes.per.cluster=5

# 绘图
pdf(paste0("walktrap_DEG_Top_",n.genes.per.cluster,".pdf"))

plotDEheatmap(con = con, 
              groups = as.factor(con$clusters$leiden$groups), 
              de =walktrap.de, 
              n.genes.per.cluster = 5, 
              column.metadata=list(samples=con$getDatasetPerCell()), 
              row.label.font.size = 7)
dev.off()


# 6.样本组之间的差异表达 ####

##1.查看数据（每个样本的合并计数矩阵）####
str(con$getClusterCountMatrices(), 1)

##2.构建样本分组信息samplegroups ####
samplegroups <- list(
  bm = c("MantonBM1_HiSeq_1","MantonBM2_HiSeq_1"),
  cb = c("MantonCB1_HiSeq_1","MantonCB2_HiSeq_1")
)

##3.计算每个细胞类型在两个分组之间的差异分析 ####
de.info <- getPerCellTypeDE(con, groups=as.factor(new.annot), 
                            sample.groups = samplegroups, 
                            ref.level='bm', # 设置对照组
                            n.cores=1)

## 4.检查输出 
str(de.info[1:3], 2)
## 5.检查B细胞
res <- de.info[['B cells']]$res
head(res[order(res$padj,decreasing = FALSE),])

# 7.更好的聚类####
# 设置主要参数
alignment.strength = 0.3
con$buildGraph(k=15, k.self=5, 
               alignment.strength=alignment.strength, 
               space='PCA',  # CPCA;CCA;PCA
               ncomps=30, 
               n.odgenes=2000, 
               matching.method='mNN', 
               metric='angular', 
               score.component.variance=TRUE, 
               verbose=TRUE)

# 可视化
con$findCommunities()
pdf("plotClusterBarplots_2.pdf")
plotClusterBarplots(con, legend.height = 0.1)
dev.off()







