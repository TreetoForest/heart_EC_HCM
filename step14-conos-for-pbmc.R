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

# install.packages('conosPanel', repos='https://kharchenkolab.github.io/drat/', type='source')
# install.packages("conos") 

# 2.加载数据 ####

## 2.1 读取pbmc3k和5k数据集 ----
library(conosPanel)
options(stringsAsFactors = F)
load('pbmc3k.Rdata')
pbmc_3k=pbmc
pbmc_5k=readRDS('pbmc_5k_v3.rds') 
library(Seurat)
panel.preprocessed.seurat <- list(
  pbmc_3k=pbmc_3k,pbmc_5k=pbmc_5k
)

## 2.2 构建Conos对象  ----
con <- Conos$new(panel.preprocessed.seurat, n.cores=1) 
space='PCA' # 可以选择 PCA, CPCA,CCA
con$buildGraph(k=30, k.self=5, 
               space=space,  # PCA, CPCA,CCA
               ncomps=30, 
               n.odgenes=2000, 
               matching.method='mNN', 
               metric='angular', 
               score.component.variance=TRUE, 
               verbose=TRUE)
plotComponentVariance(con, space=space)  

## 2.3 聚类分群   ----
resolution = 1 # 可以适当修改分群
con$findCommunities(method=leiden.community, resolution=resolution) # 相当于Seurat包中的FindClusters函数
con$plotPanel(font.size=4) # 绘图
table(con$clusters$leiden$groups)  
con$embedGraph(method='largeVis')
con$plotGraph(clustering='leiden')


## 2.4 整合后后的效果展示 ----
library(patchwork)
con$plotGraph(clustering='leiden') + con$plotGraph(color.by='sample', mark.groups=FALSE, alpha=0.1, show.legend=TRUE)


# 可视化每个cluster不同样本的占比
plotClusterBarplots(con, legend.height = 0.1) 
con$plotGraph(clustering='leiden', size=0.1)

# 假如需要替换坐标体系
# # 设置主要参数 , 做 UMAP
# spread = 15
# min.dist = 0.01
# min.prob.lower = 1e-3
# # 计算
# con$embedGraph(method="UMAP", 
#                min.dist=min.dist, 
#                spread=spread, 
#                min.prob.lower=min.prob.lower) # 默认min.prob.lower = 1e-3
# # 再次绘图，默认坐标体系变成了 UMAP 
# con$plotGraph(clustering='leiden', size=0.1)


## 4.gene expression ####
# 可视化基因表达情况，相当于FeaturePlot

### 1.embedGraph by cluster#### 
con$plotGraph(alpha=0.1) 
### 2.embedGraph by sample ####
con$plotGraph(color.by='sample', 
              mark.groups=FALSE, 
              alpha=0.1, # 调整点的透明度
              show.legend=TRUE)
### 3.gene expression ####
# 可视化基因表达情况，相当于FeaturePlot

# T Cells (CD3D, CD3E, CD8A), 
# B cells (CD19, CD79A, MS4A1 [CD20]), 
# Plasma cells (IGHG1, MZB1, SDC1, CD79A), 
# Monocytes and macrophages (CD68, CD163, CD14),
# NK Cells (FGFBP2, FCG3RA, CX3CR1),  
library(patchwork)
gene_to_check= c('PTPRC', 'CD3D', 'CD3E','IL7R','CD4','CD8A','CD19', 'CD79A', 'MS4A1')
pl = lapply(gene_to_check,  function(gene){
  con$plotGraph(gene = gene,title=paste0(gene,' expression'),
                alpha=0.5)
})
wrap_plots(pl, byrow = T, nrow = 3)

save(con,file = 'merge_by_conos.Rdata')
table(con$clusters$leiden$groups)  


