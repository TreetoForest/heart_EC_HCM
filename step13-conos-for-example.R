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
## 2.1 测试数据 ----
library(conosPanel)
panel <- conosPanel::panel 
# panel是一个List，包含4个单细胞样本的表达量稀疏矩阵
# 而且都是3000个细胞，3万多个基因
lapply(panel, dim)
### 用 Seurat 对4个单细胞样品都进行预处理
library(Seurat)
panel.preprocessed.seurat <- lapply(panel, basicSeuratProc)

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

## 2.3 leiden.community法分群  ----
resolution = 1 # 可以适当修改分群
con$findCommunities(method=leiden.community, resolution=resolution) # 相当于Seurat包中的FindClusters函数
con$plotPanel(font.size=4) # 绘图
table(con$clusters$leiden$groups) 
con$embedGraph(method='largeVis') 

## 2.4 整合后后的效果展示 ----
library(patchwork)
con$plotGraph(clustering='leiden') + con$plotGraph(color.by='sample', mark.groups=FALSE, alpha=0.1, show.legend=TRUE)

# 可视化每个cluster不同样本的占比
plotClusterBarplots(con, legend.height = 0.1) 

