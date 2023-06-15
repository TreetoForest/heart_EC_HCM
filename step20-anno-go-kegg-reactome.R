library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final  
library(Seurat)
table(Idents(sce))
DimPlot(sce,label = T)
input_sce = sce
table(Idents(input_sce))
pro = 'cosg_seurat_clusters'
library(COSG)
marker_cosg <- cosg(
  input_sce,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=100)

save(marker_cosg,file = paste0(pro,'_marker_cosg.Rdata'))
head(marker_cosg)
 

symbols_list <-  as.list(as.data.frame(apply(marker_cosg$names,2,head,100)))
source('/share/home/linmiao/bio/2023_Ying_PAH_scRNA/PAH_scRNA/1script/com_go_kegg_ReactomePA_mice.R')
com_go_kegg_ReactomePA_mice(symbols_list, pro='pbmc' )


