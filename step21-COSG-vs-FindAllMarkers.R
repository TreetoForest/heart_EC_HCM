library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final   
library(Seurat)
DimPlot(sce,label = T,repel = T) 
as.data.frame(table(Idents(sce)))


library(future)
# check the current active plan
plan()
plan("multiprocess", workers = 4)
plan()
start = Sys.time()
sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
end = Sys.time() 
dur = end-start 
dur

DT::datatable(sce.markers)
pro='FindAllMarkers'
write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
save(sce.markers,file = paste0(pro, '_sce.markers.Rdata'))


#remotes::install_github(repo = 'genecell/COSGR')
start = Sys.time()
library(COSG)
marker_cosg <- cosg(
  sce,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=10)
end = Sys.time() 
dur = end-start 
dur


findG = split(sce.markers$gene,sce.markers$cluster)
names(findG)
tmp = unlist(  lapply(names(findG), function(x){
  sum(marker_cosg$names[,x] %in% findG[[x]])
}))
names(tmp) = names(findG)
as.data.frame(tmp)

library(stringr)  
library(ggplot2)  
g1 =head( marker_cosg$names$DC  ,10)
g2 = head(findG$DC,10)
th= theme(axis.text.x = element_text(angle = 90))
p1 = DotPlot(sce, features = g1,assay='RNA'  ) + th +NoLegend()
p2 = DotPlot(sce, features = g2,assay='RNA'  ) + th
library(patchwork)
p1+p2




