library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final  
library(Seurat)
table(Idents(sce))
DimPlot(sce,label = T)
av <-AverageExpression(sce , 
                       assays = "RNA")
av=av[[1]] 
cg=names(tail(sort(apply(av, 1, sd)),1000)) 
pheatmap::pheatmap(cor(av[cg,])) 

head(av)
dim(av)

library(GSVA) 
library(GSEABase)
library(msigdbr)
library(clusterProfiler)

all_gene_sets = msigdbr(species = "Homo sapiens",
                        category='H') 
msigdbr_collections()
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category = "C2", subcategory =   "CP:REACTOME"  ) 
# geneList = av[,1] 
gl = apply(av, 2, function(geneList){
  geneList=sort(geneList,decreasing = T)
  print(head(geneList))
  print(tail(geneList))
  geneList=geneList[geneList>0.1]
  egmt <- GSEA(geneList, TERM2GENE= all_gene_sets[,c('gs_name','gene_symbol')] , 
               minGSSize = 20, 
               pvalueCutoff = 1,
               verbose=FALSE)
  head(egmt)
  egmt@result 
  gsea_results_df <- egmt@result 
  return(gsea_results_df) 
  
})
lapply(gl, dim)
tmp=gl[[1]]
path = unique(all_gene_sets$gs_name)
es.max <- do.call(cbind,
        lapply(gl, function(x){
          x[path,'enrichmentScore']
        }))
rownames(es.max) = path
head(es.max)  
es.max=na.omit(es.max)
pheatmap::pheatmap(es.max,show_colnames =T,show_rownames = F) 

#每个单细胞亚群的特异性top5基因集的 富集分析结果
library(dplyr)
df = do.call(rbind,
             lapply(1:ncol(es.max), function(i){
               dat= data.frame(
                 path  = rownames(es.max),
                 cluster =   colnames(es.max)[i],
                 sd.1 = es.max[,i], #每一列原值
                 sd.2 = apply(es.max[,-i], 1, median)  #除当列以外每行的中位值
               )
             })) 
df$fc = df$sd.1 - df$sd.2#两值相减，变化越大说明越有意义（从中挑出top5）
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)#找出每个细胞类型的前五个通路
n=es.max[unique(top5$path),]
rownames(n)
rownames(n)[grepl('FGFR',rownames(n))]
rownames(n)=gsub('REACTOME_','',rownames(n))
rownames(n)=substring(rownames(n),1,30)
pheatmap::pheatmap(n,show_rownames = T) 

 