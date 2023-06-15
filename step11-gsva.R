library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final   
library(Seurat)
table(Idents(sce))
p1=DimPlot(sce,label = T)

sce$celltype = Idents(sce)
av <-AverageExpression(sce ,
                       group.by = "celltype",
                       assays = "RNA") 
av=av[[1]]
head(av)
write.csv(av,file = 'AverageExpression-0.8.csv')

cg=names(tail(sort(apply(av, 1, sd)),1000))
pheatmap::pheatmap(cor(av[cg,]))
pheatmap::pheatmap(cor(av[cg,]),
                   file = 'AverageExpression-0.8.pdf')
dev.off()

library(msigdbr)
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category='C7')
length(unique(table(all_gene_sets$gs_name)))
tail(sort(table(all_gene_sets$gs_name)))



library(gplots)
library(ggplot2) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA) # BiocManager::install('GSVA')
library(GSEABase)
gs=split(all_gene_sets$gene_symbol,all_gene_sets$gs_name)
gs = lapply(gs, unique)
gsc <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, gs, names(gs)))
gsc
geneset <- gsc
es.max <- gsva(as.matrix( sce@assays$RNA@counts), geneset, 
               mx.diff=FALSE, verbose=FALSE, 
               parallel.sz=8)
save(es.max,file = 'gsva_all_pbmc.Rdata')

X = av 
es.max <- gsva(X, geneset, 
               mx.diff=FALSE, verbose=FALSE, 
               parallel.sz=8)
save(es.max,file = 'gsva_celltype_pbmc.Rdata')

library(pheatmap)
pheatmap(es.max,show_rownames = F)

load('es.max.Rdata')
cg = names(tail(sort(apply(es.max, 1, sd) ),10))
library(pheatmap)
pheatmap(es.max[cg,])
pheatmap(es.max[cg,],show_rownames = F)


library(dplyr) 
df = do.call(rbind,
        lapply(1:ncol(es.max), function(i){
          dat= data.frame(
            path  = rownames(es.max),
            cluster =   colnames(es.max)[i],
            sd.1 = es.max[,i],
            sd.2 = apply(es.max[,-i], 1, median)  
          )
        }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster)
rownames(rowcn) = paste0(1:nrow(rowcn))
n = es.max[top5$path,]
rownames(n) = paste0(1:nrow(rowcn))
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = F)
#  load('es.max.Rdata')
cg = names(tail(sort(apply(es.max, 1, sd) ),10))
library(pheatmap)
pheatmap(es.max[cg,])
pheatmap(es.max[cg,],show_rownames = F)


library(dplyr) 
df = do.call(rbind,
             lapply(1:ncol(es.max), function(i){
               dat= data.frame(
                 path  = rownames(es.max),
                 cluster =   colnames(es.max)[i],
                 sd.1 = es.max[,i],
                 sd.2 = apply(es.max[,-i], 1, median)  
               )
             }))
df$fc = df$sd.1 - df$sd.2
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
rowcn = data.frame(path = top5$cluster)
rownames(rowcn) = paste0(1:nrow(rowcn))
n = es.max[top5$path,]
rownames(n) = paste0(1:nrow(rowcn))
pheatmap(n,
         annotation_row = rowcn,
         show_rownames = F)

# top5$path 需要标记到  pheatmap上面


