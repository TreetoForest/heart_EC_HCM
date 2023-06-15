library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final  

library(Seurat)
library(ggplot2)
th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 


genes_to_check= rownames(sce)[grepl('^MS[0-9][AB]',rownames(sce))]
genes_to_check
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_KRT_genes.pdf',width = 10,height = 8)