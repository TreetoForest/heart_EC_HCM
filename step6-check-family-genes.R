library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final  
library(Seurat)

genes_to_check = c("CD4","CD3E","IL7R", "KLF2", "CCR7","TCF7",
                   "SELL", "CCL4", "CCL5", "PRF1",  "GZMB",
                   "GZMK", "FGFBP2", "CX3CR1", "RORC","CXCL13",
                   "CXCR5", "FOXP3", "IL2RA","IL5", "IL1RL1",
                   "GATA3", "PTGDR2")

library(stringr)  
library(ggplot2)  
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  ) +theme(axis.text.x = element_text(angle = 90))

p

th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 
 

DimPlot(sce,reduction = "umap",label=T )  
library(patchwork)
library(Seurat) 
FeaturePlot(sce,features = c('EPCAM','ESR1',
                             'KRT5','KRT14','ACTA2',
                             'MKI67'),
            reduction = 'umap',pt.size = 1,
            min.cutoff = "q9")

library(ggplot2) 
genes_to_check = c(  'PROM1', 'CD44' , 'THY1','ACTA2','MKI67', 'ESR1',
                     'EPCAM', 'KRT19',  'ALDH1A1', 'CD24' )
p1 <- DotPlot(sce, features = genes_to_check, 
             assay='RNA'  )  + coord_flip() +th

p1
p2=DimPlot(sce,label = T,repel = T)
p1+p2
ggsave('umap_dotplot_epi_all_markers.pdf',width = 10,height = 5)

library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'C1QA',  'C1QB',  # mac
                   'S100A9', 'S100A8', 'MMP19',# monocyte
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'KLRB1','NCR1', # NK 
                   'FGF7','MME', 'ACTA2',
                   'PECAM1', 'VWF', 
                   'MKI67','TOP2A',
                   'EPCAM' , 'KRT19', 'PROM1', 'ALDH1A1' )
library(stringr)   
p_all_markers <- DotPlot(sce , features = genes_to_check,
                        assay='RNA'  )  + coord_flip() +th

p_all_markers 
ggsave('p_all_markers.pdf',width = 10,height = 8)

 
# GSE88715 (including 38 TNBC and 38 normal control)
# A total of 949 DEGs were identified in TNBC (469 up regulated genes, and 480 down regulated genes),  
# overexpressed mRNA (ADAM15, BATF, NOTCH3, ITGAX and SDC1)
# underexpressed mRNA (RPL4, EEF1G, RPL3, RBMX and ABCC2)
genes_to_check = c('ADAM15', 'BATF', 'NOTCH3', 'ITGAX' , 'SDC1',
                   'RPL4', 'EEF1G', 'RPL3', 'RBMX' , 'ABCC2')  
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_GSE88715_markers.pdf',width = 10,height = 5)


# 2018-prognostic immunity markers in breast cancer （看17个免疫基因和19个增殖基因）
proliferation_genes ='AURKA, BIRC5, CCNB1, CCNE1, CDC20, CDC6, CENPF, 
CEP55,EXO1, MKI67, KIF2C, MELK, MYBL2, NDC80, ORC6, PTTG1, RRM2, TYMS, UBE2C'
immunity_genes = 'APOBEC3G, CCL5, CCR2, CD2, CD27,CD3D, CD52, CORO1A, CXCL9, GZMA, GZMK, HLA-DMA, IL2RG, LCK, PRKCB,PTPRC, SH2D1A'
genes_to_check=trimws(strsplit(proliferation_genes,',')[[1]])
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_proliferation_genes.pdf',width = 10,height = 5)
genes_to_check=trimws(strsplit(immunity_genes,',')[[1]])
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_immunity_genes.pdf',width = 10,height = 5)

genes_to_check= rownames(sce)[grepl('^S100',rownames(sce))]
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_S100_genes.pdf',width = 10,height = 5)

genes_to_check= rownames(sce)[grepl('^KRT',rownames(sce))]
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_KRT_genes.pdf',width = 10,height = 8)

genes_to_check= rownames(sce)[grepl('^MUC',rownames(sce))]
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_MUC_genes.pdf',width = 10,height = 5)

genes_to_check= rownames(sce)[grepl('^MMP',rownames(sce))]
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_MMP_genes.pdf',width = 10,height = 6)

genes_to_check= rownames(sce)[grepl('^CCL',rownames(sce))]
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_CCL_genes.pdf',width = 10,height = 6)

genes_to_check= rownames(sce)[grepl('^CXCL',rownames(sce))]
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_CXCL_genes.pdf',width = 10,height = 6)

genes_to_check= rownames(sce)[grepl('^CCR',rownames(sce))]
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_CCR_genes.pdf',width = 10,height = 6)

genes_to_check= rownames(sce)[grepl('^CXCR',rownames(sce))]
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_CXCR_genes.pdf',width = 10,height = 6)


#mesenchymal genes (AGER, FN1, MMP2, SNAI2, VIM, ZEB2) 
#epithelial genes (CDH1, CDH3, CLDN4, EPCAM, MAL2, and ST14)
emt='AGER, FN1, MMP2, SNAI2, VIM, ZEB2,CDH1, CDH3, CLDN4, EPCAM, MAL2,   ST14'
genes_to_check=trimws(strsplit(emt,',')[[1]])
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  )  + coord_flip() +th

p  
ggsave('p_EMT_genes.pdf',width = 10,height = 8)


library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final   

library(Seurat)
DimPlot(sce,label = T,repel = T) 

library(ggplot2)

th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5)) 

cg=c('C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')
p1 <- DotPlot(sce, features = cg, 
              assay='RNA'  )  + coord_flip() +th

p1
p2=DimPlot(sce,label = T,repel = T)
p1+p2
ggsave('umap_dotplot_house_keeping.pdf',width = 10,height = 5)


 