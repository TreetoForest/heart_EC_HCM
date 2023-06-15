library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final  

library(Seurat)
cg=c('PIP','FABP4','ACKR1','TM4SF1','GNG11',
     'CLDN5','SELE','AQP1','PLVAP','RAMP2','KRT14',
     'KRT6A','KRT5','KRT16','KRT1','LYZ','40 LTB',
     'S100A9','LY6D','S100A8','DMKN','S100A7','CCL18',
     'HLA.DRA','C1QA','CXCL8','C1QB','AIF1','C1QC',
     'IL32','CD7','CXCR4','CD69','DUSP2','TRBC2',
     'NKG7','CCL5','TRBC1','CFD','TRAC','DCN','APOD',
     'COL1A2','COL1A1','SFRP2','PTGDS','LUM','COL3A1',
     'FBLN1','IGKC','MS4A1','CD79A','CD83','GPR183',
     'IGLC2','TPSB2','CTSG','HPGD','HPGDS','RGS13',
     'CMA1','GATA2','VWA5A','ANXA1','ACTA2','MYL9',
     'TAGLN','TPM2','MYH11','CCL2','CALD1','MYLK',
     'TPM1','S100B','PLP1','PMP22','CLU','RGS5',
     'GPM6B','PMEL','CAPN3','MITF','QPCT',
     'NOV','KRT19','MUCL1','SPARCL1','NDRG2',
     'GAPDH','CRYAB','PEBP1','CNN3','TYRP1',
     'DCT','KRT7','AQP5','MLANA','IGHM','TPSAB1',
     'AZGP1','DCD','HLA-DPA1','HLA-DPB1','JCHAIN',
     'CHCHD6','SCGB2A2','SCGB1B2P','SCGB1D2')
library(stringr)  
library(ggplot2)  
p <- DotPlot(sce, features = cg,
             assay='RNA'  ) +theme(axis.text.x = element_text(angle = 90))

p

library(Seurat)
genes_to_check = c("CD14",'PTPRC','CD68','FCGR3A')
FeaturePlot(sce,genes_to_check)
DimPlot(sce,label = T,repel = T)
FeaturePlot(sce,c("CD14" ,'FCGR3A'),blend = T)


library(stringr)  
library(ggplot2)  
p <- DotPlot(sce, features = genes_to_check,
             assay='RNA'  ) +theme(axis.text.x = element_text(angle = 90))

p

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


library(stringr)
cg = str_to_upper(
  c(
    'CD3D', 'CD3E', 'CD4','CD8A',
    'CD19', 'CD79A', 'MS4A1' 
  )
)
cg
FeaturePlot(sce,features = cg)
DoHeatmap(sce,features = cg,size = 3)
gl = list(
  Tcells =  cg[1:4],
  Bcells =  cg[5:7]
)
gl
T_mat = sce@assays$RNA@counts[gl[[1]],]
B_mat = sce@assays$RNA@counts[gl[[2]],]

sce$t_sum = colSums(T_mat) 
p1=FeaturePlot(sce,'t_sum') 
sce =  AddModuleScore(object = sce,features = gl[1])
colnames(sce@meta.data)
p2=FeaturePlot(sce,'Cluster1')
library(patchwork)
p1+p2

table(colSums(T_mat) > 1 ,
      colSums(B_mat) > 1 )

gs=list(
  DC1 = c( 'Clec9a', 'Xcr1',   'Wdfy4'), 
  DC2 = c('Itgax', 'Sirpa',   'Cd209a'), 
  mregDCs= c('Ccr7', 'Cd80', 'Cd200',   'Cd247') ,
  hypoxia=c('Hif1a', 'Slc2a1', 'Vegfa', 'Hmox1', 
            'Bnip3', 'Nos2', 'Mmp2', 'Sod3', 
            'Cited2', 'Ldha')
)
gs = lapply(gs, toupper)
sce =  AddModuleScore(object = sce,gs)
colnames(sce@meta.data)
FeaturePlot(sce,'Cluster1')
VlnPlot(sce,'Cluster4')
ncol(sce@meta.data)
ac=sce@meta.data[,4,drop=F]
dat= sce@meta.data[,8:11]
colnames(dat) = names(gs)
pheatmap::pheatmap(dat,
                   show_rownames = F,
                   annotation_row = ac)


p=VlnPlot(sce,'Cluster4')
library(ggpubr)
df = aggregate(p$data$Cluster4,list(p$data$ident),median)
ggbarplot(df,'Group.1','x') + coord_flip()


th=theme(axis.text.x = element_text(angle = 45, 
                                    vjust = 0.5, hjust=0.5))  
myeloids = list(
  Mac=c("C1QA","C1QB","C1QC","SELENOP","RNASE1","DAB2","LGMN","PLTP","MAF","SLCO2B1"),
  mono=c("VCAN","FCN1","CD300E","S100A12","EREG","APOBEC3A","STXBP2","ASGR1","CCR2","NRG1"),
  neutrophils =  c("FCGR3B","CXCR2","SLC25A37","G0S2","CXCR1","ADGRG3","PROK2","STEAP4","CMTM2" ),
  pDC = c("GZMB","SCT","CLIC3","LRRC26","LILRA4","PACSIN1","CLEC4C","MAP1A","PTCRA","C12orf75"),
  DC1 = c("CLEC9A","XCR1","CLNK","CADM1","ENPP1","SNX22","NCALD","DBN1","HLA-DOB","PPY"),
  DC2=c( "CD1C","FCER1A","CD1E","AL138899.1","CD2","GPAT3","CCND2","ENHO","PKIB","CD1B"),
  DC3 =  c("HMSD","ANKRD33B","LAD1","CCR7","LAMP3","CCL19","CCL22","INSM1","TNNT2","TUBB2B")
)
p <- DotPlot(sce , features = myeloids,
             assay='RNA'  )  +th

p
ggsave(plot=p, filename="check_myeloids_marker_by_celltype.pdf") 

