library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final   
library(Seurat)
DimPlot(sce,label = T,repel = T) 
as.data.frame(table(Idents(sce)))

pbmc = pbmc3k.final
# PROGENy （Pathway RespOnsive GENes for activity inference）
# 是2018年发表在Nature Communication的R包
# 参考：https://mp.weixin.qq.com/s/5QSOaHvz__xoMBRuljKuSg

library(progeny)
pbmc <- progeny(pbmc, scale=FALSE, organism="Human", top=500, perm=1, 
                return_assay = TRUE)

## We can now directly apply Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
pbmc <- Seurat::ScaleData(pbmc, assay = "progeny") 

## We transform Progeny scores into a data frame to better handling the results
library(tidyverse)
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(pbmc, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 
head(progeny_scores_df)

## We match Progeny scores with the cell clusters.
CellsClusters <- data.frame(Cell = names(Idents(pbmc)), 
                            CellType = as.character(Idents(pbmc)),
                            stringsAsFactors = FALSE)
head(CellsClusters)
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)

## We summarize the Progeny scores by cellpopulation
summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))

## We prepare the data for the plot
summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 

paletteLength = 100
myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))
progeny_hmap = pheatmap::pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)
dev.off()
progeny_hmap


av <-AverageExpression(sce , 
                       assays = "RNA")
av=av[[1]] 
cg=names(tail(sort(apply(av, 1, sd)),1000)) 
pheatmap::pheatmap(cor(av[cg,])) 
head(av)
dim(av)
gene_expression= av
# 计算通路活性
pathways <- progeny(gene_expression, scale=TRUE,
                    organism="Human",  
                    top = 100, perm = 1)
head(pathways)[1:5,1:5]  
library(pheatmap)
myColor = colorRampPalette(c("Darkblue", "white","red"))(100)
pheatmap::pheatmap(t(pathways),fontsize=14, show_rownames = T,
         color=myColor, main = "PROGENy", angle_col = 45, treeheight_col = 0,  
         border_color = NA)








