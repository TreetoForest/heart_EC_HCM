library(SeuratData) #加载seurat数据集  
getOption('timeout')
options(timeout=10000)
#InstallData("pbmc3k")  
data("pbmc3k")  
sce <- pbmc3k.final  
library(Seurat)
table(Idents(sce))
DimPlot(sce,label = T)

library(future)
# check the current active plan
plan()
plan("multiprocess", workers = 4)
plan()

sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)
DT::datatable(sce.markers)
pro='markers'
write.csv(sce.markers,file=paste0(pro,'_sce.markers.csv'))
save(sce.markers,file = paste0(pro, 'sce.markers.Rdata'))

pro='markers'
load(file = paste0(pro, 'sce.markers.Rdata'))
sce$celltype = Idents(sce)
library(dplyr) 
top3 <- sce.markers %>% group_by(cluster) %>% top_n(3, avg_log2FC)
DoHeatmap(sce ,top3$gene,size=3)
 
sce.all <- ScaleData(sce,features =  top3$gene)  
library(paletteer) 
color <- c(paletteer_d("awtools::bpalette"),
           paletteer_d("awtools::a_palette"),
           paletteer_d("awtools::mpalette"))
unique(sce.all$celltype)
ord = c('Naive CD4 T' ,'Memory CD4 T', 'CD8 T', 'NK', 
        'CD14+ Mono', 'FCGR3A+ Mono' ,'DC',  'B','Platelet')
sce.all$celltype = factor(sce.all$celltype ,levels = ord)
ll = split(top10$gene,top10$cluster)
ll = ll[ord]
rmg=names(table(unlist(ll))[table(unlist(ll))>1])
ll = lapply(ll, function(x) x[!x %in% rmg])
library(ggplot2)
DoHeatmap(sce.all,
          features = unlist(ll),
          group.by = "celltype",
          assay = 'RNA',
          group.colors = color,label = F)+
   scale_fill_gradientn(colors = c("white","grey","firebrick3"))

ggsave(filename = "marker_pheatmap.pdf",units = "cm",width = 36,height = 42)



df = as.data.frame(table(sce.markers$cluster))
colnames(df) = c('celltype','NofDEGs')
sce$NofDEGs = df[match(sce$celltype,df$celltype),2]
library(patchwork)
DimPlot(sce,label = T,repel = T)+FeaturePlot(sce,'NofDEGs')

dat = as.data.frame(sce@reductions$umap@cell.embeddings)
dat$cluster=as.character(Idents(sce))
dat$NofDEGs = sce$NofDEGs
head(dat)
library(ggplot2)
library(cowplot)
library(ggrepel)

save(dat,file = 'dat.Rdata')

my_FeaturePlot = ggplot()+
   geom_point(data=dat,mapping=aes(x=UMAP_1,
                                       y=UMAP_2,
                                       col= NofDEGs),
              size=AutoPointSize(data=dat))+
   scale_colour_gradient2(low="lightgrey", high = "blue",guide = 'colourbar')+
   theme_cowplot() +  
   theme(
      plot.title = element_text(hjust = 0.5),
      legend.title=element_blank()
   ) +ggtitle("NofDEGs")  
my_FeaturePlot

library(patchwork)
my_FeaturePlot
FeaturePlot(sce,'NofDEGs') + my_FeaturePlot

table(dat$cluster)
my_FeaturePlot2 = my_FeaturePlot  +
   stat_ellipse(data=dat,mapping=aes(x= UMAP_1,
                                         y= UMAP_2, group= cluster),
                colour ='black',
                geom = "path",
                linetype = 2,        ###圆圈线的类型
                size=1,              ###圆圈线的粗细
                alpha=0.5)
my_FeaturePlot2

library(tidyverse)
median_df <- dat %>% group_by(cluster) %>%
   summarise(median.1 = median(UMAP_1),
             median.2 = median(UMAP_2)) 
head(median_df) 
my_FeaturePlot2
my_FeaturePlot2 + ggplot2::geom_text(data = median_df,
                                     aes(median.1,
                                         median.2, label =cluster) )

my_FeaturePlot2 + ggrepel::geom_text_repel(data = median_df,
                                           aes(median.1,
                                               median.2, label =cluster),
                                           min.segment.length = 0)

library(dplyr) 
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(sce,top10$gene,size=3)
ggplot2::ggsave(filename=paste0(pro,'_sce.markers_heatmap.pdf'),height = 15)


as.data.frame(table(sce.markers$cluster))
deg_list=split(sce.markers$gene,
               sce.markers$cluster)
library(UpSetR)
data <- fromList(deg_list)
upset(data,nsets = 9)

pdf('upset_total.pdf',width = 10,height = 8)
upset(data, 
      #sets = c("Action", "Adventure"),#查看特定的几个集合
      mb.ratio = c(0.75, 0.25),#控制上方条形图以及下方点图的比例
      order.by = "freq", #如何排序,这里freq表示从大到小排序展示
      keep.order = TRUE, #keep.order按照sets参数的顺序排序
      #number.angles = 30, #调整柱形图上数字角度
      matrix.color="red", #交集点颜色
      point.size = 2, line.size = 1, #点和线的大小
      mainbar.y.label = "Genre Intersections", 
      sets.x.label = "Set Size", #坐标轴名称
      text.scale = c(1.3, 1.3, 1, 1, 1.5, 1)) #六个数字,分别控制c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
dev.off()


