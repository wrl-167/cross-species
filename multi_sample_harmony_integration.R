setwd("your/project/directory")

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(devtools)



data_dir <- "data/"
dir_name <- list.files(data_dir)
# 读取数据
scRNAlist <- list()
for(i in 1:length(dir_name)){
  counts <- Read10X(data.dir = paste("data/", dir_name[i], sep = ""))
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = dir_name[i],min.cells = 3, min.features = 300)
}



####计算线粒体比例####
for(i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]]
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")#计算以"MT-"开头的基因比例，具体视样本而定
  scRNAlist[[i]] <- sc
  # 删除sc
  rm(sc)
}



####批量绘制质控前小提琴图####
violin_before <- list()
for(i in 1:length(scRNAlist)){
  violin_before[[i]] <- VlnPlot(scRNAlist[[i]],
                                features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                                pt.size = 0.01, 
                                ncol = 4) 
}
# 合并图片
violin_before_merge <- CombinePlots(plots = violin_before,nrow=length(scRNAlist),legend='none')
# 将图片输出到画板上
violin_before_merge
# 保存图片
ggsave("dir_name", plot = violin_before_merge, width = 15, height =7)

violin_before

####批量过滤####
scRNAlist <- lapply(X = scRNAlist, FUN = function(x){
  x <- subset(x, 
              subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & 
                mt_percent < 10 & 
                nCount_RNA < quantile(nCount_RNA,0.97) & 
                nCount_RNA > 1000)})

View(scRNAlist[[1]]@meta.data)
#没有固定的阈值标准


####merge合并样本####
scRNAlist <- merge(x=scRNAlist[[1]],y=scRNAlist[-1])
## 统计细胞数
table(scRNAlist[[]]$orig.ident)

##过滤后可视化
violin_after <- VlnPlot(scRNAlist,
                        features = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent"), 
                        pt.size = 0.01,
                        ncol = 4)
# 将图片输出到画板上 
violin_after
# 保存图片
ggsave("dir_name", plot = violin_after, width = 15, height =7) 


####数据归一化、筛选高变基因与PCA降维####
scRNAlist <- NormalizeData(scRNAlist) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30, verbose = T)#npcs：计算和存储的PC数（默认为 50）
a=DimPlot(scRNAlist,reduction = "pca",group.by = "orig.ident")
a


# 提取前15个高变基因ID
top15 <- head(VariableFeatures(scRNAlist), 15) 
plot1 <- VariableFeaturePlot(scRNAlist) 
plot2 <- LabelPoints(plot = plot1, points = top15, repel = TRUE, size=3) 
# 合并图片
feat_15 <- CombinePlots(plots = list(plot1,plot2),legend = "bottom")
feat_15
# 保存图片
ggsave(file = "dir_name",plot = feat_15,he = 10,wi = 15 )


# Seurat是5.0版在要运行下这个代码才能进行差异分析(Seurat V4可忽略此代码)：
scRNAlist <- JoinLayers(scRNAlist)#连接layers层的count数据




####RunHarmony去批次####
# 整合需要指定Seurat对象和metadata中需要整合的变量名。
scRNA_harmony <- RunHarmony(scRNAlist, group.by.vars = "orig.ident")
scRNA_harmony@reductions[["harmony"]][[1:5,1:5]]
b=DimPlot(scRNA_harmony,reduction = "harmony",group.by = "orig.ident")
b
# 合并图片
pca_harmony_integrated <- CombinePlots(list(a,b),ncol=1)
pca_harmony_integrated



####聚类、umap/tsne降维降维####
ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony")
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:10) %>% FindClusters(resolution = 0.05)#分辨率可以自己调
##umap/tsne降维
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:10)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:10)

# 绘图
umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")
umap_integrated2 <- DimPlot(scRNA_harmony, reduction = "umap", label = TRUE)
tsne_integrated1 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = "orig.ident") 
tsne_integrated2 <- DimPlot(scRNA_harmony, reduction = "tsne", label = TRUE)
# 合并图片
umap_tsne_integrated <- CombinePlots(list(tsne_integrated1,tsne_integrated2,umap_integrated1,umap_integrated2),ncol=2)
# 将图片输出到画板
umap_tsne_integrated
# 保存图片
ggsave("dir_name",umap_tsne_integrated,wi=25,he=15)
#保存数据
save(scRNA_harmony,scRNAlist,file = "dir_name")
load("dir_name")

table(scRNA_harmony@meta.data$seurat_clusters)



####细胞注释####


# 差异分析
markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)



all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val_adj<0.05)

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
View(top10)
write.csv(top10,"dir_name",row.names = T)


# 可视化maker
DoHeatmap(scRNA_harmony, features = top10$gene, slot="data") + NoLegend()#slot默认使用scaledata里边只有2k个高变gene表达 这里要使用data数据，不然有些gene会找不到表达量
VlnPlot(scRNA_harmony,features = top10$gene[1:10])#选前10个makergene看看,这里选的都是cluster0的差异基因

VlnPlot(scRNA_harmony,features = c()





View(scRNA_harmony@meta.data)
table(scRNA_harmony@meta.data$seurat_clusters)


####找细胞marker####
library(ggplot2) 
p <- DotPlot(scRNA_harmony, features = top10$gene[1:10],
             assay='RNA' ,group.by = 'seurat_clusters' ) + coord_flip()+ggtitle("")

genes_to_check = c()














