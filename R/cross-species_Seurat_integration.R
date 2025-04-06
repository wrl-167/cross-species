library(Seurat)  
library(tidyverse)  

# 读取同源基因表格  
homolog_genes <- read.csv("path/to/your/homologTable.csv")  # 修改为相对路径  
species_A_genes <- homolog_genes$Species_A_symbol  
species_B_genes <- homolog_genes$Species_B_symbol  

# 读取Seurat对象  
species_A_seurat <- readRDS("path/to/species_A_seurat_object.rds")  # 修改为相对路径  
species_B_seurat <- readRDS("path/to/species_B_seurat_object.rds")  # 修改为相对路径  

# 获取人类、物种 A 和物种 B共有基因  
human_genes <- homolog_genes$Human_symbol  # 确保数据集中有一列表示人类基因符号  
common_genes <- intersect(toupper(human_genes), toupper(species_B_genes))  

# 筛选Seurat对象，仅保留交集基因  
seurat_list <- list(species_A_seurat, species_B_seurat)  
seurat_list_filtered <- lapply(seurat_list, function(seurat_obj) {  
  subset(seurat_obj, features = intersect(rownames(seurat_obj), common_genes))  
})  

#去除重复基因  
remove_duplicate_genes <- function(seurat_obj) {  
  gene_names <- rownames(seurat_obj)  
  expression_data <- GetAssayData(seurat_obj, slot = "data")  
  
  duplicated_indices <- which(duplicated(gene_names) | duplicated(gene_names, fromLast = TRUE))  
  
  if (length(duplicated_indices) > 0) {  
    df <- data.frame(gene = gene_names, mean_expr = rowMeans(expression_data), stringsAsFactors = FALSE)  
    duplicate_genes <- unique(gene_names[duplicated_indices])  
    genes_to_keep <- unique(duplicate_genes[match(duplicate_genes, gene_names)])  
    non_duplicate_genes <- df$gene[!df$gene %in% duplicate_genes]  
    final_genes_to_keep <- c(genes_to_keep, non_duplicate_genes)  
  } else {  
    final_genes_to_keep <- gene_names  
  }  
  
  seurat_obj_filtered <- subset(seurat_obj, features = final_genes_to_keep)  
  
  cat("原始基因数：", length(gene_names), "\n")  
  cat("保留的基因数：", length(final_genes_to_keep), "\n")  
  cat("移除的重复基因数：", length(gene_names) - length(final_genes_to_keep), "\n")  
  
  return(seurat_obj_filtered)  
}  

# 去重并清理Seurat对象  
seurat_list_filtered_unique_cleaned <- lapply(seurat_list_filtered, remove_duplicate_genes)  

# 打印清理后的对象信息  
lapply(seurat_list_filtered_unique_cleaned, function(obj) {  
  cat("对象信息：\n")  
  cat("细胞数：", ncol(obj), "\n")  
  cat("基因数：", nrow(obj), "\n")  
  cat("行名（基因名）：\n")  
  print(head(rownames(obj), 10))  
  cat("行名是否唯一：", length(rownames(obj)) == length(unique(rownames(obj))), "\n")  
  cat("\n")  
})  

# 合并清理后的Seurat对象  
unintegrated <- merge(  
  seurat_list_filtered_unique_cleaned[[1]],   
  y = seurat_list_filtered_unique_cleaned[[2]],   
  add.cell.ids = c("Species_A", "Species_B")  
)  

# 打印合并结果信息  
cat("合并后的细胞数：", ncol(unintegrated), "\n")  
cat("合并后的基因数：", nrow(unintegrated), "\n")  

# 对数据进行标准化处理  
unintegrated <- NormalizeData(unintegrated)  

# 找到高可变特征  
unintegrated <- FindVariableFeatures(  
  unintegrated,   
  selection.method = "vst",   
  nfeatures = length(features)  # 使用所有有效基因  
)  

# 对数据进行缩放 (归一化)  
unintegrated <- ScaleData(  
  unintegrated,   
  features = features  # 使用所有有效基因  
)  

# 运行主成分分析 (PCA)  
unintegrated <- RunPCA(  
  unintegrated,   
  features = features,  # 使用所有有效基因  
  npcs = 50,  # 提取前50个主成分  
  verbose = TRUE  # 打印信息  
)  

# 可视化 PCA 结果  
pca_plot <- DimPlot(unintegrated, reduction = "pca", group.by = "orig.ident")  
print(pca_plot)  # 打印 PCA 结果  

# 使用 UMAP 进行降维和可视化  
unintegrated <- RunUMAP(  
  unintegrated,   
  dims = 1:50  # 基于前50个主成分计算 UMAP，数量可根据需要调整  
)  

# 可视化 UMAP 结果  
umap_plot <- DimPlot(unintegrated, reduction = "umap", group.by = "orig.ident")

#################### 应用不同整合方法
######----------(1) Seurat CCA 整合--------######  

runtime_cca <- system.time({  
  integrated_cca<- IntegrateLayers(  
    object = unintegrated,  
    method = CCAIntegration,  # 利用 CCA 方法进行整合  
    orig.reduction = "pca",  
    new.reduction = "integrated.cca"  
  ) %>%  
  RunUMAP(reduction = "integrated.cca", dims = 1:50)  
})  

# 记录运行时间  
integrated_cca@misc$runtime <- runtime_cca['elapsed']  

# 可视化运行时间  
print(paste("运行时间:", integrated_cca@misc$runtime, "秒"))  
saveRDS(integrated_cca, file = "integrated_cca_result.rds")

p1 <- DimPlot(integrated_cca, reduction = "umap", group.by = "species") +  
    ggtitle("Species_A vs Species_B") +  
    scale_color_manual(values = c("Species_A" = "darkred", "Species_B" = "steelblue"))  
print(p1)

integrated_cca <- FindNeighbors(integrated_cca, reduction = "integrated.cca", dims = 1:30) 
integrated_cca<- FindClusters(integrated_cca, resolution = 0.2)  
p2 <- DimPlot(integrated_cca, reduction = "umap", label = TRUE) +   
  ggtitle("CCA Integration Clusters")  
print(p2)


######----------(2) Seurat RPCA 整合--------######  
runtime_rpca <- system.time({  
  integrated_rpca <- IntegrateLayers(  
    object = unintegrated,  
    method = RPCAIntegration,  
    orig.reduction = "pca",  
    new.reduction = "integrated.rpca",  
    verbose = TRUE  
  ) %>%  
    RunUMAP(reduction = "integrated.rpca", dims = 1:30)  
})    

integrated_rpca@misc$runtime <- runtime_rpca['elapsed']  

# 可视化运行时间  
print(paste("运行时间:", integrated_rpca@misc$runtime, "秒"))   
saveRDS(integrated_rpca, file = "integrated_rpca_result.rds")
p1 <- DimPlot(integrated_rpca, reduction = "umap", group.by = "species") +  
    ggtitle("Species_A vs Species_B") +  
    scale_color_manual(values = c("Species_A" = "darkred", "Species_B" = "steelblue"))  
print(p1)  
 
integrated_rpca <- FindNeighbors(integrated_rpca,   
                                reduction = "integrated.rpca",   
                                dims = 1:30)  
 
integrated_rpca <- FindClusters(integrated_rpca,   
                              resolution = 0.2,  
                              verbose = TRUE)  

p2 <- DimPlot(integrated_rpca, reduction = "umap", label = TRUE) +   
  ggtitle("RPCA Integration Clusters")  
print(p2)
 
 
  ######--------(3) Harmony 整合--------###### 
# 使用 Harmony 进行整合并记录运行时间  
runtime_harmony <- system.time({  
  integrated_harmony <- IntegrateLayers(  
    object = unintegrated,  
    method = HarmonyIntegration,  # 改用 Harmony 方法进行整合  
    orig.reduction = "pca",  
    new.reduction = "harmony",    # 更改降维结果名称  
    verbose = TRUE                # 显示进度信息  
  ) %>%  
  RunUMAP(reduction = "harmony", dims = 1:30)  # 注意这里使用 harmony 作为降维结果  
})  

# 记录运行时间  
integrated_harmony@misc$runtime <- runtime_harmony['elapsed']  

# 可视化运行时间  
print(paste("Harmony 整合运行时间:", integrated_harmony@misc$runtime, "秒")) 
saveRDS(integrated_harmony, file = "integrated_harmony_result.rds")
p1 <- DimPlot(integrated_harmony, reduction = "umap", group.by = "species") +  
    ggtitle("Species_A vs Species_B") +  
    scale_color_manual(values = c("Species_A" = "darkred", "Species_B" = "steelblue"))  
print(p1) 

integrated_harmony <- FindNeighbors(integrated_harmony,   
                                 reduction = "harmony",   
                                 dims = 1:30)  

integrated_harmony <- FindClusters(integrated_harmony,   
                                resolution = 0.4,  
                                verbose = TRUE)
p2 <- DimPlot(integrated_harmony,   
             reduction = "umap",   
             label = TRUE) +   
  ggtitle("Harmony Integration Clusters")
p2  
