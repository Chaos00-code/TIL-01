ibrary(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
#-----------------------------加载数据-----------------------------------------------------------#
load("sce_all.RData")
sce.all
sce.all@meta.data[1:3,]
table(sce.all$antigen)
table(sce.all$mhc_allele)
table(substr(sce.all$barcode,1,5))

sce.all$sample = substr(sce.all$barcode,1,5)

sce.all
table(sce.all$sample)
table(sce.)
dev.off()
gc()
#doublet finder
install.packages("devtools")
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

rmdoubletmat <- function(input.mat, min.features = 50){
  library(DoubletFinder)
  
  obj_tmp <- CreateSeuratObject(input.mat, min.features = min.features)
  obj_tmp <- NormalizeData(obj_tmp, verbose = F)
  obj_tmp <- FindVariableFeatures(obj_tmp, selection.method = "vst", nfeatures = 2000, verbose = F)
  obj_tmp <- ScaleData(obj_tmp, verbose = F)
  obj_tmp <- RunPCA(obj_tmp, npcs = 20, verbose = F)
  obj_tmp <- FindNeighbors(obj_tmp, dims = 1:20, verbose = F)
  obj_tmp <- FindClusters(obj_tmp, resolution = 0.5, verbose = F)
  obj_tmp <- RunUMAP(obj_tmp, dims = 1:20, verbose = F)
  library(DoubletFinder)
  sweep.res.list_kidney <- paramSweep(obj_tmp, PCs = 1:10, sct = FALSE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- Idents(obj_tmp)
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(obj_tmp@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  
  obj_tmp <- doubletFinder(obj_tmp, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  pANN.col.name<-colnames(obj_tmp@meta.data)[6]
  obj_tmp <- doubletFinder(obj_tmp, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = pANN.col.name, sct = FALSE)
  
  
  ## Plot results --------------------------------------------------------------------------------------------------------------
  obj_tmp@meta.data[,"DF_hi.lo"] <- obj_tmp@meta.data[,7]
  obj_tmp@meta.data$DF_hi.lo[which(obj_tmp@meta.data$DF_hi.lo == "Doublet" & obj_tmp@meta.data[,8] == "Singlet")] <- "Doublet_lo"
  obj_tmp@meta.data$DF_hi.lo[which(obj_tmp@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
  # DimPlot(obj_tmp, reduction = "umap", label = TRUE, pt.size = 1, group.by = "DF_hi.lo")
  print(table(obj_tmp$DF_hi.lo))
  
  cells_singlet = row.names(subset(obj_tmp@meta.data, obj_tmp@meta.data$DF_hi.lo == "Singlet"))
  
  output.mat = input.mat[,cells_singlet]
  return(output.mat)
}

exp1 = subset(sce.all, sample == "Z2__1")@assays$RNA@counts
exp2 = subset(sce.all, sample == "Z2__2")@assays$RNA@counts


exp1 = rmdoubletmat(exp1)
exp2 = rmdoubletmat(exp2)

exp_merge = cbind(exp1, exp2)

meta_merge = sce.all@meta.data[colnames(exp_merge),]

#checksample
library(harmony)
library(data.table)
library(dplyr)
library(ggplot2)
obj = CreateSeuratObject(exp_merge, meta.data = meta_merge)
# 安装和加载必要的包
library(Seurat)
library(dplyr)
library(patchwork)
#------------------------------------------------------整合scRNA-seq和CITE-seq----------------------------------#
# 确保细胞名唯一，避免冲突
obj <- RenameCells(obj, add.cell.id = "scRNA")
sce.cite.all <- RenameCells(sce.cite.all, add.cell.id = "CITEseq")

# 将 ADT 数据从 sce.cite.all 添加到 obj 中
adt.data <- GetAssayData(sce.cite.all, assay = "RNA", layer = "counts")
# 检查两个对象的细胞名
head(Cells(obj))         # obj 的细胞名
head(Cells(sce.cite.all)) # sce.cite.all 的细胞名
# 去除细胞名前缀，统一 obj 和 sce.cite.all 的细胞名
obj <- RenameCells(obj, new.names = sub("^scRNA_", "", Cells(obj)))
sce.cite.all <- RenameCells(sce.cite.all, new.names = sub("^CITEseq_", "", Cells(sce.cite.all)))

# 再次检查细胞名是否一致
common_cells <- intersect(Cells(obj), Cells(sce.cite.all))
length(common_cells)  # 应该大于 0

# 根据公共细胞名子集化
obj <- subset(obj, cells = common_cells)
sce.cite.all <- subset(sce.cite.all, cells = common_cells)

# 获取 ADT 数据
adt.data <- GetAssayData(sce.cite.all, assay = "RNA", layer = "counts")

# 将 ADT 数据作为新的 Assay 添加到 obj
obj[["ADT"]] <- CreateAssayObject(counts = adt.data)

print(Assays(obj)) 
length(Cells(obj))
length(Cells(sce.cite.all))

# RNA 分析流程：标准化、特征选择、PCA
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# ADT 分析流程：标准化、特征选择、PCA
DefaultAssay(obj) <- "ADT"
VariableFeatures(obj) <- rownames(obj[["ADT"]]) # 使用所有 ADT 特征
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2) %>%
  ScaleData() %>%
  RunPCA(reduction.name = "apca")

# 计算加权最近邻（WNN）邻居图
# 检查 RNA PCA 的 elbow point
#ElbowPlot(obj, reduction = "pca")
# 检查 ADT PCA 的 elbow point
#ElbowPlot(obj, reduction = "apca")
obj <- FindMultiModalNeighbors(
  obj, 
  reduction.list = list("pca", "apca"),  # 使用 RNA 和 ADT 的 PCA
  dims.list = list(1:30, 1:18),         # 分别选择 PCA 的前 20 和 18 个成分
  modality.weight.name = "RNA.weight"  # RNA 模态权重
)

# 使用 WNN 进行降维 (UMAP)
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

# 使用 WNN 图聚类
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

# 可视化 WNN UMAP 的聚类结果
p1 <- DimPlot(obj, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5,cols = cluster_cols,pt.size = 0.8) + NoLegend()
p1

# 确保 Seurat 对象中存在细胞分群信息（通常存储在 meta.data 中的某列）
# 假设分群信息存储在 "seurat_clusters" 列中
obj_filtered <- subset(obj, idents = c("28", "29"), invert = TRUE)

# 再次绘制 DimPlot，不包含 cluster28 和 cluster29
p1 <- DimPlot(obj_filtered, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5,cols = cluster_cols,pt.size = 0.8) + NoLegend()
p1

table(obj_filtered$seurat_clusters)
dev.off()
gc()
# 对 RNA 和 ADT 分别运行 UMAP 降维
obj <- RunUMAP(obj_filtered, reduction = "pca", dims = 1:30, assay = "RNA", reduction.name = "rna.umap", reduction.key = "rnaUMAP_")
obj <- RunUMAP(obj_filtered, reduction = "apca", dims = 1:18, assay = "ADT", reduction.name = "adt.umap", reduction.key = "adtUMAP_")

# 可视化 RNA 和 ADT 的降维结果
p2 <- DimPlot(obj, reduction = "wnn.umap", label = TRUE, repel = TRUE, label.size = 5,cols = cluster_cols,pt.size = 0.5)# + NoLegend()
p3 <- DimPlot(obj, reduction = "adt.umap", label = TRUE, repel = TRUE, label.size = 5,cols = cluster_cols) + NoLegend()
p2 + p3
pdf('5_DimPlot_wnnUMAP_seurat_clusters.pdf', width = 10, height = 8)
p2
dev.off()
gc()
# 特征表达可视化：WNN UMAP 上的 RNA 和 ADT 特征
p4 <- FeaturePlot(obj, features = c("CD4.1", "CD8A.1"),
                  reduction = "wnn.umap", max.cutoff = 2, cols = c("lightgrey", '#23452F'), ncol = 2,label = F,pt.size = 0.8)


p5 <- FeaturePlot(obj, features = c("rna_CD4", "rna_CD8A"), 
                  reduction = "wnn.umap", max.cutoff = 3,cols = c("lightgrey",'#CC3333'), ncol = 2,label = F,pt.size = 0.8)
pdf('5_Dimplot_UMAP_cd4_cd8.pdf',width = 14,height = 12)
p4 / p5
dev.off()
gc()

#--------------------------对scRNA-seq数据降维分析并CD4/CD8分群可视乎------------------------------------------#     

DefaultAssay(obj) <- "RNA"
all.markers  <- FindAllMarkers(obj, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.75)
significant.markers  <- all.markers [all.markers $p_val_adj < 0.2, ]
write.csv(significant.markers, file = "5_obj_significant.markers.csv")


DoHeatmap(obj, 
          features = markers,
          group.by = "seurat_clusters",
          assay = "RNA")
table(obj$seurat_clusters)
new.cluster.ids <- c("0"="CD8 T", 
                     "1"="CD8 T", 
                     "2"="CD8 T", 
                     "3"="CD8 T", 
                     "4"="CD8 T", 
                     "5"="CD8 T", 
                     "6"="CD4 T", 
                     "7"="CD8 T", 
                     "8"="CD8 T", 
                     "9"="CD4 T", 
                     "10"="CD4 T", 
                     "11"="CD8 T", 
                     "12"="CD8 T", 
                     "13"="CD8 T", 
                     "14"="CD8 T", 
                     "15"="CD8 T", 
                     "16"="CD8 T", 
                     "17"="CD8 T", 
                     "18"="CD8 T", 
                     "19"="CD4 T", 
                     "20"="CD8 T", 
                     "21"="CD8 T",
                     "22"="CD8 T",
                     "23"="CD4 T", 
                     "24"="CD8 T", 
                     "25"="CD8 T",
                     "26"="CD8 T",
                     "27"="CD4 T" 
                     )
obj <- RenameIdents(obj, new.cluster.ids)     
# save general cell type names to metadata and relevel based on order in figures
levels(obj) <- c("CD4 T","CD8 T")
obj$celltype <- obj@active.ident
pdf('5_Dimplot_UMAP_rename_CD4_CD8.pdf',width = 10,height = 8)
DimPlot(obj, group.by = "celltype",cols =c('#23452F','#CC3333'),pt.size = 0.8 )
dev.off()
gc()
save(obj, file = "5_obj_CD4_CD8_rename.RData")
dev.off()
gc()

#------------------------------细胞亚群分析----------------------------------------
#method1传统
# 提取 CD8 T 亚群
DefaultAssay(obj) <- "RNA"
CD8_T <- subset(obj, subset = celltype == "CD8 T")

CD8_T <- ScaleData(CD8_T, verbose = FALSE)
CD8_T <- FindVariableFeatures(CD8_T, selection.method = "vst", nfeatures = 4000)
CD8_T <- RunPCA(CD8_T, npcs = 50, verbose = FALSE)
CD8_T <- FindNeighbors(CD8_T, dims = 1:50)
CD8_T <- FindClusters(CD8_T, resolution = 0.9)
CD8_T <- RunUMAP(CD8_T, dims = 1:50)
pdf('6_Dimplot_UMAP_CD8_2.pdf',width = 6,height = 6)
DimPlot(CD8_T, reduction = "umap", label = TRUE, pt.size = 1.2,cols = my36colors)+NoLegend()
dev.off()
gc()




CD8_T_cellmarker <- c("MKI67",'STMN1','TUBB',#CD8_T_STNM1
                       'LAG3','PDCD1','FASLG',#CD8_T_PDCD1
                       'ARG2',"ALDOC",'HILPDA',#CD8_T_ARG2
                       'CTLA4','IL2RA','WDR74',#CD8_T_IL2RA
                       'CCT2','ILF2',#CD8_T_ILF2
                       'IFIT3','IFIT1','RSAD2',#CD8_T_IFIT1
                       'SELL','ANKRD55','KLF2'#CD8_T_SELL
)

CD8_T_cellmarker_24_common_genes<-c('CXCL13','HMOX1','PDCD1','LAYN','CD27','HAVCR2','TNFRSF9','MIR155HG','BATF','TIGIT',
                      'GZMH','CD70','TMEM121','LRRN3','NHS','TTN','LINC01281','ASB2','SIRPG','ANKS1B',
                      'HLA-DR','LGALS3','KLRG1','PTPRC-CD45RO','PTPRC-CD45RA','ENTPD1','TIGIT.1','CD39_FBC','IL7R') 
#CD39=ENTPD1,
cols = c("gray", "#23452F")
pdf('6_featureplot_CD8_T_cellmarker_24_common_genes.pdf',width = 30,height = 30)
FeaturePlot(CD8_T,CD8_T_cellmarker_24_common_genes,ncol = 4,pt.size = 1,cols = cols)
dev.off()
gc()


library(ggplot2)
pdf('6_Dotplot_CD8_T_cellmarker_xu.pdf',width = 12,height = 17)
DotPlot(CD8_T, features = CD8_T_cellmarker)+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=0.5,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))
dev.off()
gc()

cols = c("gray", "#23452F")
pdf('6_featureplot_CD8_T_cellmarker_SCIENCE_genes.pdf',width = 12,height = 9)
FeaturePlot(CD8_T,c("PRF1","ENTPD1","PDCD1","TIGIT","ITGAE","TOX","CXCR6","IFNG","LINC01871"),cols = cols)
dev.off()
gc()

# 加载 ComplexHeatmap
library(ComplexHeatmap)
library(circlize)

# 绘制热图
pdf('6_Heatmap_CD8_T_cellmarker_xu.pdf',width = 12,height = 17)
Heatmap(
  heatmap_matrix,
  name = "Expression", 
  cluster_rows = FALSE, 
  cluster_columns = TRUE,
  row_names_side = "left",
  column_names_side = "top",
  show_heatmap_legend = TRUE
)
dev.off()
gc()

pdf(file="6_Vlnplot_CD8_T_cellmarker_xu.pdf",width=20,height=16)
VlnPlot(CD8_T, features = CD8_T_cellmarker,pt.size=0,ncol=5)
dev.off()
#FeaturePlot(CD8_T,"TIGIT",pt.size = 1)



DefaultAssay(CD8_T) <- "RNA"
CD8_T_all.markers  <- FindAllMarkers(CD8_T, 
                                     only.pos = TRUE, 
                                     min.pct = 0.25, 
                                     logfc.threshold = 0.75)
write.csv(CD8_T_all.markers, file = "5_CD8_T_all.markers.csv")

#CD8_T <- subset(CD8_T, idents = c("1","8","9"), invert = TRUE)
#new.cluster.ids <- c("0"="Macrophage", 
#                     "2"="T cell", 
#                     "3"="Macrophage", 
#                     "4"="mDC", 
#                     "5"="Neutrophil", 
#                     "6"="Macrophage", 
#                     "7"="Macrophage", 
#                     "10"="Mast")
#CD8_T <- RenameIdents(CD8_T, new.cluster.ids)                        
#CD8_T$celltype <- CD8_T@active.ident
#DimPlot(CD8_T, label = T,pt.size = 1,group.by = "celltype")
##----------------------------------------ssGSEA注释CD8------------------------------##

# 加载必要的库
install.packages('readr')
BiocManager::install('GSVA',force = TRUE)
library(GSVA)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(pheatmap)
library(tibble)
# 加载数据
signature_genes <- read.csv("9_publication signature_roseberge.csv", header = TRUE)
cd8_markers <- CD8_T_all.markers


# 加载必要的 R 包
library(GSVA)
library(dplyr)
library(tidyr)
library(readr)

# 查看数据结构
head(signature_genes)
head(cd8_markers)

# 转换为基因表达矩阵
expr_matrix <- cd8_markers %>%
  select(gene, cluster, avg_log2FC) %>%  # 提取基因、cluster 和表达值
  pivot_wider(names_from = cluster, values_from = avg_log2FC) %>%
  column_to_rownames(var = "gene")

# 填充 NA 值为 0
expr_matrix[is.na(expr_matrix)] <- 0
head(expr_matrix)
class(expr_matrix)
# 确保矩阵格式
expr_matrix <- as.matrix(expr_matrix)
# 将 signature_genes 转换为列表
signature_list <- signature_genes %>%
  gather(key = "signature_name", value = "gene") %>%  # 宽转长格式
  group_by(signature_name) %>%                       # 按 signature_name 分组
  summarise(genes = list(gene)) %>%                  # 汇总为基因列表
  deframe()                                          # 转换为命名列表

gsvaPar<-gsvaParam(expr_matrix,
                   signature_list,
                   kcdf = "Gaussian")
gsvaPar
gsva.es<-gsva(gsvaPar,verbose=FALSE)
dim(gsva.es)
gsva.es[1:5,1:5]
write.csv(gsva.es,'9_gsvapar_signature_CD8t_expr.csv')

gsvassGSEAPar<-ssgseaParam(exprData = as.matrix(expr_matrix),
                           geneSets = signature_list,
                           normalize = TRUE)
ssgsea_matrix<-gsva(gsvassGSEAPar,verbose=FALSE)
dim(ssgsea_matrix)
ssgsea_matrix[1:5,1:5]
write.csv(ssgsea_matrix,'9_ssgsea_matrix_signature_CD8t_expr.csv')
#-------------------------------------NMF方法分析 CD8+T BEAM_Positive-scRNASEQ-----------------------------------------#              
# Step 4: 区分 Positive/Negative Antigen
table(CD8_T$antigen)
CD8_T$log_antigen_umi <- log(CD8_T$antigen_umi + 1)
quantile(CD8_T$antigen_specificity_score)
FeaturePlot(CD8_T,features = c('antigen_specificity_score'),pt.size = 0.5,reduction = 'umap',raster = F,cols = c("grey97","red")) & theme(aspect.ratio=1)
dev.off()
gc()

# 1. 创建数据框
dfdf <- data.frame(
  antigen_specificity_score = CD8_T$antigen_specificity_score,
  antigen = CD8_T$antigen
)

cluster_antigen <- function(df) {
  print(df)  # 输出当前分组的数据，便于调试
  if (nrow(df) < 2) {
    df$cluster <- ifelse(df$antigen_specificity_score > 0, "positive", "negative")
    df$boundary <- NA
    return(df)
  }
  
  if (length(unique(df$antigen_specificity_score)) < 2) {
    df$cluster <- "positive"
    df$boundary <- unique(df$antigen_specificity_score)
    return(df)
  }
  
  if (nrow(df) == 2) {
    threshold <- mean(df$antigen_specificity_score)
    df$cluster <- ifelse(df$antigen_specificity_score > threshold, "positive", "negative")
    df$boundary <- threshold
    return(df)
  }
  
  kmeans_result <- kmeans(df$antigen_specificity_score, centers = 2)
  df$cluster <- ifelse(kmeans_result$cluster == which.max(kmeans_result$centers), "positive", "negative")
  df$boundary <- mean(kmeans_result$centers)
  return(df)
}


# 3. 分组聚类
library(dplyr)
dfdf_clustered <- dfdf %>%
  group_by(antigen) %>%
  group_modify(~ cluster_antigen(.x)) %>%
  ungroup()

# 4. 将聚类结果添加到原始数据中
CD8_T$antigen_positive <- dfdf_clustered$cluster

#方法 2：改用其他聚类方法（如阈值分割）
#如果部分分组数据较少，kmeans 可能不适合。可以尝试基于阈值的分割方法。例如：
#cluster_antigen <- function(df) {
#  threshold <- mean(df$antigen_specificity_score)  # 使用均值作为分界点
#  df$cluster <- ifelse(df$antigen_specificity_score > threshold, "positive", "negative")
#  df$boundary <- threshold
#  return(df)
#}

# Step 5: 可视化
# UMAP 聚类图
pdf('7_DimPlot_Positive_BEAM_UMAP_boundary_NMF.pdf', width = 15, height = 12)
DimPlot(CD8_T, reduction = "umap", label = F, pt.size = 0.5, cols = c("grey90", "red"),
        split.by = "antigen",
        group.by = "antigen_positive", ncol = 5, raster = F) & theme(aspect.ratio = 1)
dev.off()

boundaries <- dfdf_clustered %>%
  group_by(antigen) %>%
  summarise(boundary = first(boundary))

pdf('7_Histogram_beam_antigen_specificity_score——boundary_NMF.pdf', width = 12, height = 8)
ggplot(dfdf_clustered, aes(x = antigen_specificity_score, fill = cluster)) +
  geom_histogram(binwidth = 0.2, color = NA, alpha = 0.6, position = "identity") +
  geom_vline(data = boundaries, aes(xintercept = boundary), color = "red", linetype = "dashed", size = 1) +
  facet_wrap(~ antigen, scales = "free") +
  scale_y_continuous(trans = "pseudo_log", breaks = c(0, 1, 10, 100, 1000)) +
  scale_fill_manual(values = c("negative" = "skyblue", "positive" = "red")) +
  labs(
    title = "Antigen Specificity Score Clustering (Positive vs Negative) by Antigen",
    x = "Antigen Specificity Score",
    y = "Frequency",
    fill = "Cluster"
  ) +
  theme_minimal()
dev.off()
gc()


library(dplyr)
library(ggplot2)

# 确保 "cluster" 和 "antigen_positive" 列存在于 Seurat 对象的 meta.data
CD8_T$seurat_clusters <- Idents(CD8_T)  # 如果 cluster 信息存储为 Idents

# 统计每个 cluster 中 Positive 的百分比
cluster_stats <- CD8_T@meta.data %>%
  group_by(seurat_clusters, antigen_positive) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(percent = count / sum(count) * 100)

# 打印统计结果
print(cluster_stats)

# 绘制百分比柱状图
pdf('7_Histogram_cluster_beam_positive_percent_CD8_T.pdf', width = 12, height = 6)
ggplot(cluster_stats, aes(x = seurat_clusters, y = percent, fill = antigen_positive)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = c("positive" = "red", "negative" = "gray")) +
  labs(
    title = "Percentage of Positive and Negative Cells in Each Cluster",
    x = "Cluster",
    y = "Percentage",
    fill = "Antigen Positive"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
gc()
pdf('7_Histogram_cluster_beam_positive_percent25_CD8_T.pdf', width = 12, height = 6)
ggplot(cluster_stats, aes(x = seurat_clusters, y = percent, fill = antigen_positive)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = c("positive" = "red", "negative" = "gray")) +
  labs(
    title = "Percentage of Positive and Negative Cells in Each Cluster",
    x = "Cluster",
    y = "Percentage",
    fill = "Antigen Positive"
  ) +
  scale_y_continuous(limits = c(0, 25)) +  # 设置纵坐标范围
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
gc()
# 绘制分布图按样本分面
# 确保 seurat_clusters 和 antigen_positive 列存在于元数据
cluster_stats <- CD8_T@meta.data %>%
  group_by(seurat_clusters, antigen_positive, antigen) %>%  # 包括 antigen 列
  summarise(count = n(), .groups = "drop") %>%
  group_by(seurat_clusters, antigen) %>%
  mutate(percent = count / sum(count) * 100)

pdf('7_Histogram_cluster_beam_positive_percent_split_by_seurat_CD8_T.pdf', width = 24, height = 9)
ggplot(cluster_stats, aes(x = seurat_clusters, y = percent, fill = antigen_positive)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = c("positive" = "red", "negative" = "gray")) +
  labs(
    title = "Percentage of Positive and Negative Cells in Each Cluster by Antigen",
    x = "Cluster",
    y = "Percentage",
    fill = "Antigen Positive"
  ) +
  facet_wrap(~ antigen, scales = "free_y",nrow = 3) +  # 按 antigen 分面，并允许 Y 轴独立缩放
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10, face = "bold")  # 调整分面标题样式
  )
dev.off()
gc()

pdf('7_Histogram_cluster_beam_positive_percent_split_by_seurat_1_CD8_T.pdf', width = 24, height = 9)
ggplot(cluster_stats, aes(x = factor(seurat_clusters), y = percent, fill = antigen_positive)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = c("positive" = "red", "negative" = "gray")) +
  labs(
    title = "Percentage of Positive and Negative Cells in Each Cluster by Antigen",
    x = "Cluster",
    y = "Percentage",
    fill = "Antigen Positive"
  ) +
  facet_wrap(~ antigen, scales = "free",nrow = 3) +  # 按 antigen 分面，并允许 X 和 Y 轴独立缩放
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 调整 X 轴文字角度
    strip.text = element_text(size = 10, face = "bold"), # 调整分面标题样式
    panel.spacing = unit(1, "lines")  # 增加分面之间的间距
  )
dev.off()
gc()

library(dplyr)
library(ggplot2)

# 过滤出 Positive 数据并重新计算百分比
positive_stats <- cluster_stats %>%
  filter(antigen_positive == "positive") %>%  # 只保留 Positive 数据
  group_by(seurat_clusters, antigen) %>%
  summarise(percent = sum(percent), .groups = "drop")  # 计算每个 cluster 的 Positive 百分比

# 绘制水平柱状图
pdf('7_Histogram_cluster_beam_positive_CD8_T.pdf', width = 20, height = 20)
ggplot(positive_stats, aes(x = percent, y = factor(seurat_clusters), fill = antigen)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = my36colors[1:length(unique(positive_stats$antigen))]) +  # 自动匹配 antigen 数量
  labs(
    title = "Positive Proportion of Each Cluster by Antigen",
    x = "Percentage of Positive Cells (%)",
    y = "Cluster",
    fill = "Antigen"
  ) +
  facet_wrap(~ antigen, scales = "free_x") +  # 按 antigen 分面，允许横坐标独立缩放
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 15, face = "bold"),  # 调整纵坐标文本样式
    strip.text = element_text(size = 15, face = "bold"),  # 调整分面标题样式
    panel.spacing = unit(1, "lines")  # 增加分面之间的间距
  )

dev.off()
gc()

library(ggplot2)

# 创建颜色映射，每个 cluster 分配一个固定颜色
cluster_colors_1 <- setNames(my36colors[1:length(unique(positive_stats$seurat_clusters))], 
                             unique(positive_stats$seurat_clusters))
pdf('7_Histogram_beam_positive_percent_cluster_split_by_antigen_CD8_T.pdf', width = 18, height = 18)
ggplot(positive_stats, aes(x = percent, y = factor(seurat_clusters), fill = factor(seurat_clusters))) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = cluster_colors_1) +  # 使用一致的 cluster 颜色
  labs(
    title = "Positive Proportion of Each Cluster by Antigen",
    x = "Percentage of Positive Cells (%)",
    y = "Cluster",
    fill = "Cluster"
  ) +
  facet_wrap(~ antigen, scales = "free_x") +  # 每个 antigen 一个子图
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),  # 调整纵坐标文字样式
    strip.text = element_text(size = 12, face = "bold"),  # 调整分面标题样式
    panel.spacing = unit(1, "lines")  # 增加分面之间的间距
  )
dev.off()
gc()


library(dplyr)
library(tidyr)
library(pheatmap)

# 1. 提取 antigen_positive == "positive" 的行
positive_data <- CD8_T@meta.data %>%
  filter(antigen_positive == "positive")

# 2. 构建矩阵
heatmap_data <- positive_data %>%
  select(antigen, clonotype_id, antigen_specificity_score) %>%
  group_by(antigen, clonotype_id) %>%
  summarise(antigen_specificity_score = mean(antigen_specificity_score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = clonotype_id, values_from = antigen_specificity_score, values_fill = 0)

# 将 antigen 设置为行名
heatmap_matrix <- as.matrix(heatmap_data[,-1])  # 去掉 antigen 列
rownames(heatmap_matrix) <- heatmap_data$antigen

# 3. 筛选列（保留总和大于 25 的 clonotype_id）
col_sums <- colSums(heatmap_matrix)  # 计算每列的总和
heatmap_matrix_filtered <- heatmap_matrix[, col_sums > 25]  # 筛选出总和大于 25 的列

# 4. 确保矩阵中没有缺失值
heatmap_matrix_filtered[is.na(heatmap_matrix_filtered)] <- 0

# 5. 绘制热图
pdf('4_heatmap_antigen_specificity_clustermap.pdf', width = 20, height = 9)
pheatmap(
  heatmap_matrix_filtered,
  cluster_rows = TRUE,            # 行聚类
  cluster_cols = TRUE,            # 列聚类
  scale = "none",                 # 不进行数据缩放
  color = colorRampPalette(c("white", "blue"))(50),  # 颜色映射
  main = "Filtered Antigen Specificity Score Heatmap",
  fontsize_row = 10,
  fontsize_col = 8
)

dev.off()
gc()



#---------------------------NMF方法分析 CD8+T BEAM_Positive-scRNASEQ，boundary设置为90----------------------------------#    

###boundary设置为90
# 创建数据框
dfdf <- data.frame(
  antigen_specificity_score = CD8_T$antigen_specificity_score,
  antigen = CD8_T$antigen
)

# 聚类函数
cluster_antigen <- function(df) {
  df$boundary <- 90  # 将 boundary 设置为固定值 0.9
  df$cluster <- ifelse(df$antigen_specificity_score > df$boundary, "positive", "negative")  # 使用固定阈值分类
  return(df)
}

# 分组聚类
library(dplyr)
dfdf_clustered <- dfdf %>%
  group_by(antigen) %>%
  group_modify(~ cluster_antigen(.x)) %>%
  ungroup()

# 将聚类结果添加到原始数据中
CD8_T$antigen_positive <- dfdf_clustered$cluster
pdf('7_DimPlot_Positive_BEAM_UMAP_boundary90.pdf', width = 18, height = 12)
DimPlot(CD8_T, reduction = "umap", label = F, pt.size = 0.5, cols = c("grey90", "red"),
        split.by = "antigen",
        group.by = "antigen_positive", ncol = 5, raster = F) & theme(aspect.ratio = 1)
dev.off()
#-------------------------------------NMF方法分析 CD8+T BEAM_Positive-scRNASEQ-----------------------------------------#    

#-------------------------------------scRNA+scTCR联合分析-------------------------------------------------#
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
install.packages('ggpubr')
library(ggpubr)
library(ComplexHeatmap)
# 加载必要库
library(dplyr)
library(Seurat)

# 读取样本 Z2_1 和 Z2_2 的 filtered_contig_annotations.csv 和 clonotypes.csv
Z2_1_annotations <- read.csv("抗原特异性T细胞 单细胞转录组+TCR/2024_analysis_reports_zhaohaichao_NEO-SP-164_zhaohaichao_NEO-SP-164_10xBEAM_human_20240914/zhaohaichao_NEO-SP-164_10xBEAM_human_20240914/Z2-1/vdj_t/filtered_contig_annotations.csv")
Z2_1_clonotypes <- read.csv("抗原特异性T细胞 单细胞转录组+TCR/2024_analysis_reports_zhaohaichao_NEO-SP-164_zhaohaichao_NEO-SP-164_10xBEAM_human_20240914/zhaohaichao_NEO-SP-164_10xBEAM_human_20240914/Z2-1/vdj_t/clonotypes.csv")
Z2_2_annotations <- read.csv("抗原特异性T细胞 单细胞转录组+TCR/2024_analysis_reports_zhaohaichao_NEO-SP-164_zhaohaichao_NEO-SP-164_10xBEAM_human_20240914/zhaohaichao_NEO-SP-164_10xBEAM_human_20240914/Z2-2/vdj_t/filtered_contig_annotations.csv")
Z2_2_clonotypes <- read.csv("抗原特异性T细胞 单细胞转录组+TCR/2024_analysis_reports_zhaohaichao_NEO-SP-164_zhaohaichao_NEO-SP-164_10xBEAM_human_20240914/zhaohaichao_NEO-SP-164_10xBEAM_human_20240914/Z2-2/vdj_t/clonotypes.csv")

# 添加样本标识
Z2_1_annotations$sample <- "Z2_1"
Z2_2_annotations$sample <- "Z2_2"

# 合并两个样本的 annotations
combined_annotations <- bind_rows(Z2_1_annotations, Z2_2_annotations)


# 确保 barcode 信息一致性，并直接使用 raw_clonotype_id 作为 clonotype 列
combined_annotations <- combined_annotations %>%
  filter(!is.na(raw_clonotype_id)) %>% # 过滤掉 raw_clonotype_id 中的缺失值
  mutate(clonotype = raw_clonotype_id) %>% # 将 clonotype 直接赋值为 raw_clonotype_id
  select(barcode, sample, clonotype) # 保留 barcode, sample, 和 clonotype 列


# 检查克隆型信息
combined_clonotypes <- bind_rows(Z2_1_clonotypes, Z2_2_clonotypes)

# 加载 dplyr 库
library(dplyr)

# 移除不需要的列
obj@meta.data <- obj@meta.data %>%
  select(-sample.y, -clonotype.x, -sampleTCR, -clonotype.y, -sample, -clonotype)

# 确保条形码的一致性
obj@meta.data$barcode <- rownames(obj@meta.data)


# 合并前后的行数对比
# 查看 scRNA-seq 数据的 meta.data 列
colnames(obj@meta.data)
table(obj@meta.data$sample)
nrow_before <- nrow(obj@meta.data) # 合并前
#obj@meta.data <- obj@meta.data %>%
  left_join(combined_annotations, by = c("barcode" = "barcode"))
#nrow_after <- nrow(obj@meta.data) # 合并后
#这里不做合并，使用原来metadata中的raw_clonotype_id列名
# 打印行数对比
cat("行数变化: 合并前 =", nrow_before, ", 合并后 =", nrow_after, "\n")


# 检查是否存在 NA 值
sum(is.na(obj@meta.data$raw_clonotype_id)) # 检查 clonotype 列的 NA 数量
sum(is.na(obj@meta.data$sample.x))    # 检查 sample 列的 NA 数量

# 检查 sample 列的唯一值
unique(obj@meta.data$sample.x)

# 检查 clonotype 和 sample 的长度是否一致
length(obj@meta.data$raw_clonotype_id)
length(obj@meta.data$sample.x)

table(obj@meta.data$raw_clonotype_id, obj@meta.data$sample.x)

DimPlot(obj, group.by = "raw_clonotype_id")

table(obj@meta.data$raw_clonotype_id, obj@meta.data$sample.x)

obj

colnames(obj@meta.data)
table(obj@meta.data$seurat_clusters)


CD8_T
colnames(CD8_T@meta.data)
table(CD8_T@meta.data$seurat_clusters)

library(dplyr)

# 提取 metadata 信息
metadata <- CD8_T@meta.data

# 按 cluster 和 clonotype_id 统计 unique clone 的数量
clone_stats <- metadata %>%
  group_by(seurat_clusters) %>%
  summarise(
    total_cells = n(),  # 每个 cluster 的细胞总数
    unique_clones = n_distinct(clonotype_id),  # 每个 cluster 的 unique clone 数量
    clone_proportion = n_distinct(clonotype_id) / n()  # unique clone 占比
  )

# 查看结果
print(clone_stats)

library(ggplot2)
cluster_colors_1 <- setNames(my36colors[1:length(unique(positive_stats$seurat_clusters))], 
                             unique(positive_stats$seurat_clusters))
cluster_colors_2 <- setNames(cluster_cols[1:length(unique(CD8_T$seurat_clusters))], 
                             unique(CD8_T$seurat_clusters))
# 绘制每个 cluster 的 unique clone 数量柱状图
ggplot(clone_stats, aes(x = seurat_clusters, y = unique_clones),) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Unique Clones per Cluster", x = "Cluster", y = "Unique Clones")

# 绘制每个 cluster 的 unique clone 比例柱状图，带颜色
#ggplot(clone_stats, aes(x = seurat_clusters, y = clone_proportion, fill = seurat_clusters)) +
#  geom_bar(stat = "identity") +
#  scale_fill_manual(values = cluster_colors_1) + # 使用自定义颜色
#  theme_minimal() +
#  labs(title = "Proportion of Unique Clones per Cluster", x = "Cluster", y = "Proportion of Unique Clones")
# 去掉图注
pdf("10_Histogram_Proportion_of_Uniqu_Clones_per_Cluster.pdf", width = 8, height = 4)
ggplot(clone_stats, aes(x = seurat_clusters, y = clone_proportion, fill = seurat_clusters)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors_1) + # 使用自定义颜色
  theme_minimal() +
  labs(title = "Proportion of Unique Clones per Cluster", x = "Cluster", y = "Proportion of Unique Clones") +
  theme(legend.position = "none") # 隐藏图注
dev.off()
gc()


library(dplyr)
install.packages('ineq')
library(ineq) # 用于计算基尼系数

# 提取所需的元数据
clonotype_data <- CD8_T@meta.data %>%
  select(seurat_clusters, clonotype_id) %>%
  filter(!is.na(clonotype_id)) # 移除 NA 克隆型

# 计算每个 cluster 中的克隆分布
cluster_clonotypes <- clonotype_data %>%
  group_by(seurat_clusters, clonotype_id) %>%
  summarise(clone_count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(clone_proportion = clone_count / sum(clone_count))

# 按 cluster 计算基尼系数
gini_results <- cluster_clonotypes %>%
  group_by(seurat_clusters) %>%
  summarise(
    gini_coefficient = ineq(clone_proportion, type = "Gini"),
    mean_clone_count = mean(clone_count),
    clone_diversity = n() # unique clones in the cluster
  )

# 查看基尼系数结果
print(gini_results)
library(ggplot2)

# 绘制基尼系数分布
pdf("10_Histogram_Gini_coefficients_throughout_clusters.pdf", width = 8, height = 4)

# 创建柱状图并映射颜色
ggplot(gini_results, aes(x = as.factor(seurat_clusters), y = gini_coefficient, fill = as.factor(seurat_clusters))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors_1) + # 使用自定义颜色
  theme_minimal() +
  labs(title = "Gini Coefficient of Clonal Expansion",
       x = "Cluster",
       y = "Gini Coefficient",
       fill = "Cluster") # 为图例添加标题

ggplot(gini_results, aes(x = as.factor(seurat_clusters), y = gini_coefficient, fill = as.factor(seurat_clusters))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors_1) + # 使用自定义颜色
  theme_minimal() +
  theme(legend.position = "none") + # 去掉图注
  labs(title = "Gini Coefficient of Clonal Expansion",
       x = "Cluster",
       y = "Gini Coefficient")
dev.off()
gc()


# 按 seurat_clusters 分组计算克隆数量
library(dplyr)

# 统计每个 cluster 的总克隆数量
clonotype_summary <- CD8_T@meta.data %>%
  group_by(seurat_clusters) %>%
  summarize(total_clones = n_distinct(raw_clonotype_id)) # 计算 unique 克隆数

print(clonotype_summary)
library(ggplot2)

library(ggplot2)

# 可视化每个 cluster 的总克隆数量
# 设置每个 cluster 的颜色
cluster_colors_1 <- setNames(my36colors[1:length(unique(clonotype_summary$seurat_clusters))], 
                             unique(clonotype_summary$seurat_clusters))

# 绘制柱状图并修改颜色和去掉图例
pdf("10_Histogram_clonotype_summary_throughout_clusters.pdf", width = 8, height = 4)
ggplot(clonotype_summary, aes(x = seurat_clusters, y = total_clones, fill = seurat_clusters)) +
  geom_bar(stat = "identity") + # 使用柱状图
  scale_fill_manual(values = cluster_colors_1) + # 使用指定的颜色
  labs(
    title = "Total Clones per Cluster",
    x = "Cluster",
    y = "Total Clones"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # 调整 X 轴文字角度
    plot.title = element_text(hjust = 0.5), # 标题居中
    legend.position = "none" # 去掉图例
  )
dev.off()
gc()

library(immunarch)
library(ggplot2)

# 准备 immunarch 格式数据
# 假设 CD8_T@meta.data 包含 clonotype_id 和 cluster 信息
clonotype_data <- CD8_T@meta.data %>%
  filter(!is.na(clonotype_id)) %>%  # 过滤掉缺失的克隆型
  select(clonotype_id, seurat_clusters)  # 提取必要的列

# 将数据加载为 immunarch 格式
immdata <- list(
  meta = clonotype_data, # 元数据
  data = clonotype_data  # 克隆型数据
)

# 使用 quantContig 探索每个 cluster 的 unique 克隆型数量
p2 <- quantContig(
  immdata, 
  cloneCall = "gene+nt",  # 使用基因和核苷酸序列定义克隆型
  scale = FALSE           # 不进行比例缩放
)

# 自定义图形样式
p2 <- p2 + 
  theme(
    axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# 显示图形
print(p2)


library(immunarch)
library(ggplot2)

library(dplyr)

# 构建 immdata 数据
clonotype_data <- CD8_T@meta.data %>%
  filter(!is.na(cdr3s_aa)) %>%  # 仅保留具有 CDR3 氨基酸序列的条目
  select(
    barcode,          # 细胞的唯一条形码
    cdr3aa = cdr3s_aa, # 重命名 cdr3s_aa 为 cdr3aa
    seurat_clusters,   # cluster 信息
    clonotype_id,      # 克隆型 ID
    v_gene,            # V 基因
    j_gene             # J 基因
  )

# 创建 immunarch 格式的数据
immdata <- list(
  meta = clonotype_data, # 元数据
  data = clonotype_data  # 克隆型数据
)

# 查看 immdata 数据
head(immdata$data)
colnames(immdata$data)

library(ggplot2)

# 使用 repExplore 分析 unique 克隆型数量
unique_clones <- repExplore(
  immdata$data,
  .method = "len",  # 按克隆型长度计算
  .col = "seurat_clusters"  # 按 cluster 分组
)
# 检查数据中的 cdr3aa 列是否为空
sum(is.na(immdata$data$cdr3aa))  # 应返回 0
# Load necessary libraries
library(immunarch)
library(ggplot2)
library(dplyr)

# Load the metadata file
CD8_metadata <-CD8_T@meta.data

# Prepare the data for immunarch analysis
# Extract the relevant columns and rename them
immdata <- list(
  meta = CD8_metadata,
  data = CD8_metadata %>%
    select(barcode, seurat_clusters = sample, cdr3aa = cdr3s_nt) %>%
    filter(!is.na(cdr3aa))  # Remove rows without CDR3 sequences
)

# Ensure cluster is a factor
immdata$data$seurat_clusters <- as.factor(immdata$data$seurat_clusters)

# Explore unique clones per cluster using quantContig
p2 <- quantContig(
  immdata$data,
  cloneCall = "gene+nt",  # Define clones by gene and nucleotide sequence
  scale = FALSE           # Disable scaling
) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75))  # Customize x-axis text

# Display the plot
print(p2)
?quantContig
library(scRepertoire)
?clonalQuant
library(BiocStyle)
library(scater)
BiocManager::install('scater',force = T)
install.packages('scater')
## -----------------------------------------------------------------------------
suppressMessages(library(scRepertoire))


# 可视化
plot(unique_clones) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.85, hjust = 0.75))

###------香浓指数
# 加载 dplyr 以便数据操作
library(dplyr)

# 提取克隆型数据和分群信息
clonotype_data <- CD8_T@meta.data %>%
  group_by(seurat_clusters, clonotype_id) %>%
  summarise(clone_count = n()) %>%   # 每个 cluster 内的克隆型计数
  ungroup() %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = clone_count / sum(clone_count)) %>% # 计算比例
  summarise(shannon_index = -sum(proportion * log(proportion), na.rm = TRUE)) # 香浓指数

# 查看每个 cluster 的香浓指数
print(clonotype_data)

library(ggplot2)

ggplot(clonotype_data, aes(x = seurat_clusters, y = shannon_index)) +
  geom_bar(stat = "identity") +
  labs(title = "Shannon Index by Cluster", x = "Cluster", y = "Shannon Index")

install.packages("dplyr")
if (!requireNamespace("immunarch", quietly = TRUE)) {
  remotes::install_github("immunomind/immunarch")
}

library(dplyr)

# 提取元数据和克隆型信息
clonotype_data <- CD8_T@meta.data %>%
  select(barcode, clonotype_id, seurat_clusters) %>%
  filter(!is.na(clonotype_id)) # 排除没有克隆型信息的条目

# 按 cluster 计算每个 cluster 内的 Shannon 指数
shannon_index <- clonotype_data %>%
  group_by(seurat_clusters) %>%
  summarise(
    Shannon = -sum((prop.table(table(clonotype_id)) * log(prop.table(table(clonotype_id)))))
  )

# 查看结果
print(shannon_index)

# 加载 immunarch 包
library(immunarch)

# 准备数据：每个样本的克隆型信息
quant_contig_data <- clonotype_data %>%
  group_by(seurat_clusters) %>%
  summarise(
    Unique_Clones = n_distinct(clonotype_id), # Unique 克隆型个数
    Total_Cells = n(), # 总细胞数
    Clone_Ratio = Unique_Clones / Total_Cells # Unique 克隆型比例
  )

# 查看结果
print(quant_contig_data)

library(ggplot2)
pdf("10_Histogram_Shannon_Index_by_Cluster.pdf", width = 8, height = 4)
ggplot(shannon_index, aes(x = seurat_clusters, y = Shannon)) +
  geom_bar(stat = "identity") +
  labs(title = "Shannon Index by Cluster", x = "Cluster", y = "Shannon Index") +
  theme_minimal()

pdf("10_Histogram_Shannon_coefficients_throughout_clusters.pdf", width = 8, height = 4)
ggplot(quant_contig_data, aes(x = seurat_clusters, y = Clone_Ratio)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Unique Clone Ratio by Cluster", x = "Cluster", y = "Clone Ratio") +
  theme_minimal()
dev.off()
gc()


obj
table(obj@meta.data$celltype)

library(dplyr)

# 提取样本和细胞类型信息
cell_proportions <- obj@meta.data %>%
  filter(celltype %in% c("CD4 T", "CD8 T")) %>%  # 筛选 CD4 和 CD8 T 细胞
  group_by(sample, celltype) %>%
  summarise(CellCount = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(Proportion = CellCount / sum(CellCount) * 100) %>%  # 计算每种细胞类型的比例
  ungroup()

# 查看整理后的数据
print(cell_proportions)


# 自定义颜色，基于示例图的配色
custom_colors <- c(
  "CD4 T" = '#CC3333',  # 示例图的 CD4 T 类似颜色
  "CD8 T" = "#33A02C"   # 示例图的 CD8 T 类似颜色
)

# 绘制水平方向的堆积柱状图
pdf("10_CD4_CD8_T_Proportions_by_Sample.pdf", width = 8, height = 4)
ggplot(cell_proportions, aes(x = sample, y = Proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +  # 转换为水平方向
  theme_minimal() +
  labs(
    title = "CD4 and CD8 T Cell Proportions by Sample",
    x = "Sample",
    y = "Proportion (%)",
    fill = "Cell Type"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # 格式化 y 轴为百分比
  scale_fill_manual(values = custom_colors)  # 使用自定义颜色
dev.off()
gc()

#-------------------------------------scRNA+scTCR联合分析-------------------------------------------------#

# UMAP 可视化克隆分布
pdf("UMAP_Clonotype_Distribution.pdf", width = 20, height = 20)
DimPlot(obj, reduction = "wnn.umap", group.by = "clonotype", label = FALSE, pt.size = 0.5) +
  ggtitle("Clonotype Distribution on UMAP") +
  theme_minimal()
dev.off()

# 按克隆频率绘制 UMAP
pdf("UMAP_Clonotype_Frequency.pdf", width = 10, height = 8)
FeaturePlot(CD8_T, features = "frequency", reduction = "umap", cols = c("lightgrey", "blue")) +
  ggtitle("Clonotype Frequency on UMAP") +
  theme_minimal()
dev.off()
gc()

FeaturePlot(CD8_T,"TIGIT")

# 在 UMAP 上按克隆型分布绘图
pdf("UMAP_Clonotype_Distribution.pdf", width = 20, height = 20)
DimPlot(
  obj, 
  reduction = "wnn.umap", 
  group.by = "clonotype", 
  label = FALSE, 
  pt.size = 0.5
) + ggtitle("Clonotype Distribution on UMAP") +
  theme_minimal()
dev.off()

# 按克隆频率绘制 UMAP
pdf("UMAP_Clonotype_Frequency.pdf", width = 10, height = 8)
FeaturePlot(
  CD8_T, 
  features = "frequency", 
  reduction = "umap", 
  cols = c("lightgrey", "blue")
) + ggtitle("Clonotype Frequency on UMAP") +
  theme_minimal()
dev.off()

# 分类规则
clone_categories <- function(freq) {
  if (freq <= 1) {
    return("Single")
  } else if (freq <= 5) {
    return("Small")
  } else if (freq <= 20) {
    return("Medium")
  } else if (freq <= 100) {
    return("Large")
  } else {
    return("Hyperexpanded")
  }
}

# 添加分类到元数据
obj$cloneType <- sapply(obj$frequency, clone_categories)
CD8_T$cloneType <- sapply(CD8_T$frequency, clone_categories)
# 检查分类分布
table(obj$cloneType)

# 定义颜色
colorblind_vector <- colorRampPalette(rev(c(
  "#0D0887FF", "#FDC926FF","#DC050C","#9C179EFF")))

# 绘制按 cloneType 的 UMAP
pdf("UMAP_CloneType_Distribution.pdf", width = 20, height = 20)
p6<-DimPlot(
  CD8_T, 
  reduction = "umap", 
  group.by = "cloneType", 
  label = TRUE, 
  pt.size = 0.8
) + scale_color_manual(values = colorblind_vector(4)) +
  ggtitle("Clone Type Distribution on UMAP") +
  theme_minimal()

dev.off()
pdf("UMAP_CloneType_Distribution.pdf", width = 20, height = 20)
p2/p6
p4/p5
# 统计每种细胞类型的克隆型分布
clonotype_celltype_stats <- obj@meta.data %>%
  group_by(seurat_clusters, clonotype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(percent = count / sum(count) * 100)

# 保存统计结果
write.csv(clonotype_celltype_stats, "Clonotype_Celltype_Stats.csv")

# 绘制按细胞类型的克隆分布柱状图
pdf("Clonotype_Distribution_By_Celltype.pdf", width = 12, height = 6)
ggplot(clonotype_celltype_stats, aes(x = seurat_clusters, y = percent, fill = clonotype)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(
    title = "Clonotype Distribution Across Cell Types",
    x = "Cell Type",
    y = "Percentage",
    fill = "Clonotype"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# 构建聚类和克隆型的矩阵
heatmap_data <- clonotype_cluster_stats %>%
  pivot_wider(names_from = clonotype, values_from = percent, values_fill = 0) %>%
  column_to_rownames("seurat_clusters")

# 绘制热图
library(pheatmap)
pdf("Clonotype_Cluster_Heatmap.pdf", width = 12, height = 8)
pheatmap(
  as.matrix(heatmap_data),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "blue"))(50),
  main = "Clonotype Distribution Across Clusters"
)
dev.off()
colnames(obj@meta.data)

FeaturePlot(CD8_T,"clonotype3073")




# 模态权重（RNA.weight）的小提琴图
VlnPlot(obj, features = "RNA.weight", group.by = "seurat_clusters", sort = TRUE, pt.size = 0.1) + NoLegend()
# 查看 Seurat 对象中抗体标签（ADT）的所有膜蛋白名称
rownames(obj[["ADT"]])
rownames(obj[['RNA']])
# 提取 ADT 数据矩阵
adt_data <- GetAssayData(obj, assay = "ADT", slot = "data")
write.csv(adt_data,"2_adt_data.csv")
# 查看某个膜蛋白的表达情况，例如 CD4 和 CD8
adt_data[c("CD4.1", "CD8A.1"), ]
# 查找所有包含 "CD" 的膜蛋白
grep("CD", rownames(obj[["ADT"]]), value = TRUE)

# 查找具体膜蛋白，比如含有 "HLA" 的蛋白
grep("HLA", rownames(obj[["ADT"]]), value = TRUE)
# 确保 sce.cite.all 包含 ADT Assay
Assays(sce.cite.all)  # 检查 Assay 列表

FeaturePlot(obj, features = "RNA.weight", reduction = "wnn.umap")
FeaturePlot(obj, features = "ADT.weight", reduction = "wnn.umap")
