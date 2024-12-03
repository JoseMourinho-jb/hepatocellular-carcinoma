########## Part 1: Setup and Initialization ##########

# Clear the workspace
rm(list = ls())

# Load required libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(DoubletFinder)
library(data.table)
library(ggrepel)
library(hdf5r)

# Define custom color palettes
mycol <- c('#ca6354','#338fc4','#fec645','#709ed1','#68c46e','#ab9fb6','#91ba90',
           '#f8a34b','#d5dca0','#c38876','#91d8f6','#dea4ca',
           '#f58888','#64ccd3','#fee8e6','#f2edba','#d6cde5',
           '#4c9e72','#b9dda4','#f29987','#90c3e8','#f8a289',
           '#7aca97','#fedbc4','#8fbac8','#baabab','#fbd94c',
           "#efb3af",'#9aa8b2',"#c0b39e","#d6e7f7")

# Additional qualitative color palette
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
mycol1 <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Set working directory

########## Part 2: Data Preprocessing and Quality Control ##########

# Load raw data files
data_directory <- list.files()
project_name <- c("1NT_P21", "1T_C21", "2NT_P24", "2T_C24", "3NT_P25", "3T_C25", 
                  "4NT_P29", "4T_C29", "5NT_P36", "5T_C36")
samples <- project_name

# Process the first sample and initialize Seurat object
sample1 <- make_seurat_object_and_doublet_removal(data_directory[1], samples[1])
seu_list <- sample1

# Process and merge all samples
for (i in 2:length(samples)) {
  sc.i <- make_seurat_object_and_doublet_removal(data_directory[i], samples[i])
  seu_list <- merge(seu_list, sc.i)
}

# Normalize, find variable features, and scale data
scRNA_harmony <- seu_list
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE)

# Perform integration using Harmony
library(harmony)
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:40) %>%
  FindClusters(resolution = 1.0)
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:40)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony", dims = 1:40)

# Save cluster identities
Idents(scRNA_harmony) <- 'seurat_clusters'
table(scRNA_harmony$seurat_clusters)


########## Part 3: Visualization ##########

# Save UMAP and TSNE plots
setwd('/Users/luocheng/Desktop/单细胞分析/0322-HCC-55w-paper/')

# UMAP visualizations
p1 <- DimPlot(scRNA_harmony, reduction = "umap", label = FALSE, cols = mycol)
p2 <- DimPlot(scRNA_harmony, reduction = "umap", split.by = 'orig.ident', cols = mycol)
p3 <- DimPlot(scRNA_harmony, reduction = "umap", group.by = 'orig.ident', cols = mycol)
ggsave(filename = "Figure1_UMAP_seurat_cluster.pdf", plot = p1, width = 5, height = 4)
ggsave(filename = "Figure2_UMAP_seurat_cluster_bySample.pdf", plot = p2, width = 30, height = 4)
ggsave(filename = "Figure3_UMAP_seurat_cluster_bySample_Overlay.pdf", plot = p3, width = 5, height = 4)

# TSNE visualizations
p4 <- DimPlot(scRNA_harmony, reduction = "tsne", label = FALSE, cols = mycol)
p5 <- DimPlot(scRNA_harmony, reduction = "tsne", split.by = 'orig.ident', cols = mycol)
p6 <- DimPlot(scRNA_harmony, reduction = "tsne", group.by = 'orig.ident', cols = mycol)
ggsave(filename = "Figure4_TSNE_seurat_cluster.pdf", plot = p4, width = 5, height = 4)
ggsave(filename = "Figure5_TSNE_seurat_cluster_bySample.pdf", plot = p5, width = 30, height = 4)
ggsave(filename = "Figure6_TSNE_seurat_cluster_bySample_Overlay.pdf", plot = p6, width = 5, height = 4)


########## Part 4: Doublet Removal and Quality Control ##########

# Function to create a Seurat object and remove doublets
make_seurat_object_and_doublet_removal <- function(data_directory, project_name) {
  # Load raw data
  colon.data <- Read10X(data.dir = data_directory)
  
  # Create Seurat object
  currentSample <- CreateSeuratObject(counts = colon.data, 
                                      project = project_name, 
                                      min.cells = 3, 
                                      min.features = 40)
  currentSample[["percent.mt"]] <- PercentageFeatureSet(currentSample, pattern = "^MT-")
  
  # Generate QC plots before filtering
  pdf(paste0("./result/qc_plots_", project_name, "_prefiltered.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3, pt.size = 0.05))
  dev.off()
  
  # Filter cells based on QC metrics
  currentSample <- subset(currentSample, 
                          subset = nFeature_RNA > 400 & 
                            nFeature_RNA < 4000 & 
                            percent.mt < 20 & 
                            nCount_RNA > 1000)
  
  # Normalize and scale the data
  currentSample <- seurat_standard_normalize_and_scale(currentSample, FALSE)
  
  # Run DoubletFinder to detect and remove doublets
  nExp_poi <- round(0.08 * ncol(currentSample) / 10000)
  seu_colon <- doubletFinder_v3(currentSample, PCs = 1:20, 
                                pN = 0.25, pK = 0.09, nExp = nExp_poi)
  
  # Filter out doublets
  seu_colon <- subset(seu_colon, subset = doublet.class != "Doublet")
  
  # Return filtered Seurat object
  return(seu_colon)
}

# Function for normalization and scaling
seurat_standard_normalize_and_scale <- function(seurat_obj, cluster, cluster_resolution = 1.0) {
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
  if (cluster) {
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
    seurat_obj <- FindClusters(seurat_obj, resolution = cluster_resolution)
  }
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
  return(seurat_obj)
}


########## Part 5: Cell Type Annotation and Visualization ##########

# Annotate cell types manually
Idents(scRNA_harmony) <- "seurat_clusters"
scRNA_harmony <- RenameIdents(scRNA_harmony, 
                              "0" = "Myeloid_cells", "1" = "Malignant_cells", 
                              "2" = "T_cells", "3" = "B_cells", "4" = "NK_cells", 
                              "5" = "Endothelial", "6" = "Mesenchymal_cells", 
                              "7" = "HPCs", "8" = "MAST")
scRNA_harmony@meta.data$celltype <- Idents(scRNA_harmony)

# Save annotated Seurat object
save(scRNA_harmony, file = "scRNA_harmony_Anno.Rda")

# Visualize annotated cell types using UMAP
p1 <- DimPlot(scRNA_harmony, group.by = "celltype", label = TRUE, cols = mycol)
ggsave(filename = "CellType_UMAP.pdf", plot = p1, width = 6, height = 4)


########## Part 6: Proportion Analysis ##########

library(tidyverse)
library(reshape2)

# Check the number of cells in each sample and cell type
table(scRNA_harmony$orig.ident)
table(scRNA_harmony$celltype)

# 1. Calculate proportions by group
pB2_df <- table(scRNA_harmony@meta.data$celltype, scRNA_harmony@meta.data$group) %>% melt()
colnames(pB2_df) <- c("Cluster", "Sample", "Number")

# Plot cell type proportions across groups
pB4group <- ggplot(data = pB2_df, aes(x = Sample, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.8, position = "fill") +
  scale_fill_manual(values = mycol) +
  theme_bw() +
  labs(x = "", y = "Proportion") +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1))

ggsave(filename = "Figure1_CellType_Proportion_Group.pdf", plot = pB4group, width = 6, height = 4)

# 2. Calculate proportions by sample
pB2_df <- table(scRNA_harmony@meta.data$celltype, scRNA_harmony@meta.data$orig.ident) %>% melt()
colnames(pB2_df) <- c("Cluster", "Sample", "Number")

# Plot cell type proportions across samples
pB4sample <- ggplot(data = pB2_df, aes(x = Sample, y = Number, fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.8, position = "fill") +
  scale_fill_manual(values = mycol) +
  theme_bw() +
  labs(x = "", y = "Proportion") +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black", angle = 45, hjust = 1))

ggsave(filename = "Figure2_CellType_Proportion_Sample.pdf", plot = pB4sample, width = 6, height = 4)



########## Part 7: Differential Gene Expression (DGE) Analysis ##########

# Set identity to "group" for comparison
Idents(scRNA_harmony) <- "group"

# Create a folder to save results
dir.create("DGE_Results/")
setwd("DGE_Results/")

# Initialize an empty data frame for DGE results
r.deg <- data.frame()

# Perform DGE analysis for each cell type
type <- unique(scRNA_harmony@meta.data$celltype)
for (i in 1:length(type)) {
  Idents(scRNA_harmony) <- "celltype"
  deg <- FindMarkers(scRNA_harmony, ident.1 = "Tumor", ident.2 = "NT", 
                     group.by = "group", subset.ident = type[i])
  
  # Save results to CSV
  write.csv(deg, file = paste0(type[i], "_Tumor_vs_NT.csv"))
  
  # Annotate results for visualization
  deg$gene <- rownames(deg)
  deg$celltype <- type[i]
  r.deg <- rbind(r.deg, deg)
}

# Filter significant DEGs
r.deg <- subset(r.deg, p_val_adj < 0.1 & abs(avg_log2FC) > 0.3)
r.deg$threshold <- ifelse(r.deg$avg_log2FC > 0, "Up", "Down")

########## Part 8: Marker Gene Visualization ##########

# UMAP heatmap visualization of marker genes
marker_genes <- c("CD79A", "IGHA2", "CLDN5", "PECAM1", "ACTA2", "KIT", "TF", "CD3D", "NKG7")

# Generate violin plots for marker genes
p <- VlnPlot(scRNA_harmony, features = marker_genes, group.by = "celltype", pt.size = 0) + NoLegend()
ggsave(filename = "ViolinPlot_MarkerGenes.pdf", plot = p, width = 10, height = 6)

# Dot plot for marker genes
p <- DotPlot(scRNA_harmony, features = marker_genes, scale = FALSE, cols = c("#E6E0B0", "#D20A13")) + RotatedAxis()
ggsave(filename = "DotPlot_MarkerGenes.pdf", plot = p, width = 8, height = 6)


########## Part 9: Subtype Analysis ##########

# Subset malignant cells
sc.subtype <- subset(scRNA_harmony, ident = "Malignant_cells")

# Normalize and integrate malignant cells
sc.subtype <- NormalizeData(sc.subtype) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose = FALSE)
sc.subtype <- RunHarmony(sc.subtype, group.by.vars = "orig.ident")
sc.subtype <- FindNeighbors(sc.subtype, reduction = "harmony", dims = 1:15) %>% 
  FindClusters(resolution = 0.3)
sc.subtype <- RunUMAP(sc.subtype, reduction = "harmony", dims = 1:15)

# Save results
save(sc.subtype, file = "MalignantCell_Subtype.Rda")

# Visualize UMAP for subtypes
p1 <- DimPlot(sc.subtype, reduction = "umap", label = TRUE, cols = mycol)
ggsave(filename = "UMAP_Malignant_Subtypes.pdf", plot = p1, width = 6, height = 4)


########## Part 10: GO Enrichment Analysis ##########

library(clusterProfiler)
library(org.Hs.eg.db)

# Select DEGs for GO analysis
degs <- subset(r.deg, threshold == "Up" & p_val_adj < 0.05)
genes <- degs$gene

# Perform GO enrichment
go_results <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL", 
                       pvalueCutoff = 0.05, qvalueCutoff = 0.05)

# Extract results and save
go_df <- as.data.frame(go_results@result)
write.csv(go_df, "GO_Enrichment_Results.csv")

# Visualize GO terms
top_terms <- go_df[1:10, ]
p <- ggplot(top_terms, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top GO Enrichment Terms", x = "GO Term", y = "-log10(p.adjust)")

ggsave(filename = "GO_Enrichment_Barplot.pdf", plot = p, width = 8, height = 6)


########## Part4. MVI  ##########
Idents(scRNA_harmony) <- c('group')

scRNA_harmony_T <- subset(scRNA_harmony, ident = 'Tumor')

Idents(scRNA_harmony_T) <- "orig.ident"

table(scRNA_harmony_T$orig.ident)
# 1T_C21 2T_C24 3T_C25 4T_C29 5T_C36 
# 4053   3566   3136   2596   5040 

# 在active.ident里直接改名
scRNA_harmony_T <- RenameIdents(scRNA_harmony_T,"1T_C21"="MVI_present","2T_C24"="MVI_present","3T_C25"="MVI_present",
                                "4T_C29"="MVI_absent","5T_C36"="MVI_absent")

table(scRNA_harmony_T@active.ident)
# MVI_present  MVI_absent 
# 10755        7636

# 将当前的active.ident写入metadata
scRNA_harmony_T@meta.data$MVI <- scRNA_harmony_T@active.ident

p5 <- DimPlot(scRNA_harmony_T, group.by ="celltype", split.by ='MVI', label=T, label.size=3, cols = mycol)
p6 <- DimPlot(scRNA_harmony_T, reduction = "tsne", group.by ="celltype", split.by ='MVI', label=T, label.size=3, cols = mycol)

ggsave(filename = "MVI_present_UMAP_celltype_byGroup.pdf", plot = p5, width = 7.5, height = 4)
ggsave(filename = "MVI_present_TSNE_celltype_byGroup.pdf", plot = p6, width = 7.5, height = 4)

save(scRNA_harmony_T, file = 'scRNA_harmony_Anno_Tumor_MVI.Rda')

#-------------------------------------------------------------------------------


library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(reshape2)

# scRNA_harmony <- scRNA_harmony_T

# 加载数据
load("scRNA_harmony_Anno.Rda")

# 查看有几个样本
table(scRNA_harmony$orig.ident)

# 查看每个细胞类型多少个细胞
table(scRNA_harmony$celltype)

# 01. 总体比例图
# Group比例
pB2_df <- table(scRNA_harmony@meta.data$celltype,scRNA_harmony@meta.data$group) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")

# 【可视化】画比例图 4*6
pB4group <- ggplot(data = pB2_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Figure1_Proportion_Group.pdf", plot = pB4group, width = 6, height = 4)

p1 <- ggplot(data = pB2_df, aes(x =Cluster, y = Number, fill =  Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Figure2_Proportion_Group_2.pdf", plot = p1, width = 6, height = 4)


# Orig.ident比例
pB2_df <- table(scRNA_harmony@meta.data$celltype,scRNA_harmony@meta.data$orig.ident) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
# 【可视化】画比例图 4*6
pB4sample <- ggplot(data = pB2_df, aes(x =Number, y = Sample, fill =Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Figure3_Proportion_Samples.pdf", plot = pB4sample, width = 6, height = 5)


p1 <- ggplot(data = pB2_df, aes(x =Cluster, y = Number, fill =  Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Figure4_Proportion_Samples_2.pdf", plot = p1, width = 6, height = 4)


# 其他metadata连图
pB2_df <- table(scRNA_harmony@meta.data$celltype,scRNA_harmony@meta.data$orig.ident) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")

# 【可视化1】画比例图 4*6
pB4sample <- ggplot(data = pB2_df, aes(x =Number, y = Cluster, fill =Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Proportion_orig_ident.pdf", plot = pB4sample, width = 5, height = 4)


pB2_df <- table(scRNA_harmony@meta.data$celltype,scRNA_harmony@meta.data$group) %>% melt()
colnames(pB2_df) <- c("Cluster","Group","Number")

# 【可视化2】画比例图 4*6
pB4sample <- ggplot(data = pB2_df, aes(x =Number, y = Cluster, fill =Group)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol[c(13,14)]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Proportion_group.pdf", plot = pB4sample, width = 5, height = 4)


pB2_df <- table(scRNA_harmony@meta.data$celltype) %>% melt()
colnames(pB2_df) <- c("Cluster","Cells")

# 【可视化3】画比例图 4*6
pB4sample <- ggplot(data = pB2_df, aes(x =Cells, y = Cluster, fill = Cluster)) +
  stat_summary(geom = "bar",fun = "mean",position = position_dodge(0.9))+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Proportion_cells.pdf", plot = pB4sample, width = 5, height = 4)

#-------------------------------------------------------------------------------



library(ggplot2)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(pheatmap)

# Marker基因展示1: 热图
# 计算差异基因
Idents(scRNA_harmony) <- 'celltype'
degs <- FindAllMarkers(scRNA_harmony, logfc.threshold = 0.5, test.use = "roc", return.thresh = 0.25, min.pct = 0.3, only.pos = T) 

# 筛选显著的差异基因 
degs_sig <- degs %>% filter(pct.1 > 0.3 & power > 0.25) %>% filter(cluster != "other") %>% arrange(cluster, -power)  

# 选择用来画热图的基因
degs_top50 <- degs_sig %>% group_by(cluster) %>% top_n(50, power) %>% top_n(50, avg_diff) %>% arrange(cluster, -power)

write.csv(degs_top50, file = 'Top50genes_all_celltype_forheatmap.csv')

# 针对data数据，将50个基因取平均值：得到每个基因，在每个细胞类型中的平均表达
avgData <- scRNA_harmony@assays$RNA@data[degs_top50$gene,] %>% 
  apply(1, function(x){
    tapply(x, scRNA_harmony$celltype, mean) # ExpMean
  }) %>% t

# 计算z-score（每次画热图都要用这一步！！！）超过2也赋值2
phData <- MinMax(scale(avgData), -2, 2) # z-score
rownames(phData) <- 1:nrow(phData)

celltype <- colnames(phData)

table(Idents(scRNA_harmony))

# [1] "Myeloid_cells"     "Malignant_cells"   "T_cells"           "B_cells"           "NK_cells"          "Endothelial"      
# [7] "Mesenchymal_cells" "HPCs"              "MAST"  

ann_colors <- list(cluster = c("Myeloid_cells" = mycol[1],"Malignant_cells"= mycol[2],"T_cells"= mycol[3],
                               "B_cells" = mycol[4],"NK_cells"= mycol[5],"Endothelial"= mycol[6],
                               "Mesenchymal_cells" = mycol[7],"HPCs" = mycol[8],"MAST" = mycol[9]))

colnames(phData)

phres <- pheatmap(
  phData, 
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#225ea8", "white", "#d7301f"))(99),
  scale = "row",
  cluster_rows = F, #不按行聚类
  cluster_cols = F, #按列聚类
  clustering_method = "complete",
  show_rownames = F, #显示cluster名
  annotation_row = data.frame(cluster = degs_top50$cluster)
) 

ggsave(filename = "Figure1_Topgene_heatmap_plot_order.pdf", plot = phres, width = 5, height = 5)


table(Idents(scRNA_harmony))

Myeloid_cells <- c('LYZ')
Malignant_cells <- c('TF')
MAST_cells <- c('KIT','MS4A2','GATA2')
B_lymphocytes <- c('CD79A')
NK_cells <- c('NKG7','KLRD1')
Endothelial <- c('PECAM1','CLDN5','FLT1','VWF')
Mesenchymal_cells <- c('ACTA2')
HPCs <- c('EPCAM')
T_cells <- c('CD3D','CD3G')


# # all features
features <- c(Myeloid_cells,Malignant_cells,T_cells,B_lymphocytes,NK_cells,Endothelial,Mesenchymal_cells,HPCs,MAST_cells)

features <- unique(features)

Idents(scRNA_harmony) <- c('celltype')

p1 <- VlnPlot(scRNA_harmony, features=features, stack = T, fill.by='ident', cols = mycol, pt.size=0) + NoLegend() 
ggsave(filename = "Figure2_Marker_violin_plot_order.pdf", plot = p1, width = 10, height = 5)

#-------------------------------------------------------------------------------


########## Part4. Augur  ##########
sc.augur <- scRNA_harmony_T

rm(scRNA_harmony_T)

library(Augur)
augur <- calculate_auc(sc.augur, cell_type_col = 'celltype', label_col = 'MVI')
augur$AUC

augur.auc <- as.data.frame(augur$AUC)

ordered_types <- unique(augur.auc$cell_type)

# 【可视化】Augur的AUC结果
augur.auc$cell_type_status <- ifelse(augur.auc$auc > 0.8, 'sig', 'nonsig')

p <- ggplot(augur.auc, aes(x=cell_type, y=auc)) + 
  geom_segment(aes(x=cell_type, xend=cell_type, y=0, yend=auc, color=cell_type_status), linewidth=1.5) +
  geom_point(aes(color=cell_type_status), size=2.3) +
  scale_x_discrete(limits=augur.auc$cell_type) +
  scale_color_manual(values = c('sig' = "#D20A13", 'nonsig' = "#223D6C")) +
  geom_hline(yintercept=0.6, linetype = "dashed") +
  coord_flip() +
  theme_bw() +
  labs(title='Augur')

ggsave(filename = "Augur_AllCellType_MVI.pdf", plot = p, width = 5, height = 5)
write.csv(augur.auc, file = 'Augur_result_celltype.csv')


## 计算扰动最大的cluster
augur <- calculate_auc(sc.augur, cell_type_col = 'seurat_clusters', label_col = 'MVI')

augur$AUC

augur.auc <- as.data.frame(augur$AUC)

ordered_types <- unique(augur.auc$cell_type)

# 【可视化】Augur的AUC结果
augur.auc$cell_type_status <- ifelse(augur.auc$auc > 0.8, 'sig', 'nonsig')

p <- ggplot(augur.auc, aes(x=cell_type, y=auc)) + 
  geom_segment(aes(x=cell_type, xend=cell_type, y=0, yend=auc, color=cell_type_status), linewidth=1.5) +
  geom_point(aes(color=cell_type_status), size=2.3) +
  scale_x_discrete(limits=augur.auc$cell_type) +
  scale_color_manual(values = c('sig' = "#D20A13", 'nonsig' = "#223D6C")) +
  geom_hline(yintercept=0.6, linetype = "dashed") +
  coord_flip() +
  theme_bw() +
  labs(title='Augur')

ggsave(filename = "Augur_AllCellType_2.pdf", plot = p, width = 4, height = 5)
write.csv(augur.auc, file = 'Augur_result_seurat_clusters.csv')

#-------------------------------------------------------------------------------


########## Part6. 差异散点图  ##########
Idents(scRNA_harmony) <- "group"

table(Idents(scRNA_harmony))
# NT Tumor 
# 23679 18391


# scRNA_harmony <- scRNA_harmony_T
# Idents(scRNA_harmony) <- "MVI"
# table(Idents(scRNA_harmony))
# MVI_present  MVI_absent 
# 10755        7636 


# 01. 创建文件夹
dir.create("deg_MVI_present_vs_MVI_absent/")
setwd("deg_MVI_present_vs_MVI_absent/")

r.deg=data.frame()

table(scRNA_harmony@meta.data$orig.ident)
type <- unique(scRNA_harmony@meta.data$celltype)

# 02. 批量计算样本A vs. 样本B中，每一个细胞类型的差异基因(Wilcox)
for (i in 1:length(type)) {
  Idents(scRNA_harmony)="celltype"
  deg=FindMarkers(scRNA_harmony,ident.1 = "MVI_present",ident.2 = "MVI_absent",
                  group.by = "MVI",subset.ident =type[i])
  write.csv(deg,file = paste0(type[i],'deg.csv') )
  deg$gene=rownames(deg) # 当前细胞类型的差异基因
  deg$celltype=type[i]   # 当前细胞类型
  deg$unm=i-1            # 细胞类型的代号
  r.deg=rbind(deg,r.deg)
}

table(r.deg$unm) 

# 备份计算好的差异表达矩阵，防止出错（r.deg2为备份）
r.deg2 <- r.deg

r.deg <- subset(r.deg, p_val_adj < 0.1 & abs(avg_log2FC) > 0.3)
r.deg$threshold <- as.factor(ifelse(r.deg$avg_log2FC > 0 , 'Up', 'Down'))

# p.adj判定该基因是Highly/Lowly
r.deg$adj_p_signi <- as.factor(ifelse(r.deg$p_val_adj < 0.05 , 'Highly', 'Lowly'))
r.deg$thr_signi <- paste0(r.deg$threshold, "_", r.deg$adj_p_signi)
r.deg$unm %<>% as.vector(.) %>% as.numeric(.)

type <- unique(scRNA_harmony@meta.data$celltype)
r.deg <- r.deg[r.deg$celltype %in% type,]


# 03. 挑选图里展示名称的基因数目（自定义显示想要展示的基因名）
# 这里挑选log2FC为top5的基因进行展示
top_up_label <- r.deg %>% 
  subset(., threshold%in%"Up") %>% 
  group_by(unm) %>% 
  top_n(n = 5, wt = avg_log2FC) %>% 
  as.data.frame()

top_down_label <- r.deg %>% 
  subset(., threshold %in% "Down") %>% 
  group_by(unm) %>% 
  top_n(n = -5, wt = avg_log2FC) %>% 
  as.data.frame()

top_label <- rbind(top_up_label,top_down_label)

top_label$thr_signi %>% 
  factor(., levels = c("Up_Highly","Down_Highly","Up_Lowly","Down_Lowly"))

colnames(r.deg)

# 04. 绘制灰色背景和色带
background_position <- r.deg %>%
  dplyr::group_by(unm) %>%
  dplyr::summarise(Min = min(avg_log2FC) - 0.2, Max = max(avg_log2FC) + 0.2) %>%
  # 灰色背景最上面比最高的点高，最低的要低（0.2）
  as.data.frame()

## `summarise()` ungrouping output (override with `.groups` argument)
background_position$unm %<>% as.vector(.) %>% as.numeric(.)
background_position$start <- background_position$unm - 0.4
background_position$end <- background_position$unm + 0.4

### 准备绘制中间区域cluster彩色bar所需数据
cluster_bar_position <- background_position
cluster_bar_position$start <- cluster_bar_position$unm - 0.5
cluster_bar_position$end <- cluster_bar_position$unm + 0.5
cluster_bar_position$unm %>% 
  factor(., levels = c(0:max(as.vector(.))))

# 05. 画图
# 设置填充颜色
cols_thr_signi <- c("Up_Highly" = "#d7301f","Down_Highly" = "#225ea8","Up_Lowly" = "black", "Down_Lowly" = "black")

cols_cluster <- c("0" = mycol[1],"1" = mycol[2],"2" = mycol[7],"3" = mycol[6],
                  "4" = mycol[5],"5" = mycol[4],"6" = mycol[8],"7" = mycol[9],
                  "8" = mycol[2])

# ann_colors <- list(cluster = c('Stromal' = mycol[1],'Epithelial' = mycol[2],'T_lymphocytes' = mycol[3],
#                                'NK_cells' = mycol[4],'Endothelial' = mycol[5],'Myeloid_cells' = mycol[6],
#                                'Neutrophil' = mycol[7],'Others' = mycol[8]))


p = ggplot() +
  geom_rect(data = background_position, aes(xmin = start, xmax = end, ymin = Min,
                                            ymax = Max),
            fill = "#525252", alpha = 0.1) + ###添加灰色背景色
  geom_jitter(data = r.deg, aes(x =unm, y = avg_log2FC, colour = thr_signi),
              size = 1,position = position_jitter(seed = 1)) +
  scale_color_manual(values = cols_thr_signi) +
  scale_x_continuous(limits = c(-0.5, max(r.deg$unm) + 0.5),
                     breaks = seq(0, max(r.deg$unm), 1),
                     label = seq(0, max(r.deg$unm),1)) + #修改坐标轴显示刻度
  # 根据top_label标注基因名
  geom_text_repel(data = top_label, aes(x =unm, y = avg_log2FC, label = gene),
                  position = position_jitter(seed = 1), show.legend = F, size = 2.5,
                  box.padding = unit(0, "lines")) +
  geom_rect(data = cluster_bar_position, aes(xmin = start, xmax = end, ymin = -0.2,
                                             ymax = 0.2, fill = factor(unm)), color = "black", 
            alpha = 1, show.legend = F) +
  scale_fill_manual(values = cols_cluster) +
  labs(x = "CellType", y = "average log2FC") +
  theme_bw()

plot1 <- p + theme(panel.grid.minor = element_blank(), ##去除网格线
                   panel.grid.major = element_blank(),
                   axis.text.y = element_text(colour = 'black', size = 14),
                   axis.text.x = element_text(colour = 'black', size = 14, vjust =55), 
                   panel.border = element_blank(), ## 去掉坐标轴
                   axis.ticks.x = element_blank(), ## 去掉的坐标刻度线
                   axis.line.y = element_line(colour = "black")) #添加y轴坐标轴

# plot1
ggsave(filename = "Diff_expr_gene.pdf", plot = plot1, width = 8, height = 6)



deg=FindMarkers(scRNA_harmony,ident.1 = "AERD",ident.2 = "CRSwNP",group.by = "group", subset.ident ='Tcells', logfc.threshold = 0.1, min.pct = 0.1 )
write.csv(deg, file = 'AERD_vs_CRSwNP_Tcells_degs.csv') 


deg=FindMarkers(scRNA_harmony,ident.1 = "Tumor",ident.2 = "NT",group.by = "group", subset.ident ='Endothelial')
write.csv(deg, file = 'Endothelial_Tumor_vs_NT_deg.csv') 




#-------------------------------------------------------------------------------


########## Part7. Malignant细分 ##########
library('harmony')

sc.subtype <- scRNA_harmony_T
Idents(sc.subtype) <- 'celltype'
table(Idents(sc.subtype))

# 提取感兴趣的细胞类型
sc.subtype.1 <- subset(sc.subtype,ident=c("Malignant_cells"))

# 细胞亚群1
# 由于有多个样本，提取后是针对data，要再次harmony
library(harmony)
sc.subtype.1 <- sc.subtype.1
sc.subtype.1 <- NormalizeData(sc.subtype.1) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
sc.subtype.1 <- RunHarmony(sc.subtype.1, group.by.vars = "orig.ident")
sc.subtype.1 <- FindNeighbors(sc.subtype.1, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 0.3)
sc.subtype.1 <- RunUMAP(sc.subtype.1, reduction = "harmony", dims = 1:15)
sc.subtype.1 <- RunTSNE(sc.subtype.1, reduction = "harmony", dims = 1:15)

setwd('/Users/luocheng/Desktop/单细胞分析/0322-HCC-55w-paper')

# 差异基因计算
# 计算各Cluster markers
markers <- FindAllMarkers(sc.subtype.1, only.pos = TRUE, min.pct = 0.15)
write.csv(markers,"Malignant_all_clusters_specific_expression.csv",row.names = FALSE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,"Malignant_all_clusters_specific_top10_expression.csv",row.names = FALSE)

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20,"Malignant_all_clusters_specific_top20_expression.csv",row.names = FALSE)


# 【可视化】样本整合后的图
p1 <- DimPlot(sc.subtype.1, reduction = "umap", label = T, cols = mycol)  # 4*4 Landscape
p2 <- DimPlot(sc.subtype.1, reduction = "umap", split.by ='MVI', cols = mycol) # 6*4 Landscape
p3 <- DimPlot(sc.subtype.1, reduction = "umap", group.by ='MVI', cols = mycol) # 4*4 Landscape
p4 <- DimPlot(sc.subtype.1, reduction = "tsne", label = T, cols = mycol)  # 4*4 Landscape
p5 <- DimPlot(sc.subtype.1, reduction = "tsne", split.by ='MVI', cols = mycol) # 6*4 Landscape
p6 <- DimPlot(sc.subtype.1, reduction = "tsne", group.by ='MVI', cols = mycol) # 4*4 Landscape

ggsave(filename = "Malignant_UMAP_by_cluster.pdf", plot = p1, width = 5, height = 4)
ggsave(filename = "Malignant_UMAP_by_cluster_split_MVI.pdf", plot = p2, width = 8, height = 4)
ggsave(filename = "Malignant_UMAP_by_cluster_merge_MVI.pdf", plot = p3, width = 5, height = 4)
ggsave(filename = "Malignant_TSNE_by_cluster.pdf", plot = p4, width = 5, height = 4)
ggsave(filename = "Malignant_TSNE_by_cluster_split_MVI.pdf", plot = p5, width = 8, height = 4)
ggsave(filename = "Malignant_TSNE_by_cluster_merge_MVI.pdf", plot = p6, width = 5, height = 4)


sc.subtype.1 <- RenameIdents(sc.subtype.1,"0"="MCs_1","1"="MCs_2","2"="MCs_3","3"="MCs_4",
                             '4'='MCs_5',"5"="MCs_6","6"="MCs_7", "7"="MCs_8","8"="MCs_9")

sc.subtype.1@meta.data$celltype = Idents(sc.subtype.1)

save(sc.subtype.1, file = 'scRNA_harmony_Malignant.Rda')

# 【可视化】样本整合后的图
p1 <- DimPlot(sc.subtype.1, reduction = "umap", label = T, cols = mycol)  # 4*4 Landscape
p2 <- DimPlot(sc.subtype.1, reduction = "umap", split.by ='MVI', cols = mycol) # 6*4 Landscape
p3 <- DimPlot(sc.subtype.1, reduction = "umap", group.by ='MVI', cols = mycol) # 4*4 Landscape
p4 <- DimPlot(sc.subtype.1, reduction = "tsne", label = T, cols = mycol)  # 4*4 Landscape
p5 <- DimPlot(sc.subtype.1, reduction = "tsne", split.by ='MVI', cols = mycol) # 6*4 Landscape
p6 <- DimPlot(sc.subtype.1, reduction = "tsne", group.by ='MVI', cols = mycol) # 4*4 Landscape

ggsave(filename = "Malignant_Anno_UMAP_by_cluster.pdf", plot = p1, width = 5, height = 4)
ggsave(filename = "Malignant_Anno_UMAP_by_cluster_split_sample.pdf", plot = p2, width = 8, height = 4)
ggsave(filename = "Malignant_Anno_UMAP_by_cluster_merge_group.pdf", plot = p3, width = 5, height = 4)
ggsave(filename = "Malignant_Anno_TSNE_by_cluster.pdf", plot = p4, width = 5, height = 4)
ggsave(filename = "Malignant_Anno_TSNE_by_cluster_split_sample.pdf", plot = p5, width = 8, height = 4)
ggsave(filename = "Malignant_Anno_TSNE_by_cluster_merge_group.pdf", plot = p6, width = 5, height = 4)


# 比例
# Group比例
pB2_df <- table(sc.subtype.1@meta.data$celltype,sc.subtype.1@meta.data$MVI) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")

# 【可视化】画比例图 4*6
pB4group <- ggplot(data = pB2_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Malignant_Proportion_MVI.pdf", plot = pB4group, width = 6, height = 4)

p1 <- ggplot(data = pB2_df, aes(x =Cluster, y = Number, fill =  Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Malignant_Proportion_MVI_2.pdf", plot = p1, width = 6, height = 4)


# Orig.ident比例
pB2_df <- table(sc.subtype.1@meta.data$celltype,sc.subtype.1@meta.data$orig.ident) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")

# 【可视化】画比例图 4*6
pB4sample <- ggplot(data = pB2_df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Malignant_Proportion_Samples.pdf", plot = pB4sample, width = 6, height = 4)


p1 <- ggplot(data = pB2_df, aes(x =Cluster, y = Number, fill =  Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Malignant_Proportion_Samples_2.pdf", plot = p1, width = 6, height = 4)


# inferCNV对象制作
Idents(scRNA_harmony_T) <- 'celltype'
# Myeloid_cells    Malignant_cells     T_cells     B_cells   NK_cells    Endothelial   Mesenchymal_cells   HPCs   MAST

scRNA_harmony_immu <- subset(scRNA_harmony_T, ident=c('T_cells','B_cells','NK_cells'))
scRNA_harmony_inferCNV <- merge(scRNA_harmony_immu, sc.subtype.1)

save(scRNA_harmony_inferCNV, file = 'scRNA_harmony_Anno_CNV.Rda')


#-------------------------------------------------------------------------------



########## Part8. T细胞及髓系细胞细分 ##########
library('harmony')

setwd('/Users/luocheng/Desktop/单细胞分析/0322-HCC-55w-paper/10_Myeloid_subtype')

sc.subtype <- scRNA_harmony_T
Idents(sc.subtype) <- 'celltype'
table(Idents(sc.subtype))

# 提取感兴趣的细胞类型
sc.subtype.1 <- subset(sc.subtype,ident=c("Myeloid_cells"))

# 细胞亚群1
# 由于有多个样本，提取后是针对data，要再次harmony
library(harmony)
sc.subtype.1 <- sc.subtype.1
sc.subtype.1 <- NormalizeData(sc.subtype.1) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
sc.subtype.1 <- RunHarmony(sc.subtype.1, group.by.vars = "orig.ident")
sc.subtype.1 <- FindNeighbors(sc.subtype.1, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 0.3)
sc.subtype.1 <- RunUMAP(sc.subtype.1, reduction = "harmony", dims = 1:15)
sc.subtype.1 <- RunTSNE(sc.subtype.1, reduction = "harmony", dims = 1:15)


Idents(sc.subtype.1) <- 'seurat_clusters'

# Monocytes	IL1B、FPR1
Monocytes <- c('S100A8','S100A9')
p <- DotPlot(sc.subtype.1, features = Monocytes, scale=F, cols = mycol) + RotatedAxis()
ggsave(filename = "Marker_plot_Monocytes.pdf", plot = p, width = 4, height = 6)

# DCs	CLU、AREG
Dendritic_cells <- c('CD1C','CD1E')
p <- DotPlot(sc.subtype.1, features = Dendritic_cells, scale=F, cols = mycol) + RotatedAxis()
ggsave(filename = "Marker_plot_Dendritic_cells.pdf", plot = p, width = 4, height = 6)

# Macrophages	CD163、CD86
Macrophages <- c('APOE','C1QA','C1QB')
p <- DotPlot(sc.subtype.1, features = Macrophages, scale=F, cols = mycol) + RotatedAxis()
ggsave(filename = "Marker_plot_Macrophages.pdf", plot = p, width = 4, height = 6)


Idents(sc.subtype.1) <- 'seurat_clusters'
# Monocytes
# Dendritic_cells
# Macrophages

sc.subtype.1 <- RenameIdents(sc.subtype.1,"0"="Macrophages","1"="DCs","2"="Macrophages","3"="Monocytes",
                             '4'='Macrophages','5'='Macrophages','6'='Macrophages','7'='DCs','8'='Macrophages')

sc.subtype.1@meta.data$celltype = Idents(sc.subtype.1)

save(sc.subtype.1, file = 'scRNA_harmony_Myeloid_Anno.Rda')


p1 <- DimPlot(sc.subtype.1, reduction = "umap", label = T, cols = mycol)  # 4*4 Landscape
p2 <- DimPlot(sc.subtype.1, reduction = "umap", split.by ='MVI', cols = mycol) # 6*4 Landscape
p3 <- DimPlot(sc.subtype.1, reduction = "umap", group.by ='MVI', cols = mycol) # 4*4 Landscape
p4 <- DimPlot(sc.subtype.1, reduction = "tsne", label = T, cols = mycol)  # 4*4 Landscape
p5 <- DimPlot(sc.subtype.1, reduction = "tsne", split.by ='MVI', cols = mycol) # 6*4 Landscape
p6 <- DimPlot(sc.subtype.1, reduction = "tsne", group.by ='MVI', cols = mycol) # 4*4 Landscape

ggsave(filename = "Myeloid_Anno_UMAP_by_anno.pdf", plot = p1, width = 5, height = 4)
ggsave(filename = "Myeloid_Anno_UMAP_by_anno_split_sample.pdf", plot = p2, width = 8, height = 4)
ggsave(filename = "Myeloid_Anno_UMAP_by_anno_merge_group.pdf", plot = p3, width = 5, height = 4)
ggsave(filename = "Myeloid_Anno_TSNE_by_anno.pdf", plot = p4, width = 5, height = 4)
ggsave(filename = "Myeloid_Anno_TSNE_by_anno_split_sample.pdf", plot = p5, width = 8, height = 4)
ggsave(filename = "Myeloid_Anno_TSNE_by_anno_merge_group.pdf", plot = p6, width = 5, height = 4)


sc.subtype.2 <- subset(sc.subtype,ident=c("T_cells","B_cells","NK_cells","Endothelial"))
save(sc.subtype.2, file = 'scRNA_harmony_Tumor_Immune_Anno.Rda')


#-------------------------------------------------------------------------------

library(Seurat)
library(SeuratDisk)

scRNA_harmony

obj = DietSeurat(object = scRNA_harmony, counts = T, data = T, scale.data = F, assays = "RNA", dimreducs =c('pca','tsne','umap','harmony'))

SaveH5Seurat(obj, filename = "scRNA_harmony_2.h5Seurat")

Convert("scRNA_harmony_2.h5Seurat", dest = "h5ad")


Idents(sc.subtype.1) <- 'orig.ident'

sc.subtype.1 <- RenameIdents(sc.subtype.1,"1NT_P21"="NT","1T_C21"="Tumor","2NT_P24"="NT","2T_C24"="Tumor","3NT_P25"="NT",
                             "3T_C25"="Tumor","4NT_P29"="NT","4T_C29"="Tumor","5NT_P36"="NT","5T_C36"="Tumor")

table(sc.subtype.1@active.ident)
# NT MVI_T Tumor 
# 9608  1333   406 

# 将当前的active.ident写入metadata
sc.subtype.1@meta.data$group <- sc.subtype.1@active.ident



# 淋巴系细胞细分
library('harmony')

Idents(scRNA_harmony) <- 'celltype'

table(Idents(scRNA_harmony))

# 提取感兴趣的细胞类型
sc.subtype.1 <- subset(scRNA_harmony,ident=c("T_cells",'B_cells','NK_cells'))


# 细胞亚群1
# 由于有多个样本，提取后是针对data，要再次harmony
library(harmony)
sc.subtype.1 <- sc.subtype.1
sc.subtype.1 <- NormalizeData(sc.subtype.1) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
sc.subtype.1 <- RunHarmony(sc.subtype.1, group.by.vars = "orig.ident")
sc.subtype.1 <- FindNeighbors(sc.subtype.1, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.6)
sc.subtype.1 <- RunUMAP(sc.subtype.1, reduction = "harmony", dims = 1:30)
sc.subtype.1 <- RunTSNE(sc.subtype.1, reduction = "harmony", dims = 1:30)


DimPlot(sc.subtype.1, reduction = "umap", label = T, cols = mycol)

save(sc.subtype.1, file = 'scRNA_harmony_Lym.Rda')

# 差异基因计算
# 计算各Cluster markers
markers <- FindAllMarkers(sc.subtype.1, only.pos = TRUE, min.pct = 0.15)
write.csv(markers,"Lym_all_clusters_specific_expression.csv",row.names = FALSE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,"Lym_all_clusters_specific_top10_expression.csv",row.names = FALSE)

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20,"Lym_all_clusters_specific_top20_expression.csv",row.names = FALSE)

Idents(sc.subtype.1) <- 'seurat_clusters'

sc.subtype.1 <- subset(sc.subtype.1, ident = c('0','1','2','3','4','5','6','8','10','11'))

sc.subtype.1 <- RenameIdents(sc.subtype.1,"0"="CD4T","1"="Plasma","2"="aNK","3"="ProL",
                             '4'='nB','5'='CTL','6'='IgL-Plasma','8'='CD8T','10'='Treg','11'='secB')

sc.subtype.1@meta.data$celltype = Idents(sc.subtype.1)

save(sc.subtype.1, file = 'scRNA_harmony_Lym.Rda')

p1 <- DimPlot(sc.subtype.1, reduction = "umap", label = F, cols = mycol)  # 4*4 Landscape
p2 <- DimPlot(sc.subtype.1, reduction = "umap", split.by ='group', cols = mycol) # 6*4 Landscape
p3 <- DimPlot(sc.subtype.1, reduction = "umap", group.by ='group', cols = mycol) # 4*4 Landscape
p4 <- DimPlot(sc.subtype.1, reduction = "tsne", label = F, cols = mycol)  # 4*4 Landscape
p5 <- DimPlot(sc.subtype.1, reduction = "tsne", split.by ='group', cols = mycol) # 6*4 Landscape
p6 <- DimPlot(sc.subtype.1, reduction = "tsne", group.by ='group', cols = mycol) # 4*4 Landscape
p7 <- DimPlot(sc.subtype.1, reduction = "umap", label = T, cols = mycol)  # 4*4 Landscape
p8 <- DimPlot(sc.subtype.1, reduction = "tsne", label = T, cols = mycol)  # 4*4 Landscape

ggsave(filename = "Lym_UMAP_by_cluster.pdf", plot = p1, width = 5.5, height = 4)
ggsave(filename = "Lym_UMAP_by_cluster_split_sample.pdf", plot = p2, width = 6, height = 3)
ggsave(filename = "Lym_UMAP_by_cluster_merge_group.pdf", plot = p3, width = 4.5, height = 4)
ggsave(filename = "Lym_TSNE_by_cluster.pdf", plot = p4, width = 5.5, height = 4)
ggsave(filename = "Lym_TSNE_by_cluster_split_sample.pdf", plot = p5, width = 6, height = 3)
ggsave(filename = "Lym_TSNE_by_cluster_merge_group.pdf", plot = p6, width = 4.5, height = 4)
ggsave(filename = "Lym_UMAP_Label.pdf", plot = p7, width = 5.5, height = 4)
ggsave(filename = "Lym_TSNE_Label.pdf", plot = p8, width = 5.5, height = 4)

Idents(sc.subtype.1) <- 'celltype'

# 比例
# Group比例
pB2_df <- table(sc.subtype.1@meta.data$celltype,sc.subtype.1@meta.data$MVIgroup) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")

# 【可视化】画比例图 4*6
pB4group <- ggplot(data = pB2_df, aes(x =Sample, y = Number, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Proportion_group.pdf", plot = pB4group, width = 4, height = 6)



library(Augur)

# 分组进行Augur
Idents(sc.subtype.1) <- 'MVIgroup'
# NT  MVI_T  Tumor 


# 01 NT - Tumor
augur <- subset(sc.subtype.1, ident = c('NT','Tumor'))

augur <- subset(sc.subtype.1, ident = c('MVI_T','Tumor'))

# augur
augur <- calculate_auc(augur, cell_type_col = 'celltype', label_col = 'MVIgroup')
augur$AUC
augur.auc <- as.data.frame(augur$AUC)
ordered_types <- unique(augur.auc$cell_type)

# 【可视化】Augur的AUC结果
augur.auc$cell_type_status <- ifelse(augur.auc$auc > 0.75, 'sig', 'nonsig')

p <- ggplot(augur.auc, aes(x=cell_type, y=auc)) + 
  geom_segment(aes(x=cell_type, xend=cell_type, y=0, yend=auc, color=cell_type_status), linewidth=1.5) +
  geom_point(aes(color=cell_type_status), size=2.3) +
  scale_x_discrete(limits=augur.auc$cell_type) +
  scale_color_manual(values = c('sig' = "#e41a1b", 'nonsig' = "#3182bd")) +
  geom_hline(yintercept=0.65, linetype = "dashed") +
  coord_flip() +
  theme_bw() +
  labs(title='Augur MVI vs Tumor')

ggsave(filename = "Augur_Lym_MVI_vs_Tumor.pdf", plot = p, width = 4, height = 3)

write.csv(augur.auc, file = 'Augur_Lym_MVI_vs_Tumor.csv')


library(Seurat)
library(SeuratDisk)
obj = DietSeurat(object = sc.subtype.1, counts = T, data = T, scale.data = F, assays = "RNA", dimreducs =c('pca','tsne','umap','harmony'))
SaveH5Seurat(obj, filename = "scRNA_harmony_Lym.h5Seurat")
Convert("scRNA_harmony_Lym.h5Seurat", dest = "h5ad")

Idents(sc.subtype.1) <- 'celltype'

# CD8T Tumor NT group
# CD4T MVI_T Tumor MVIgroup
deg <- FindMarkers(sc.subtype.1,ident.1 = "MVI_T",ident.2 = "Tumor",group.by = "MVIgroup",subset.ident ='CD4T')
write.csv(deg, file = 'CD4T_MVI_vs_Tumor_MVIgroup_DEGs.csv')

library(Seurat)
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)

# GO
degs.list <- deg
degs.list$SYMBOL <- rownames(degs.list)

# degs.list <- degs.list[degs.list$p_val < 0.05 & degs.list$avg_log2FC >0,]
degs.list <- degs.list[degs.list$p_val < 0.05 & degs.list$avg_log2FC <0,]
degs.list <- degs.list$SYMBOL

# GO
erich.go <- enrichGO(gene =degs.list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL", pvalueCutoff = 0.5, qvalueCutoff = 0.5)
erich.go.results <- erich.go@result
erich.go.results$log10p <- -log10(erich.go.results$p.adjust)


library(ggplot2)
library(dplyr)
library(ggtext)

erich.go.results <- erich.go.results %>% arrange(desc(log10p))


erich.go.results.10 <- erich.go.results[c(1:5),]


p <-ggplot(data = erich.go.results.10) +
  geom_bar(mapping = aes(x = reorder(Description, log10p), y = log10p, fill = Description), stat = 'identity') +
  coord_flip() +
  theme(aspect.ratio = 1/1) +
  scale_fill_manual(values = mycol) +
  labs(x = "Description", y = "-log10(p value)", title = "GO Term Enrichment")+
  theme(legend.position = "none") +
  theme(
    legend.position = "none",  # 删除图注
    panel.background = element_blank(),  # 去除面板背景
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank(),  # 去除次网格线
    panel.border = element_blank(),  # 去除面板边框
    axis.line = element_line(color = "black"),  # 添加轴线
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_markdown(size = 10)  # 使用element_markdown来渲染高亮文本
  )

ggsave(filename = "CD4T_Down-DEG_MVI_vs_Tumor.pdf", plot = p, width = 6, height = 6)






# 髓系细胞细分
library('harmony')

Idents(scRNA_harmony) <- 'celltype'

table(Idents(scRNA_harmony))

# 提取感兴趣的细胞类型
sc.subtype.1 <- subset(scRNA_harmony,ident=c("Myeloid_cells"))

# 细胞亚群1
# 由于有多个样本，提取后是针对data，要再次harmony
library(harmony)
sc.subtype.1 <- sc.subtype.1
sc.subtype.1 <- NormalizeData(sc.subtype.1) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
sc.subtype.1 <- RunHarmony(sc.subtype.1, group.by.vars = "orig.ident")
sc.subtype.1 <- FindNeighbors(sc.subtype.1, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 0.3)
sc.subtype.1 <- RunUMAP(sc.subtype.1, reduction = "harmony", dims = 1:15)
sc.subtype.1 <- RunTSNE(sc.subtype.1, reduction = "harmony", dims = 1:15)


DimPlot(sc.subtype.1, reduction = "umap", label = T, cols = mycol)

save(sc.subtype.1, file = 'scRNA_harmony_Lym.Rda')

# 差异基因计算
# 计算各Cluster markers
markers <- FindAllMarkers(sc.subtype.1, only.pos = TRUE, min.pct = 0.15)
write.csv(markers,"Mye_all_clusters_specific_expression.csv",row.names = FALSE)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,"Mye_all_clusters_specific_top10_expression.csv",row.names = FALSE)

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20,"Mye_all_clusters_specific_top20_expression.csv",row.names = FALSE)

Idents(sc.subtype.1) <- 'seurat_clusters'

sc.subtype.1 <- subset(sc.subtype.1, ident = c('0','1','2','3','5','6','7'))

Idents(sc.subtype.1) <- 'seurat_clusters'

sc.subtype.1 <- RenameIdents(sc.subtype.1,"0"="TRM","1"="cDC","2"="actMac","3"="Mono",
                             '5'='cDC1','6'='actDC','7'='pDC')

sc.subtype.1@meta.data$celltype = Idents(sc.subtype.1)

save(sc.subtype.1, file = 'scRNA_harmony_Mye.Rda')

p1 <- DimPlot(sc.subtype.1, reduction = "umap", label = F, cols = mycol)  # 4*4 Landscape
p2 <- DimPlot(sc.subtype.1, reduction = "umap", split.by ='group', cols = mycol) # 6*4 Landscape
p3 <- DimPlot(sc.subtype.1, reduction = "umap", group.by ='group', cols = mycol) # 4*4 Landscape
p4 <- DimPlot(sc.subtype.1, reduction = "tsne", label = F, cols = mycol)  # 4*4 Landscape
p5 <- DimPlot(sc.subtype.1, reduction = "tsne", split.by ='group', cols = mycol) # 6*4 Landscape
p6 <- DimPlot(sc.subtype.1, reduction = "tsne", group.by ='group', cols = mycol) # 4*4 Landscape
p7 <- DimPlot(sc.subtype.1, reduction = "umap", label = T, cols = mycol)  # 4*4 Landscape
p8 <- DimPlot(sc.subtype.1, reduction = "tsne", label = T, cols = mycol)  # 4*4 Landscape

ggsave(filename = "Mye_UMAP_by_cluster.pdf", plot = p1, width = 5.5, height = 4)
ggsave(filename = "Mye_UMAP_by_cluster_split_sample.pdf", plot = p2, width = 6, height = 3)
ggsave(filename = "Mye_UMAP_by_cluster_merge_group.pdf", plot = p3, width = 4.5, height = 4)
ggsave(filename = "Mye_TSNE_by_cluster.pdf", plot = p4, width = 5.5, height = 4)
ggsave(filename = "Mye_TSNE_by_cluster_split_sample.pdf", plot = p5, width = 6, height = 3)
ggsave(filename = "Mye_TSNE_by_cluster_merge_group.pdf", plot = p6, width = 4.5, height = 4)
ggsave(filename = "Mye_UMAP_Label.pdf", plot = p7, width = 5.5, height = 4)
ggsave(filename = "Mye_TSNE_Label.pdf", plot = p8, width = 5.5, height = 4)


# 比例
# Group比例
pB2_df <- table(sc.subtype.1@meta.data$celltype,sc.subtype.1@meta.data$MVIgroup) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")

# 【可视化】画比例图 4*6
pB4group <- ggplot(data = pB2_df, aes(x =Sample, y = Number, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=mycol) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)) 

ggsave(filename = "Proportion_group.pdf", plot = pB4group, width = 4, height = 6)


library(Seurat)
library(SeuratDisk)
obj = DietSeurat(object = sc.subtype.1, counts = T, data = T, scale.data = F, assays = "RNA", dimreducs =c('pca','tsne','umap','harmony'))
SaveH5Seurat(obj, filename = "scRNA_harmony_Mye.h5Seurat")
Convert("scRNA_harmony_Mye.h5Seurat", dest = "h5ad")

library(Augur)

# 分组进行Augur
Idents(sc.subtype.1) <- 'group'


# 01 NT - Tumor
augur <- subset(sc.subtype.1, ident = c('NT','Tumor'))

# augur
augur <- calculate_auc(augur, cell_type_col = 'celltype', label_col = 'group')
augur$AUC
augur.auc <- as.data.frame(augur$AUC)
ordered_types <- unique(augur.auc$cell_type)

# 【可视化】Augur的AUC结果
augur.auc$cell_type_status <- ifelse(augur.auc$auc > 0.70, 'sig', 'nonsig')

p <- ggplot(augur.auc, aes(x=cell_type, y=auc)) + 
  geom_segment(aes(x=cell_type, xend=cell_type, y=0, yend=auc, color=cell_type_status), linewidth=1.5) +
  geom_point(aes(color=cell_type_status), size=2.3) +
  scale_x_discrete(limits=augur.auc$cell_type) +
  scale_color_manual(values = c('sig' = "#e41a1b", 'nonsig' = "#3182bd")) +
  geom_hline(yintercept=0.65, linetype = "dashed") +
  coord_flip() +
  theme_bw() +
  labs(title='Augur Tumor vs Normal')

ggsave(filename = "Augur_Mye_Tumor vs Normal.pdf", plot = p, width = 4, height = 4)

write.csv(augur.auc, file = 'Augur_Mye_Tumor_vs_Normal.csv')


Idents(sc.subtype.1) <- 'celltype'
deg <- FindMarkers(sc.subtype.1,ident.1 = "Tumor",ident.2 = "NT",group.by = "group",subset.ident ='TRM')
write.csv(deg, file = 'TRM_Tumor_vs_NT_group_DEGs.csv')





library(ggplot2)
library(ggpubr)

genes <- c('APOA1', 'APOA2', 'APOC1', 'APOC3', 'APOE',
           'RBP4','FABP5','GC','HLA-A','HLA-B','IGKC',
           'S100A8','HMOX1','IFI6', 'IFI27', 'IFI44L')

for (i in 1:length(genes)) {
  pdf(paste('Genes_expr_Vln_',genes[i], ".pdf", sep=""), width = 4.5, height = 3.5)
  p <- VlnPlot(sc.subtype.1, genes[i], group.by = "group", cols = mycol[2:1]) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    ) +
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.7, outlier.shape = NA) +
    labs(title = paste0(genes[i],' Expression by Group'), x = "Gender", y = paste0(genes[i], ' Expression')) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    stat_compare_means(method = 't.test')
  print(p)
  dev.off()
}


write.csv(scRNA_harmony@meta.data, file = 'Supp_Data1.csv')



