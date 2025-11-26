source("00_Setup_and_Functions.R")

# ==============================================================================
# 1. Preprocessing & Integration
# ==============================================================================
# 1.1 Load Data (e.g., subset from AllCells or raw objects)
# Combining Cohort 1 and Cohort 2 tumor samples
epi_c1 <- LoadH5Seurat("data/raw/epi_c1_raw.h5seurat")
epi_c2 <- LoadH5Seurat("data/raw/epi_c2_mets.h5seurat")

# Merge and Integration
epi_merge <- merge_intersect_genes(epi_c1, epi_c2) %>%
             quick_fastmnn(k = 20, reduction.dims = 30, findcluster.res = 0.6)

# 1.2 Quality Control (Filtering Low Quality Clusters)
# Based on your notes: Remove clusters with low nFeature or high MALAT1 (dead cells)
# Example: Cluster 7, 10, 11 in your raw analysis
low_quality_clusters <- c("7", "10", "11") 
epi_clean <- subset(epi_merge, idents = low_quality_clusters, invert = TRUE) %>%
             quick_fastmnn(k = 20, reduction.dims = 20)

# ==============================================================================
# 2. Annotation & Malignancy Definition
# ==============================================================================
# 2.1 Marker Visualization
epi_markers <- list(
  Basal = c("KRT5", "KRT14", "TP63"),
  Luminal_Normal = c("KRT8", "KRT18", "MSMB"),
  Luminal_Malignant = c("KLK3", "TMPRSS2", "AR", "AMACR"),
  Epi_EMT = c("VIM", "FN1") 
)
dotp(epi_clean, unlist(epi_markers), group.by = "seurat_clusters", save.name = "figures/dot_epi_clusters")

# 2.2 Assign Cell Types (Normal vs Malignant)
# Logic: Basal/Ductal are usually normal; Luminal can be normal or malignant
epi_clean <- RenameIdents_my(epi_clean, c(
  "0" = "Malignant_Luminal", "1" = "Normal_Luminal", 
  "2" = "Basal", "3" = "Ductal" 
  # ... add specific cluster mappings
))
epi_clean$celltype3 <- Idents(epi_clean)

# 2.3 InferCNV Analysis (Validation of Malignancy)
# Prepare Annotation File
write.table(data.frame(cell = colnames(epi_clean), group = epi_clean$celltype3), 
            "processed_data/infercnv_anno.txt", sep = "\t", col.names = F, quote = F)

# Run InferCNV (Wrapper for the standard pipeline)
# Reference: Basal and Normal_Luminal cells
runinferCNV(
  counts_matrix = GetAssayData(epi_clean, slot = "counts"),
  anno_file = "processed_data/infercnv_anno.txt",
  gene_order_file = "data/resources/hg38_gencode_v27.txt",
  ref_group_names = c("Basal", "Normal_Luminal"),
  output_dir = "processed_data/infercnv_output"
)

# ==============================================================================
# 3. Transcriptional Regulation (SCENIC)
# ==============================================================================
# 3.1 Load SCENIC Results (AUC Matrix)
# Assuming pySCENIC has been run externally on the loom file
auc_mtx <- read.csv("data/scenic/epi/AUC_results.csv", row.names = 1)
# Format: Rows = Regulons, Cols = Cells
epi_clean[['scenic']] <- CreateAssayObject(data = t(as.matrix(auc_mtx)))

# 3.2 Identify Regulon Markers for Subtypes
Idents(epi_clean) <- "celltype3"
tf_markers <- FindAllMarkers(epi_clean, assay = "scenic", only.pos = T, logfc.threshold = 0)

# 3.3 Visualization (Heatmap of TFs)
top_tfs <- tf_markers %>% group_by(cluster) %>% top_n(5, wt = avg_log2FC) %>% pull(gene)
p1 <- seurat_mean_heatmap(epi_clean, assay = "scenic", features = top_tfs, 
                          cluster = "celltype3", scale = "row", colors = modulescore_colors)
ggsave("figures/heat_epi_scenic_regulons.pdf", p1$plot)

# Specific TF Feature Plots (e.g., AR, FOXA1, HOXB13)
FeaturePlot(epi_clean, features = c("AR(+)", "HOXB13(+)", "FOXA1(+)"), reduction = "umap")

# ==============================================================================
# 4. Trajectory & Metabolism (scVelo & COMPASS)
# ==============================================================================
# 4.1 scVelo (Focus on Malignant Differentiation)
malignant_cells <- subset(epi_clean, celltype3 == "Malignant_Luminal")
# Export for Python scVelo
# adata <- scvelo_pip(malignant_cells, "epi_malignant", "data/raw/merged.loom", ".")

# 4.2 COMPASS (Metabolic Reaction Scores)
# Load COMPASS results (Reaction consistency scores)
compass_scores <- read.csv("data/compass/epi_reaction_scores.csv", row.names = 1)
epi_clean[['compass']] <- CreateAssayObject(data = as.matrix(compass_scores))

# Compare Metabolism between Epi Subtypes (e.g., Epi1 vs Epi2)
metabolic_degs <- FindMarkers(epi_clean, ident.1 = "Epi2", ident.2 = "Epi1", assay = "compass")
# Enrichment of Metabolic Pathways
# ... (Standard GSEA or ORA on significant reactions)

# Save Final Object
SaveH5Seurat(epi_clean, "processed_data/Epithelial_Final.h5seurat")