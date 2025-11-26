source("00_Setup_and_Functions.R")

# ==============================================================================
# 1. Preprocessing & Integration
# ==============================================================================
# 1.1 Subset B cells from All Cells
# Assuming 'allcells_final' exists or load raw B cell data
bcell_obj <- subset(LoadH5Seurat("data/processed/allcells_clean.h5seurat"), celltype_l1 == 'Bcell')

# 1.2 Integration (FastMNN)
bcell_clean <- quick_fastmnn(bcell_obj, k = 20, reduction.dims = 20, findcluster.res = 0.6)

# ==============================================================================
# 2. Annotation (Developmental Stages)
# ==============================================================================
# 2.1 Define Markers
b_markers <- list(
  Pan_B = c('CD19', 'MS4A1', 'CD79A'),
  Pro_Pre_B = c('CD34', 'RAG1', 'RAG2', 'VPREB1', 'EBF1'), # Early development
  Immature_B = c('IL7R', 'MME'), # CD10
  Naive_B = c('FCER2', 'TCL1A', 'IL4R', 'IGHD'), # IgD+
  Memory_B = c('CD27', 'AIM2', 'TNFRSF13B'), # CD27+
  Plasma = c('XBP1', 'SDC1', 'MZB1', 'IGHG1'), # CD138
  Cycling = c('MKI67', 'TOP2A')
)

# 2.2 Visualization
dotp(bcell_clean, unlist(b_markers), group.by = "seurat_clusters", save.name = "figures/dot_bcell_clusters")

# 2.3 Rename Idents
# Mapping based on your analysis:
# preB (Pre-B), immB (Immature), Bn (Naive), Bm (Memory), cycleB (Cycling)
new_ids <- c(
  "0" = "Bn", "1" = "Bm", "2" = "Bn", 
  "3" = "immB", "4" = "cycleB", "5" = "preB"
)
bcell_clean <- RenameIdents_my(bcell_clean, new_ids)
bcell_clean$celltype_l2 <- Idents(bcell_clean)
# Order levels for plotting
bcell_clean$celltype_l2 <- factor(bcell_clean$celltype_l2, levels = c('preB','immB','Bn','Bm','cycleB'))

# 2.4 UMAP Visualization
bcell_colors <- c(preB='#0fa3b1', immB='#b5e2fa', Bn='#eddea4', Bm='#f7a072', cycleB='#1e96fc')
DimPlot(bcell_clean, group.by = "celltype_l2", cols = bcell_colors, pt.size = 0.8) + theme_void()
ggsave("figures/umap_bcell_final.pdf", width = 5, height = 5)

# ==============================================================================
# 3. Tissue Preference Analysis (Odds Ratio)
# ==============================================================================
# Calculate distribution preference across sample types (e.g., Healthy vs Mets)
# Using custom function 'calculate_odds_ratios' (logic preserved from your scripts)

bcell_or <- calculate_odds_ratios(
  seurat_obj = bcell_clean, 
  cluster_column = 'celltype_l2', 
  tissue_column = 'sample_type' # Ensure this meta.data column exists (Healthy, Benign, preMets, Mets)
)

# Heatmap of Odds Ratios
# High OR (>1) indicates enrichment, Low OR (<1) indicates depletion
p1 <- plotORheat2(
  bcell_or, 
  colors = paletteer_c("grDevices::cm.colors", 100), 
  cellwidth = 20, cellheight = 20, 
  cluster_cols = FALSE, cluster_rows = TRUE
)
ggsave("figures/heat_OR_bcell_tissue_distribution.pdf", p1, width = 6, height = 4)

# ==============================================================================
# 4. Save Final Object
# ==============================================================================
SaveH5Seurat(bcell_clean, "processed_data/Bcell_Final.h5seurat")
