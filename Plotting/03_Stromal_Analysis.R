source("00_Setup_and_Functions.R")

# ==============================================================================
# 1. Osteoblasts (Osteo) Analysis
# ==============================================================================
# 1.1 Load & Integrate
osteo_obj <- LoadH5Seurat("data/raw/osteo_clean.h5seurat")
# Integration with Public (GSE147390, GSE253355)
public_osteo <- LoadH5Seurat("data/public/osteo_public_merge.h5seurat")
osteo_merge <- merge_intersect_genes(osteo_obj, public_osteo) %>% quick_fastmnn()

# 1.2 Annotation
# Differentiation trajectory: MSC -> Pre-OB -> OB1 -> OB2 -> OB3
osteo_markers <- list(
  MSC_like = c('LEPR','CXCL12'),
  Pre_OB = c('PRSS23','SP7'),
  OB1 = c('ALPL','RUNX2','IBSP'), # Immature/Matrix
  OB2 = c('SPP1','MMP13'),        # Mature
  OB3 = c('DMP1','SOX9')          # Osteocyte-like
)
dotp(osteo_merge, unlist(osteo_markers), group.by="celltype")

# 1.3 Trajectory Analysis (CytoTRACE2 & Monocle)
# CytoTRACE2 for differentiation potential
ct_res <- CytoTRACE2::cytotrace2(GetAssayData(osteo_merge, slot="data"), is_seurat=F)
osteo_merge$CytoTRACE <- ct_res$CytoTRACE2_Score
# Monocle 2
cds <- monocle_seurat(osteo_merge, order_genes = FindAllMarkers(osteo_merge)$gene)
plot_cell_trajectory(cds, color_by = "CytoTRACE")

# 1.4 SCENIC (Transcription Factors)
# Load AUC matrix from pySCENIC output
auc_mtx <- read.csv("data/scenic/osteo/AUC_results.csv")
osteo_merge[['scenic']] <- CreateAssayObject(data = t(auc_mtx))
# Correlate TFs with Pseudotime
tf_cor <- genes_cor_metadata_seurat(osteo_merge, assay="scenic", meta.data="Pseudotime")
# Plot key regulators (e.g., MEF2C, RUNX2)
FeaturePlot(osteo_merge, features = c("MEF2C(+)","RUNX2(+)"), reduction = "umap")

# 1.5 NicheNet: Cancer -> Osteoblast
# Determine Ligands from Cancer cells (Epi) affecting OB1/Pre-OB
receiver_cells <- c("OB1","Pre_OB")
sender_cells <- subset(allcells, celltype == "Epithelial" & sample_type == "Mets")
geneset_abnormal <- FindMarkers(osteo_merge, ident.1="Mets", ident.2="Healthy")$gene
# Run NicheNet wrapper
# nichenet_res <- run_nichenet(receiver=receiver_cells, sender=sender_cells, geneset=geneset_abnormal)
# Visualize Ligand-Target Heatmap
# nichenet_res$heatmap


# ==============================================================================
# 2. Endothelial & Fibroblasts
# ==============================================================================
# 2.1 Endothelial
endo_obj <- LoadH5Seurat("data/processed/endo_clean.h5seurat")
# Clusters: Artery, Vein, Capillary (CapEC), Tip cells
endo_markers <- list(Art=c('GJA5'), Ven=c('ACKR1'), Cap=c('CA4','RGCC'), Tip=c('ESM1','CXCR4'))
# NicheNet: Osteo -> Endo (Angiogenesis signals like VEGF)

# 2.2 Fibroblasts (CAFs)
fibro_obj <- LoadH5Seurat("data/processed/fibro_clean.h5seurat")
# Marker checks: POSTN, ASPN (Cancer Associated)
FeaturePlot(fibro_obj, features = c("POSTN","ASPN"), split.by = "sample_type")