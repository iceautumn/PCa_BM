source("00_Setup_and_Functions.R")

# ==============================================================================
# 1. CD8+ T Cells
# ==============================================================================
# 1.1 Load & Integrate (Cohorts 1 & 2 + Healthy + Benign)
cd8_c1 <- LoadH5Seurat("data/processed/cd8_c1_clean.h5seurat")
cd8_c2 <- LoadH5Seurat("data/processed/cd8_c2_bm.h5seurat")
# Merge logic
cd8_merge <- merge_intersect_genes(cd8_c1, cd8_c2) %>% 
             quick_fastmnn(k=20, reduction.dims=30, findcluster.res=0.6)

# 1.2 Annotation
cd8_markers <- list(
  Tn = c('CCR7','LEF1'), Tcm = c('IL7R','CD28'),
  Teff = c('CX3CR1','GZMA'), Tem = c('GZMK','CCL4'),
  Trm = c('ZNF683','ITGAE'), Tex = c('HAVCR2','PDCD1','TOX')
)
# Plot Markers
dotp(cd8_merge, features = unlist(cd8_markers), group.by = "seurat_clusters")
# Assign types
cd8_merge$celltype <- Idents(cd8_merge)

# 1.3 Functional Scoring (Exhaustion & Hypoxia)
exhaustion_genes <- c('PDCD1','TOX','TIGIT','LAG3','CTLA4')
hypoxia_genes <- msigdb_hallmarks[['HALLMARK_HYPOXIA']]
cd8_merge <- AddModuleScore_muti2(cd8_merge, list(Exhaustion=exhaustion_genes, Hypoxia=hypoxia_genes))

# 1.4 Trajectory (Monocle2) & Pseudo-time Correlation
# Order cells based on DEGs
diff_genes <- FindAllMarkers(cd8_merge, only.pos=T)$gene
cds <- monocle_seurat(cd8_merge, assay="soupx", order_genes = diff_genes)
plot_cell_trajectory(cds, color_by = "celltype")
# Correlate Pseudo-time with Exhaustion Score
plot_seurat_feature_lm_meta(cd8_merge, meta="Pseudotime", feature="Exhaustion", group="celltype")

# 1.5 Survival Analysis (SU2C)
# Use CD8 Teff signature to stratify SU2C patients
su2c_data <- readRDS("data/public/su2c_bulk.rds")
teff_sig <- cd8_markers$Teff
ssgsea_res <- GSVA::gsva(su2c_data$exp, list(Teff=teff_sig), method="ssgsea")
# Survival plot code (ggsurvplot) using su2c_data$clinical


# ==============================================================================
# 2. CD4+ T Cells
# ==============================================================================
# 2.1 Load & Integrate
cd4_merge <- LoadH5Seurat("data/processed/cd4_merge_final.h5seurat")

# 2.2 Annotation (Treg focus)
cd4_markers <- list(Tn=c('TCF7'), Tcm=c('FOSL2'), Tem=c('GZMK'), Treg=c('FOXP3','IL2RA'))
cd4_merge <- RenameIdents_my(cd4_merge, c("0"="Tn", "1"="Treg", "2"="Tem", ...))

# 2.3 Treg Correlation with Myeloid IL1B
# Extract Treg proportions per sample
treg_prop <- prop.table(table(cd4_merge$orig.ident, cd4_merge$celltype), 1)[,'Treg']
# Extract IL1B expression from Myeloid object (from 01 script)
il1b_exp <- AverageExpression(myeloid_obj, features="IL1B", group.by="orig.ident")$soupx
# Correlate
cor_data <- data.frame(Treg=treg_prop, IL1B=il1b_exp[1, names(treg_prop)])
ggplot(cor_data, aes(x=IL1B, y=Treg)) + geom_point() + geom_smooth(method='lm') + stat_cor()


# ==============================================================================
# 3. NK & DNT Cells
# ==============================================================================
nk_merge <- LoadH5Seurat("data/processed/nk_merge_final.h5seurat")
# Annotation: CD56bright (rNK), CD56dim (mNK), Adaptive (aNK)
nk_markers <- list(rNK=c('NCAM1','XCL1'), mNK=c('FCGR3A','CX3CR1'), aNK=c('KLRC2','CD3E'))
# Check Dysfunction/Checkpoints
checkpoints <- c('TIGIT','LAG3','LILRB1')
compare_gene_pos_ratio(nk_merge, gene=checkpoints, group="sample_type", comparisons=list(c("Healthy","Mets")))