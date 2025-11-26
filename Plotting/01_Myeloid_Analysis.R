source("00_Setup_and_Functions.R")

# ==============================================================================
# Part 1: Neutrophils (Neu) & TANs
# ==============================================================================
# 1.1 Load & Preprocessing
neu_obj <- LoadH5Seurat("data/raw/neu_raw.h5seurat")
# Filter: remove BS1 sample and high MT cells
neu_clean <- subset(neu_obj, sample_type != 'CRPC_BS1' & percent.mt_soupx <= 30) %>% 
             quick_fastmnn(k = 20, reduction.dims = 20, findcluster.res = 0.6)

# 1.2 Annotation
neu_markers <- list(
  NeuP = c('PRTN3','CTSG','KIT','CD38'),
  preNeu = c('FUT4','MS4A3','DEFA3'),
  immNeu = c('CAMP','LTF','LCN2','OLR1'),
  mNeu = c('S100A4','MMP9','FPR1'),
  TAN = c('CCR3','CXCL8','LGALS9','ICAM1','IL1B','CD274')
)
dotp(neu_clean, features = unlist(neu_markers), group.by = "seurat_clusters", save.name = "figures/neu_markers_dotplot")
# Rename clusters (example mapping based on your code)
neu_clean <- RenameIdents_my(neu_clean, c("0"="mNeu", "1"="TAN", "2"="preNeu", ...)) # Replace with actual ID mapping
neu_clean$celltype <- Idents(neu_clean)

# 1.3 Validation with Public Data (Integration)
# Integrating with GSE201333, GSE142786, and Pan-cancer datasets
public_neu <- LoadH5Seurat("data/public/neu_merge_frac_donor.h5seurat")
neu_integrated <- merge_intersect_genes(neu_clean, public_neu) %>% quick_fastmnn()
FeaturePlot(neu_integrated, features = c('CCR3','ICAM1','IL1B'), split.by = 'dataset')

# 1.4 Functional Analysis: AUCell (IFN & IGF1 Signatures)
# Load gene sets from DEGs
ifn_genes <- read.csv("data/signatures/GSE218535_IFNg_deg.csv") %>% filter(padj<0.05) %>% pull(symbol)
igf1_genes <- readRDS("data/signatures/mouse_neu_IGF1_deg.rds") # Converted to human
# Run AUCell
cells_rankings <- AUCell_buildRankings(GetAssayData(neu_clean, slot="data"))
cells_auc <- AUCell_calcAUC(list(IFN=ifn_genes, IGF1=igf1_genes), cells_rankings)
neu_clean$IGF1_Score <- getAUC(cells_auc)["IGF1",]

# 1.5 Wet-lab Validation (IGF1 & TAN)
# 1.5.1 IGF1 Bioavailability (Blood/IHC)
igf1_blood <- readxl::read_xlsx("data/wet_lab/IGF1_clinical.xlsx") %>% melt()
p1 <- barplot_wet(igf1_blood, x="variable", y="value", fill="variable", 
                  comparisons=list(c('Control','Bone Mets')), colors=sample_type_colors)
# 1.5.2 Mouse Model (BLI/Survival)
bli_data <- readxl::read_xlsx("data/wet_lab/IGF1_induced_BLI.xlsx") 
# PlotBLI logic (geom_line + geom_point)
ggplot(bli_data, aes(x=Days, y=log(Mean+1), color=Group)) + geom_line() + theme_classic()


# ==============================================================================
# Part 2: Monocytes & Macrophages
# ==============================================================================
# 2.1 Load & Merge Cohorts (C1, C2, Healthy, Benign)
mono_c1 <- LoadH5Seurat("data/processed/mono_subset_c1.h5seurat")
mono_c2 <- LoadH5Seurat("data/processed/mono_subset_c2_bm.h5seurat")
kf_benign <- LoadH5Seurat("data/public/kf_mono_benign.h5seurat")
healthy <- LoadH5Seurat("data/public/healthy_mono.h5seurat")

mono_merge <- merge_intersect_genes(mono_c1, mono_c2) %>% 
              merge_intersect_genes(kf_benign) %>% 
              merge_intersect_genes(healthy) %>%
              quick_fastmnn(k=20, reduction.dims=30, findcluster.res=0.6)

# 2.2 Annotation
# Monocytes: cMono1/2/3, ncMono
# Macrophages: TAM1 (Ag-presenting), TAM2 (Inflammatory), TAM3 (Lipid/SPP1), RTM (Resident)
mono_markers <- list(
  cMono = c('S100A9','CD14'), ncMono = c('FCGR3A','CX3CR1'),
  TAM1_HLA = c('HLA-DRA','CD74'), TAM2_Inflam = c('IL1B','CXCL8'), TAM3_SPP1 = c('SPP1','FABP4')
)
mono_merge <- RenameIdents_my(mono_merge, c("cMono1", "cMono2", "ncMono", "TAM1", "TAM2", "TAM3"))

# 2.3 Trajectory Analysis (scVelo)
# Remove dissociation genes before velocity
core_genes <- read.csv("data/resources/coregene_df.csv")[,7]
clean_feats <- setdiff(rownames(mono_merge), core_genes)
mono_velo <- subset(mono_merge, features = clean_feats)
# Run scVelo wrapper (Python integration)
# adata <- scvelo_pip(mono_velo, "mono_velo", "data/raw/loom/merged.loom", ".")

# 2.4 M1/M2 Scoring
m1_m2_sigs <- list(
  M1 = c("NOS2","IL1B","TNF","CD86"),
  M2 = c("ARG1","CD163","MRC1","IL10")
)
mono_merge <- AddModuleScore_muti2(mono_merge, m1_m2_sigs)
FeatureScatter(mono_merge, feature1 = "M1", feature2 = "M2", group.by = "celltype")

# 2.5 Cell-Cell Communication (CPDB - Macrophage to Stromal)
# Focus on interactions: Macro (Ligand) -> Epi/Endo (Receptor)
# Load CPDB results
cpdb_res <- load_cpdb_res("data/cpdb_output/macro_stromal/")
# Plot specific interactions
plot_cpdb_cellsubset_point(cpdb_res, celltype1=c("TAM2","TAM3"), celltype2=c("CapEC","Osteo1"),
                           interaction_names = c("VEGFA_KDR", "SPP1_ITGAV"))

# 2.6 LGALS9 Expression Check
# Compare LGALS9 in PreMets vs Benign
compare_gene_pos_ratio(mono_merge, gene="LGALS9", group="sample_type", 
                       comparisons=list(c("Benign","preMets")))