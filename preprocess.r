# ---
# Single-Cell RNA Sequencing Analysis Pipeline
#
# This script processes 17 prostate cancer samples, performing ambient RNA correction,
# quality control, doublet detection, clustering, and cell type annotation.
#
# Author: [Your Name]
# Date: [Date]
# ---

# 1. SETUP AND CONFIGURATION
# ========================

# Load required libraries
library(magrittr)
library(Seurat)
library(dplyr)
library(DoubletFinder)

# Set the working directory
setwd('/home01/tangyk/scrna/Prostate_Cancer/17sample_final_matrix/0_rawdata')


# 2. DATA LOADING AND INITIAL OBJECT CREATION
# ==========================================

# --- 2.1. Load Raw and Filtered Data ---

# Get the paths to the raw feature-barcode matrices
raw_list <- list.dirs(path = './0_rawdata//raw_feature_bc_matrix', full.names = TRUE)[-1]
names(raw_list) <- sapply(strsplit(raw_list, split = '/'), '[', 5)

# Get the paths to the filtered feature-barcode matrices
filtered_list <- list.dirs('0_rawdata//filtered_feature_bc_matrix', full.names = TRUE, recursive = TRUE)[-1]
names(filtered_list) <- sapply(strsplit(filtered_list, split = '/'), '[', 3)

# Read the 10X data
raw.matrix <- Read10X(raw_list)
filtered.matrix <- Read10X(filtered_list)

# --- 2.2. Create Initial Seurat Objects ---

# Create Seurat objects from the raw and filtered matrices
raw <- CreateSeuratObject(counts = raw.matrix)
fil <- CreateSeuratObject(counts = filtered.matrix)

# Save the initial Seurat objects
saveRDS(raw, 'raw.cellranger.17samples.rds')
saveRDS(fil, 'filtered.cellranger.17samples.rds')


# 3. AMBIENT RNA CORRECTION (SOUPX)
# ================================

#' @title Run SoupX in Batches
#' @description This function applies the SoupX algorithm to correct for ambient RNA
#'              contamination in a batch-wise manner.
#'
#' @param raw.obj A Seurat object containing the raw count matrix.
#' @param fil.obj A Seurat object containing the filtered count matrix.
#' @param batch The metadata column specifying the batch for each cell.
#' @param Rho The expected contamination fraction.
#' @param cluster The metadata column containing cell cluster information.
#' @param nonExpressedGeneList A list of genes known to be unexpressed in certain cell types.
#' @param assay_name The name of the new assay to store the corrected counts.
#' @param nonExpressedGeneList_based_on_celltype A boolean indicating if the nonExpressedGeneList is cell-type specific.
#' @param celltype_ratio_threshold The proportion threshold for a cell type to be considered for contamination calculation.
#'
#' @return A Seurat object with a new assay containing the SoupX-corrected counts.

run_soupx_batch2 <- function(raw.obj, fil.obj, batch, Rho = 0.2, cluster, nonExpressedGeneList = NULL, assay_name = 'SoupX', nonExpressedGeneList_based_on_celltype = FALSE, celltype_ratio_threshold = 0.001) {

    # The nonExpressedGeneList is based on the cell types present in each batch and
    # is used to estimate the contamination fraction.
    # For example, if the nonExpressedGeneList includes genes specific to certain
    # cell types, and a batch does not contain those cell types, the corresponding
    # genes will not be used to calculate the contamination for that batch.

    # Split Seurat objects by batch
    raw.list <- SplitObject(raw.obj, split.by = batch)
    fil.list <- SplitObject(fil.obj, split.by = batch)

    # Run SoupX for each batch
    if (nonExpressedGeneList_based_on_celltype) {
        print('The names in nonExpressedGeneList should correspond to the clusters in fil.obj')
        fil.soupx.list <- list()
        for (i in names(fil.list)) {
            check_type <- names(nonExpressedGeneList)[names(nonExpressedGeneList) %in% fil.list[[i]]@meta.data[, cluster]]
            check_type_propotion <- prop.table(table(fil.list[[i]]@meta.data[, cluster])) %>%
                as.data.frame() %>%
                .[.$Var1 %in% check_type & .$Freq > celltype_ratio_threshold, 'Var1']

            nonExpressedGeneList_i <- nonExpressedGeneList[which(names(nonExpressedGeneList) %in% check_type_propotion)]
            fil.soupx.list[[i]] = run_SoupX(raw.obj = raw.list[[i]], fil.obj = fil.list[[i]], cluster = cluster, Rho = Rho, nonExpressedGeneList = nonExpressedGeneList_i)
        }
    } else {
        fil.soupx.list <- list()
        for (i in names(raw.list)) {
            fil.soupx.list[[i]] = run_SoupX(raw.obj = raw.list[[i]], fil.obj = fil.list[[i]], cluster = cluster, Rho = Rho, nonExpressedGeneList = nonExpressedGeneList)
        }
    }

    # Merge the SoupX results
    fil.soupx <- fil.soupx.list[[1]]
    for (i in 2:length(fil.soupx.list)) {
        fil.soupx <- merge(fil.soupx, fil.soupx.list[[i]])
    }

    # Create a new assay with the corrected counts
    fil.obj[[assay_name]] <- CreateAssayObject(counts = fil.soupx@assays$SoupX@counts)
    return(fil.obj)
}



celltype <- c('Neu', 'Neu', 'T', 'Osteo', 'Endothelial', 'Neu', 'Neu', 'Neu', 'Neu', 'Neu', 'Tumor', 'Ery', 'Tumor', 'GMP', 'Monolytic', 'Tumor', 'Monolytic', 'Tumor', 'Tumor', 'Tumor',
              'Osteo', 'Mast', 'Tumor', 'Eos', 'B', 'Monolytic', 'Osteo', 'Endothelial', 'Plasma', 'Tumor', 'Monolytic', 'Osteo', 'Osteo')
names(celltype) <- 0:33
sample_filtered <- RenameIdents(sample_filtered, celltype)

sample_filtered$celltype <- Idents(sample_filtered)
saveRDS(sample_filtered, '17samples_nfeature_filtered_neu300_other800_percentMT50_fastmnn_k10_dims50.rds')

# --- 5.3. Apply SoupX Correction ---

# Define a list of non-expressed genes for different cell types
nonExpressedGeneList = list(
    B = c('IGKC', "IGHA1", "IGHA2", "IGHG1"),
    Ery = c('HBB', 'HBD', 'HBA2'),
    Neu = c('MPO', 'ELANE', 'LTF'),
    Osteo = c('COL1A1', 'COL1A2', 'COL3A1')
)

# Run SoupX on the filtered data
sample_soupx <- run_soupx_batch2(
    raw.obj = raw,
    fil.obj = sample_filtered,
    batch = 'orig.ident',
    cluster = 'celltype',
    nonExpressedGeneList = nonExpressedGeneList,
    assay_name = 'SoupX',
    nonExpressedGeneList_based_on_celltype = TRUE
)


# 6. DOUBLET DETECTION (DOUBLETFINDER)
# ====================================

# --- 6.1. DoubletFinder Function ---

#' @title Run DoubletFinder in Batches
#' @description This function runs the DoubletFinder algorithm on a Seurat object
#'              to identify potential doublets.
#'
#' @param seurat_obj A Seurat object.
#' @param batch The metadata column specifying the batch for each cell.
#' @param dims The number of principal components to use.
#' @param GT A boolean indicating if ground truth doublet information is available.
#' @param cluster The metadata column containing cell cluster information.
#' @param DoubletRate The estimated doublet rate.
#' @param pN The proportion of artificial doublets to generate.
#' @param sct A boolean indicating if SCTransform has been used for normalization.
#' @param reuse.pANN A boolean indicating whether to reuse existing pANN values.
#'
#' @return A list of data frames containing the DoubletFinder results for each batch.

run_doubletfinder <- function(seurat_obj, batch, dims, GT = FALSE, cluster, DoubletRate = 0.02, pN = 0.25, sct = FALSE, reuse.pANN = FALSE, ...) {
    library(DoubletFinder)
    library(dplyr)
    library(Seurat)

    obj.list <- SplitObject(seurat_obj, split.by = batch)
    result.list <- list()

    for (i in names(obj.list)) {
        tmp <- obj.list[[i]]
        sweep.res.list <- paramSweep_v3(tmp, PCs = 1:dims, sct = sct)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = GT)
        bcmvn <- find.pK(sweep.stats)
        pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
        homotypic.prop <- modelHomotypic(tmp@meta.data[, cluster])

        nExp_poi <- round(DoubletRate * ncol(tmp))
        nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
        tmp <- doubletFinder_v3(
            tmp, PCs = 1:dims, pN = pN, pK = pK_bcmvn,
            nExp = nExp_poi.adj, reuse.pANN = reuse.pANN, sct = sct
        )
        result.list[[i]] <- tmp@meta.data
    }

    return(result.list)
}


# 7. FINAL FILTERING AND ANALYSIS
# ===============================

# --- 7.1. Apply DoubletFinder and Filter ---

# Final workflow:
# 1. Run DoubletFinder on the initial object -> [doublet]
# 2. Filter cells (nFeature >= 800) and run SoupX -> [sample_soupx] (already done)
# 3. Remove the identified doublets from the SoupX-corrected object

# Read the previously saved Seurat objects
sample_soupx <- readRDS('/home01/tangyk/scrna/Prostate_Cancer/17sample_final_matrix/1_preprocess/3_final_pipline/17samples_nfeature_filtered_neu300_other800_percentMT50_fastmnn_k10_dims50.rds')
sample <- readRDS('17samples_nfeature300_percentMT50_fastmnn_k10_dims50.rds')
sample_filtered <- readRDS('17samples_nfeature_filtered_neu300_other800_percentMT50_fastmnn_k10_dims50.rds')

# Run DoubletFinder
sample_df_results <- run_doubletfinder(seurat_obj = sample, batch = 'orig.ident', dims = 50, cluster = 'RNA_snn_res.0.6', DoubletRate = 0.06)

# Extract doublet cell names
doublet <- sample_df_results[[1]] %>% {rownames(.)[which(.[, 10] != 'Singlet')]}

# Add doublet information to the Seurat objects
sample$doublet <- sapply(colnames(sample), FUN = function(x) {ifelse(x %in% doublet, '1', '0')})
sample_soupx$doublet <- sapply(colnames(sample_soupx), FUN = function(x) {ifelse(x %in% doublet_filtered, '1', '0')})

# Subset to remove doublets
sample_soupx_filtered <- subset(sample_soupx, doublet == '0')

# --- 7.2. Alternative Workflow: Doublet Detection before Filtering ---

# Check an alternative workflow: nFeature_RNA > 300 -> DoubletFinder -> nFeature_RNA > 800 -> SoupX

# Filter for singlets and re-cluster
sample_singlet <- subset(sample, doublet == '0') %>%
    quick_fastmnn(k = 10, reduction.dims = 50, findcluster.res = 0.6)

# Identify neutrophils in the singlet object
# Note: Cluster 28 here is different from the previous neutrophil definition
# as it includes a part of cluster 41 from the original 'sample' object's res 1.2,
# which was difficult to separate.
sample_singlet$check_neutrophil[sample_singlet$RNA_snn_res.1 == '28'] <- 'neu'

# Subset and filter singlets
sample_singlet_neu <- subset(sample_singlet, check_neutrophil == 'neu')
sample_singlet_other <- subset(sample_singlet, check_neutrophil != 'neu' & nFeature_RNA >= 800)

# Merge and re-cluster
sample_singlet_filtered <- merge(sample_singlet_neu, sample_singlet_other) %>%
    quick_fastmnn(k = 10, reduction.dims = 50, findcluster.res = 0.6)

# Re-annotate cell types
new.id <- c("Neu", "Neu", "T", "Osteo", "Endothelial", "Neu", "Neu", "Neu", "Neu", "Neu", "Tumor", "Monolytic", "Neu", "GMP", "Tumor", "Monolytic", "Osteo", "Ery", "Tumor", "Tumor", "Tumor", "Tumor", "Mast", "Tumor", "Eos", "Ery", "Osteo", "Endothelial", "Tumor", "B", "B", "Osteo")
names(new.id) <- 0:31
sample_singlet_filtered <- RenameIdents(sample_singlet_filtered, new.id)
sample_singlet_filtered$celltype <- Idents(sample_singlet_filtered) %>% as.character()

# Run SoupX on this filtered singlet object
sample_singlet_filtered_soupx <- run_soupx_batch2(
    raw.obj = raw,
    fil.obj = sample_singlet_filtered,
    batch = 'orig.ident',
    cluster = 'celltype',
    nonExpressedGeneList = nonExpressedGeneList,
    assay_name = 'SoupX',
    nonExpressedGeneList_based_on_celltype = TRUE
)
