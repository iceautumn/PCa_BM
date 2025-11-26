# ==============================================================================
# 00_Setup_and_Functions.R
# Purpose: Load libraries, set paths, and define custom wrapper functions.
# ==============================================================================

# 1. Libraries
library(Seurat); library(SeuratDisk); library(SeuratWrappers)
library(dplyr); library(ggplot2); library(patchwork); library(reshape2)
library(ComplexHeatmap); library(paletteer); library(ggpubr); library(ggthemes)
library(nichenetr); library(monocle); library(GSVA)
library(reticulate) # For scVelo interaction

# 2. Global Settings
setwd("./Prostate_Cancer_Project") # Set your root directory
options(stringsAsFactors = FALSE)
dir.create("figures", showWarnings = F)
dir.create("processed_data", showWarnings = F)

# 3. Color Palettes (Consolidated from your scripts)
my_pal <- paletteer_d("ggsci::nrc_npg")
sample_type_colors <- c(Healthy='#40A2E3', Benign='#B7E5B4', preMets='#CCD7CF', Mets='#EA8579')
modulescore_colors <- colorRampPalette(c('#2B307A','#77C2F3','#F7EEF6','#D8A0C7','#A96FB0'))(50)

# 4. Custom Functions Definitions

#' Wrapper for FastMNN Integration pipeline
#' Splits object by batch, normalizes, finds HVGs, runs FastMNN, UMAP, and Clustering.
quick_fastmnn <- function(seurat_obj, k = 20, reduction.dims = 30, findcluster.res = 0.6, 
                          batch_col = "orig.ident", assay = "soupx") {
  obj_list <- SplitObject(seurat_obj, split.by = batch_col)
  obj_list <- lapply(obj_list, function(x) {
    x <- NormalizeData(x, assay = assay)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, assay = assay)
    return(x)
  })
  seurat_obj <- RunFastMNN(object.list = obj_list, k = k, d = reduction.dims)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "mnn", dims = 1:reduction.dims)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "mnn", dims = 1:reduction.dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = findcluster.res)
  return(seurat_obj)
}

#' Custom DotPlot Wrapper
dotp <- function(object, features, group.by, angle.x = 90, save.name = NULL) {
  p <- DotPlot(object, features = features, group.by = group.by) + 
       coord_flip() + 
       theme(axis.text.x = element_text(angle = angle.x, hjust = 1),
             panel.background = element_rect(fill = NA, colour = "black"),
             panel.grid.major = element_blank()) +
       scale_color_distiller(palette = "RdBu")
  if(!is.null(save.name)) ggsave(paste0(save.name, ".pdf"), p, width = 5, height = length(features)*0.3)
  return(p)
}

#' Rename Idents Wrapper
RenameIdents_my <- function(object, new_names_vector) {
  # new_names_vector should be named: c("OldName" = "NewName") or matching levels order
  current.ids <- levels(object)
  if(is.null(names(new_names_vector))) names(new_names_vector) <- current.ids
  object <- RenameIdents(object, new_names_vector)
  return(object)
}

#' Barplot for Wet-lab Data (qPCR, Flow, BLI)
barplot_wet <- function(data, x, y, fill, comparisons, method='wilcox.test', colors) {
  ggplot(data, aes_string(x=x, y=y)) +
    stat_summary(geom = 'bar', fun = 'mean', width = 0.6, aes_string(fill=fill)) + 
    geom_jitter(width = 0.2, size = 2, color='black', shape=21) +
    stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="errorbar", width=0.2) +
    scale_fill_manual(values = colors) +
    stat_compare_means(comparisons = comparisons, method = method) +
    theme_classic() + NoLegend()
}

#' Helper to merge and intersect genes (preventing errors during merge)
merge_intersect_genes <- function(obj1, obj2) {
  genes <- intersect(rownames(obj1), rownames(obj2))
  obj1 <- subset(obj1, features = genes)
  obj2 <- subset(obj2, features = genes)
  return(merge(obj1, obj2))
}