setwd('/home01/tangyk/scrna/Prostate_Cancer/17sample_final_matrix/2_subcells/1_tumor/final_figures/compass/compass_results_in_final_figures')
#load data
microp <- read.table('/home01/tangyk/scrna/Prostate_Cancer/17sample_final_matrix/2_subcells/1_tumor/compass/231019_microcluster_30//results//micropools.tsv')
microp$microcluster <- paste0('cluster_',microp$microcluster)
microp_celltype <- microp
microp_celltype[colnames(cancer_scenic),'celltype'] <- cancer_scenic$celltype2
microp_celltype <- na.omit(microp_celltype) %>% as.data.frame

#micropool construction
for(i in unique(microp_celltype$microcluster)){
     check_celltype <- microp_celltype %>% {.[.$microcluster==i,2]} %>% table %>% as.matrix
     max_ratio = max(check_celltype[,1])/sum(check_celltype[,1])
     microp_celltype[microp_celltype$microcluster==i,'microcluster_celltype_ratio'] <- max_ratio
}

#filtered the cells whose purity is less than 90%
microp_celltype[ !duplicated(microp_celltype$microcluster),] %>% .[.$microcluster_celltype_ratio>0.9,]
microcluster_keep <- microp_celltype[ !duplicated(microp_celltype$microcluster),] %>% subset(microcluster_celltype_ratio>0.9)
rownames(microcluster_keep) <- microcluster_keep$microcluster

#calculate scpre and create object
compass_bm <- read.table('/home01/tangyk/scrna/Prostate_Cancer/17sample_final_matrix/2_subcells/1_tumor/compass/231019_microcluster_30//results//reactions.tsv')
get_reaction_consistencies <- function(compass_reaction_penalties, min_range=1e-3){
  #https://yoseflab.github.io/Compass/notebooks/Demo.html
  # Apply transformation to the entries in the dataframe
  df <- -log(compass_reaction_penalties + 1)
  # Filter rows where max value - min value is gte to min_range
  df <- df[apply(df, 1, function(x) max(x) - min(x)) >= min_range, ]
  # Subtract minimum value from the whole dataframe
  df <- df - min(apply(df, 2, min)) 
  return(df)
}
compass_score <- get_reaction_consistencies(compass_bm[,rownames(microcluster_keep)])
epi_compass <- CreateSeuratObject(counts = compass_score,data = compass_score,project = 'epi_compass',meta.data = microcluster_keep)
SeuratDisk::SaveH5Seurat(epi_compass,'epi_compass.h5seurat')
#Find markers, merge metadata
Idents(epi_compass) <- 'celltype'
marker_epi_compass_celltype2 <- FindAllMarkers(epi_compass,logfc.threshold = 0,min.pct = 0,only.pos = F)
