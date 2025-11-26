#
setwd('/home01/tangyk/scrna/Prostate_Cancer/other_data/GSE126030_Tcells_activate')
matrix.files = list.files('1_matrix/',full.names = T)
matrix.list <- lapply(matrix.files,FUN = function(x) read.table(x,row.names = 1,header = T))

#uniforem genes
rownames_matrix = rownames(matrix.list[[1]])
#check all matrix contain the same genes
sapply(matrix.list,function(x){
        all(rownames_matrix%in%rownames(x))
    })
matrix.list <- lapply(matrix.list,FUN = function(x){x <- x[rownames_matrix,];return(x)})
matrix_all <- rlist::list.cbind(matrix.list2)#merge
saveRDS(matrix_all,'matrix_all.rds')
matrix_all <- as.data.frame(matrix_all[,-1],row.names=matrix_all[,1])

#duplicated barcodes
matrix_all %>% colnames %>% duplicated %>% table
#remove duplicated
#duplicated_barcodes <- colnames(matrix_all) %>% {.[duplicated(.)]}
#matrix_all_clean <- matrix_all %>% {.[,-which(colnames(.)%in%duplicated_barcodes)]}

meta.data <- read.csv('0_rawdata/41467_2019_12464_MOESM9_ESM_metadata_fig6.csv')
matrix_all_filtered <- matrix_all[,meta.data$barcode[meta.data$barcode%in%colnames(matrix_all)]]#intersect with meta.data provided by authors
matrix_all_filtered <- matrix_all_filtered %>% {.[-which(rownames(.) %in% rownames(.)[which(duplicated(rownames(.)))]),]}#remove duplicated genes

#rownames(meta.data) <- meta.data$barcode
meta.data_filtered <- meta.data %>% {.[-which(.$barcode%in%.$barcode[duplicated(.$barcode)]),]}#in meta.data, there existed duplicated barcodes, remove them
rownames(meta.data_filtered) <- meta.data_filtered$barcode

sample <- CreateSeuratObject(counts = matrix_all_filtered[,rownames(meta.data_filtered)])


sample@meta.data <- meta.data_filtered[colnames(sample),]
SeuratDisk::SaveH5Seurat(sample,filename = 'uniuqe_barcode_object2.h5seurat')

#
#tcell_activate <- quick_fastmnn(tcell_activate,k = 20,batch_remove_based_on = 'donor',reduction.dims = 30,findcluster.res = 0.4)
tcell_activate <- RenameAssays(tcell_activate,RNA = 'soupx')
tcell_activate <- SeuratDisk::LoadH5Seurat('~/scrna/Prostate_Cancer//other_data//GSE126030_Tcells_activate//uniuqe_barcode_object2.h5seurat')

tcell_activate_cd8 <- subset(tcell_activate,cd4cd8_status=='CD8')

 Anchors_tcell_activate <- FindTransferAnchors(
        reference = cd8_3,
        query = tcell_activate_cd8,
        normalization.method = "LogNormalize",
        k.anchor = 5,
        features = intersect(rownames(cd8_3),rownames(tcell_activate_cd8)),
        reference.reduction = "mnn",
        dims = 1:30)

query_tcell_activated <- MapQuery(
    anchorset = anchors,
    query = cd8_3,
    reference = pbmc,
    refdata = list(l1 = "celltype2"),
    reference.reduction = "mnn",
    reduction.model = "wnn.umap") 
###

pbmc <- readRDS('/home01/tangyk/scrna/Prostate_Cancer/other_data/Tcells/tmp1/PBMC_vaccine_CITE.rds')

nchors <- FindTransferAnchors(
    reference = pbmc,
    query = cd8_3,
    normalization.method = "SCT",
    k.anchor = 5,
    reference.reduction = "spca",
    dims = 1:50)


query_scRNA <- MapQuery(
    anchorset = anchors,
    query = cd8_3,
    reference = pbmc,
    refdata = list(
      l1 = "celltypel1",
      l2 = "celltypel2",
      l3 = "celltypel3"),
    reference.reduction = "spca",
    reduction.model = "wnn.umap") 
  
