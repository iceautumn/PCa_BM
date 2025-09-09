#batch correction and clustering
    library(Seurat)
    library(SeuratWrappers)
    library(magrittr)
quick_fastmnn <- function(seurat_obj,k=5, normalizeddata = T,batch_remove_based_on='orig.ident',reduction.dims=30,features = 2000,findcluster.res=c(0.4,0.6,0.8),...){

    seurat_obj <- DietSeurat(seurat_obj)
    if(normalizeddata){seurat_obj <- NormalizeData(seurat_obj)}
    seurat_obj <- seurat_obj %>% FindVariableFeatures %>% {RunFastMNN(object.list = SplitObject(., split.by = batch_remove_based_on),k=k, features = features)} %>% 
                                RunUMAP(reduction = "mnn", dims = 1:reduction.dims,...) %>% FindNeighbors(reduction = "mnn", dims = 1:reduction.dims) 

                                for(i in findcluster.res){seurat_obj <- FindClusters(seurat_obj,resolution=i)}#find cluster
    return( seurat_obj)
}