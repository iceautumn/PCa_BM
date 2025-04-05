#
setwd('/home01/tangyk/scrna/Prostate_Cancer/other_data/GSE120221_Heathy_bone_marrow/0_rawdata/matrix')
list.files('./matrix/',pattern = '_barcodes_') %>% gsub(pattern = 'GSM[0-9][0-9][0-9][0-9][0-9][0-9][0-9]_|barcodes_|.tsv.gz','',x = .) -> samples_name
for(i in samples_name){
    dir.create(i)
    pattern_i = paste0('_',i,'.')
    files_tmp = list.files(path = '.',pattern = pattern_i)
    files_new_name = gsub(pattern = 'GSM[0-9][0-9][0-9][0-9][0-9][0-9][0-9]|_','',x = files_tmp) %>% gsub(pattern = i,'',x = .) %>% gsub(pattern = 'genes', 'features', x = .) %>% {paste0(i,'/',x = .)}
    for(j in 1:length(files_tmp)){
        file.rename(from = files_tmp[i],to = files_new_name[i])
    }
}

raw.matrix.dir = samples_name
names(raw.matrix.dir) = samples_name
raw.matrix = Read10X(raw.matrix.dir)

raw.obj[['percent.mt']] <- PercentageFeatureSet(raw.obj,pattern = '^MT-')
p1 <- VlnPlot(raw.obj,features = c('nFeature_RNA','nCount_RNA','percent.mt'),stack = T,flip = T)+NoLegend()
#r$> SeuratDisk::SaveH5Seurat(raw.obj, filename = 'raw.integrated.h5Seurat')
filtered = subset(raw.obj,nFeature_RNA>=800&percent.mt<=10)

#Doublet finder:
r$> sample_df_results_0.06 <- run_doubletfinder(seurat_obj = filtered,batch = 'orig.ident',assay = 'RNA',dims = 20,cluster_resolution = 0.4,DoubletRate = 0.06)
r$> saveRDS(sample_df_results_0.06,'sample_df_dims20_louvain0.4_results_0.06.rds')

r$> doublet_0.06_list = lapply(sample_df_results_0.06,function(x){colnames(x)[7:8] = c('pANN','df');return(x)}) 
doublet_0.06_list = lapply(doublet_0.06_list, function(x){x[,'barcode'] = rownames(x);return(x)})
r$> doublet_0.06_res = rlist::list.rbind(doublet_0.06_list)
r$> rownames(doublet_0.06_res)=doublet_0.06_res$barcode


#soupx
r$> filered_singlet <- subset(filtered,df == 'Singlet') %>% quick_fastmnn(k = 20,findcluster.res = 0.4,reduction.dims = 20)
##anno


new.id <- c("T", "Mono", "T", "T", "Ery", "B", "HSCs", "Ery", "Mono", "Mono", "B", "Ery", "Ery", "pDC", "Mono", "B", "B", "T", "Mono", "Mono")


nonExpressedGeneList = list(B = c('IGKC',"IGHG3", "IGHA1", "IGLC2"),Ery=c('HBB','HBD','HBA1'),Mono = c('CD14','FCGR3A'),Tcell = c('CD3E','CD3D','CD3G'))
sample_singlet_filtered_soupx <- run_soupx_batch2(raw.obj = raw.obj,fil.obj = filtered_singlet,batch = 'orig.ident',cluster = 'celltype',nonExpressedGeneList = nonExpressedGeneList,assay_name = 'soupx',nonExpressedGeneList_based_on_celltype = T)
DefaultAssay(sample_singlet_filtered_soupx) = 'soupx'
#SeuratDisk::SaveH5Seurat(sample_singlet_filtered_soupx, filename = 'sample_singlet_filtered_soupx_notReclustered.h5Seurat')

##
sample_singlet_filtered_soupx <- quick_fastmnn(sample_singlet_filtered_soupx,k = 20,reduction.dims = 20,findcluster.res = 0.4)
#SeuratDisk::SaveH5Seurat(sample_singlet_filtered_soupx_recluster, filename = 'sample_singlet_filtered_soupx_Reclustered_based_on_soupx.h5Seurat')

#erythroid
r$> ery_healthy = subset(healthy,celltype == 'Ery')
r$> ery_healthy <- quick_fastmnn(ery_healthy,k = 20,reduction.dims = 30,findcluster.res = 0.6)







#tcell
r$> tcell <- subset(sample_singlet_filtered_soupx_recluster,celltype == 'Tcell')
r$> SeuratDisk::SaveH5Seurat(tcell,'tcell.h5Seurat')

tcell_health = tcell
tcell_health <- FindClusters(tcell_health,algorithm = 4, resolution  =1 )
tcell_health$celltype_tmp = sapply(tcell_health$soupx_snn_res.1,function(x){
        if(x  %in% c(8,14)){x <- 'NK/T'}
        else if(x %in% c(1,4,5,9,10)){x <- 'CD4'}
        else if(x %in% c(2,3,7,11,12)){x <- 'CD8'}
        else{x <- 'mix'}
    })

r$> tcell_healthy_mix <- subset(tcell_health,celltype_tmp == 'mix') %>% quick_fastmnn(k = 20,reduction.dims = 20,findcluster.res = 0.4)
r$> tcell_healthy_mix <- FindClusters(tcell_healthy_mix,resolution = 1)
tcell_healthy_mix$celltype <- sapply(tcell_healthy_mix$soupx_snn_res.1,function(x){
        if(x %in% c(1,2,6)){x <- 'CD8_T'}
        else{x <- 'CD4_T'}
    })
SeuratDisk::SaveH5Seurat(tcell_healthy_mix, '/home01/tangyk/scrna/Prostate_Cancer/other_data/GSE120221_Heathy_bone_marrow/1_tcell/tcell_healthy_mix.h5seurat')

r$> tcell_health_cd8 <- subset(tcell_health,celltype_tmp == 'CD8') %>% quick_fastmnn(k = 20,reduction.dims = 20,findcluster.res = 0.4)
r$> fp(tcell_health_cd8,obj_name = 'tcell_heath_cd8',features = c('SELL','CCR7','IL7R','LEF1','CD28','CD27','PECAM1','CX3CR1','GZMA','GZMB','IFNG'))
r$> fp(tcell_health_cd8,obj_name = 'tcell_heath_cd8',features = c('CD4','CD8A','CD8B'))

tcell_health_cd8$celltype = sapply(tcell_health_cd8$soupx_snn_res.0.4,function(x) ifelse(x == 0, x <- 'CD4_T',x <- 'CD8_T'))



r$> tcell_health_cd4 <- subset(tcell_health,celltype_tmp == 'CD4') %>% quick_fastmnn(k = 20,reduction.dims = 20,findcluster.res = 0.4)
#celltypist
r$> tcell_health_cd4_matrix <- GetAssayData(tcell_health_cd4,slot = 'counts',assay = 'soupx')
write.table(tcell_health_cd4_matrix,'/home01/tangyk/scrna/Prostate_Cancer/other_data/GSE120221_Heathy_bone_marrow/1_tcell/tcell_health_cd4_counts_soupx.txt',sep = '\t',quote = F,row.names = T)
#bsah
celltypist --indata /home01/tangyk/scrna/Prostate_Cancer/other_data/GSE120221_Heathy_bone_marrow/1_tcell/tcell_health_cd4_counts_soupx.txt \
--outdir /home01/tangyk/scrna/Prostate_Cancer/other_data/GSE120221_Heathy_bone_marrow/1_tcell \
--transpose-input


r$> tcell_health_nkt <- subset(tcell_health,celltype_tmp=='NK/T') %>% quick_fastmnn(k = 20,reduction.dims = 20,findcluster.res = 0.4)
tcell_health_nkt_matrix <- GetAssayData(tcell_health_nkt,slot = 'counts',assay = 'soupx')
write.table(tcell_health_nkt_matrix,'/home01/tangyk/scrna/Prostate_Cancer/other_data/GSE120221_Heathy_bone_marrow/1_tcell/tcell_health_nkt_counts_soupx.txt',sep = '\t',quote = F,row.names = T)

r$> celltypist_nkt_health = read.csv('/home01/tangyk/scrna/Prostate_Cancer/other_data/GSE120221_Heathy_bone_marrow/1_tcell/predicted_labels.csv',row.names = 1)

#monocytic 
r$> healthy <- SeuratDisk::LoadH5Seurat('/home01/tangyk/scrna/Prostate_Cancer/other_data/GSE120221_Heathy_bone_marrow/0_rawdata/sample_singlet_filtered_soupx_Reclustered_based_on_soupx.h5Seurat',assay = 'soupx')
monocytic <- subset(healthy, celltype == 'Mono')

r$> healthy <- SeuratDisk::LoadH5Seurat('/home01/tangyk/scrna/Prostate_Cancer/other_data/GSE120221_Heathy_bone_marrow/0_rawdata/sample_singlet_filtered_soupx_Reclustered_based_on_soupx.h5Seurat')
monocytic <- subset(healthy, celltype == 'Mono')