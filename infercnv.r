#inferCNV========================================================================================================
library(infercnv)
library(Seurat)
setwd('/home01/tangyk/scrna/Prostate_Cancer/new_PCa_data_Ren241013/4_rmCSL/1_subcells/2_Epi/induvidual_analysis/FYQ_WGS_WWZ_XJG/infercnv/241225')


epi_mets8_test <- SeuratDisk::LoadH5Seurat('/home01/tangyk/scrna/Prostate_Cancer/new_PCa_data_Ren241013/4_rmCSL/1_subcells/2_Epi/induvidual_analysis/FYQ_WGS_WWZ_XJG/epi_mets8.FYQ_WGS_WWZ_XJG.h5seurat')
obj <- epi_mets8_test
#not separate


get_infercnv_input_data(obj,split_group=NULL,anno_group='RNA_snn_res.0.6',save_dir='./data',assay='RNA')

output_dir = '/home01/tangyk/scrna/Prostate_Cancer/new_PCa_data_Ren241013/4_rmCSL/1_subcells/2_Epi/induvidual_analysis/FYQ_WGS_WWZ_XJG/infercnv/241225/output'

runinferCNV(anno_group  = '/home01/tangyk/scrna/Prostate_Cancer/new_PCa_data_Ren241013/4_rmCSL/1_subcells/2_Epi/induvidual_analysis/FYQ_WGS_WWZ_XJG/infercnv/241225/data/infercnv_anno_cluster.txt',
            counts = '/home01/tangyk/scrna/Prostate_Cancer/new_PCa_data_Ren241013/4_rmCSL/1_subcells/2_Epi/induvidual_analysis/FYQ_WGS_WWZ_XJG/infercnv/241225/data/infercnv_count.txt',
               gene_location = '/home01/tangyk/scrna/data/Genome_reference/hg38_gencode_v27.txt',
               ref_group_names = c("4","9","0"),output_dir = output_dir,HMM = F)