#scenic
conda activate scenic
python ~/scrna/data/MyFunction.py/pySCENIC.script.py \
--loom /home01/tangyk/scrna/Prostate_Cancer/new_PCa_data_Ren241013/4_rmCSL/1_subcells/2_Epi/induvidual_analysis/FYQ_WGS_WWZ_XJG/scenic/epi_mets8_test.loom \
--out /home01/tangyk/scrna/Prostate_Cancer/new_PCa_data_Ren241013/4_rmCSL/1_subcells/2_Epi/induvidual_analysis/FYQ_WGS_WWZ_XJG/scenic/res \
--tfs_path /home01/tangyk/scrna/data/SCENIC_database/allTFs_hg38.txt \
--db_feather /home01/tangyk/scrna/data/SCENIC_database/motif_database/hg38/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--motif_path /home01/tangyk/scrna/data/SCENIC_database/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
