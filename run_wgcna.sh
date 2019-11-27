#!/usr/bin/bash 

timestamp="191127a"

time Rscript ./wgcna.R  --raw.data_path /projects/jonatan/pub-perslab/18-mousebrain/RObjects/mousebrain_L5_subset_raw.csv.gz  --meta.data_path /projects/jonatan/pub-perslab/18-mousebrain/RObjects/mousebrain_L5_subset_metadata.csv   --project_dir /projects/jonatan/pub-perslab/18-mousebrain/18-mousebrain_7/  --scratch_dir /scratch/rkm916/tmp-wgcna/ --data_prefix ClusterName_prior --run_prefix ${timestamp}  --metadata_subset_col ClusterName  --data_organism mmusculus --min.cells 20  --minCellClusterSize 50 --regress_out 'c(NULL)' --genes_use PCA --n_genes_use 5000 --pca_genes var.genes  --corFnc cor  --networkType 'c("signed hybrid")'  --hclustMethod average  --minClusterSize "c(15)"  --deepSplit "c(2)"  --pamStage "c(TRUE)"  --moduleMergeCutHeight "c(0.15)"  --jackstrawnReplicate 0  --TOMnReplicate 0 --autosave T --fuzzyModMembership kME --kM_reassign F --kM_signif_filter T --PPI_filter F --RAM_Gb_max 400 &> run_wgcna_log${timestamp}.txt

#time Rscript ./wgcna.R  --seurat_path /projects/jonatan/pub-perslab/18-mousebrain/RObjects/L5_Neurons_sub_20190212_seurat2.RDS.gz  --project_dir /projects/jonatan/pub-perslab/18-mousebrain/18-mousebrain_7/  --data_prefix Neurons_sub_ClusterName_7.2  --run_prefix repr_3  --metadata_subset_col ClusterName  --data_organism mmusculus --min.cells 20  --minCellClusterSize 50 --regress_out 'c(NULL)' --genes_use PCA  --pca_genes var.genes  --corFnc cor  --networkType 'c("signed hybrid")'  --hclustMethod average  --minClusterSize "c(15)"  --deepSplit "c(2)"  --pamStage "c(TRUE)"  --moduleMergeCutHeight "c(0.15)"  --jackstrawnReplicate 0  --TOMnReplicate 0 --autosave T --fuzzyModMembership kME --kM_reassign F --kM_signif_filter T --PPI_filter F &> run_wgcna_log.txt
