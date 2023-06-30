#PBS -N Seurat_Analysis_R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B
#PBS -o /mnt/hpcscratch/vfajardo/R24/jobs_scripts/seurat_analysis//R24_Cancer_Myeloid/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B_11-08-2022qc-mye-spc_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B_general_seurat_analysis.out.txt
#PBS -e /mnt/hpcscratch/vfajardo/R24/jobs_scripts/seurat_analysis//R24_Cancer_Myeloid/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B_11-08-2022qc-mye-spc_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B_general_seurat_analysis.err.txt
#PBS -m abe
#PBS -M vfajardo@lji.org
#PBS -q default
#PBS -l nodes=1:ppn=8
#PBS -l mem=130gb
#PBS -l walltime=200:00:00

echo -e "\n\n######### ------------- Job to Run Seurat Analysis ------------ #########\n\n"
echo -e "### --------------------------- Arguments --------------------------- ###"
echo -e "### --------------------- Job general arguments --------------------- ###"
echo -e "Project name: R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B"
echo -e "Job output directory: /mnt/hpcscratch/vfajardo/R24/jobs_scripts/seurat_analysis//R24_Cancer_Myeloid/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B_11-08-2022qc-mye-spc_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL"
echo -e "### ------------------- Seurat analysis arguments ------------------- ###"
echo -e "Seurat script version (absolute path): /home/vfajardo/scripts/seurat_analysis/general_seurat_analysis.2.3.5.R"
echo -e "Project name: R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B"
echo -e "Absolute path to count matrix file: "
echo -e "Feature ID: ensembl"
echo -e "Minimum and maximum thresholds for counts number per cell: 1500 and 13000"
echo -e "Minimum and maximum thresholds for genes number per cell: 350 and 4000"
echo -e "Maximum mitochondrial genes percentage per cell: 15"
echo -e "Normalization method: LogNormalize"
echo -e "Variable features selection method: vst"
echo -e "Number of most variable genes to be considered: 30"
echo -e "Resolution argument for cluster analysis: 'c(0.1, 0.2, 0.3, 0.4, 0.5)'"
echo -e "PCs to be taken into account for dimensionality reduction: 'c(30)'\n\n"
echo -e "Mean cutoff: 0.01"
echo -e "Should likely low quality cells be filtered out? TRUE"
echo -e "Absolute path to markers file (or NULL value): /home/vfajardo/scripts/seurat_analysis/general_data/Myeloid-cell_markers_ensembl.1.0.RData"
echo -e "### --------------------- Annotations arguments --------------------- ###"
echo -e "Should annotations be added? TRUE (If FALSE, this chunk's variables' values don't really matter)."
echo -e "Absolute path to annotations table: /mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/05-14-2022/aggr/data/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_aggr_table_annotated.1.0.csv"
echo -e "Lane ID: 'species.tag;cell.type.tag;donors.sort.batch.tag;donors.no.tag;tissue.tag'"
echo -e "10X Chromium Batch ID: chrom.batch.tag"
echo -e "Sequencing batch ID: seq.batch.tag"
echo -e "HTO ID: hashtag.tag"
echo -e "Demuxlet ID: "
echo -e "Overlay ID: donor.tag"
echo -e "### ---------------------- Filtering arguments ---------------------- ###"
echo -e "Should there be any kind of preprocess filetring? TRUE (If FALSE, this chunk's variables' values don't really matter)."
echo -e "Tags-related criteria: --TagsCriteria /mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/data_for_subsetting/subset_hto-all/TagsSubsetCriteria.csv"
echo -e "Features-related criteria: "

echo -e "Here starts analysis with R package seurat...\n\n"

Rscript /home/vfajardo/scripts/seurat_analysis/general_seurat_analysis.2.3.5.R --ReportsPath /mnt/hpcscratch/vfajardo/R24/seurat_analysis//R24_Cancer_Myeloid/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B_11-08-2022qc-mye-spc_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL --PrjName R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-B --DataFile /mnt/hpcscratch/vfajardo/R24/seurat_analysis/R24_Cancer_Myeloid/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-A/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-A_11-08-2022qc-std_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL/seurat_objects/SeuratObjectForPrjR24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-A_WithArgs_NoPCs_30.RDS --Raw10x /mnt/hpcscratch/vfajardo/sequencing_data/05-14-2022/aggr/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung/outs/filtered_feature_bc_matrix --InputType seurat  --DoPreSubset TRUE --PreSubsetCriteria /mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/data_for_subsetting/custom_subsets//R24_Cancer_Myeloid//R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-A/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_Subset-A_11-08-2022qc-std_var-30_pc-30_hto-all_harmony-seq.batch.tag_regresscc-NULL/Subset-B/TagsSubsetCriteria.csv --FeatureID ensembl --minCounts 1500 --maxCounts 13000 --minFeatures 350 --maxFeatures 4000 --maxMP 15 --GenNormMethod LogTransform --NormMethod LogNormalize --FVFsMethod vst --FeatsForDSA 30 --PCs 'c(30)' --ForHarmony 'c("seq.batch.tag")' --Resolution 'c(0.1, 0.2, 0.3, 0.4, 0.5)' --IsFivePrime FALSE --MeanCutoff 0.01 --PropCutoff 0.001 --VarsToRegress 'c("nCount_RNA", "percent.mt")'  --FilterOut TRUE --MarkersFile /home/vfajardo/scripts/seurat_analysis/general_data/Myeloid-cell_markers_ensembl.1.0.RData --DEA FALSE --DimReduction 'c("umap")' --RAM 130 --DoAnnotate TRUE --AggrTable /mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/05-14-2022/aggr/data/R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung_aggr_table_annotated.1.0.csv --DonorsMetaData /mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/donors_metadata/GeneralMetaDataProcessedTableForR24CancerBatches-1-to-20.1.0.csv  --MergeBy donor.tag --LaneID 'species.tag;cell.type.tag;donors.sort.batch.tag;donors.no.tag;tissue.tag' --MultAnnPerLane TRUE --ChromBatchID chrom.batch.tag --SeqBatchID seq.batch.tag --HTOID hashtag.tag --OverlayID donor.tag  --DoSubset TRUE --TagsCriteria /mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/data_for_subsetting/subset_hto-all/TagsSubsetCriteria.csv 

echo -e "Job completed!\nCheck for errors if any."
