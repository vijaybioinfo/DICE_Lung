**------ CTWAS Wrapper scripts: ------**
These scripts can be run in command line or using snakemake rules.

I have provided both options. 
Please refer to ./shell/MANUAL_SLURM_ARRAY directory for individual bash scripts that set up for a herman cluster.

Alternatively please configure configs for snakemake in ./config directory.


**1_prepare_referece**
Constructs the reference data required for cTWAS by extracting variants from the LD reference panel, 
harmonizing variant identifiers, and computing or subsetting LD matrices for the variants used in the analysis. 
It ensures that SNPIDs, allele orientation, and genomic coordinates are consistent across reference, weights, and GWAS inputs.

**Script arguments**
- **`--region_path`**  
  Path to the genomic region definition file used to partition the genome into LD regions for cTWAS analysis.

- **`--ld_dir`**  
  Directory containing LD reference data (e.g., SNP covariance/correlation matrices) used to compute variant correlation structure.

- **`--genome_version`** *(default: `b38`)*  
  Genome build used for variant coordinates and region definitions (e.g., `b37` or `b38`). (This is a place-holder parameter. Not needed.)

- **`--chromosome`** *(optional)*  
  Integer specifying a chromosome to subset the analysis to a single chromosome instead of processing all regions.

- **`--outpath`**  
  Output directory where processed reference objects (`.RDS` files) will be written for downstream cTWAS steps.

**2_Process_GWAS**
Pre-precess GWAS summary statistics to produce / extract GWAS Z-score from GWAS summary statistics.
This step harmonises GWAS SNPIDs with the LD reference panel, checks allele orientation, removes problematic variants
and formats the GWAS summary statistics for running cTWAS.

**Script arguments**
- **`--GWAS_Path`**  
  Path to the GWAS summary statistics file that will be processed and harmonized for use in the cTWAS pipeline.

- **`--SNP_map_path`**  
  Path to the SNP reference map used to align GWAS variants with the reference panel (e.g., mapping rsIDs or variant IDs to reference coordinates).

- **`--genome_version`** *(default: `b38`)*  
  Genome build used for variant coordinates (e.g., `b37` or `b38`) to ensure consistency between GWAS data and the reference panel.

- **`--chromosome`** *(optional)*  
  Integer specifying a chromosome to subset the analysis to a single chromosome rather than processing the full GWAS dataset.

- **`--outpath`**  
  Directory where the processed GWAS objects (typically `.RDS` files) will be written for downstream cTWAS steps.

- **`--disease`** *(optional)*  
  Label or identifier for the trait/disease being analyzed, typically used for naming output files or organizing results.


**3_process_weights**
Pre-process prediction model weights (e.g. PrediXcan - elastic net or mashR models) by harmonising SNPIDs and ensure that 
LD information is available for each SNP. SNPs are then harmonised with LD panels and GWAS.

**Script arguments**
- **`--weight_files`**  
  Path to gene expression prediction model weight files (e.g., PrediXcan/mashR models) used to construct genetically predicted expression variables.

- **`--region_path`**  
  Path to the genomic region definition file used to assign variants and gene models to LD regions.

- **`--chromosome`** *(optional)*  
  Integer specifying a chromosome to subset the analysis to a single chromosome instead of processing all regions.

- **`--snp_map_path`**  
  Path to the SNP reference map used to harmonize variant identifiers between weight models, GWAS data, and the LD reference.

- **`--gwas_snpid_path`**  
  Path to the file containing SNP identifiers from the processed GWAS dataset, used to ensure consistency between GWAS variants and prediction model SNPs.

- **`--prediction_model`** *(optional)*  
  Identifier for the prediction model being used (e.g., specific tissue or model set), typically used for labeling outputs.

- **`--type`** *(optional)*  
  QTL type associated with the prediction model (e.g., `eQTL`, `sQTL`, `pQTL`).

- **`--context`** *(optional)*  
  Biological context for the prediction model, commonly indicating the tissue or cell type.

- **`--LD_map_path`** *(optional)*  
  Path to a table containing filenames for LD matrices and associated SNP information for each genomic region. Required when LD is loaded externally rather than from PredictDB resources.

- **`--outpath`**  
  Output directory where processed weight objects and associated metadata (typically `.RDS` files) will be written.

- **`--ncores`** *(default: `1`)*  
  Number of CPU cores used for parallel processing during weight preparation.


**4_cTWAS_Runner**
Runs the core cTWAS analysis, integrating GWAS summary statistics, LD reference data, and gene expression prediction weights and applying a fine-mapping utilising a Bayesian statistical approach. 
**This method jointly models SNP and gene effects** to estimate posterior inclusion probabilities (PIPs) and effect size distributions for both variant and gene features.


**Script arguments**
- **`--harmonised_gwas_score`**  
  Path to the harmonized GWAS summary statistics containing standardized association statistics (e.g., Z-scores) aligned to the reference SNP set.

- **`--harmonised_weights`**  
  One or more paths to processed prediction model weight files that have been harmonized with the reference SNP map. Multiple files can be supplied (e.g., for multiple tissues or conditions).

- **`--region_path`**  
  Path to the genomic region definition file used to assign variants and gene models to LD regions for cTWAS analysis.

- **`--snp_map`**  
  Path to the SNP reference map used to align variant identifiers between GWAS data, prediction models, and LD reference panels.

- **`--LD_map`**  
  Path to a table containing filenames and metadata for LD matrices corresponding to each genomic region.

- **`--maxSNP`** *(default: `20000`)*  
  Maximum number of SNPs allowed within a genomic region during analysis to control computational complexity.

- **`--min_group_size`** *(default: `100`)*  
  Minimum number of variants required within a group or region for it to be included in the cTWAS analysis.

- **`--ncore`** *(default: `6`)*  
  Number of CPU cores used for parallel processing.

- **`--genome_version`** *(default: `b38`)*  
  Genome build used for variant coordinates (e.g., `b37` or `b38`).

- **`--chromosome`** *(optional)*  
  Integer specifying a chromosome to subset the analysis to a single chromosome.

- **`--outpath`**  
  Directory where cTWAS result objects (typically `.RDS` files) will be written.

- **`--fname`**  
  Optional filename prefix used when saving output results.


**5_summarise**
Aggregates and formats cTWAS outputs across regions or tissues into analysis-ready tables.
This step typically combines PIPs, effect size summaries, and variance explained metrics (e.g., PVE), producing final result files for downstream interpretation, comparison, and visualization.

**Script arguments**
- **`--finemap_path`**  
  One or more paths to cTWAS fine-mapping result files (e.g., SuSiE output tables) that will be aggregated during the summarization step.

- **`--param_path`**  
  One or more paths to parameter files containing estimated model parameters from the cTWAS runs (e.g., priors, variance estimates).

- **`--snp_map_path`** *(optional)*  
  Path to the SNP reference map used to annotate or harmonize SNP identifiers when generating the final summary tables.

- **`--gwas_n`** *(optional)*  
  GWAS sample size used to compute or rescale downstream metrics (e.g., variance explained or effect size summaries).

- **`--outpath`**  
  Directory where the aggregated summary outputs (e.g., combined tables or `.RDS` files) will be written.
