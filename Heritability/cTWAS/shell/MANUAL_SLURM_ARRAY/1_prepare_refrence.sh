#!/bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=PR_PE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/multi_thread_in_script.out
#SBATCH --error=/multi_thread_in_script.err
#SBATCH --mem=30GB
#SBATCH --time=05:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=###############

########################################################
# Prepare.reference.R -> region info, snp_map, LD_map
########################################################
#  Here to consider - per population and per chormosome 
# region path and LD matrix on a per population basis.

# Explanation - 
# LD structure and region definitions depend on population (EUR, vs AFR, vs ASN)
# Genome build must match the PredictDB weights and GWAS variant positions.
# These files are independent of GWAS or tissue - they are static reference objects.

########################################################
########################################################

/usr/bin/time -v \
Rscript "../../scripts/1_prepare_referece.R" \
    --region_path "/mnt/biohome/jottensmeier/miniforge3/envs/cTWAS/lib/R/library/ctwas/extdata/ldetect/EUR.b38.ldetect.regions.RDS" \
    --genome_version "b38" \
    --ld_dir  "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/data/LD_reference/LD_BioBank/" \
    --outpath "/../../intermediary/"

########################################################
########################################################

