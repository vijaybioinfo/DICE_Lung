#!/bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=PGWAS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logging/harmonised_gwas_zsnps/PGWAS.out
#SBATCH --error=logging/harmonised_gwas_zsnps/PGWAS.err
#SBATCH --mem=30GB
#SBATCH --time=05:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=###############

########################################################
# Process.GWAS.R -> GWAS Z-scores
########################################################
# Here to consider - to be run per GWAS study (and potentiall population - 
# otherwise , the union of all SNPs included per pipulations so that this is run only once per study.

# Explanation
# Each GWAS has its own summary statistics and sample size.
# Allele harmonisation depends on GWAS SNP IDs vs reference SNPs.
# Independent of tissues.

# Output -> z_snp (Z-scores harmonised with reference).
# Keep a unique output per GWAS so multuple GWAS can resuse the reference.

# Process per population - process per

GWAS_Path="../../data/INPUT/VCF/ebi-a-GCST004131.vcf.gz"
/usr/bin/time -v \
Rscript "../../scripts/2_Process_GWAS.R" \
        --GWAS_Path $GWAS_Path \
        --SNP_map_path "/../../intermediary/snp_map/ukb_b38_0.1_chrom_all.rds" \
        --genome_version "b38" \
        --outpath "/../../intermediary/" \

########################################################
########################################################
