#!/bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=P_We
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=/logging/4_process_weights.out
#SBATCH --error=/logging/4_process_weights.err
#SBATCH --mem=64GB
#SBATCH --time=05:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=###############

########################################################
# Preprocess_weights.R -> extract weights 
########################################################
# Here to consider - to be run per cell-type (and potentially per population)
# otherwise , the union of all SNPs included per populations so that this is run only once per cell-type.

# Explanation
# Tissue specific weights (liver vs adipose) are static.
# Weights do not depend on phenotype or GWAS, only onwhich SNPs exist in the reference.
# Important note
# Should provide gwas_snp_ids as the union of all GWAS SNPs that are planned to be analysed,
# not only a single GWAS study, so that the weight preprocessing is future-proof.
# This avoids recomputing weights for every new GWAS.

# Output -> A processed weight object per tissue (or QTL type).

########################################################
########################################################

Rscript "../../scripts/4_process_weights.R" \
    --weight_files "../../scripts/eqtl_weights_no_model.yaml" \
    --region_path "../../data/INPUT/HAPLOBLOCK_REF/ukb_b38_0.1_chrom_all.rds" \
    --snp_map_path "/../../intermediary/snp_map/ukb_b38_0.1_chrom_all.rds" \
    --gwas_snpid_path "../../intermediary/VCF/ebi-a-GCST004131.RDS" \
    --outpath "../../intermediary/" \
    --LD_map_path "../../LD_map/ukb_b38_0.1_chrom_all.rds" \
    --ncores 2 

# Or just provide one model manually 

Rscript "../../scripts/4_process_weights.R" \
    --weight_files "../../scripts/eqtl_weights_no_model.yaml" \
    --region_path "../../data/INPUT/HAPLOBLOCK_REF/ukb_b38_0.1_chrom_all.rds" \
    --snp_map_path "/../../intermediary/snp_map/ukb_b38_0.1_chrom_all.rds" \
    --gwas_snpid_path "../../intermediary/VCF/ebi-a-GCST004131.RDS" \
    --outpath "../../intermediary/" \
    --LD_map_path "../../LD_map/ukb_b38_0.1_chrom_all.rds" \
    --context "${context}" \
    --prediction_model "${prediction_model}" \
    --type "${qtl_type}" \
    --ncores 2

########################################################
########################################################
