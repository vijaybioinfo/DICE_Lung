#!/bin/bash

#SBATCH --job-name=cTWAS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=3:00:00
#SBATCH --output=/logging/ctwas_%A_%a.out
#SBATCH --error=/logging/ctwas_ar_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=###############

############# ARRAY CONFIG ###############
#SBATCH --array=0-5%6

start=`date +%s`
##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

echo "------------------------------------------------------"
echo -n "Job is running on node "; $SLURM_JOB_NODELIST
echo "------------------------------------------------------"
echo SLURM: qsub is running on $SLURM_SUBMIT_HOST
echo SLURM: originating queue is $SLURM_JOB_PARTITION
echo SLURM: executing queue is $SLURM_JOB_PARTITION
echo SLURM: working directory is $SLURM_SUBMIT_DIR
echo SLURM: job identifier is $SLURM_JOB_ID
echo SLURM: job name is $SLURM_JOB_NAME
echo SLURM: node file is $SLURM_JOB_NODELIS
echo "------------------------------------------------------"

# The working directory for the job is inside the scratch directory

echo "------------------------------------------------------"
echo -n "Job is running on node "; $SLURM_JOB_NODELIST
echo "------------------------------------------------------"
echo " "
echo " "

###############################################################
#                                                             #
#                        Set Variables                        #
#                                                             #
###############################################################
# Set bash option to nullglob -s to ensure empty glob expantion results
# in an array of length 0 rather than 1.
shopt -s nullglob
weight_files=(/../../intermediary/harmonised_weights/ukb-d-30780_irnt/*.rds)
echo ${#weight_files[@]}
# Unset previous parameter to avoid unexpected behaviours down-stream
shopt -u nullglob
###############################################################
###############################################################
# Check whether the bash array is empty or not.
# If empty -> exit code.
if (( ${#weight_files[@]} == 0 )); then
  echo "No weight files found" >&2
  exit 1
fi
###############################################################
###############################################################
# Fileter files in above dir for files and remove any directories.
files=()
for p in "${weight_files[@]}"; do
  [[ -f "$p" ]] && files+=("$p")
  done

weight_filess=("${files[@]}")
###############################################################
###############################################################
echo "Found ${#weight_files[@]} weight files"
###############################################################
###############################################################

###############################################################
#                                                             #
#                      Main Program                           #
#                                                             #
###############################################################
# weight_files=(*)
# Extract weight files from array based on jobid (0-n)
FILE_PATH="${weight_files[$SLURM_ARRAY_TASK_ID]}"
# FILE_PATH="${weight_files[0]}"

# Slice the file path with / as sep and extract last element (fname).
filename="${FILE_PATH##*/}"
# Remove appendix e.g. .tsv to give clean file name.
FILE_NAME="${filename%.*}"

echo "harmonised_weights file $FILE_PATH"
echo "fname $FILE_NAME"
harmonised_gwas_score_path="../../intermediary/harmonised_weights/ukb-d-30780_irnt/predixcan_mashr_eqtl_Adipose.rds"
harmonised_gwas_score_sname="${harmonised_gwas_score_path##*/}"
harmonised_gwas_score_sname="${harmonised_gwas_score_sname%.*}"
# --snp_map_path "/home/jottensmeier/BioAdHoc/Projects/DICE-LUNG/cTWAS/results/snp_map/ukb_b38_0.1_chrom_all.rds"

###############################################################
###############################################################
outpath="/../../data/OUT/single/raw/$harmonised_gwas_score_sname/"
###############################################################
###############################################################
if [[ ! -d $outpath ]]; then
  mkdir -p "$outpath"
fi
###############################################################
###############################################################


/usr/bin/time -v \
Rscript ../../scripts/5_cTWAS_Runner.R \
  --harmonised_gwas_score "$harmonised_gwas_score_path" \
  --harmonised_weights "$FILE_PATH" \
  --region_path "../../data/INPUT/HAPLOBLOCK_REF/ukb_b38_0.1_chrom_all.rds" \
  --snp_map "/../../intermediary/snp_map/ukb_b38_0.1_chrom_all.rds" \
  --LD_map "../../LD_map/ukb_b38_0.1_chrom_all.rds" \
  --outpath "$outpath" \
  --fname "$FILE_NAME"

########################################################
########################################################
