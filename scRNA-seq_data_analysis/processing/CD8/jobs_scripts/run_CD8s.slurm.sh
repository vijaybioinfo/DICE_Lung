#!/bin/bash
#SBATCH --job-name=DICET_CD8s
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5g
#SBATCH --time=120:00:00
#SBATCH --output=/mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/github_repo/DICE_Lung/scRNA-seq_data_analysis/processing/CD8/jobs_scripts/run_CD8s.out.txt
#SBATCH --error=/mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/github_repo/DICE_Lung/scRNA-seq_data_analysis/processing/CD8/jobs_scripts/run_CD8s.err.txt

source activate snakemake
cd /mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/github_repo/DICE_Lung/scRNA-seq_data_analysis/processing/CD8

WORKDIR=/mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/github_repo/DICE_Lung/scRNA-seq_data_analysis/processing/CD8
LOGS_PATH=${WORKDIR}/logs/

snakemake -p --jobs 10 --latency-wait 60 --snakefile ${WORKDIR}/jobs_scripts/Snakefile --configfile ${WORKDIR}/configs/config.yaml --cluster-config ${WORKDIR}/configs/cluster.json --cluster "sbatch --time={cluster.walltime} --nodes=1 --ntasks=1 --cpus-per-task={cluster.cores} --mem={cluster.memory} -e ${LOGS_PATH}/{rule}.{jobid}.err.txt -o ${LOGS_PATH}/{rule}.{jobid}.out.txt --export ALL --parsable" --stats ${LOGS_PATH}/snakemake.stats >& ${LOGS_PATH}/snakemake.log --rerun-incomplete --rerun-triggers mtime --cluster-status /mnt/bioadhoc-temp/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/pipelines/pipeline_DICE/smk-slurm.sh