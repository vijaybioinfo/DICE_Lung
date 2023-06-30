#!/bin/bash
#SBATCH --job-name=DICET_CD4s
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5g
#SBATCH --time=120:00:00
#SBATCH --output=/mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/github_repo/DICE_Lung/scRNA-seq_data_analysis/processing/CD4/jobs_scripts/run_CD4s.out.txt
#SBATCH --error=/mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/github_repo/DICE_Lung/scRNA-seq_data_analysis/processing/CD4/jobs_scripts/run_CD4s.err.txt

conda activate snakemake
cd /mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/github_repo/DICE_Lung/scRNA-seq_data_analysis/processing/CD4

WORKDIR=/mnt/bioadhoc-temp/Groups/vd-vijay/vfajardo/R24/paper_developments/R24_Cancer/github_repo/DICE_Lung/scRNA-seq_data_analysis/processing/CD4
LOGS_PATH=${WORKDIR}/logs/

snakemake --jobs 10 --latency-wait 60 --snakefile ./jobs_scripts/Snakefile --configfile ./configs/config.yaml --cluster-config ./configs/cluster.json --cluster "sbatch --time={cluster.walltime} --nodes=1 --ntasks=1 --cpus-per-task={threads} --mem={cluster.memory} -e ${LOGS_PATH}/logs/{rule}.{jobid}.{wildcards}.err.txt -o ${LOGS_PATH}/{rule}.{jobid}.{wildcards}.out.txt --export ALL --parsable" --stats ${LOGS_PATH}/snakemake.stats >& ${LOGS_PATH}/snakemake.log  --rerun-incomplete --use-conda --cluster-status /mnt/bioadhoc-temp/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/pipelines/pipeline_DICE/smk-slurm.sh --touch --rerun-triggers mtime 1> gen_output.txt