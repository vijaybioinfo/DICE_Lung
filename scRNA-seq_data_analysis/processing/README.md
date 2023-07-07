Single-cell RNA-seq analysis per predefined immune subset
===========

We have performed the standard scRNA-seq analysis on a predefined immune subset basis and by applying the [Seurat toolkit](https://satijalab.org/seurat/) 
(shoutout to the Seurat developers and maintainers!).<br>
For this section to be completed, we have developed a workflow that is described in detail in [another repo](https://github.com/VicenteFR/Seurat-based_scRNA-seq_Analysis_v2.3). If 
you wish to reproduce the results from this section, please follow the instructions in that repo to download the scripts and set up the workflow in your own machine.<br>

### **Important disclaimers**
- The results will not be fully reproducible if the versions of the dependencies that were used during the analysis are not the same that you have installed to your machine.
  - Indeed, we don't expect the UMAP plots nor clustering results to be fully reproducible unless the exact same versions are used. Nonetheless, pretty similar results can be
obtained that, importantly, reflect the biology that were captured during our analysis and we present in our manuscript.
  - For you to obtain similar results, you might need to tweak the contamination clusters that are removed during the process (see below for more information).
  - For you to reproduce the exact same results (which we have done 3 times), you must make sure to have the specific versions we detail for the dependencies in the file <code>Dependencies.txt</code>
within this folder.
- The data we report on in our manuscript are huge! Therefore, we required tons of resources to be spent to process these data. Specifically, we ran these analyses in a cluster that
allows for up to 300 GB in RAM (with as many processors as there might be available since it helps speed up the process). We haven't attempted to run this analysis in a local machine
(and you shouldn't either).

---
# System requirements and installation
* R (version $\geq$ v3.6) and Snakemake (we have successfully done it with v7.22.0).
* Please follow the instructions in [this repo](https://github.com/VicenteFR/Seurat-based_scRNA-seq_Analysis_v2.3) to install and set up the workflow. Do install the recommended dependencies.

---
# Data requirements
To obtain the raw data for the results presented in our study, you must download it from [dbGaP]([https://www.ncbi.nlm.nih.gov/gap/](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001703.v4.p1)). 
You must request access authorization and detailed instructions on how to do this are provided [somewhere else](https://dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?page=login). The input for
this section is the aggregated gene by cell UMI count matrix for each of the predefined immune subsets we sorted for, namely: Myeloid cells, NK cells, B cells, CD4 T cells and CD8 T cells.
Please notice that we've pooled NK and B cells in the same libraries, so actually only four aggregated matrices have been provided. The cell type deconvolution for these is applied within this workflow.

---
# Usage
After your workflow is up and running and you have obtained the raw data as detailed in the two sections above, you're ready to run this section!🥳 Also, we assumed that you've already cloned this repo to your own server.
<br>
Let us mark this again. **You must run the analysis for each predefined immune subset.** A snakefile is provided for each immune subset (within their specific folders) and you don't have to worry about modifying it at all. 
Instead, you only have two modify the config files within the <code>configs</code> subfolder, as described below.
#### File <code>config.yaml</code>
* Do replace the entry with key <code>seurat_script</code> by pasting the absolute path to the location where you've saved the master script from the [seurat analysis workflow repo](https://github.com/VicenteFR/Seurat-based_scRNA-seq_Analysis_v2.3).
Notably, at this point, you should have already modified the locations to the module files as explained in there.
* Do replace the entry with key <code>seurat_script</code>. Here you must specify the absolute path to the UMI count matrix for the corresponding cell type.
#### File <code>cluster.json</code>
We don't recommend you make any changes to this file. The process ran without major issues for us with the resources we specify in that file. However, if you'd like to give it a try with
less intense memory or processors, you're welcome to. Then again, we don't guarantee it will work since we have attempted to use as little resources as possible.
<br>
#### Final step
Finally, you can run the process as follows as an example for the CD8 T cell subset:
**Example to run the process is a slurm cluster**
```bash
source activate snakemake
WORKDIR=/path/to/this/repo/DICE_Lung/scRNA-seq_data_analysis/processing/CD8
cd ${WORKDIR}
LOGS_PATH=${WORKDIR}/logs/
snakemake -p --jobs 10 --latency-wait 60 --snakefile ${WORKDIR}/jobs_scripts/Snakefile --configfile ${WORKDIR}/configs/config.yaml --cluster-config ${WORKDIR}/configs/cluster.json --cluster "sbatch --time={cluster.walltime} --nodes=1 --ntasks=1 --cpus-per-task={cluster.cores} --mem={cluster.memory} -e ${LOGS_PATH}/{rule}.{jobid}.err.txt -o ${LOGS_PATH}/{rule}.{jobid}.out.txt --export ALL --parsable" --stats ${LOGS_PATH}/snakemake.stats >& ${LOGS_PATH}/snakemake.log --rerun-incomplete --rerun-triggers mtime --cluster-status /mnt/bioadhoc-temp/Groups/vd-vijay/Cristian/DICE_LungCancer/eQTL_pipeline/pipelines/pipeline_DICE/smk-slurm.sh
```

**Example to run the process is a torque cluster**<br>
This is a room for improvement.
