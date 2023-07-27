DICE Lung
===========

Repo for our manuscript presenting the DICE tissue project. Check out our manuscript:
- [*Single-cell eQTL analysis of activated T cell subsets reveals activation and cell type–dependent effects of disease-risk variants*](https://www.science.org/doi/10.1126/sciimmunol.abm2508#) (PLACEHOLDER; TO BE REPLACED)

We have splitted the process in seven main sections:
- Variant calling (based on WGS data)
- scRNA-seq data analysis
- scATAC-seq data analysis
- eQTL mapping
- GWAS compilation
- Overlap and colocalization analyses
- GARFIELD

Below we broadly describe these sections and within the directory corresponding to each of them we describe them in extensive detail.

---
# Variant calling (based on WGS data)

Cristian to define.


---
# scRNA-seq data analysis

This process was applied to each predefined immune subset independently.
This section comprises the following major steps:
- Preprocessing of reads by mapping them to the reference genome and collapsing them into UMI count matrices (10x's cellranger)
- Donor deconvolution based on hashtag oligos (HTOs) and based on genetic variants (Demuxlet; depends on the output from the section above).
- Clustering and dimensionality reduction.
- Population definition
- Sex-biased transcript calling

The first step is described in the methods section of our manuscript in detail. Here we provide the code for the rest of the steps.

---
# scATAC-seq data analysis

- Preprocessing of raw data was done entirely with 10x's cellranger atac.
- Clustering, dimensionalty reduction, population definition, pseudobulk aggregation and peak calling was done in R with ArchR package and custom scripts.

---
# eQTL mapping

Cristian to define

---
# GWAS compilation

The process of GWAS retrieval is divided into two main stages:

## Stage 1 : Downloading and pre-processing

The main script for this step is called `preProcess_GWAS.py`. It takes as input a metadata table for a specific disease, which lists all of its GWAS datasets. The script then downloads the raw data from the source database and pre-processes it into a standardized file with relevant columns for further analysis.

The required fields from the sample table are as follows:
- **`Sum stats`**: A value of 0 or 1 to indicate whether is a non-summary or summary statistics file, respectively.
- **`Source`**: The source of the GWAS data. The three possible values for this specific project are *Pan-UK Biobank*, *COVID-19 HGI*, or *NHGRI-EBI Catalog*.
- **`Sample name`:** A personalized ID for the GWAS study, which will be used as a prefix in the output files.
- **`Genome version`:** The corresponding genome build.
- **`Link`**: The URL address for the summary statistics file, or the local download path for non-summary statistics files.
- **`Messy dataset:`** A flag column, 1 to jump specific datasets.
- **`populations`**: Listing ancestry populations for each GWAS cohort.

The output of the `preProcess_GWAS.py` script is two files:
- *Raw file*: raw download from the source database.
- *Pre-processed file*: standardized file with relevant columns for further analysis.

This script can be run with the following command:
```bash
python3 preProcess_GWAS.py --disease [name_of_disease(comma separated)]
```

In addition to the pre-processing step, a quality control (QC) script `QC_preProcess_GWAS.py`, is also included to double-check the integrity of the GWAS data. This script will check the following:
- The consistency of the data types.
- The presence of missing values.
- The consistency of values for summary statistics metrics.

The QC script will output a report that summarizes the results of the checks.

This script can be run with the following command:
```bash
python3 QC_preProcess_GWAS.py --disease [name_of_disease(comma separated)]
```

## Stage 2 : Standardization

Job to define.

---
# Overlap and colocalization analyses

Job to define.

---
# GARFIELD

Eli to define.
