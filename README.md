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

Angel to define.

---
# eQTL mapping

Cristian to define

---
# GWAS compilation

The process of GWAS retrieval supposed two main stages:

## Stage 1 : Downloading and pre-processing

The main script for this step is called `preProcess_GWAS.py`, receiving as input the metadata table for a specific disease listing all its GWAS datasets and giving as output two files:

- *Raw file*: raw download from source database
- *Pre-processed file*: standardized file with relevant columns for further analysis

```bash
python3 preProcess_GWAS.py --disease [name_of_disease(comma separated)]
```

In addition, a QC script is included to double check GWAS integrity `QC_preProcess_GWAS.py`

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
