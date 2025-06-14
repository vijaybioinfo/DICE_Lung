DICE Lung
===========

Repo for our manuscript presenting the DICE tissue project. Check out our manuscript (available soon):
- *Tissue-resident immune cells drive genetic risk in autoimmune and infectious diseases*

We have splitted the process into seven main sections:
- Variant calling (based on WGS data)
- scRNA-seq data analysis
- scATAC-seq data analysis
- eQTL mapping
- GWAS compilation
- Colocalization analysis
- GARFIELD

Below we broadly describe these sections and within the directory corresponding to each of them we describe them in extensive detail.

---
# Variant calling (based on WGS data)

This pipeline takes the fastq files as input and retrieves a VCF file with the variants called for all the donors in our cohort filtered by MAF, HWE, DP, and QC metrics.
- Map reads to the reference genome.
- Mark duplicates and recalibrates base quality scores.
- Called variants and merged all donor's files together.
- Recalibrates and filter variant calls.
- Imputate missing genotype
- Generate VCF file, dosage file, and matrix eQTL-like files.

Check the submodule and methods section in our manuscript for more information. 

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
- Clustering, dimensionality reduction, population definition, pseudobulk aggregation and peak calling was done in R with ArchR package and custom scripts.

---
# eQTL mapping

- Normalize the pseudo-bulk expression data.
- Calculate PCA from genotype.
- Calculate PEER factors from expression data.
- eQTL analysis by MatrixeQTL
- Multiple correction test ([eigenMT](https://github.com/joed3/eigenMT) and permutation-based)

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

This workflow attempts to standardize GWAS data format, such as column names and data types, and aligns reference and alternative alleles based on a reference panel.

For detailed information, go to GWAS/GWAS_Associations and check the readme file.

---
# Colocalization analyses

Colocalization analysis infers whether two phenotypes are likely to be influenced by the same causal genetic variant in a given region. Coloc is one of the available tools for performing colocalization analysis.

For detailed information, go to GWAS/GWAS_Associations and check the readme file.

---
# GARFIELD

GARFIELD-v2.0 software was executed using the default parameters, following the [documentation](https://www.ebi.ac.uk/birney-srv/GARFIELD/documentation-v2/GARFIELD-v2.pdf) from *Valentina Iotchkova, Graham R.S. Ritchie, Matthias Geihs, Sandro Morganella, Josine L. Min, Klaudia Walter, Nicholas J. Timpson, UK10K Consortium, Ian Dunham, Ewan Birney and Nicole Soranzo. GARFIELD - GWAS Analysis of Regulatory or Functional Information Enrichment with LD correction. doi: https://doi.org/10.1101/085738*.

After preparing the inputs, the following lines were executed for each study from each disease:
```bash
./garfield-prep-chr -ptags prunetags -ctags clumptags -maftss maftssd -pval pvalue -ann annot -o output -excl 895,975,976,977,978,979,980
```
```bash
Rscript garfield-Meff-Padj.R -i prepfile -o outfile
```
```bash
Rscript garfield-test.R -i prepfile -o outfile -l linkfile -pt pthreshs -b binning -c condition -s subset
```
```bash
Rscript garfield-plot.R -i test.out -o output_path_prefix -t plot_title -f min \
-padj thresh
```
