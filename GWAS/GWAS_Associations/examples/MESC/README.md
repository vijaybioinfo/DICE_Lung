Heritability Analysis
-------------------------------------
MESC (mediated expression score regression) estimates the proportion 
of trait heritability that is mediated through gene expression. This 
method breaks down SNP-based heritability to show how much is driven 
by gene expression, helping pinpoint relevant cell types in disease.


For further reading see:
https://www.nature.com/articles/s41588-020-0625-2


Input Files
-------------------
Let's say we have data for 2 traits Asthma GWAS data and eQTL data
and one dataset for each. The analysis requires tsv files with the 
columns listed below. Please check exaple data in data/ folder.

*GWAS Data*
1.- Chromosome
2.- Position
3.- SNPID
4.- P
5.- BETA
6.- SE
7.- AF
8.- Z

*eQTL Data*
1.- Chromosome
2.- GENEID
3.- Position
4.- SNPID
5.- P
6.- BETA
7.- SE
8.- AF
9.- Z

*Sample tables*
This is a csv file listing the sample/datasets to use in the
analysis (for a given trait) along with some metadata 
such as sample_size. Don't forget to provide the path to this 
file in your PEP files.

Configuration
----------------
To start running you'll need to configure 3 files.

1.- config.yaml: This file contains parameters 
such as input data and output directories.
2.- PEP_GWAS.yaml: Here you specified the samples 
to use as well as the location of each sample.
3.- PEP_QTL.yaml: Same as PEP_GWAS.yaml but for
your QTL data.

How to Run
-------------------
It is highly recomended to run this on a HPC for maximum performance.
Please check cluster_profile/config.yaml as an example for slurm 
configuration.

Once you have created your config.yaml, PEP_GWAS.yaml
and PEP_QTL.yaml you can test run the workflow as 
follows.

#enviroment activation
mamba activate R24;
#dry-run
snakemake -np --configfile config.yaml --snakefile ../../R24/MESC/workflow/Snakefile
#run
snakemake --profile ../cluster_profile --configfile config.yaml --snakefile ../../R24/MESC/workflow/Snakefile


Output
---------------
Your output folder will contain 3 folders:

1.- heritability: This folder contains a
single file for every gwas-cell type 
comparison listing the heritability 
estimates from mesc.

2.- filtered: These are filtered results based on
heritability estimates between 0 and 1

3.- summary: Contains all gwas-cell type pairs 
listing the heritability estimates from MESC 
per disease.


