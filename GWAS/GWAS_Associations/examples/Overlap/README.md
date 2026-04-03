Overlap Analysis
-------------------------------------
Overlap analysis identifies snps that are associated
to two different traits. For example Asthma GWAS snps
and eQTLs.

Input Files
-------------------
Let's say we have data for 2 traits Asthma GWAS data and eQTL data
and one dataset for each. The analysis requires tsv files with the 
columns listed below. Please check exaple data in data/ folder.

*Data*
1.- SNPID
2.- P
#For eQTLs only
3.- GENEID 
4.- GENENAME
5.- BETA

*Sample tables*
This is a csv file listing the sample/datasets to use in the
analysis (for a given trait) along with some metadata. Don't 
forget to provide the path to this file in your PEP files.

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
snakemake -np --configfile config.yaml --snakefile ../../R24/Overlap/workflow/Snakefile
#run
snakemake --profile ../cluster_profile --configfile config.yaml --snakefile ../../R24/Overlap/workflow/Snakefile


Output
---------------
Your output folder will contain 3 folders:

1.- raw: Here you can find all snps shared between
your eQTLs and GWAS snps. Example: all tests between GWAS
loci that overlap with eQTL loci.

2.- filtered: Top significant snps for each Feature (GENE).

3.- summary: Contains matrix like files listing 
all top associations.


