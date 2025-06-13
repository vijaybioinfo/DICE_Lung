GWAS
-------------------------------------
This workflow normalizes GWAS data format
and aligns snps to a reference panel of
choice


Input Files
-------------------
Each file must tab separated and must contain
the required columns listed below. Input format
can be a GWAS summary statistic file or GWAS 
Catalog format.

*Data*
#Required
1.- Chromosome
2.- Position
3.- EA
4.- NEA
5.- P
#Optional
5.- BETA
6.- SE
7.- Z
8.- AF
9.- RSID

*Sample tables*
This is a csv file listing the sample/datasets to use in the
analysis along with some metadata. Don't forget to provide 
the path to this file in your PEP file.

*SNP Database*
This is a binary db file containing one single table
that should be called "Variants" with columns described below.

1.- Chromosome: str
2.- Position: int
3.- REF: str
4.- ALT: str
5,- SNPID: str    #format: 1:100:REF:ALT
6.- RSID: str

*ChainFiles*
Chain files to perform liftover for each genome version
you'd like as output.

Configuration
----------------
To start running you'll need to configure 2 files.

1.- config.yaml: This file contains parameters 
such as input data and output directories.
2.- PEP.yaml: Here you specified the samples 
to use as well as the location of each sample.

How to Run
-------------------
It is highly recomended to run this on a HPC for maximum performance.
Please check cluster_profile/config.yaml as an example for slurm 
configuration.

Once you have created your config.yaml and PEP.yaml 
you can test run the workflow as follows.

#enviroment activation
mamba activate R24;
#dry-run
snakemake -np --configfile config.yaml --snakefile ../../R24/GWAS/workflow/Snakefile
#run
snakemake --profile ../cluster_profile --configfile config.yaml --snakefile ../../R24/GWAS/workflow/Snakefile


Output
---------------
Your output folder will contain 2 folders:

1.- data: Here you can find all input files 
in different genome versions as requested.

2.- sample_tables: folder with original
sample tables. These tables contain meta-information
for each sample. (You can use these for Coloc analysis ).


