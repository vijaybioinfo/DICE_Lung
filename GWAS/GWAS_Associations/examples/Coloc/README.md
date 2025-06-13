Colocalization Analysis
-------------------------------------
Colocalization analysis infers whether
two phenotypes are likely to be under
the influence of the same causal genetic
variant in a given region. Coloc is one
of the available tools to perform
colocalization analysis.


For further reading see:
https://pubmed.ncbi.nlm.nih.gov/19039033/
https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383
https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008720


Input Files
-------------------
Let's say we have data for 2 traits Asthma GWAS data and eQTL data
and one dataset for each. The analysis requires tsv files with the 
columns listed below. Please check exaple data in data/ folder.

*Data*
1.- Chromosome
2.- Position
3.- SNPID
4.- P
5.- BETA
6.- SE
7.- AF

*Sample tables*
This is a csv file listing the sample/datasets to use in the
analysis (for a given trait) along with some metadata 
such as cc_ratio. Don't forget to provide the path to this 
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
snakemake -np --configfile config.yaml --snakefile ../../R24/Coloc/workflow/Snakefile
#run
snakemake --profile ../cluster_profile --configfile config.yaml --snakefile ../../R24/Coloc/workflow/Snakefile


Output
---------------
Your output folder will contain 3 folders:

1.- raw: Here you can find all tests performed by
the coloc package. Example: all tests between GWAS
loci that overlap with eQTL loci.

2.- filtered: These are filtered results based on
your PPH4 threshold

3.- summary: Contains matrix like files listing 
all filtered associations


