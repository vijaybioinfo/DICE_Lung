Identifying potential causal genes relevant in disease.
--------------------------------------------------------
Developer: Job H. Rocha
La Jolla Institute for Immunology

La Jolla, CA 92037, USA
*************************
This repository contains Snakemake-based gene-prioritization 
analyses to identify potential causal genes relevant in disease.

Workflows:

1.- GWAS: This workflow attempts to clean GWAS data to keep a 
consistent data format and aligns alleles to a reference panel.

2.- Coloc: Powered by Coloc package, this workflow performs 
Colocalization analysis between two traits.


Description
----------------
In this repository, you'll find 2 folders and 1 yaml file:

1.- examples: contains example data and configuration files
needed to run each workflow.

2.- R24: In this folder, you can find all scripts and snakefiles
for all three workflows.

3.- envs/R24.yaml: yaml file needed to set up the main environment 
needed by the workflows.

4.- envs/coloc.yaml: yaml file needed to set up the coloc environment 
required by the colocalization workflow.


*Colocalization analysis*
--------------------------
This analysis identifies colocalizing signals by asking 
whether two different traits (e.g GWAS trait and eQTLs) 
share causal genetic variants in the  user-defined genomic 
regions. The tool employed for this analysis is Coloc 
(described in https://github.com/chr1swallace/coloc).

Steps:
1.- Run coloc.abf function using user-defined GWAS loci 
and QTL information.

2.- Summarize and filter results using user-defined PPH4.

For a detailed description, go to the examples/Coloc folder
and check the README.md file.


*GWAS*
-------
This workflow attempts to standardize GWAS data format 
such as column names and data types, and aligns 
reference and alternative alleles based on a reference 
panel.

Steps:

1.- Normalize column names.

2.- Retrieve missing coordinates when possible.

3.- Align alleles to the  reference panel.

4.- Adjust statistics such as BETA and Z if SNP was 
adjusted to match the reference panel.

For a detailed description, go to the examples/GWAS folder and 
Check the README.md file.


Set up Environment
---------------------
You'll need to install mamba to create the R24 environment.
See https://mamba.readthedocs.io/en/latest/installation.html.

Once mamba is installed on your system, follow the next steps.

1.- Clone this repository

2.- cd GWAS_Associations;

3.- mamba env create -f envs/R24.yaml

4.- mamba env create -f envs/coloc.yaml

5.- pip install -e ./


Now everything is ready to run the workflows.


Contact
--------------
Please email jrocha@lji.org for any questions.
