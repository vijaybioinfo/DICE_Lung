Identifying potential causal genes relevant in disease.
--------------------------------------------------------
Developer: Job H. Rocha
La Jolla Institute for Immunology

La Jolla, CA 92037, USA
*************************
This repository contains snakemake-based gene-prioritization 
analyses to identify potential causal genes relevant in disease.

Workflows:
1.- GWAS: This workflow attempts to clean GWAS data to keep a 
consistent data format and aligns alleles to a reference panel.
2.- Overlap: Perfoms overlap between GWAS snps and QTL data.
3.- Coloc: Powered by coloc packages this workflow performs 
colocalization analysis between two traits.


Description
----------------
In this repository you'll find 2 folders and 1 yaml file:

1.- examples: contains example data and configuration files
needed to run each workflow.

2.- R24: in this folder you can find all scripts and snakefiles
for all three workflows.

3.- R24.yaml: yaml file needed to set up the main envioroment 
needed by the workflows.


*Overlap analysis*
This analysis consists in finding overlaping snps that
are both QTLs and GWAS-associated variants or in LD with
a GWAS-associated variant.

Steps:
1.- Find significant GWAS snps (P<5e-08).
2.- Get LD snps for every GWAS significant snps.
3.- Overlap GWAS significant + LD snps with QTL.

To see detailed information for this analysis go to
R24/Overlap folder and check README.md.


*Colocalization analysis*
This analysis identifies colocalizing signals by asking 
wheter two different traits (e.g GWAS trait and eQTLs) 
share causal genetic variants in user-defined genomic 
regions. The tool employed for this analysis is coloc 
(described in https://github.com/chr1swallace/coloc).

Steps:
1.- Run coloc.abf function using user-defined GWAS loci 
and QTL information.
2.- Summarize and filter results using user-defined PPH4.

For detailed description and go to R24/Coloc folder
and check README.md file.


*GWAS*
This workflow attempts to standarize GWAS data format 
such as column names and data types, and aligns 
reference and alternativealleles based on a reference 
panel.

Steps:
1.- Normalize column names.
2.- Retrive missing coordinate when possible.
3.- Align alleles to reference panel.
4.- Adjust statistics such as BETA and Z if snp was 
adjusted to match reference panel.

For detailed description and go to R24/GWAS folder and 
check README.md file.


Set up Enviroment
---------------------
You'll need to install mamba to create the R24 enviroment.
See https://mamba.readthedocs.io/en/latest/installation.html.

Once mamba is install in your system follow the next steps.

1.- Clone this repository

2.- cd GWAS_Associations;

3.- mamba env create -f R24.yaml

4.- pip install -e ./


Now everything is ready to run the workflows.


Contact
--------------
Please email to jrocha@lji.org for any questions.
