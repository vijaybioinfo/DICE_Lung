import os
import numpy as np 
import pandas as pd

from enum import Enum
from collections import namedtuple
from R24.GWAS.module.utils.snpdb import ReferenceSNP

class Cols(Enum):
    SNP_COORD = ["Position", "SNP_COORD"]
    GENE_COORD = ["Gene_Position", "GENE_COORD"]

############################################
############# Helper Functions #############
############################################
def read_file(file: str):
    return pd.read_csv(file, sep="\t", chunksize=10_000_000)

def read_annot(file: str):
    read = pd.read_csv if file.endswith(".csv") else pd.read_table
    annot = read(file)

    return annot.set_index("Geneid").to_dict()

def gene_coord(start, end):
    return int(start + abs((end - start) / 2))

def column_exists(data: pd.DataFrame, options: list[str]):
    columns = data.columns

    for opt in options.value:
        
        if opt in columns:
            data.rename(columns={opt: options.name}, inplace=True)
            return True

    return False


##############################################
############# Prepare MESC input #############
##############################################
def prepare_mesc(chunk: pd.DataFrame, annot: dict[str, str], db, sample_size):
    #### Check SNP and GENE Coordinates ####
    if not column_exists(chunk, Cols.SNP_COORD):
        chunk["SNP_COORD"] = chunk.SNPID.str.split(":").str[1]

    if not column_exists(chunk, Cols.GENE_COORD):
        chunk = chunk[chunk.GENEID.isin(list(annot["gene_name"].keys()))]
        chunk["GENE_COORD"] = [gene_coord(annot["Start"][gene], annot["End"][gene]) for gene in chunk["GENEID"]]

    #### Sample Size and Chromosome ####
    chunk["N"] = sample_size
    chunk["Chromosome"] = chunk["Chromosome"].replace("X", 23)

    #### Translate SNPID to RSID ####
    reference = db.search_snpid(chunk["SNPID"].unique().tolist())
    translate = reference.set_index("SNPID").to_dict()["RSID"]
    chunk["SNP"] = [translate.get(snp, snp) for snp in chunk.SNPID]

    #### Rearrange Columns ####
    chunk = chunk[["GENEID", "GENE_COORD", "SNP", "Chromosome", "SNP_COORD", "N", "Z"]]
    chunk.columns = ["GENE", "GENE_COORD", "SNP", "CHR", "SNP_COORD", "N", "Z"]
    
    return chunk

#########################################
############# MAIN Function #############
#########################################
def main(args):
    chunks = read_file(args.file)
    annot = read_annot(args.annot)
    db = ReferenceSNP(args.db)

    chromosomes = {}
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    for chunk in chunks:
        chunk = prepare_mesc(chunk, annot, db, args.N)
        
        for chrom in chunk.CHR.unique():
            if not int(chrom) in chromosomes:
                chromosomes[int(chrom)] = []
            chromosomes[int(chrom)].append(chunk[chunk.CHR == chrom])

    for chrom in sorted(chromosomes.keys()):
        if chrom == 23: continue
        data = pd.concat(chromosomes[chrom])
        data.sort_values(["GENE_COORD", "GENE", "SNP_COORD"], inplace=True)
        data.to_csv(f"{args.out}.{chrom}.tsv", index=False, sep="\t")


########################################################
################### Argument Parsing ###################
########################################################
def cli():
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--file", required=True, type=str, help="eQTL summary statistics file")
    parser.add_argument("--annot", required=True, type=str, help="Gene annotation file")
    parser.add_argument("--db", required=True, type=str, help="DBSNP database file")
    parser.add_argument("--out", required=True, type=str, help="Output prefix")

    args = parser.parse_args()

    return args

def from_snakemake(snake):
    args = namedtuple("args", ["file", "annot", "db", "N", "out"])
    
    args.file = snake.input.file
    args.annot = snake.input.annot
    args.db = snake.input.db
    args.N = snake.params.N
    args.out = snake.params.out

    return args



if __name__ == "__main__":
    if "snakemake" in locals():
        args = from_snakemake(snakemake)
    else:
        args = cli()
    
    main(args)