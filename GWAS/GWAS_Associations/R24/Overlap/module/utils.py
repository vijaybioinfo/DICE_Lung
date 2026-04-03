import os
import re
import numpy as np
import pandas as pd
import datatable as dt

from subprocess import run
from typing import Any, Callable
from functools import partial

from R24.LD import LDFactory
from R24.GWAS.module.utils.columns import Columns as GWASColumns

class Columns(GWASColumns):
    GENENAME = "GENENAME"
    GENEID = "GENEID"
    FDR = "FDR"
    CELL = "CELL"

class LDMethods:
    @staticmethod
    def plink(args: dict):
        return LDFactory.build("PLINK", args)
    
shell: Callable = partial(run, shell=True)

def read_file(file: str):
    data = dt.fread(file).to_pandas()

    if data.empty:
        return pd.DataFrame(columns=data.columns)

    if Columns.CHR in data.columns:    
        data[Columns.CHR] = data[Columns.CHR].astype(str).str.replace("chr", "")

        data[Columns.CHR].replace({"True": "1", 
                                   True: "1", 
                                   "23": "X", 
                                   23: "X", 
                                   "24": "Y", 
                                   24: "Y"},
                                   inplace=True
                                )

        data[Columns.CHR] = data[Columns.CHR].str.replace("\.0", "")
        data[Columns.CHR] = data[Columns.CHR].replace("", np.nan)

    data = data.replace({"True": 1, True: 1})
    
    if Columns.POS in data.columns:
        data[Columns.POS] = pd.to_numeric(data[Columns.POS])

    return data

def merge_files(files: list[str]):
    combined = []
    for file in files:
        #data = read_file(file)
        data = pd.read_table(file, dtype={Columns.CELL: "str"})
        
        data["dataset"] = trim_file_extension(get_filename(file))
        combined.append(data)

    combined = pd.concat(combined).reset_index(drop=True)
    return combined

def make_matrix(data: pd.DataFrame, metadata: pd.DataFrame = None):
    top = data.loc[data.groupby([Columns.GENEID, Columns.GENENAME, Columns.CELL])[Columns.PVAL].idxmin()]

    matrix = pd.pivot_table(top, 
                            index=[Columns.GENEID, Columns.GENENAME], 
                            columns=Columns.CELL, 
                            values=[Columns.PVAL])
    
    matrix.columns = matrix.columns.get_level_values(1)

    if not metadata is None:
        for sample in metadata["sample_name"]:
            if not sample in matrix.columns:
                matrix[sample] = np.nan

        order = metadata.sort_values("order").sample_name.tolist()
        matrix = matrix[order]

    return matrix

def most_significant(over: pd.DataFrame, by: str = Columns.GENEID) -> pd.DataFrame:
    return over.loc[over.groupby(by)[Columns.PVAL].idxmin()]

def parent_directory(path: str) -> str:
    return os.path.abspath(os.path.join(path, os.pardir))

def get_filename(path: str) -> str:
    return os.path.basename(path)

def trim_file_extension(path: str) -> str:
    return os.path.splitext(path)[0]

def mkdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def get_rsids(snps: list[str], db) -> dict[str, str]:
    ds = {}

    for s in snps:
        C = s.split(":")[0]

        if not C in ds:
            ds[C] = []

        ds[C].append(s)

    ids = pd.concat( [ db.search_snpid( list( set(ls) ), [Columns.SNPID, Columns.RSID] ) for ls in ds.values() ] )

    return ids.set_index(Columns.SNPID).to_dict()[Columns.RSID]

################### INTERNAL USE ONLY ######################
def _parse_file_name(file: str) -> tuple[str, Any]:
    name = os.path.basename(file)
    disease = os.path.splitext(name)[0]

    if category := re.search("\d+$", disease):
        category = category.group()
        disease = re.sub("_\d+$", "", disease)
        return disease, category
        
    return disease, np.nan

def make_master_summary(files: list[str], metadata: pd.DataFrame) -> pd.DataFrame:
    '''This function takes a list of matrix files (output of Overlap workflow) and returns a master table combining all diseases'''
    def process_file(file: str):
        data = pd.read_table(file)
        disease, category = _parse_file_name(file)

        data.insert(0, "Category", category)
        data.insert(0, "Disease", disease)

        n = metadata[(metadata.disease == disease) & (metadata.Category == int(category))].shape[0]
        data.insert(6, "Number_of_GWAS_studies", n)

        best = data[data.columns[7:]].max(axis=1)
        n_asso = (data[data.columns[7:]] > 0).sum(axis=1)

        data["Max_P"] = best
        data["Number_of_celltypes"] = n_asso

        return data

    return pd.concat([process_file(file) for file in files]).reset_index(drop=True)

        

def get_files(paths: list[str], metadata: pd.DataFrame):
    files = []
    for disease, dataset in zip(metadata["disease"], metadata["sample_name"]):
        locations = [f"{folder}/{disease}/{dataset}.tsv" for folder in paths if os.path.exists(f"{folder}/{disease}/{dataset}.tsv")]
        
        if not locations:
            raise Exception(f"Dataset {dataset} of Disease {disease} couldn't be found in folder list")
        elif len(locations) > 1:
            raise Exception(f"Dataset {dataset} of Disease {disease} is present in more than one path")
        
        files.append(locations[0])

    return files

    
def create_replication_matrix(data: pd.DataFrame):
    ##### Replication Function #####
    replication = lambda x: pd.DataFrame({"GWAS_datasets": ", ".join(x.dataset.tolist()), "Number_of_GWAS_datasets": x.shape[0]}, index=[0])

    ##### Make Matrix #####
    matrix = make_matrix(data)
    matrix = -np.log10(matrix).fillna(0)
    
    core = matrix.index.to_frame().reset_index(drop=True)
    reps = data.drop_duplicates([Columns.GENEID, "dataset"])

    #### Replication Counts ####
    if reps.shape[0] <= 1:
        reps = pd.DataFrame({Columns.GENEID: reps[Columns.GENEID], "GWAS_datasets": reps["dataset"], "Number_of_GWAS_datasets": 1}, index=[0])
    else:
        reps = reps.groupby(Columns.GENEID).apply(replication).reset_index(level=0)

    ##### Merge matrix and replication counts #####
    joined = pd.merge(core, reps)
    joined = pd.merge(joined, matrix.reset_index())

    return joined


def make_master_table_from_files(files: list[str], samples: list = None):       #type: ignore
    merged = merge_files(files)
    
    #if samples:
    #    merged = merged[merged["CELL"].isin(samples)]
    
    matrix = create_replication_matrix(merged)

    if samples:
        for sample in samples:
            if not sample in matrix.columns:
                matrix[sample] = np.nan
    
        matrix = matrix[list(matrix.columns[:4]) + samples]

    n = len(files)
    matrix.insert(4, "Number_of_GWAS_studies", n)

    best = matrix[matrix.columns[5:]].max(axis=1)
    n_asso = (matrix[matrix.columns[5:]] > 0).sum(axis=1)

    matrix["Max_P"] = best
    matrix["Number_of_celltypes"] = n_asso

    return matrix