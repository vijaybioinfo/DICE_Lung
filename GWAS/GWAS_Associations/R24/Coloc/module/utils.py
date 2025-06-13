import os
import numpy as np
import pandas as pd

from typing import Union
from R24.GWAS.module.utils.columns import Columns
from R24.GWAS.module.utils.misc import summarize_loci


def parent_directory(path: str) -> str:
    return os.path.abspath(os.path.join(path, os.pardir))

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

def make_coloc_matrix(data: pd.DataFrame, 
                      metadata: Union[pd.DataFrame, None] = None, 
                      rows: Union[list[str], str] = "Feature2",
                      cols: Union[list[str], str] = "dataset2") -> pd.DataFrame:

    if isinstance(rows, str):
        rows = [rows]

    if isinstance(cols, str):
        cols = [cols]

    top = max_pp4(data, rows + cols)
    matrix = pd.pivot_table(top, index=rows, columns=cols, values="PPH4")

    if not metadata is None and not metadata.empty:
        for sample in metadata["sample_name"]:
            if not sample in matrix.columns:
                matrix[sample] = np.nan

        order = metadata.sort_values("order").sample_name.tolist()
        return matrix[order]
    
    return matrix

def max_pp4(data: pd.DataFrame, groups: Union[list[str], str] = "Feature2"):
    data = data.reset_index(drop=True)
    return data.loc[data.groupby(groups).PPH4.idxmax()]


def merge_coloc_files(files: list[str], column: str):
    merged = []
    for file in files:
        data = pd.read_table(file)
        data[column] = os.path.split(file)[-1][:-4]
        merged.append(data)

    return pd.concat(merged)


################### INTERNAL USE ONLY ######################
def _merge_by_group(path: str, datasets: list[str]):
    _all = []

    for dat in datasets:
        files = []

        for file in os.listdir(f"{path}/{dat}"):
            if file.endswith(".tsv") and not "snps" in file:
                files.append(f"{path}/{dat}/{file}")
        
        merged = merge_coloc_files(files, "dataset2")
        merged["dataset1"] = dat

        _all.append(merged)

    return pd.concat(_all).reset_index(drop=True)

def _matrix_by_group(path: str, datasets: list[str], threshold: float = 0.8):
    merged = _merge_by_group(path, datasets)
    matrix = make_coloc_matrix(merged).reset_index().fillna(0)

    keep = merged.loc[merged.PPH4 >= threshold]
    
    if keep.empty:
        return pd.DataFrame()
    elif keep.shape[0] == 1:
        keep = pd.DataFrame({"Feature2": keep["Feature2"], "GWAS_datasets": keep["dataset1"], "Number_of_GWAS_datasets": [1]})
    else:
        keep = keep.groupby("Feature2").apply(lambda x: pd.DataFrame({
                                                                    "GWAS_datasets": ", ".join(x.dataset1.unique()),
                                                                    "Number_of_GWAS_datasets": x.dataset1.nunique()
                                                                    },
                                                                    index=[0]
                                                                    )
                                            ).reset_index(level=0)
    
    matrix = pd.merge(keep, matrix)
    matrix.insert(3, "Number_of_GWAS_studies", len(datasets))

    best = matrix[matrix.columns[4:]].max(axis=1)
    n_asso = (matrix[matrix.columns[4:]] >= threshold).sum(axis=1)

    matrix["Max_PPH4"] = best
    matrix["Number_of_celltypes"] = n_asso

    return matrix


def make_master_summary(path: str, metadata: pd.DataFrame, genes: dict[str, str] = None, threshold: float = 0.8) -> pd.DataFrame:
    collection = []
    groups = metadata[["disease", "Category"]].drop_duplicates()

    for d, c in zip(groups.disease, groups.Category):
        datasets = metadata[(metadata.disease == d) & (metadata.Category == c)]["Sample name"].tolist()
        
        matrix = _matrix_by_group(f"{path}/{d}", datasets, threshold)
        matrix.insert(0, "Category", c)
        matrix.insert(0, "Disease", d)

        collection.append(matrix)

    collection = pd.concat(collection).rename(columns={"Feature2": "GENEID"}).reset_index(drop=True)

    if not genes is None:
        collection.insert(3, "GENENAME", collection["GENEID"].replace(genes))

    return collection


############# Associate Loci to features ###############
def associate_coloc_to_haploblocks(haploblocks: pd.DataFrame, features: pd.DataFrame, summary_haploblock: bool = True):
    if summary_haploblock:
        haploblocks = summarize_haploblocks(haploblocks)
    
    features = features[["Feature1", "Feature2", "PPH4", "dataset2"]]

    features = pd.pivot_table(features, index=["Feature1", "Feature2"], columns="dataset2", values="PPH4").reset_index().rename(columns={"Feature1": Columns.LOCI})

    return pd.merge(haploblocks, features, how="left")