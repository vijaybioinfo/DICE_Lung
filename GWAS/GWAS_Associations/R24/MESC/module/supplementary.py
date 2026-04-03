import os
import re
import numpy as np
import pandas as pd

from enum import Enum
from utils import parent_directory, mkdir

####################################################
#### Define columns for the supplementary table ####
####################################################
class Supplementary_Columns(Enum): 
    disease = "Disease"
    disease_category = "Level of disease association"
    dataset = "GWAS dataset name"
    Cell = "Predifined immune cell population"
    Cluster = "Cell cluster"
    h2 = "Total SNP heritability"
    h2_SE = "SE of total SNP heritability"
    h2med = "Gene expression-mediated heritability"
    h2med_SE = "SE of gene expression-mediated heritability"
    h2med_h2 = "Proportion of heritability mediated by gene expression"
    h2med_h2_SE = "SE of proportion of heritability mediated by gene expression"

#####################################################################
######################### Read Data #################################
#####################################################################
class Read_Data:
    @staticmethod
    def from_path(path: str):
        filename = path.split("/")[-1][:-4]
        disease, category = filename.rsplit("_", 1)

        data = pd.read_table(path)
        data["disease"] = disease
        data["disease_category"] = int(float(category.split(".")[0]))

        return data
    
    @staticmethod
    def from_list(l: list[str]):
        return [Read_Data.from_path(file) for file in l]
    
    @staticmethod
    def from_metadata(meta: pd.DataFrame, path: str):
        collect = []
        mapper = meta[["Disease", "disease"]].drop_duplicates().set_index("disease").to_dict()["Disease"]

        for disease in meta.disease.unique():
            folder = f"{path}/{disease}"
            if not os.path.exists(folder): continue

            ### Collect all files for the disease ###
            files = [f"{folder}/{file}" for file in os.listdir(folder) if file.endswith(".all.h2med.ss") or file.endswith(".all.h2med")]
            data = pd.concat(Read_Data.from_list(files))
            
            ### Subset to the datasets in the metadata ###
            subset = meta[meta.disease == disease]
            data = data[data.dataset.isin(subset.sample_name.tolist())]
            
            collect.append(data)

        data = pd.concat(collect).drop_duplicates().reset_index(drop=True)
        data["disease"] = data["disease"].replace(mapper)

        return data

#################################################################
######################### Order #################################
#################################################################
class NoOrder:
    def order(self, data):
        return data

class CategoryOrder:
    def __init__(self, categories: list[str]):
        self.categories = categories

    def order(self, data: pd.DataFrame):
        raise NotImplementedError("Category order not implemented yet...")

class CellTypeOrder(CategoryOrder):
    def order(self, data: pd.DataFrame):
        data["CELL"] = data["CELL"].astype("category")
        data["CELL"] = data['CELL'].cat.reorder_categories(self.categories, ordered=True)
        data.sort_values(["CELL"], ascending=True, inplace=True)

        return data


def load_sort_method(choice: str, project_metadata):
    if choice == "cell":
        return CellTypeOrder(project_metadata.sample_name.tolist())

    return NoOrder()


def match_celltype_name(name: str):
    if re.search(".+_\d+", name):
        return ( 
                re.sub("_\d", "",  name), 
                "_".join(re.findall("_(\d)", name)) 
                )
    
    return (name, np.nan)
###############################################################################
######################### Supplementary Table #################################
###############################################################################
class Supplementary_Table:
    COLUMNS = Supplementary_Columns
    
    def __init__(self, order_method = None):
        if order_method is None:
            self.order_method = NoOrder()
        else:
            self.order_method = order_method
    
    def make_table(self, data: pd.DataFrame):
        cols_1 = ["disease", "disease_category", "dataset", "CELL", "Estimate", "SE(Estimate)"]
        cols_2 = ["disease", "disease_category", "dataset", "CELL", "Estimate_over_h2", "SE(Estimate_over_h2)"]

        h2 = data.loc[data.Quantity == "h2", cols_1].rename(columns={"Estimate": "h2", "SE(Estimate)": "h2_SE"})
        h2med = data.loc[data.Quantity == "h2med", cols_1].rename(columns={"Estimate": "h2med", "SE(Estimate)": "h2med_SE"})
        proportion = data.loc[data.Quantity == "h2med", cols_2].rename(columns={"Estimate_over_h2": "h2med_h2", "SE(Estimate_over_h2)": "h2med_h2_SE"})
        
        results = pd.merge(h2, h2med).merge(proportion)
        results = self.order_method.order(results)
        
        assert results.shape[0] == h2.shape[0], "number of rows doesn't match..."
        
        ###Add cell type info
        cell_name = list(map(match_celltype_name, results.CELL))
        results["Cell"] = [i[0] for i in cell_name]
        results["Cluster"] = [i[1] for i in cell_name]
        
        ### Sort rows by disease, cell type, and dataset ###
        
        results.sort_values(["disease", "disease_category", "dataset", "CELL"], ascending=True, inplace=True)
        
        ###Rename and adjust columns
        results = results[[col.name for col in self.COLUMNS]]
        results.columns = [col.value for col in self.COLUMNS]

        return results





def supplementary_table_wrapper(metadata: pd.DataFrame, path: str, project_metadata: pd.DataFrame, method: str = None):
    '''Wrapper function to generate the supplementary table for MESC results.'''
    
    ### Load and process data ###
    project_metadata.sort_values("order", ascending=True, inplace=True)
    method = load_sort_method(method, project_metadata)
    
    ### Generate supplementary table ###
    engine = Supplementary_Table(method)
    data = Read_Data.from_metadata(metadata, path)
    supp = engine.make_table(data)
    
    return supp

def main(args):
    ### Load metadata ###
    metadata = pd.read_csv(args.metadata)
    metadata = metadata[metadata["MESC"]]
    project_metadata = pd.read_csv(args.project_metadata)
    
    ### Generate supplementary table ###
    supplementary = supplementary_table_wrapper(metadata, args.path, project_metadata, args.method)
    
    ### Save supplementary table ###
    mkdir(parent_directory(args.out))
    supplementary.to_csv(f"{args.out}.csv", index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--metadata", required=True, help="GWAS metadata file")
    parser.add_argument("--project-metadata", required=True, help="Project metadata file")
    parser.add_argument("--path", required=True, help="Summary results path")
    parser.add_argument("--method", default="cell", choices=["noorder", "cell", "metric"])
    parser.add_argument("--out", help="Output prefix")

    args = parser.parse_args()

    main(args)