import re
import os
import numpy as np
import pandas as pd

from enum import Enum
from R24.GWAS.module.utils.snpdb import ReferenceSNP
from utils import  get_filename, get_rsids, trim_file_extension, parent_directory, mkdir

class Supplementary_Columns(Enum): 
    disease = "Disease"
    disease_category = "Level of disease association"
    dataset = "GWAS dataset name"
    ld = "LD analysis applied to GWAS dataset"
    Cell = "Predifined immune cell population"
    Cluster = "Cell cluster"
    GENEID = "eGene Ensmble ID"
    GENENAME = "eGene ID"
    SNPID = "eQTL SNP ID"
    Position = "Postion of eQTL SNP"
    REF = "REF"
    ALT = "ALT"
    P = "Adj. association P value (FDR)"
    BETA = "Effect size (slope)"
    GWAS_P = "GWAS P value of eQTL SNP"
    GWAS_BETA = "GWAS effect size (slope) of eQTL SNP"
    Lead = "Associated GWAS lead SNP ID"
    GWAS_Lead_Position = "Position of GWAS lead SNP"
    GWAS_Lead_ref = "REF"
    GWAS_Lead_alt = "ALT"
    GWAS_Lead_P = "GWAS P value of lead SNP"
    GWAS_Lead_BETA = "GWAS effect size (slope) of lead SNP"


def match_celltype_name(name: str):
    if re.search(".+_\d+", name):
        return ( 
                re.sub("_\d", "",  name), 
                "_".join(re.findall("_(\d)", name)) 
                )
    
    return (name, np.nan)

class Supplementary_Table:
    COLUMNS = Supplementary_Columns

    def is_ready(self):
        return hasattr(self, "snp_db")
    
    def check(self):
        if not self.is_ready:
            raise Exception("snp database not specified")

    def set_snp_db(self, db: str):
        self.snp_db = ReferenceSNP(db)

    def make_table(self, data: pd.DataFrame):
        self.check()

        ### Check if LD was performed
        data["ld"] = "ld_category" in data.columns

        ### LD associated columns
        if not "ld_category" in data.columns:
            data["ld"] = False
            data["Lead"] = data["SNPID"].tolist()
            data["GWAS_P"] = data["GWAS_Lead_P"]
            data["GWAS_BETA"] = data["GWAS_Lead_BETA"]

        ### Add cell type
        cell_name = list(map(match_celltype_name, data.CELL))
        data["Cell"] = [i[0] for i in cell_name]
        data["Cluster"] = [i[1] for i in cell_name]

        ### Add SNPID info
        SNP = data["SNPID"].str.split(":")
        data["Position"] = "Chr" + SNP.str[0] + ":" + SNP.str[1].map(lambda x: format(int(x), ","))
        data["SNPID"] = data["SNPID"].replace( get_rsids(data["SNPID"].tolist(), self.snp_db) )

        ### Add GWAS SNP info
        LEAD = data["Lead"].str.split(":")
        data["GWAS_Lead_ref"] = LEAD.str[2]
        data["GWAS_Lead_alt"] = LEAD.str[3]

        data["GWAS_Lead_Position"] = "Chr" + LEAD.str[0] + ":" + LEAD.str[1].map(lambda x: format(int(x), ","))
        data["Lead"] = data["Lead"].replace( get_rsids(data["Lead"].tolist(), self.snp_db))                

        ### Rename and adjust columns
        data = data[[col.name for col in Supplementary_Columns]]
        data.columns = [col.value for col in Supplementary_Columns]

        return data

class Read_Data:
    @staticmethod
    def from_path(path: str):
        filename = trim_file_extension(get_filename(path))
        disease, category = filename.rsplit("_", 1)

        data = pd.read_table(path)
        data["disease"] = disease
        data["disease_category"] = int(float(category))
        data.rename(columns={"Category": "ld_category"}, inplace=True)

        return data

    @staticmethod
    def from_metadata(meta: pd.DataFrame, path: str):
        files  = []
        mapper = meta[["Disease", "disease"]].drop_duplicates().set_index("disease").to_dict()["Disease"]

        for D in mapper.keys():
            if os.path.exists(f"{path}/{D}"):
                files.extend([f"{path}/{D}/{i}" for i in os.listdir(f"{path}/{D}") if re.match(f"{D}_\d+.tsv", i) or re.match(f"{D}_\d+.0.tsv", i)])

        data = pd.concat(Read_Data.from_list(files))
        data["disease"] = data["disease"].replace(mapper)

        return data

    @staticmethod
    def from_list(l: list[str]):
        return [Read_Data.from_path(file) for file in l]


def supplementary_table_wrapper(metadata: pd.DataFrame, paths: list[str], db: str):

    if isinstance(paths, str):
        paths = [paths]
    
    data = [Read_Data.from_metadata(metadata, p) for p in paths]
    
    engine = Supplementary_Table()
    engine.set_snp_db(db)
    
    supp = pd.concat([engine.make_table(item) for item in data])

    return supp.sort_values(["Disease", "Level of disease association", "GWAS dataset name", "Predifined immune cell population", "Cell cluster", "eGene Ensmble ID"])


def main(args):
    metadata = pd.read_csv(args.metadata)
    supp = supplementary_table_wrapper(metadata, args.path, args.db)
    
    mkdir(parent_directory(args.out))
    supp.to_csv(f"{args.out}.csv", index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--metadata", help="GWAS metadata file")
    parser.add_argument("--path", nargs="+", help="Overlap results folder(s)")
    parser.add_argument("--db", help="SNP database file")
    parser.add_argument("--out", help="Output prefix")

    args = parser.parse_args()

    main(args)