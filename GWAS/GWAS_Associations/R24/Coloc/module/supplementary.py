import os
import re
import numpy as np
import pandas as pd

from enum import Enum
from utils import get_rsids, parent_directory, mkdir
from R24.GWAS.module.utils.snpdb import ReferenceSNP

class Supplementary_Columns(Enum): 
    disease = "Disease"
    disease_category = "Level of disease association"
    dataset = "GWAS dataset name"
    Cell = "Predifined immune cell population"
    Cluster = "Cell cluster"
    locus = "COLOC locus"
    GENEID = "eGene Ensmble ID"
    GENENAME = "eGene ID"
    PPH0 = "PPH0"
    PPH1 = "PPH1"
    PPH2 = "PPH2"
    PPH3 = "PPH3"
    PPH4 = "PPH4"
    PPH4_PPH3 = "Ratio PPH4 / PPH3"
    SNPID = "Colocalized SNP ID"
    Position = "Position of colocalized SNP"
    REF = "REF"
    ALT = "ALT"
    SNP_PPH4 = "SNP PPH4"

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

    def has_annot(self): 
        return hasattr(self, "annotation")
    
    def check(self):
        if not self.is_ready():
            raise Exception("snp database not specified")

    def set_snp_db(self, db: str):
        self.snp_db = ReferenceSNP(db)

    def set_annotation(self, annot: str):
        self.annotation = pd.read_csv(annot, sep="\t", usecols=["Geneid", "gene_name"]).rename(columns={"Geneid": "geneid"})
        self.annotation["geneid"] = self.annotation.geneid.str.replace(r"\.d+", "")
        self.annotation = self.annotation.set_index("geneid").to_dict()["gene_name"]

    def make_table(self, data: pd.DataFrame):
        self.check()

        ###Add cell type info
        cell_name = list(map(match_celltype_name, data.dataset2))
        data["Cell"] = [i[0] for i in cell_name]
        data["Cluster"] = [i[1] for i in cell_name]

        ###Add feature info
        ##### Add comma format to locus ####
        f1 = data["Feature1"].str.split(":")
        f2 = data["Feature2"].str.split("_").str[0]

        chromosome, region = f1.str[0], f1.str[1].str.split("-")
        start, end = region.str[0], region.str[1]

        data["GENEID"] = f2
        data["locus"] = "Chr" + chromosome.str[:] + ":" + start.map(lambda x: format(int(x), ",")) + "-" + end.map(lambda x: format(int(x), ","))
       
        if self.has_annot():
            data["GENENAME"] = [self.annotation.get(g, g) for g in data["GENEID"]]

        ###Add SNP info
        tmp = data["SNPID"].str.split(":")
        data["Position"] = "Chr" + tmp.str[0] + ":" + tmp.str[1].map(lambda x: format(int(x), ","))
        data["REF"] = tmp.str[2]
        data["ALT"] = tmp.str[3]

        snps = get_rsids(data["SNPID"].tolist(), self.snp_db)
        data["SNPID"] = [snps.get(s, s) for s in data["SNPID"]]

        ###rename some columns
        data.rename(columns = {"PPH4/PPH3": "PPH4_PPH3", "SNP.PP.H4": "SNP_PPH4"}, inplace=True)

        ###Rename and adjust columns
        data = data[[col.name for col in Supplementary_Columns]]
        data.columns = [col.value for col in Supplementary_Columns]

        return data


class Read_Data:
    @staticmethod
    def from_path(path: str):
        filename = path.split("/")[-1][:-4]
        disease, category = filename.rsplit("_", 1)

        data = pd.read_table(path)
        data["disease"] = disease
        data["disease_category"] = int(float(category))
        data.rename(columns={"Category": "ld_category", "dataset1": "dataset"}, inplace=True)

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


def supplementary_table_wrapper(metadata: pd.DataFrame, path: str, db: str, annot: str):
    data = Read_Data.from_metadata(metadata, path)
    
    engine = Supplementary_Table()
    engine.set_snp_db(db)

    if not annot is None:
        engine.set_annotation(annot)
    
    supp = engine.make_table(data).reset_index(drop=True)
    supp = supp.loc[supp.groupby(["Disease", "Level of disease association", "GWAS dataset name", "Predifined immune cell population", "Cell cluster", "eGene Ensmble ID"], dropna=False).PPH4.idxmax()]
    
    return supp.sort_values(["Disease", "Level of disease association", "GWAS dataset name", "Predifined immune cell population", "Cell cluster", "eGene Ensmble ID"])


def main(args):
    metadata = pd.read_csv(args.metadata)
    supp = supplementary_table_wrapper(metadata, args.path, args.db, args.annot)
    
    mkdir(parent_directory(args.out))
    supp.to_csv(f"{args.out}.csv", index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--metadata", help="GWAS metadata file")
    parser.add_argument("--annot", required=False, help="Gene annotation file")
    parser.add_argument("--path", help="Coloc results folder(s)")
    parser.add_argument("--db", help="SNP database file")
    parser.add_argument("--out", help="Output prefix")

    args = parser.parse_args()

    main(args)
