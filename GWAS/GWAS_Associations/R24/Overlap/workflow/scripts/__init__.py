import os
import pandas as pd
from peppy import Project
from subprocess import run
from tempfile import gettempdir
from R24.Overlap.module.formatting import QTL_TYPES

############################ Interfaces ##############################
class DATA:
    PROJECT = "project"
    DATASET = "dataset"

    def __init__(self, pep: Project) -> None:
        self.pep = pep

    def _make_sample_table(self):
        table = pd.concat([pd.DataFrame(i.to_dict(), index=range(len(i.file))) for i in self.pep.samples])
        table = table.drop_duplicates().reset_index(drop=False)
        self._sample_table = table

    @property
    def sample_table(self):
        if not hasattr(self, "_sample_table"):
            self._make_sample_table()
        return self._sample_table

    def subset(self, column: str, value: str):
        return self.sample_table[self.sample_table[column] == value]

    def get_column(self, column: str, unique=False):
        if unique:
            return self.sample_table[column].drop_duplicates()
        return self.sample_table[column]

    def Projects(self, unique=False):
        return self.get_column(self.PROJECT, unique)

    def Datasets(self, unique: bool = False) -> list[str]:
        return self.get_column("sample_name", unique)
    
    def get_datasets(self, project: str):
        return self.subset(self.PROJECT, project)["sample_name"].tolist()

    def file(self, project: str, dataset: str):
        selected = self.subset(self.PROJECT, project)
        return selected.loc[selected["sample_name"] == dataset, "file"]

    @property
    def get_metadata(self) -> str:
        return self.pep._config.sample_table


class GWAS(DATA):
    PROJECT = "disease"
    CATEGORY = "Category"

    def __init__(self, pep: Project):
        super().__init__(pep)

        if "OVERLAP" in self.sample_table.columns:
            self.sample_table["OVERLAP"].replace({"True": True, "False": False}, inplace=True)
            self._sample_table = self.sample_table[self.sample_table["OVERLAP"]]

    def _has_categories(self):
        return self.CATEGORY in self.sample_table.columns

    def file(self, wildcards):
        return super().file(wildcards.disease, wildcards.dataset)

    def Diseases(self, unique: bool=False) -> list[str]:
        return self.Projects(unique)

    def get_datasets(self, wildcards) -> list[str]:
        return super().get_datasets(wildcards.disease)

    def Categories(self, unique=False):
        return self.get_column(self.CATEGORY, unique) if self._has_categories() else ["1"]
    
    def get_by_category(self, wildcards):
        selected = self.subset(self.PROJECT, wildcards.disease)
        if self._has_categories():
            selected = selected[selected[self.CATEGORY] == wildcards.category]
        return selected["sample_name"]


class QTL(DATA):
    def get_project_datasets(self, project: str) -> list[str]:
        return self.get_datasets(project)

    def file(self, wildcards):
        return super().file(wildcards.project, wildcards.cell)



class CONFIGURATION:
    def __init__(self, config: dict):
        self.config = config

    @property
    def prefix(self):
        return os.path.join(self.config["OutDir"], "Overlap")

    @property
    def LD(self): return self.config.get("LD")
    
    def get_feature(self) -> str:
        return self.config["QTL_PEP"]["Feature"]

    @property
    def type(self): return self.config["QTL_PEP"].get("type", "eQTL")


################ Functions #################
def clean_logs(path: str):
    if os.path.exists(path) and os.listdir(path):
        run(f"rm -r {path}/*", shell=True)

def temp_dir(path: str):
    if path is None:
        user = run("echo `whoami`", capture_output=True, shell=True).stdout
        user = user.decode("UTF-8").strip()
        path = os.path.join(gettempdir(), user)
    
    return os.path.join(path, "Overlap")

def make_pep_project(config: dict, key: str):
    return Project(config[key]["pepfile"], amendments=config[key].get("amend"))

def initialize(config: dict):
    qtl = QTL(make_pep_project(config, "QTL_PEP"))
    GS = GWAS(make_pep_project(config, "GWAS_PEP"))
    CFG = CONFIGURATION(config)

    ### Chekc QTL type ####
    if not CFG.type in QTL_TYPES:
        raise Exception(f"QTL type {CFG.type} is not supported. Supported QTL are [{' ,'.join(QT_TYPES)}]")

    return CFG, GS, qtl