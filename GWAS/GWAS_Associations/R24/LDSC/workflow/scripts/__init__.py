import os
import pandas as pd

from peppy import Project
from subprocess import run
from tempfile import gettempdir

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

    def Datasets(self, unique: bool = False):
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

    def file(self, wildcards):
        return super().file(wildcards.disease, wildcards.dataset)

    def Diseases(self, unique: bool=False):
        return self.Projects(unique)

    def get_datasets(self, wildcards) -> list[str]:
        return super().get_datasets(wildcards.disease)

    def Categories(self, unique=False): return self.get_column("Category", unique)
    
    def get_by_category(self, wildcards):
        selected = self.subset(self.PROJECT, wildcards.disease)
        return selected[selected["Category"] == wildcards.category]["sample_name"]
    
    def N(self, wildcards):
        selected = self.subset(self.PROJECT, wildcards.disease)
        return selected.loc[selected["sample_name"] == wildcards.dataset, "sample_size"].values[0]


class InternalProject(DATA):
    def get_project_datasets(self, project: str) -> list[str]:
        return self.get_datasets(project)

    def file(self, wildcards):
        return super().file(wildcards.project, wildcards.cell)


class BASE:
    def __init__(self, pep: Project) -> None:
        self.pep = pep

    def _make_sample_table(self):
        table = pd.concat([pd.DataFrame(i.to_dict(), index=[0]) for i in self.pep.samples])
        table = table.drop_duplicates().reset_index(drop=False)
        self._sample_table = table

    @property
    def sample_table(self):
        if not hasattr(self, "_sample_table"):
            self._make_sample_table()
        return self._sample_table

    def subset(self, chromosome):
        return self.sample_table[self.sample_table.sample_name == chromosome]

    @property
    def Chromosomes(self):
        return self.sample_table["sample_name"].tolist()

    def bfile(self, wildcards): return self.subset(wildcards.chr)["bfile"].values[0]

    def get_plink_files(self, wildcards):
        prefix = self.bfile(wildcards)
        return [f"{prefix}.{ext}" for ext in ["bed", "bim", "fam"]]

    def get_bim_file(self, wildcards): return f"{self.bfile(wildcards)}.bim"

    @property
    def weights(self): 
        return [self.subset(chr)["weights"].values[0] for chr in self.Chromosomes]

    @property
    def frequencies(self): 
        return [self.subset(chr)["frequencies"].values[0] for chr in self.Chromosomes]

    @property
    def has_base_annotation(self): return "model" in self.sample_table.columns

    def annotation(self, wildcards):
        if self.has_base_annotation: 
            return self.subset(wildcards.chr)["model"].values[0]

    def snps(self, wildcards): 
        if "snps" in self.sample_table.columns:
            return self.subset(wildcards.chr)["snps"].values[0]
        return None

################ Functions #################
def clean_logs(path: str):
    if os.path.exists(path) and os.listdir(path):
        run(f"rm -r {path}/*", shell=True)

def temp_dir(path: str):
    if path is None:
        user = run("echo `whoami`", capture_output=True, shell=True).stdout
        user = user.decode("UTF-8").strip()
        path = os.path.join(gettempdir(), user)
    
    return os.path.join(path, "LDSC")

def make_pep_project(config: dict, key: str):
    return Project(config[key]["pepfile"], amendments=config[key].get("amend"))

def initialize(config: dict):
    PROJECT = InternalProject(make_pep_project(config, "PROJECT_PEP"))
    GS = GWAS(make_pep_project(config, "GWAS_PEP"))
    BS = BASE(make_pep_project(config, "BASE_PEP"))

    return BS, GS, PROJECT