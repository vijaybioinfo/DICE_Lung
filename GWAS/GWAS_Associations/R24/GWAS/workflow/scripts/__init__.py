import os
import pandas as pd

from peppy import Project
from subprocess import run
from tempfile import gettempdir
from R24.LiftOver.chain_map import ChainFiles
from R24.GWAS.module.utils.misc import get_filename

class GWAS:
    def __init__(self, pep: Project):
        self.pep = pep
        
        table = pd.concat([pd.DataFrame(i.to_dict(), index=range(len(i.file))) for i in self.pep.samples])
        table = table.drop_duplicates().reset_index(drop=False)

        if "Messy dataset" in table.columns:
            table = table[~table["Messy dataset"].astype(str).str.replace("\.0", "").astype(int).astype(bool)]

        self._sample_table = table

    def Diseases(self, unique: bool=False):
        if unique:
            return self._sample_table["disease"].unique().tolist()
        return self._sample_table["disease"].tolist()

    def Datasets(self, unique: bool=False):
        if unique:
            return self._sample_table["sample_name"].unique().tolist()
        return self._sample_table["sample_name"].tolist()
    
    def subset(self, column: str, value) -> pd.DataFrame:
        return self._sample_table[self._sample_table[column] == value]

    def get_datasets(self, wildcards):
        return self.subset("disease", wildcards.disease)["sample_name"].unique().tolist()          #type: ignore

    def _select_row(self, wildcards):
        selected = self.subset("disease", wildcards.disease)
        return selected[selected["sample_name"] == wildcards.dataset]

    def data_genome_versions(self, wildcards):
        return list(self._select_row(wildcards)["Genome version"].values)

    def file(self, wildcards):
        return list(self._select_row(wildcards)["file"].values)

    def get_source_files(self, wildcards):
        files = []
        selected = self.subset("disease", wildcards.disease)
        
        for item in self.get_datasets(wildcards):
            tmp = selected[selected["sample_name"] == item]
            versions = tmp["Genome version"].tolist()

            if wildcards.version in versions:
                file = list(tmp.loc[tmp["Genome version"] == wildcards.version, "file"].values)
            else:
                file = list(tmp["file"].values)

            if len(file) != 1:
                raise Exception("None of more than one source file found. This is a problem with sample table specified...")

            files.extend(file)
        
        return files


    def sumstats(self, wildcards):
        return list(map(float, self._select_row(wildcards)["Sum stats"].values))
    
    @property
    def sample_table_file(self) -> str:
        return self.pep._config.sample_table

        
class CONFIGURATION:
    def __init__(self, config: dict):
        self.config = config

    def output_genome_versions(self) -> list[str]:
        return self.config["Genome-Versions"]

    @property
    def snp_database(self):
        return self.config["SNP-Database"]

    @property
    def chain_maps(self):
        return self.config["ChainFiles"]
    
    @property
    def OutDir(self):
        return self.config["OutputDirectory"]

    def get_alternative_versions(self, wildcards):
        maps = ChainFiles.from_names(self.chain_maps)
        alternatives = maps.find_all_paths_to(wildcards.version)
        return alternatives
        

######### Functions ###########
def get_file_size(path: str):
    return round(os.stat(path).st_size / 1024**3, 1)

def get_mem_gb(path, attempt, rate: float = 30.0):
    mem = rate * get_file_size(path) + (10 * (attempt - 1))
    mem = int(max(20, mem))
    return mem

def clean_logs(path: str):
    if os.path.exists(path) and os.listdir(path):
        run(f"rm -r {path}/*", shell=True)

def temp_dir(path: str):
    if path is None:
        user = run("echo `whoami`", capture_output=True, shell=True).stdout
        user = user.decode("UTF-8").strip()
        path = os.path.join(gettempdir(), user)
    
    return os.path.join(path, "GWAS")

def make_pep_project(config: dict):
    return Project(config["PEP"]["pepfile"], amendments=config["PEP"].get("amend"))

def initialize(config: dict):
    gwas = GWAS(make_pep_project(config))
    cfg = CONFIGURATION(config)
    return cfg, gwas

def file_exists(file: str) -> bool:
    return os.path.exists(file)